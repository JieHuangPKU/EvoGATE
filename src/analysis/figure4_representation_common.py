import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import yaml
from sklearn.metrics import davies_bouldin_score, silhouette_score
from sklearn.preprocessing import StandardScaler

from src.data.frozen_protocol_loader import load_protocol_dataset
from src.train.run_frozen_protocol_feature_combo_model import resolve_graph_model_config
from src.train.run_frozen_protocol_model import normalize_model_name


rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["font.family"] = "DejaVu Sans"

UMAP_PARAMS = {
    "n_components": 2,
    "n_neighbors": 15,
    "min_dist": 0.1,
    "metric": "euclidean",
    "random_state": 1029,
}

LABEL_COLORS = {
    0: "#4C78A8",
    1: "#E45756",
}

TRANSITION_COLORS = {
    "TP_stable": "#3B7A57",
    "TN_stable": "#9AA5B1",
    "FN_to_TP_rescued": "#D62728",
    "FP_to_TN_corrected": "#9467BD",
    "FN_persistent": "#FF9896",
    "FP_persistent": "#C5B0D5",
    "TP_to_FN_regressed": "#8C564B",
    "TN_to_FP_regressed": "#BCBD22",
}

TRANSITION_ORDER = [
    "TP_stable",
    "TN_stable",
    "FN_to_TP_rescued",
    "FP_to_TN_corrected",
    "FN_persistent",
    "FP_persistent",
    "TP_to_FN_regressed",
    "TN_to_FP_regressed",
]


def species_title(protocol):
    if protocol == "fgraminearum_newlabel":
        return "Fusarium graminearum (new label)"
    if protocol == "scerevisiae":
        return "Saccharomyces cerevisiae"
    return protocol


def label_regime(protocol):
    return "newlabel" if protocol == "fgraminearum_newlabel" else "standard"


def locate_run_dir(upstream_root, protocol, model, feature_setting, seed):
    candidate = Path(upstream_root) / protocol / model / feature_setting / "run_{0}".format(seed)
    if not (candidate / "best_model.pt").exists():
        raise FileNotFoundError("Missing checkpoint: {0}".format(candidate / "best_model.pt"))
    return candidate.resolve()


def subset_indices(bundle, subset):
    if subset == "test":
        return np.asarray(bundle["test_idx"], dtype=np.int64)
    if subset == "val":
        return np.asarray(bundle["val_idx"], dtype=np.int64)
    if subset == "train":
        return np.asarray(bundle["train_idx"], dtype=np.int64)
    split_manifest = bundle["split_manifest"]
    mapping = bundle["mapping"]
    return np.asarray([mapping[item] for item in split_manifest["graph_gene_id"].astype(str)], dtype=np.int64)


def _mean_neighbor_aggregate(x, edge_index):
    src = edge_index[0]
    dst = edge_index[1]
    num_nodes = x.shape[0]
    aggregated = x.new_zeros(x.shape)
    aggregated.index_add_(0, dst, x[src])
    degree = x.new_zeros((num_nodes, 1))
    degree.index_add_(0, dst, torch.ones((dst.shape[0], 1), dtype=x.dtype, device=x.device))
    degree = degree.clamp(min=1.0)
    return aggregated / degree


def graphsage_penultimate_and_logits(x, edge_index, state_dict, aggregator_type):
    if str(aggregator_type).strip().lower() != "mean":
        raise ValueError("Figure4 utilities currently support mainline mean GraphSAGE only")
    hidden = x
    layer_index = 0
    while "layers.{0}.lin.weight".format(layer_index) in state_dict:
        prefix = "layers.{0}.lin".format(layer_index)
        neighbor_agg = _mean_neighbor_aggregate(hidden, edge_index)
        combined = torch.cat([hidden, neighbor_agg], dim=1)
        hidden = F.linear(combined, state_dict[prefix + ".weight"], state_dict[prefix + ".bias"])
        hidden = torch.relu(hidden)
        layer_index += 1
    neighbor_agg = _mean_neighbor_aggregate(hidden, edge_index)
    combined = torch.cat([hidden, neighbor_agg], dim=1)
    logits = F.linear(combined, state_dict["out_layer.lin.weight"], state_dict["out_layer.lin.bias"])
    return hidden, logits


def standardize(matrix):
    return StandardScaler().fit_transform(matrix)


def compute_umap(matrix, params=None):
    params = dict(UMAP_PARAMS if params is None else params)
    import umap.umap_ as umap

    return umap.UMAP(**params).fit_transform(standardize(matrix))


def centroid_distance(matrix, labels):
    pos = matrix[labels == 1]
    neg = matrix[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float(np.linalg.norm(pos.mean(axis=0) - neg.mean(axis=0)))


def separation_metrics(matrix, labels):
    metrics = {
        "centroid_distance": centroid_distance(matrix, labels),
        "silhouette_score": float("nan"),
        "davies_bouldin_index": float("nan"),
    }
    if len(np.unique(labels)) >= 2 and len(labels) >= 3:
        metrics["silhouette_score"] = float(silhouette_score(matrix, labels))
        metrics["davies_bouldin_index"] = float(davies_bouldin_score(matrix, labels))
    return metrics


def error_transition(gold_label, baseline_pred, esm2_pred):
    gold_label = int(gold_label)
    baseline_pred = int(baseline_pred)
    esm2_pred = int(esm2_pred)
    if gold_label == 1 and baseline_pred == 1 and esm2_pred == 1:
        return "TP_stable"
    if gold_label == 0 and baseline_pred == 0 and esm2_pred == 0:
        return "TN_stable"
    if gold_label == 1 and baseline_pred == 0 and esm2_pred == 1:
        return "FN_to_TP_rescued"
    if gold_label == 0 and baseline_pred == 1 and esm2_pred == 0:
        return "FP_to_TN_corrected"
    if gold_label == 1 and baseline_pred == 0 and esm2_pred == 0:
        return "FN_persistent"
    if gold_label == 0 and baseline_pred == 1 and esm2_pred == 1:
        return "FP_persistent"
    if gold_label == 1 and baseline_pred == 1 and esm2_pred == 0:
        return "TP_to_FN_regressed"
    if gold_label == 0 and baseline_pred == 0 and esm2_pred == 1:
        return "TN_to_FP_regressed"
    return "unclassified"


def fine_scatter(ax, coords, labels=None, transitions=None, title="", by="label", legend=True):
    ax.set_facecolor("white")
    if by == "transition":
        order = [name for name in TRANSITION_ORDER if name in set(transitions)]
        for name in order:
            subset = coords[np.asarray(transitions) == name]
            size = 5.0 if name == "FN_to_TP_rescued" else 2.0
            alpha = 0.9 if name == "FN_to_TP_rescued" else 0.55
            ax.scatter(
                subset[:, 0],
                subset[:, 1],
                s=size,
                alpha=alpha,
                c=TRANSITION_COLORS.get(name, "#7F7F7F"),
                edgecolors="none",
                label=name,
            )
    else:
        for label in [0, 1]:
            subset = coords[np.asarray(labels) == label]
            ax.scatter(
                subset[:, 0],
                subset[:, 1],
                s=1.5 if label == 0 else 3.0,
                alpha=0.55,
                c=LABEL_COLORS[label],
                edgecolors="none",
                label="non-essential" if label == 0 else "essential",
            )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    if legend:
        ax.legend(frameon=False, loc="best", fontsize=7)


def load_predictions(run_dir):
    df = pd.read_csv(run_dir / "predictions.tsv", sep="\t")
    df["graph_gene_id"] = df["graph_gene_id"].astype(str)
    return df


def load_hidden_case(runtime_config_path, upstream_root, protocol, feature_setting, model, seed, subset):
    bundle = load_protocol_dataset(runtime_config_path, protocol, feature_setting)
    run_dir = locate_run_dir(upstream_root, protocol, model, feature_setting, seed)
    base_config = yaml.safe_load(Path(runtime_config_path).read_text(encoding="utf-8"))
    model_cfg = resolve_graph_model_config(base_config, normalize_model_name(model))
    checkpoint = torch.load(str(run_dir / "best_model.pt"), map_location="cpu")
    x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
    edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
    hidden, logits = graphsage_penultimate_and_logits(
        x, edge_index, checkpoint["state_dict"], model_cfg["aggregator_type"]
    )
    selected_idx = subset_indices(bundle, subset)
    subset_df = bundle["node_manifest"].iloc[selected_idx].copy().reset_index(drop=True)
    subset_df["graph_gene_id"] = subset_df["graph_gene_id"].astype(str)
    subset_df["label"] = subset_df["label"].astype(int)
    prediction_df = load_predictions(run_dir)
    prediction_df = prediction_df[prediction_df["split"].astype(str) == subset].copy().reset_index(drop=True)
    merged = subset_df.merge(
        prediction_df[["graph_gene_id", "pred_score", "pred_label"]],
        on="graph_gene_id",
        how="left",
        validate="one_to_one",
    )
    if merged["pred_label"].isna().any():
        raise ValueError("Missing predictions after merge for {0} {1}".format(protocol, feature_setting))
    return {
        "protocol": protocol,
        "species": bundle["species"],
        "label_regime": label_regime(protocol),
        "feature_setting": feature_setting,
        "model": model,
        "seed": seed,
        "subset": subset,
        "bundle": bundle,
        "run_dir": str(run_dir),
        "checkpoint_path": str(run_dir / "best_model.pt"),
        "split_manifest_path": bundle["split_manifest_path"],
        "label_manifest_path": bundle["label_manifest_path"],
        "feature_schema": bundle["feature_schema"].copy(),
        "subset_df": merged,
        "input_matrix": np.asarray(bundle["feature_matrix"][selected_idx], dtype=np.float32),
        "hidden_matrix": hidden.detach().cpu().numpy()[selected_idx].astype(np.float32, copy=False),
        "logits": logits.detach().cpu().numpy()[selected_idx].reshape(-1).astype(np.float32, copy=False),
        "labels": merged["label"].to_numpy(dtype=int),
        "node_ids": merged["graph_gene_id"].astype(str).tolist(),
        "pred_label": merged["pred_label"].astype(int).to_numpy(),
        "pred_score": merged["pred_score"].astype(float).to_numpy(),
    }


def pair_cases(baseline_case, esm2_case):
    a = baseline_case["subset_df"].copy()
    b = esm2_case["subset_df"].copy()
    if a["graph_gene_id"].tolist() != b["graph_gene_id"].tolist():
        raise ValueError("Node order mismatch between paired cases")
    if a["label"].tolist() != b["label"].tolist():
        raise ValueError("Label order mismatch between paired cases")
    paired = a[["graph_gene_id", "label", "split", "pred_label", "pred_score"]].copy()
    paired = paired.rename(columns={"pred_label": "baseline_pred_label", "pred_score": "baseline_pred_score"})
    paired["esm2_pred_label"] = b["pred_label"].astype(int).to_numpy()
    paired["esm2_pred_score"] = b["pred_score"].astype(float).to_numpy()
    paired["transition"] = [
        error_transition(gold, base_pred, esm_pred)
        for gold, base_pred, esm_pred in zip(
            paired["label"].astype(int),
            paired["baseline_pred_label"].astype(int),
            paired["esm2_pred_label"].astype(int),
        )
    ]
    return paired


def save_pdf(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path), format="pdf", dpi=300, transparent=False, bbox_inches="tight")
    plt.close(fig)


def write_json(path, payload):
    def _convert(value):
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {str(k): _convert(v) for k, v in value.items()}
        if isinstance(value, list):
            return [_convert(v) for v in value]
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            return None
        return value
    Path(path).write_text(json.dumps(_convert(payload), indent=2, sort_keys=False) + "\n", encoding="utf-8")
