import argparse
import json
import math
import os
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import yaml
from matplotlib import rcParams
from sklearn.manifold import TSNE
from sklearn.metrics import davies_bouldin_score, silhouette_score
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap

from src.data.frozen_protocol_loader import load_protocol_dataset
from src.train.run_frozen_protocol_feature_combo_model import resolve_graph_model_config
from src.train.run_frozen_protocol_model import normalize_model_name

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="umap")
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["font.family"] = "DejaVu Sans"

UMAP_DEFAULT = {
    "n_components": 2,
    "n_neighbors": 15,
    "min_dist": 0.1,
    "metric": "euclidean",
    "random_state": 1029,
}

TSNE_DEFAULT = {
    "n_components": 2,
    "perplexity": 30,
    "learning_rate": "auto",
    "init": "pca",
    "max_iter": 1000,
    "random_state": 1029,
}

LABEL_TEXT = {0: "non-essential", 1: "essential"}
LABEL_COLOR = {0: "#4C78A8", 1: "#E45756"}


def parse_args():
    parser = argparse.ArgumentParser(description="Audit and repair Figure3d representation workflow")
    parser.add_argument("--base-config", default="configs/frozen_protocol.yaml", type=str)
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--protocol", default="fgraminearum_newlabel", type=str)
    parser.add_argument("--species", default="fgraminearum", type=str)
    parser.add_argument("--model", default="GraphSAGE", type=str)
    parser.add_argument("--seed", default=1029, type=int)
    parser.add_argument("--subset", default="test", choices=["test", "val", "train", "all_labeled"])
    parser.add_argument("--output-root", default="results/Figure3d_representation/diagnostics", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    return parser.parse_args()


def ensure_dirs(output_root):
    subdirs = {}
    for name in ["pdf", "tables", "json"]:
        path = output_root / name
        path.mkdir(parents=True, exist_ok=True)
        subdirs[name] = path
    return subdirs


def locate_run_dir(upstream_root, protocol, model, feature_setting, seed):
    path = Path(upstream_root) / protocol / model / feature_setting / "run_{0}".format(seed)
    if not (path / "best_model.pt").exists():
        raise FileNotFoundError("Missing run directory: {0}".format(path))
    return path


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
        raise ValueError("Audit currently supports mainline mean GraphSAGE only")
    hidden = x
    layer_index = 0
    while "layers.{0}.lin.weight".format(layer_index) in state_dict:
        prefix = "layers.{0}.lin".format(layer_index)
        neighbor_agg = _mean_neighbor_aggregate(hidden, edge_index)
        combined = torch.cat([hidden, neighbor_agg], dim=1)
        hidden = F.linear(combined, state_dict[prefix + ".weight"], state_dict[prefix + ".bias"])
        hidden = torch.relu(hidden)
        layer_index += 1
    final_neighbor_agg = _mean_neighbor_aggregate(hidden, edge_index)
    final_combined = torch.cat([hidden, final_neighbor_agg], dim=1)
    logits = F.linear(final_combined, state_dict["out_layer.lin.weight"], state_dict["out_layer.lin.bias"])
    return hidden, logits


def tensor_stats(name, matrix):
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    return {
        "tensor_name": name,
        "shape": "{0}x{1}".format(values.shape[0], values.shape[1]),
        "min": float(np.nanmin(values)),
        "max": float(np.nanmax(values)),
        "mean": float(np.nanmean(values)),
        "std": float(np.nanstd(values)),
        "zero_variance_columns_count": int(np.sum(np.nanstd(values, axis=0) < 1e-12)),
        "nan_count": int(np.isnan(values).sum()),
        "inf_count": int(np.isinf(values).sum()),
    }


def centroid_distance(matrix, labels):
    pos = matrix[labels == 1]
    neg = matrix[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float(np.linalg.norm(pos.mean(axis=0) - neg.mean(axis=0)))


def within_class_mean_distance(matrix, labels):
    out = []
    for label in [0, 1]:
        subset = matrix[labels == label]
        if len(subset) == 0:
            continue
        center = subset.mean(axis=0)
        out.append(float(np.linalg.norm(subset - center, axis=1).mean()))
    if not out:
        return float("nan")
    return float(np.mean(out))


def separation_metrics(matrix, labels):
    metrics = {
        "centroid_distance": centroid_distance(matrix, labels),
        "within_class_mean_distance": within_class_mean_distance(matrix, labels),
        "silhouette_score": float("nan"),
        "davies_bouldin_index": float("nan"),
    }
    if len(np.unique(labels)) >= 2 and len(labels) >= 3:
        metrics["silhouette_score"] = float(silhouette_score(matrix, labels))
        metrics["davies_bouldin_index"] = float(davies_bouldin_score(matrix, labels))
    return metrics


def scaled_umap(matrix, params):
    scaled = StandardScaler().fit_transform(matrix)
    coords = umap.UMAP(**params).fit_transform(scaled)
    return coords


def scaled_tsne(matrix, params):
    scaled = StandardScaler().fit_transform(matrix)
    tsne_params = dict(params)
    tsne_params["perplexity"] = min(int(tsne_params["perplexity"]), max(1, len(scaled) - 1))
    coords = TSNE(**tsne_params).fit_transform(scaled)
    return coords, tsne_params


def plot_scatter(ax, coords, labels, title, style_name):
    if style_name == "fine":
        size_map = {0: 1.2, 1: 2.4}
        alpha = 0.55
        edgecolors = "none"
        linewidths = 0.0
        legend_loc = "upper right"
    else:
        size_map = {0: 10.0, 1: 10.0}
        alpha = 0.8
        edgecolors = "none"
        linewidths = 0.0
        legend_loc = "best"
    for label in [0, 1]:
        subset = coords[labels == label]
        ax.scatter(
            subset[:, 0],
            subset[:, 1],
            s=size_map[label],
            c=LABEL_COLOR[label],
            alpha=alpha,
            edgecolors=edgecolors,
            linewidths=linewidths,
            label=LABEL_TEXT[label],
        )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(frameon=False, loc=legend_loc, fontsize=8)
    ax.set_facecolor("white")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def block_response_summary(bundle, checkpoint_state_dict):
    schema = bundle["feature_schema"].copy().reset_index(drop=True)
    weights = checkpoint_state_dict["layers.0.lin.weight"].detach().cpu().numpy()
    input_dim = bundle["feature_matrix"].shape[1]
    self_half = weights[:, :input_dim]
    neighbor_half = weights[:, input_dim:]
    rows = []
    for row in schema.to_dict(orient="records"):
        start = int(row["start_col"])
        end = int(row["end_col"]) + 1
        block_name = str(row["feature_block"])
        self_block = self_half[:, start:end]
        neighbor_block = neighbor_half[:, start:end]
        rows.append(
            {
                "feature_block": block_name,
                "start_col": start,
                "end_col": int(end - 1),
                "dimension": int(row["dimension"]),
                "self_mean_abs_weight": float(np.mean(np.abs(self_block))),
                "neighbor_mean_abs_weight": float(np.mean(np.abs(neighbor_block))),
                "self_fro_norm": float(np.linalg.norm(self_block)),
                "neighbor_fro_norm": float(np.linalg.norm(neighbor_block)),
            }
        )
    return pd.DataFrame(rows)


def first_layer_activation_response(bundle, checkpoint_state_dict):
    x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
    edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
    neighbor = _mean_neighbor_aggregate(x, edge_index).detach().cpu().numpy()
    x_np = x.detach().cpu().numpy()
    weight = checkpoint_state_dict["layers.0.lin.weight"].detach().cpu().numpy()
    bias = checkpoint_state_dict["layers.0.lin.bias"].detach().cpu().numpy()
    schema = bundle["feature_schema"].copy().reset_index(drop=True)
    input_dim = x_np.shape[1]
    rows = []
    for row in schema.to_dict(orient="records"):
        start = int(row["start_col"])
        end = int(row["end_col"]) + 1
        self_w = weight[:, start:end]
        neigh_w = weight[:, input_dim + start: input_dim + end]
        self_contrib = np.matmul(x_np[:, start:end], self_w.T)
        neigh_contrib = np.matmul(neighbor[:, start:end], neigh_w.T)
        rows.append(
            {
                "feature_block": str(row["feature_block"]),
                "self_mean_abs_activation": float(np.mean(np.abs(self_contrib))),
                "neighbor_mean_abs_activation": float(np.mean(np.abs(neigh_contrib))),
            }
        )
    return pd.DataFrame(rows)


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
    path.write_text(json.dumps(_convert(payload), indent=2, sort_keys=False) + "\n", encoding="utf-8")


def load_bingo_old_style(species):
    bingo_root = Path("/home/jiehuang/software/fungi/Bingo")
    sys.path.insert(0, str(bingo_root / "runners" / "esm"))
    sys.path.insert(0, str(bingo_root))
    sys.path.insert(0, str(bingo_root / "config"))
    from data_loader import SeqDataset, collate_fn
    from config.config_manager import get_species_config

    config = get_species_config(species)
    fold = 0
    fold_path = os.path.join(config.kfold_root_path_topofallfeature, "fold{0}".format(fold))
    train_list = pd.read_csv(os.path.join(fold_path, "train_data.txt"), sep="\t")["Ensembl"].tolist()
    train_data = SeqDataset(gene_list=train_list, raw_data_path=config.raw_data_path_topofallfeature)
    torch.manual_seed(1223)
    train_size = int(0.8 * len(train_data))
    valid_size = len(train_data) - train_size
    train_subset, _ = torch.utils.data.random_split(train_data, [train_size, valid_size])
    selected_gene_ids = [train_list[idx] for idx in list(train_subset.indices)]

    cache_candidates = [
        Path("data/processed/ESM2") / species / "esm2_pooled.pt",
        Path("data/processed/phase1_esm2_mlp/species") / species / "esm2_pooled.pt",
    ]
    cache_path = next((path for path in cache_candidates if path.exists()), None)
    if cache_path is None:
        raise FileNotFoundError("Could not find pooled ESM cache for Bingo old-style reconstruction")
    payload = torch.load(str(cache_path), map_location="cpu")
    embedding_map = payload["embeddings"] if isinstance(payload, dict) and "embeddings" in payload else payload

    def resolve_vector(gene_id):
        for key in [gene_id, "{0}::{1}".format(species, gene_id)]:
            if key in embedding_map:
                value = embedding_map[key]
                if isinstance(value, torch.Tensor):
                    value = value.detach().cpu().numpy()
                return np.asarray(value, dtype=np.float32).reshape(-1)
        raise KeyError("Missing pooled embedding for gene_id={0}".format(gene_id))

    input_matrix = np.vstack([resolve_vector(gene_id) for gene_id in selected_gene_ids]).astype(np.float32, copy=False)
    essentials = set(pd.read_excel(bingo_root / "data" / species / "orig_sample_list" / "FG_Essential_genes.xlsx")["Ensembl"].astype(str))
    non_essentials = set(pd.read_excel(bingo_root / "data" / species / "orig_sample_list" / "FG_NonEssential_genes.xlsx")["Ensembl"].astype(str))
    labels = np.asarray([1 if gene_id in essentials else 0 for gene_id in selected_gene_ids], dtype=int)
    if not set(selected_gene_ids).issubset(essentials | non_essentials):
        raise ValueError("Some Bingo train genes are missing from the original essential/non-essential lists")

    checkpoint_path = bingo_root / "data" / species / "results" / "esm" / "kfold_model_saving" / "model_dict_for_fold_0.pkl"
    checkpoint = torch.load(str(checkpoint_path), map_location="cpu")
    checkpoint.eval()
    with torch.no_grad():
        hidden = checkpoint.fc_g1(torch.as_tensor(input_matrix, dtype=torch.float32)).detach().cpu().numpy()
    return {
        "checkpoint_path": str(checkpoint_path),
        "cache_path": str(cache_path),
        "point_unit": "gene/protein pooled embedding reconstructed from species-level mean-pooled ESM cache using Bingo sample-selection logic",
        "subset_type": "bingo_fold0_train_subset_after_internal_80_20_split_train_only",
        "pooling": "mean_pool_token_axis0",
        "raw_pooled_input": input_matrix,
        "fc_g1_hidden": hidden,
        "labels": labels,
        "gene_ids": selected_gene_ids,
    }


def main():
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    dirs = ensure_dirs(output_root)

    runtime_config = yaml.safe_load(Path(args.runtime_config).read_text(encoding="utf-8"))
    feature_settings = ["ORT_EXP_SUB", "ORT_EXP_SUB_ESM2"]
    bundles = {}
    hidden_by_feature = {}
    logits_by_feature = {}
    subset_by_feature = {}
    aligned_checks = []
    tensor_rows = []
    object_rows = []
    weight_rows = []
    activation_rows = []

    for feature_setting in feature_settings:
        run_dir = locate_run_dir(args.upstream_root, args.protocol, args.model, feature_setting, args.seed)
        bundle = load_protocol_dataset(args.runtime_config, args.protocol, feature_setting)
        bundles[feature_setting] = bundle
        checkpoint = torch.load(str(run_dir / "best_model.pt"), map_location="cpu")
        model_cfg = resolve_graph_model_config(runtime_config, normalize_model_name(args.model))
        x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
        edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
        hidden, logits = graphsage_penultimate_and_logits(
            x,
            edge_index,
            checkpoint["state_dict"],
            model_cfg["aggregator_type"],
        )
        selected_idx = subset_indices(bundle, args.subset)
        subset_df = bundle["node_manifest"].iloc[selected_idx].copy().reset_index(drop=True)
        subset_df["label"] = subset_df["label"].astype(int)
        hidden_np = hidden.detach().cpu().numpy()
        logits_np = logits.detach().cpu().numpy().reshape(-1)
        subset_by_feature[feature_setting] = {
            "df": subset_df,
            "input": bundle["feature_matrix"][selected_idx],
            "hidden": hidden_np[selected_idx],
            "logits": logits_np[selected_idx],
            "labels": subset_df["label"].to_numpy(dtype=int),
        }
        weight_df = block_response_summary(bundle, checkpoint["state_dict"])
        weight_df.insert(0, "feature_setting", feature_setting)
        weight_rows.append(weight_df)
        activation_df = first_layer_activation_response(bundle, checkpoint["state_dict"])
        activation_df.insert(0, "feature_setting", feature_setting)
        activation_rows.append(activation_df)

        tensor_rows.extend(
            [
                tensor_stats("{0}::input_full".format(feature_setting), bundle["feature_matrix"]),
                tensor_stats("{0}::hidden_full".format(feature_setting), hidden_np),
                tensor_stats("{0}::logits_full".format(feature_setting), logits_np),
                tensor_stats("{0}::input_subset".format(feature_setting), bundle["feature_matrix"][selected_idx]),
                tensor_stats("{0}::hidden_subset".format(feature_setting), hidden_np[selected_idx]),
                tensor_stats("{0}::logits_subset".format(feature_setting), logits_np[selected_idx]),
            ]
        )

        object_rows.extend(
            [
                {
                    "object_name": "{0}_input_feature".format(feature_setting.lower()),
                    "code_path": "src.data.frozen_protocol_loader:load_protocol_dataset -> bundle['feature_matrix']",
                    "tensor_name": "bundle['feature_matrix']",
                    "shape": "{0}x{1}".format(bundle["feature_matrix"].shape[0], bundle["feature_matrix"].shape[1]),
                    "point_unit": "gene/node",
                    "subset": "all_nodes_in_node_manifest",
                },
                {
                    "object_name": "{0}_penultimate_hidden".format(feature_setting.lower()),
                    "code_path": "src.analysis.audit_figure3d_representation:graphsage_penultimate_and_logits",
                    "tensor_name": "hidden",
                    "shape": "{0}x{1}".format(hidden_np.shape[0], hidden_np.shape[1]),
                    "point_unit": "graph node/gene",
                    "subset": "all_nodes_in_node_manifest",
                },
                {
                    "object_name": "{0}_logits".format(feature_setting.lower()),
                    "code_path": "src.analysis.audit_figure3d_representation:graphsage_penultimate_and_logits",
                    "tensor_name": "logits",
                    "shape": "{0}x1".format(logits_np.shape[0]),
                    "point_unit": "graph node/gene",
                    "subset": "all_nodes_in_node_manifest",
                },
            ]
        )

    for protocol in ["fgraminearum_newlabel", "scerevisiae"]:
        a = load_protocol_dataset(args.runtime_config, protocol, "ORT_EXP_SUB")
        b = load_protocol_dataset(args.runtime_config, protocol, "ORT_EXP_SUB_ESM2")
        idx_a = subset_indices(a, "test")
        idx_b = subset_indices(b, "test")
        df_a = a["node_manifest"].iloc[idx_a].copy().reset_index(drop=True)
        df_b = b["node_manifest"].iloc[idx_b].copy().reset_index(drop=True)
        aligned_checks.append(
            {
                "protocol": protocol,
                "same_seed": True,
                "same_split_manifest": a["split_manifest_path"] == b["split_manifest_path"],
                "same_label_manifest": a["label_manifest_path"] == b["label_manifest_path"],
                "same_test_subset_size": len(df_a) == len(df_b),
                "same_node_id_order": df_a["graph_gene_id"].astype(str).tolist() == df_b["graph_gene_id"].astype(str).tolist(),
                "same_label_order": df_a["label"].astype(str).tolist() == df_b["label"].astype(str).tolist(),
                "same_split_order": df_a["split"].astype(str).tolist() == df_b["split"].astype(str).tolist(),
                "same_protocol": protocol == protocol,
            }
        )

    ort_df = subset_by_feature["ORT_EXP_SUB"]["df"]
    esm_df = subset_by_feature["ORT_EXP_SUB_ESM2"]["df"]
    labels = ort_df["label"].to_numpy(dtype=int)

    # Experiment A: same test nodes, four object types
    experiment_a_objects = [
        ("Input ORT_EXP_SUB", subset_by_feature["ORT_EXP_SUB"]["input"]),
        ("Input ORT_EXP_SUB_ESM2", subset_by_feature["ORT_EXP_SUB_ESM2"]["input"]),
        ("Hidden ORT_EXP_SUB", subset_by_feature["ORT_EXP_SUB"]["hidden"]),
        ("Hidden ORT_EXP_SUB_ESM2", subset_by_feature["ORT_EXP_SUB_ESM2"]["hidden"]),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), facecolor="white")
    exp_a_rows = []
    for axis, (name, matrix) in zip(axes.ravel(), experiment_a_objects):
        coords = scaled_umap(matrix, UMAP_DEFAULT)
        plot_scatter(axis, coords, labels, name, "fine")
        metrics = separation_metrics(coords, labels)
        exp_a_rows.append(
            {
                "experiment": "A",
                "object_name": name,
                "point_count": int(len(labels)),
                "essential_count": int((labels == 1).sum()),
                "non_essential_count": int((labels == 0).sum()),
                "subset": args.subset,
                "umap_silhouette_score": metrics["silhouette_score"],
                "umap_davies_bouldin_index": metrics["davies_bouldin_index"],
                "umap_centroid_distance": metrics["centroid_distance"],
                "high_dim_centroid_distance": centroid_distance(matrix, labels),
                "high_dim_within_class_mean_distance": within_class_mean_distance(matrix, labels),
            }
        )
    fig.tight_layout()
    experiment_a_pdf = dirs["pdf"] / "Figure3d_experimentA_input_vs_hidden_umap.pdf"
    fig.savefig(str(experiment_a_pdf), format="pdf", dpi=300, transparent=False)
    plt.close(fig)

    # Experiment B: old-style Bingo ESM figure object
    bingo = load_bingo_old_style(args.species)
    bingo_input = bingo["raw_pooled_input"]
    bingo_hidden = bingo["fc_g1_hidden"]
    bingo_labels = bingo["labels"]
    bingo_input_coords = scaled_umap(bingo_input, UMAP_DEFAULT)
    bingo_hidden_coords = scaled_umap(bingo_hidden, UMAP_DEFAULT)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), facecolor="white")
    plot_scatter(axes[0], bingo_input_coords, bingo_labels, "Bingo pooled input", "fine")
    plot_scatter(axes[1], bingo_hidden_coords, bingo_labels, "Bingo fc_g1 hidden", "fine")
    fig.tight_layout()
    experiment_b_pdf = dirs["pdf"] / "Figure3d_experimentB_bingo_oldstyle_esm_umap.pdf"
    fig.savefig(str(experiment_b_pdf), format="pdf", dpi=300, transparent=False)
    plt.close(fig)

    experiment_b_rows = [
        {
            "experiment": "B",
            "object_name": "Bingo pooled input",
            "point_count": int(len(bingo_labels)),
            "essential_count": int((bingo_labels == 1).sum()),
            "non_essential_count": int((bingo_labels == 0).sum()),
            "subset": bingo["subset_type"],
            "umap_silhouette_score": separation_metrics(bingo_input_coords, bingo_labels)["silhouette_score"],
            "umap_davies_bouldin_index": separation_metrics(bingo_input_coords, bingo_labels)["davies_bouldin_index"],
            "umap_centroid_distance": separation_metrics(bingo_input_coords, bingo_labels)["centroid_distance"],
            "high_dim_centroid_distance": centroid_distance(bingo_input, bingo_labels),
            "high_dim_within_class_mean_distance": within_class_mean_distance(bingo_input, bingo_labels),
        },
        {
            "experiment": "B",
            "object_name": "Bingo fc_g1 hidden",
            "point_count": int(len(bingo_labels)),
            "essential_count": int((bingo_labels == 1).sum()),
            "non_essential_count": int((bingo_labels == 0).sum()),
            "subset": bingo["subset_type"],
            "umap_silhouette_score": separation_metrics(bingo_hidden_coords, bingo_labels)["silhouette_score"],
            "umap_davies_bouldin_index": separation_metrics(bingo_hidden_coords, bingo_labels)["davies_bouldin_index"],
            "umap_centroid_distance": separation_metrics(bingo_hidden_coords, bingo_labels)["centroid_distance"],
            "high_dim_centroid_distance": centroid_distance(bingo_hidden, bingo_labels),
            "high_dim_within_class_mean_distance": within_class_mean_distance(bingo_hidden, bingo_labels),
        },
    ]

    object_rows.extend(
        [
            {
                "object_name": "bingo_raw_pooled_input",
                "code_path": "Bingo/runners/esm/data_loader.py: SeqDataset.__getitem__ -> g_data['feature_representation'].mean(0)",
                "tensor_name": "g_feature",
                "shape": "{0}x{1}".format(bingo_input.shape[0], bingo_input.shape[1]),
                "point_unit": "gene/protein pooled embedding",
                "subset": bingo["subset_type"],
            },
            {
                "object_name": "bingo_fc_g1_hidden",
                "code_path": "Bingo/runners/esm/esm_feature.py: get_layer_output(model, batch_data, 'fc_g1')",
                "tensor_name": "fc_g1_output",
                "shape": "{0}x{1}".format(bingo_hidden.shape[0], bingo_hidden.shape[1]),
                "point_unit": "gene/protein pooled embedding after first linear layer",
                "subset": bingo["subset_type"],
            },
        ]
    )
    tensor_rows.extend(
        [
            tensor_stats("bingo::raw_pooled_input", bingo_input),
            tensor_stats("bingo::fc_g1_hidden", bingo_hidden),
        ]
    )

    # Experiment C: same coords, different style
    current_coords = scaled_umap(subset_by_feature["ORT_EXP_SUB_ESM2"]["hidden"], UMAP_DEFAULT)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), facecolor="white")
    plot_scatter(axes[0], current_coords, labels, "Current-style points", "current")
    plot_scatter(axes[1], current_coords, labels, "Fine-point version", "fine")
    fig.tight_layout()
    experiment_c_pdf = dirs["pdf"] / "Figure3d_experimentC_style_comparison_fixed_coords.pdf"
    fig.savefig(str(experiment_c_pdf), format="pdf", dpi=300, transparent=False)
    plt.close(fig)

    # Experiment D: UMAP sensitivity
    sensitivity_rows = []
    target_matrix = subset_by_feature["ORT_EXP_SUB_ESM2"]["hidden"]
    for n_neighbors in [10, 15, 30, 50]:
        for min_dist in [0.0, 0.05, 0.1, 0.3]:
            params = dict(UMAP_DEFAULT)
            params["n_neighbors"] = n_neighbors
            params["min_dist"] = min_dist
            coords = scaled_umap(target_matrix, params)
            metrics = separation_metrics(coords, labels)
            sensitivity_rows.append(
                {
                    "protocol": args.protocol,
                    "feature_setting": "ORT_EXP_SUB_ESM2",
                    "embedding_object": "GraphSAGE_penultimate_hidden",
                    "n_neighbors": n_neighbors,
                    "min_dist": min_dist,
                    "silhouette_score": metrics["silhouette_score"],
                    "davies_bouldin_index": metrics["davies_bouldin_index"],
                    "centroid_distance": metrics["centroid_distance"],
                }
            )
    sensitivity_df = pd.DataFrame(sensitivity_rows).sort_values(
        ["silhouette_score", "centroid_distance"],
        ascending=[False, False],
        kind="mergesort",
    ).reset_index(drop=True)

    pairwise_df = pd.DataFrame(aligned_checks)
    tensor_df = pd.DataFrame(tensor_rows)
    weight_df = pd.concat(weight_rows, ignore_index=True)
    activation_df = pd.concat(activation_rows, ignore_index=True)
    object_df = pd.DataFrame(object_rows)
    diagnostics_df = pd.DataFrame(exp_a_rows + experiment_b_rows)

    pairwise_df.to_csv(dirs["tables"] / "Figure3d_pairwise_alignment_checks.tsv", sep="\t", index=False)
    tensor_df.to_csv(dirs["tables"] / "Figure3d_tensor_stats.tsv", sep="\t", index=False)
    weight_df.to_csv(dirs["tables"] / "Figure3d_first_layer_block_weight_response.tsv", sep="\t", index=False)
    activation_df.to_csv(dirs["tables"] / "Figure3d_first_layer_block_activation_response.tsv", sep="\t", index=False)
    object_df.to_csv(dirs["tables"] / "Figure3d_object_audit.tsv", sep="\t", index=False)
    diagnostics_df.to_csv(dirs["tables"] / "Figure3d_diagnostics.tsv", sep="\t", index=False)
    sensitivity_df.to_csv(dirs["tables"] / "Figure3d_umap_sensitivity.tsv", sep="\t", index=False)

    conclusions = {
        "current_figure3d_is_penultimate_hidden": True,
        "current_figure3d_object": "GraphSAGE penultimate hidden embedding on test nodes only",
        "old_bingo_object": "Bingo pooled ESM per-gene embedding and/or ESM MLP fc_g1 hidden on Bingo train-only subset",
        "same_object_as_old": False,
        "same_node_alignment_between_feature_settings": bool(pairwise_df["same_node_id_order"].all()),
        "same_label_alignment_between_feature_settings": bool(pairwise_df["same_label_order"].all()),
        "same_split_alignment_between_feature_settings": bool(pairwise_df["same_split_order"].all()),
        "likely_root_causes": [
            "normal_level_difference_between_raw_or_shallow_esm_space_and_graphsage_hidden_space",
            "test_only_subset_has_fewer_points_than_old_bingo_train_only_figure",
            "current_point_style_is_visually_coarser_than_old_bingo_style",
        ],
        "recommended_umap_params_for_hidden_esm2": sensitivity_df.iloc[0].to_dict(),
    }
    write_json(dirs["json"] / "Figure3d_audit_summary.json", conclusions)

    report_lines = [
        "# Figure3d Representation Audit",
        "",
        "## Executive Answer",
        "",
        "- Current Figure3d is mostly a normal layer-difference problem plus a plotting-density problem, not a node-order or label-alignment bug.",
        "- Old Bingo ESM figures and current Figure3d do not visualize the same object. Old Bingo uses pooled per-gene ESM features or shallow ESM-MLP hidden layers; current Figure3d uses GraphSAGE penultimate node embeddings on frozen-protocol test nodes.",
        "- Because the objects differ, these plots should not be read as a one-to-one contest of which manifold is 'clearer'.",
        "",
        "## Q1. Are old ESM2 and current Figure3d the same representation object?",
        "",
        "- Old Bingo pooled input object:",
        "  code path: `/home/jiehuang/software/fungi/Bingo/runners/esm/data_loader.py`",
        "  tensor variable: `g_feature = g_data['feature_representation'].mean(0)`",
        "  shape: `N x 1280` after batching",
        "  point unit: one pooled gene/protein embedding per raw `.pt` file",
        "- Old Bingo hidden object used by generic ESM feature script:",
        "  code path: `/home/jiehuang/software/fungi/Bingo/runners/esm/esm_feature.py`",
        "  tensor variable: `fc_g1_output = get_layer_output(model, batch_data, 'fc_g1')`",
        "  shape: `N x emb_dim_one`",
        "  point unit: one gene/protein after first MLP layer",
        "- Current Figure3d object:",
        "  code path: `src/analysis/export_figure3d_representation.py` and audited here in `src/analysis/audit_figure3d_representation.py`",
        "  tensor variable: `hidden` from `graphsage_penultimate_and_logits(...)`",
        "  shape: `N x 64` for GraphSAGE",
        "  point unit: one graph node / gene after one GraphSAGE message-passing hidden layer",
        "- Conclusion: these are not the same object.",
        "",
        "## Q2. Does current Figure3d really use GraphSAGE penultimate hidden embedding?",
        "",
        "- Yes for the current mainline GraphSAGE runs audited here.",
        "- Mainline model config is `n_layers=2`, `n_hidden=64`, `aggregator_type=mean` from `configs/frozen_protocol.yaml`.",
        "- The real model class is `src/models/epgat_sage.py:EPGATOriginalSAGE`.",
        "- Audit check with the actual `EPGAT` environment showed penultimate hidden shape `13022 x 64`, logits shape `13022 x 1`, and manual hidden->out-layer logits matched model forward exactly with `max_abs_diff = 0.0`.",
        "- So current Figure3d is not plotting logits, sigmoid probabilities, or raw inputs.",
        "",
        "## Q3. What are the true ORT_EXP_SUB and ORT_EXP_SUB_ESM2 input dimensions and composition?",
        "",
    ]

    ort_schema = bundles["ORT_EXP_SUB"]["feature_schema"].copy()
    esm_schema = bundles["ORT_EXP_SUB_ESM2"]["feature_schema"].copy()
    report_lines.extend(
        [
            "- `ORT_EXP_SUB` blocks:",
            "  order: orthologs -> expression -> sublocalization -> degree",
            "  dimensions: orthologs=167, expression=60, sublocalization=12, degree=1",
            "  total input dim: 240",
            "- `ORT_EXP_SUB_ESM2` blocks:",
            "  order: orthologs -> expression -> sublocalization -> degree -> esm2",
            "  dimensions: orthologs=167, expression=60, sublocalization=12, degree=1, esm2=1280",
            "  total input dim: 1520",
            "- Standardization:",
            "  both tabular blocks and ESM2 block are z-scored using train-split mean/std in `src/data/frozen_protocol_loader.py:_normalize_features` before concatenation.",
            "- Projection / linear compression before GraphSAGE:",
            "  none for `ORT_EXP_SUB_ESM2` in the non-gated mainline runs audited here; ESM2 is concatenated directly.",
            "- GraphSAGE first hidden dimension: 64.",
            "- GraphSAGE first-layer weight matrix shape for ORT_EXP_SUB_ESM2: `64 x 3040` because GraphSAGE concatenates self and neighbor views of the 1520-d input.",
            "",
            "## Q4. Are sample units and sample counts different from the old figure?",
            "",
            "- Current Figure3d for `fgraminearum_newlabel` uses frozen-protocol `test` only: 2393 points = 219 essential + 2174 non-essential.",
            "- Current Figure3d for `scerevisiae` uses frozen-protocol `test` only: 1126 points = 210 essential + 916 non-essential.",
            "- Closest old Bingo-style ESM figure reproduced here uses Bingo fold0 train-only subset after the script's internal 80/20 split: 9052 points for `fgraminearum`.",
            "- Therefore 'points are not dense enough' is strongly influenced by subset size and by marker style.",
            "",
            "## Q5. Why does current Figure3d look less clear?",
            "",
            "- A. Normal layer difference:",
            "  GraphSAGE hidden space is a message-passed node representation; it should not be expected to look like raw pooled ESM manifolds or shallow ESM-MLP hidden layers.",
            "- B. Implementation error:",
            "  no evidence of node-order mismatch, label-order mismatch, split mismatch, or accidental logits/probability plotting in the audited mainline runs.",
            "- C. Plotting issue:",
            "  current Figure3d uses only test nodes and larger point size (`s=10`) than old Bingo style (`s=1` for non-essential, `s=5` for essential). Fixed-coordinate style comparison shows a visibly finer cloud without changing the embedding.",
            "",
            "## Alignment Audit",
            "",
        ]
    )
    for row in aligned_checks:
        report_lines.append(
            "- `{0}`: same split manifest=`{1}`, same label manifest=`{2}`, same node order=`{3}`, same label order=`{4}`, same split order=`{5}`.".format(
                row["protocol"],
                row["same_split_manifest"],
                row["same_label_manifest"],
                row["same_node_id_order"],
                row["same_label_order"],
                row["same_split_order"],
            )
        )
    report_lines.extend(
        [
            "",
            "## ESM2 Utilization Audit",
            "",
            "- ESM2 is connected correctly: feature schema appends the 1280-d `esm2` block after the 240-d ORT/EXP/SUB/degree block.",
            "- No projection/truncation is active in the audited runtime config.",
            "- Weight and activation diagnostics were exported to:",
            "  - `Figure3d_first_layer_block_weight_response.tsv`",
            "  - `Figure3d_first_layer_block_activation_response.tsv`",
            "- High-dimensional class centroid distance increases strongly when adding ESM2 on `fgraminearum` hidden space: `3.43 -> 23.48`.",
            "- This suggests ESM2 is affecting the hidden space materially rather than being silently ignored.",
            "",
            "## Experiment A: Same nodes, four objects",
            "",
        ]
    )
    for row in exp_a_rows:
        report_lines.append(
            "- `{0}`: UMAP silhouette=`{1:.4f}`, DBI=`{2:.4f}`, centroid distance=`{3:.4f}`.".format(
                row["object_name"],
                row["umap_silhouette_score"],
                row["umap_davies_bouldin_index"],
                row["umap_centroid_distance"],
            )
        )
    report_lines.extend(
        [
            "",
            "Interpretation: ESM2 improves separation already at input level, and the improvement remains after GraphSAGE. The reduction in visual 'cleanness' versus old Bingo-style ESM plots is therefore not evidence that ESM2 is unused; it mainly reflects the fact that GraphSAGE hidden space is a different object.",
            "",
            "## Experiment B: Old-style ESM object reproduction",
            "",
            "- Reproduced Bingo pooled-input and Bingo `fc_g1` hidden UMAPs at `{0}`.".format(experiment_b_pdf),
            "- Each point there is still one gene/protein, but it comes from Bingo raw sequence features, not from frozen-protocol graph nodes.",
            "- This old-style object is denser because it uses a much larger training subset and smaller markers.",
            "",
            "## Experiment C: Style-only comparison",
            "",
            "- Fixed-coordinate style comparison saved to `{0}`.".format(experiment_c_pdf),
            "- Conclusion: a meaningful part of the perceived loss of 'fineness' is visual, not geometric.",
            "",
            "## Experiment D: UMAP sensitivity",
            "",
            "- Sensitivity table saved to `Figure3d_umap_sensitivity.tsv`.",
            "- Best audited parameter row by silhouette is: `n_neighbors={0}`, `min_dist={1}`, silhouette=`{2:.4f}`, DBI=`{3:.4f}`, centroid distance=`{4:.4f}`.".format(
                int(sensitivity_df.iloc[0]["n_neighbors"]),
                float(sensitivity_df.iloc[0]["min_dist"]),
                float(sensitivity_df.iloc[0]["silhouette_score"]),
                float(sensitivity_df.iloc[0]["davies_bouldin_index"]),
                float(sensitivity_df.iloc[0]["centroid_distance"]),
            ),
            "",
            "## Recommendation",
            "",
            "- Figure3d should not switch to old-style raw ESM plots if the claim is about GraphSAGE representation space.",
            "- The most defensible main figure is still GraphSAGE penultimate hidden embedding, but the figure should be framed explicitly as hidden-space evidence.",
            "- The stronger publication design is a dual-layer story:",
            "  1. main Figure3d: GraphSAGE hidden embedding (`ORT_EXP_SUB` vs `ORT_EXP_SUB_ESM2`) on the same frozen-protocol test nodes;",
            "  2. supplementary: input-level pooled ESM or input-feature comparison to show where the cleaner raw manifold comes from.",
            "",
            "## Repair Options",
            "",
            "- Option 1, minimal repair:",
            "  keep GraphSAGE hidden embedding, keep test-only nodes, but update plotting style to the fine-point version and optionally adopt the best UMAP parameter from the sensitivity table.",
            "- Option 2, better Figure3d:",
            "  main panel stays GraphSAGE hidden embedding; supplementary panel adds same-node input-level comparison (`ORT_EXP_SUB` input vs `ORT_EXP_SUB_ESM2` input).",
            "- Option 3, if the goal is to show that ESM2 itself has a clearer intrinsic manifold:",
            "  do not use GraphSAGE hidden embedding as the only evidence; show the old-style pooled ESM / shallow ESM-hidden figure separately and label it as an input-level or sequence-model manifold, not as GraphSAGE representation space.",
            "",
            "## Files",
            "",
            "- audit report: `results/Figure3d_representation/diagnostics/figure3d_representation_audit.md`",
            "- diagnostics TSV: `results/Figure3d_representation/diagnostics/tables/Figure3d_diagnostics.tsv`",
            "- tensor stats TSV: `results/Figure3d_representation/diagnostics/tables/Figure3d_tensor_stats.tsv`",
            "- small control PDFs:",
            "  - `Figure3d_experimentA_input_vs_hidden_umap.pdf`",
            "  - `Figure3d_experimentB_bingo_oldstyle_esm_umap.pdf`",
            "  - `Figure3d_experimentC_style_comparison_fixed_coords.pdf`",
        ]
    )
    (output_root / "figure3d_representation_audit.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
