import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

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


rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["font.family"] = "DejaVu Sans"

DEFAULT_PROTOCOLS = ["fgraminearum_newlabel", "scerevisiae"]
DEFAULT_FEATURE_SETTINGS = ["ORT_EXP_SUB", "ORT_EXP_SUB_ESM2"]
DEFAULT_UPSTREAM_ROOTS = ["outputs/Figure3a", "outputs/Figure2a", "outputs/Figure1"]
LABEL_TO_TEXT = {0: "non-essential", 1: "essential"}
LABEL_TO_COLOR = {0: "#4C78A8", 1: "#E45756"}

UMAP_PARAMS = {
    "n_components": 2,
    "n_neighbors": 15,
    "min_dist": 0.1,
    "metric": "euclidean",
    "random_state": 1029,
}

TSNE_BASE_PARAMS = {
    "n_components": 2,
    "perplexity": 30,
    "learning_rate": "auto",
    "init": "pca",
    "max_iter": 1000,
    "random_state": 1029,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export Figure 3d GraphSAGE representations from mainline checkpoints")
    parser.add_argument("--base-config", default="configs/frozen_protocol.yaml", type=str)
    parser.add_argument("--protocols", nargs="+", default=DEFAULT_PROTOCOLS)
    parser.add_argument("--feature-settings", nargs="+", default=DEFAULT_FEATURE_SETTINGS)
    parser.add_argument("--model", default="GraphSAGE", type=str)
    parser.add_argument("--seed", default=1029, type=int)
    parser.add_argument("--subset", default="test", choices=["test", "val", "train", "all_labeled"])
    parser.add_argument("--output-root", default="results/Figure3d_representation", type=str)
    parser.add_argument("--upstream-roots", nargs="+", default=DEFAULT_UPSTREAM_ROOTS)
    return parser.parse_args()


def protocol_label_regime(protocol: str) -> str:
    return "newlabel" if protocol == "fgraminearum_newlabel" else "standard"


def feature_slug(feature_setting: str) -> str:
    return str(feature_setting).strip().lower()


def output_stem(protocol: str, feature_setting: str, subset: str, seed: int, model: str) -> str:
    return f"{protocol}_{str(model).strip().lower()}_{feature_slug(feature_setting)}_{subset}_seed{seed}"


def locate_run_dir(protocol, model, feature_setting, seed, upstream_roots):
    for root in upstream_roots:
        candidate = Path(root) / protocol / model / feature_setting / f"run_{seed}"
        if (candidate / "best_model.pt").exists():
            return candidate.resolve()
    raise FileNotFoundError(
        f"Could not find checkpoint for protocol={protocol} model={model} feature_setting={feature_setting} seed={seed}"
    )


def load_run_context(run_dir, base_config_path):
    resolved_path = run_dir / "resolved_config.yaml"
    if resolved_path.exists():
        resolved = yaml.safe_load(resolved_path.read_text(encoding="utf-8"))
        config_path = str(resolved.get("config_path") or resolved.get("metrics", {}).get("config_used") or base_config_path)
        return resolved, config_path
    return {}, base_config_path


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


def _graphsage_linear_forward(x, edge_index, weight, bias, aggregator_type):
    if aggregator_type != "mean":
        raise ValueError("Figure 3d export currently supports mainline GraphSAGE mean aggregator only")
    neighbor_agg = _mean_neighbor_aggregate(x, edge_index)
    combined = torch.cat([x, neighbor_agg], dim=1)
    return F.linear(combined, weight, bias)


def forward_penultimate_from_checkpoint(x, edge_index, state_dict, aggregator_type):
    hidden = x
    layer_index = 0
    while "layers.{0}.lin.weight".format(layer_index) in state_dict:
        prefix = "layers.{0}.lin".format(layer_index)
        hidden = _graphsage_linear_forward(
            hidden,
            edge_index,
            state_dict[prefix + ".weight"],
            state_dict[prefix + ".bias"],
            aggregator_type,
        )
        hidden = torch.relu(hidden)
        layer_index += 1
    logits = _graphsage_linear_forward(
        hidden,
        edge_index,
        state_dict["out_layer.lin.weight"],
        state_dict["out_layer.lin.bias"],
        aggregator_type,
    )
    return hidden, logits


def subset_indices(bundle, subset):
    if subset == "test":
        return np.asarray(bundle["test_idx"], dtype=np.int64)
    if subset == "val":
        return np.asarray(bundle["val_idx"], dtype=np.int64)
    if subset == "train":
        return np.asarray(bundle["train_idx"], dtype=np.int64)
    if subset == "all_labeled":
        split_manifest = bundle["split_manifest"]
        mapping = bundle["mapping"]
        return np.asarray([mapping[item] for item in split_manifest["graph_gene_id"].astype(str)], dtype=np.int64)
    raise ValueError(f"Unsupported subset '{subset}'")


def centroid_distance(matrix: np.ndarray, labels: np.ndarray) -> float:
    pos = matrix[labels == 1]
    neg = matrix[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float(np.linalg.norm(pos.mean(axis=0) - neg.mean(axis=0)))


def within_class_mean_distance(matrix: np.ndarray, labels: np.ndarray) -> float:
    values = []
    for label in [0, 1]:
        subset = matrix[labels == label]
        if len(subset) == 0:
            continue
        center = subset.mean(axis=0)
        values.append(float(np.linalg.norm(subset - center, axis=1).mean()))
    if not values:
        return float("nan")
    return float(np.mean(values))


def compute_embedding_metrics(matrix: np.ndarray, labels: np.ndarray) -> Dict[str, float]:
    metrics = {
        "centroid_distance": centroid_distance(matrix, labels),
        "within_class_mean_distance": within_class_mean_distance(matrix, labels),
    }
    if len(np.unique(labels)) >= 2 and len(labels) >= 3:
        metrics["silhouette_score"] = float(silhouette_score(matrix, labels))
        metrics["davies_bouldin_index"] = float(davies_bouldin_score(matrix, labels))
    else:
        metrics["silhouette_score"] = float("nan")
        metrics["davies_bouldin_index"] = float("nan")
    return metrics


def fit_umap_coords(embedding: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
    scaler = StandardScaler()
    scaled = scaler.fit_transform(embedding)
    reducer = umap.UMAP(**UMAP_PARAMS)
    coords = reducer.fit_transform(scaled)
    return coords, dict(UMAP_PARAMS)


def fit_tsne_coords(embedding: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
    scaler = StandardScaler()
    scaled = scaler.fit_transform(embedding)
    perplexity_limit = max(1, len(scaled) - 1)
    perplexity = min(int(TSNE_BASE_PARAMS["perplexity"]), perplexity_limit)
    params = dict(TSNE_BASE_PARAMS)
    params["perplexity"] = perplexity
    reducer = TSNE(**params)
    coords = reducer.fit_transform(scaled)
    return coords, params


def save_plot(coords_df: pd.DataFrame, x_col: str, y_col: str, projection_name: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 4.0), facecolor="white")
    ax.set_facecolor("white")

    for label_value in [0, 1]:
        subset = coords_df[coords_df["label"] == label_value].copy()
        ax.scatter(
            subset[x_col],
            subset[y_col],
            s=10,
            alpha=0.8,
            c=LABEL_TO_COLOR[label_value],
            edgecolors="none",
            label=LABEL_TO_TEXT[label_value],
        )

    ax.set_xlabel(f"{projection_name} 1")
    ax.set_ylabel(f"{projection_name} 2")
    ax.legend(frameon=False, loc="best")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", dpi=300, transparent=False)
    plt.close(fig)


def save_coords_table(
    subset_df: pd.DataFrame,
    coords: np.ndarray,
    output_path: Path,
    protocol: str,
    feature_setting: str,
    model: str,
    seed: int,
) -> pd.DataFrame:
    out = pd.DataFrame(
        {
            "node_id": subset_df["graph_gene_id"].astype(str).to_numpy(),
            "species": subset_df["species"].astype(str).to_numpy(),
            "label": subset_df["label"].astype(int).to_numpy(),
            "split": subset_df["split"].astype(str).to_numpy(),
            "feature_setting": feature_setting,
            "model": model,
            "seed": seed,
            "dim1": coords[:, 0],
            "dim2": coords[:, 1],
        }
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)
    return out


def json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_ready(val) for key, val in value.items()}
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def render_case_summary(case_rows, output_path, pairwise_rows):
    lines = [
        "# Figure3d Representation Visualization Summary",
        "",
        "- Primary figure intent: UMAP is the main Figure3d panel; t-SNE is exported as supplementary robustness.",
        "- Fixed model: `GraphSAGE`.",
        "- Fixed seed: `1029`.",
        "- Fixed feature settings: `ORT_EXP_SUB` and `ORT_EXP_SUB_ESM2`.",
        "",
    ]

    for row in case_rows:
        lines.extend(
            [
                f"## {row['protocol']} | {row['feature_setting']}",
                "",
                f"- species: `{row['species']}`",
                f"- label regime: `{row['label_regime']}`",
                f"- feature setting: `{row['feature_setting']}`",
                f"- model: `{row['model']}`",
                f"- seed: `{row['seed']}`",
                f"- split manifest path: `{row['split_manifest_path']}`",
                f"- label manifest path: `{row['label_manifest_path']}`",
                f"- checkpoint path: `{row['checkpoint_path']}`",
                f"- resolved config path: `{row['resolved_config_path']}`",
                f"- whether using test-only subset: `{row['uses_test_only_subset']}`",
                f"- number of nodes plotted: `{row['node_count']}`",
                f"- class counts: essential=`{row['class_counts']['essential']}`, non-essential=`{row['class_counts']['non_essential']}`",
                "",
                "### Projection Parameters",
                "",
                f"- UMAP: `{json.dumps(row['umap_params'], sort_keys=True)}`",
                f"- t-SNE: `{json.dumps(row['tsne_params'], sort_keys=True)}`",
                "",
                "### Separation Metrics",
                "",
                f"- UMAP silhouette score: `{row['umap_metrics']['silhouette_score']}`",
                f"- UMAP Davies-Bouldin index: `{row['umap_metrics']['davies_bouldin_index']}`",
                f"- UMAP class centroid distance: `{row['umap_metrics']['centroid_distance']}`",
                f"- t-SNE silhouette score: `{row['tsne_metrics']['silhouette_score']}`",
                f"- t-SNE Davies-Bouldin index: `{row['tsne_metrics']['davies_bouldin_index']}`",
                f"- t-SNE class centroid distance: `{row['tsne_metrics']['centroid_distance']}`",
                f"- high-dimensional centroid distance: `{row['high_dim_metrics']['centroid_distance']}`",
                f"- high-dimensional within-class mean distance: `{row['high_dim_metrics']['within_class_mean_distance']}`",
                "",
                "### Outputs",
                "",
                f"- UMAP PDF: `{row['umap_pdf']}`",
                f"- t-SNE PDF: `{row['tsne_pdf']}`",
                f"- UMAP coords: `{row['umap_coords_tsv']}`",
                f"- t-SNE coords: `{row['tsne_coords_tsv']}`",
                f"- metadata JSON: `{row['metadata_json']}`",
                "",
            ]
        )

    if pairwise_rows:
        lines.extend(["## Pairwise Validation", ""])
        for row in pairwise_rows:
            lines.append(
                f"- {row['protocol']}: same split=`{row['same_split_manifest']}`, same label manifest=`{row['same_label_manifest']}`, "
                f"same node ids=`{row['same_node_ids']}`, same labels=`{row['same_labels']}`, same subset=`{row['subset']}`."
            )
        lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines).strip() + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    pdf_dir = output_root / "pdf"
    coords_dir = output_root / "coords"
    metadata_dir = output_root / "metadata"
    summary_dir = output_root / "summary"
    for path in [pdf_dir, coords_dir, metadata_dir, summary_dir]:
        path.mkdir(parents=True, exist_ok=True)

    case_rows = []
    pairwise_inputs = {}

    for protocol in args.protocols:
        for feature_setting in args.feature_settings:
            run_dir = locate_run_dir(protocol, args.model, feature_setting, args.seed, args.upstream_roots)
            resolved_config, config_path = load_run_context(run_dir, args.base_config)
            base_config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
            bundle = load_protocol_dataset(config_path, protocol, feature_setting)
            checkpoint = torch.load(run_dir / "best_model.pt", map_location="cpu")
            model_cfg = resolve_graph_model_config(base_config, normalize_model_name(args.model))
            x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
            edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
            hidden, logits = forward_penultimate_from_checkpoint(
                x,
                edge_index,
                checkpoint["state_dict"],
                str(model_cfg.get("aggregator_type", "mean")).strip().lower(),
            )

            selected_idx = subset_indices(bundle, args.subset)
            subset_df = bundle["node_manifest"].iloc[selected_idx].copy().reset_index(drop=True)
            subset_df["species"] = bundle["species"]
            subset_df["label"] = subset_df["label"].astype(int)
            subset_df["split"] = subset_df["split"].astype(str)
            embedding = hidden.detach().cpu().numpy()[selected_idx]
            labels = subset_df["label"].to_numpy(dtype=int)

            umap_coords, umap_params = fit_umap_coords(embedding)
            tsne_coords, tsne_params = fit_tsne_coords(embedding)
            umap_metrics = compute_embedding_metrics(umap_coords, labels)
            tsne_metrics = compute_embedding_metrics(tsne_coords, labels)
            high_dim_metrics = compute_embedding_metrics(embedding, labels)

            stem = output_stem(protocol, feature_setting, args.subset, args.seed, args.model)
            umap_pdf = pdf_dir / f"{stem}_umap.pdf"
            tsne_pdf = pdf_dir / f"{stem}_tsne.pdf"
            umap_coords_tsv = coords_dir / f"{stem}_umap_coords.tsv"
            tsne_coords_tsv = coords_dir / f"{stem}_tsne_coords.tsv"
            metadata_json = metadata_dir / f"{stem}_metadata.json"

            umap_df = save_coords_table(subset_df, umap_coords, umap_coords_tsv, protocol, feature_setting, args.model, args.seed)
            tsne_df = save_coords_table(subset_df, tsne_coords, tsne_coords_tsv, protocol, feature_setting, args.model, args.seed)
            save_plot(umap_df, "dim1", "dim2", "UMAP", umap_pdf)
            save_plot(tsne_df, "dim1", "dim2", "t-SNE", tsne_pdf)

            class_counts = {
                "essential": int((labels == 1).sum()),
                "non_essential": int((labels == 0).sum()),
            }
            row = {
                "protocol": protocol,
                "species": bundle["species"],
                "label_regime": protocol_label_regime(protocol),
                "feature_setting": feature_setting,
                "model": args.model,
                "seed": int(args.seed),
                "subset": args.subset,
                "uses_test_only_subset": args.subset == "test",
                "node_count": int(len(subset_df)),
                "class_counts": class_counts,
                "split_manifest_path": bundle["split_manifest_path"],
                "label_manifest_path": bundle["label_manifest_path"],
                "checkpoint_path": str((run_dir / "best_model.pt").resolve()),
                "resolved_config_path": str((run_dir / "resolved_config.yaml").resolve()) if (run_dir / "resolved_config.yaml").exists() else "",
                "run_dir": str(run_dir),
                "umap_params": umap_params,
                "tsne_params": tsne_params,
                "umap_metrics": umap_metrics,
                "tsne_metrics": tsne_metrics,
                "high_dim_metrics": high_dim_metrics,
                "umap_pdf": str(umap_pdf),
                "tsne_pdf": str(tsne_pdf),
                "umap_coords_tsv": str(umap_coords_tsv),
                "tsne_coords_tsv": str(tsne_coords_tsv),
                "metadata_json": str(metadata_json),
            }
            metadata_json.write_text(json.dumps(json_ready(row), indent=2, sort_keys=False) + "\n", encoding="utf-8")
            case_rows.append(row)
            pairwise_inputs.setdefault(protocol, {})[feature_setting] = {
                "split_manifest_path": bundle["split_manifest_path"],
                "label_manifest_path": bundle["label_manifest_path"],
                "subset": args.subset,
                "node_ids": subset_df["graph_gene_id"].astype(str).tolist(),
                "labels": subset_df["label"].astype(int).tolist(),
            }

    pairwise_rows = []
    for protocol, feature_map in sorted(pairwise_inputs.items()):
        if len(feature_map) < 2:
            continue
        baseline = feature_map[args.feature_settings[0]]
        compare = feature_map[args.feature_settings[1]]
        pairwise_rows.append(
            {
                "protocol": protocol,
                "subset": baseline["subset"],
                "same_split_manifest": baseline["split_manifest_path"] == compare["split_manifest_path"],
                "same_label_manifest": baseline["label_manifest_path"] == compare["label_manifest_path"],
                "same_node_ids": baseline["node_ids"] == compare["node_ids"],
                "same_labels": baseline["labels"] == compare["labels"],
            }
        )

    summary_path = summary_dir / "Figure3d_representation_summary.md"
    render_case_summary(case_rows, summary_path, pairwise_rows)

    manifest_rows = []
    for row in case_rows:
        manifest_rows.append(
            {
                "protocol": row["protocol"],
                "species": row["species"],
                "label_regime": row["label_regime"],
                "feature_setting": row["feature_setting"],
                "model": row["model"],
                "seed": row["seed"],
                "subset": row["subset"],
                "uses_test_only_subset": row["uses_test_only_subset"],
                "node_count": row["node_count"],
                "essential_count": row["class_counts"]["essential"],
                "non_essential_count": row["class_counts"]["non_essential"],
                "umap_silhouette_score": row["umap_metrics"]["silhouette_score"],
                "umap_davies_bouldin_index": row["umap_metrics"]["davies_bouldin_index"],
                "umap_centroid_distance": row["umap_metrics"]["centroid_distance"],
                "tsne_silhouette_score": row["tsne_metrics"]["silhouette_score"],
                "tsne_davies_bouldin_index": row["tsne_metrics"]["davies_bouldin_index"],
                "tsne_centroid_distance": row["tsne_metrics"]["centroid_distance"],
                "high_dim_centroid_distance": row["high_dim_metrics"]["centroid_distance"],
                "high_dim_within_class_mean_distance": row["high_dim_metrics"]["within_class_mean_distance"],
                "checkpoint_path": row["checkpoint_path"],
                "split_manifest_path": row["split_manifest_path"],
                "label_manifest_path": row["label_manifest_path"],
                "umap_pdf": row["umap_pdf"],
                "tsne_pdf": row["tsne_pdf"],
                "umap_coords_tsv": row["umap_coords_tsv"],
                "tsne_coords_tsv": row["tsne_coords_tsv"],
                "metadata_json": row["metadata_json"],
            }
        )
    pd.DataFrame(manifest_rows).to_csv(metadata_dir / "Figure3d_representation_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
