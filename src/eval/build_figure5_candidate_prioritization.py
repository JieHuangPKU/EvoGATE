import argparse
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rcParams
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import torch
import yaml

from src.analysis.figure5_final_common import load_runtime_seed_list
from src.analysis.figure4_representation_common import graphsage_penultimate_and_logits
from src.data.frozen_protocol_loader import load_protocol_dataset
from src.train.run_frozen_protocol_feature_combo_model import resolve_graph_model_config
from src.train.run_frozen_protocol_model import normalize_model_name


MODEL = "GraphSAGE"
PROTOCOL = "fgraminearum_newlabel"
SPECIES = "fgraminearum"
BASELINE = "ORT_EXP_SUB"
ESM2 = "ORT_EXP_SUB_ESM2"
TOP_KS = (100, 200, 500)
PREFERRED_SEED = 1029

ROOT_DEFAULT = "results/Figure5_new_candidate_prioritization"
COLORS = {
    "other": "#cfcfcf",
    "shared": "#2b7bba",
    "esm2_specific": "#2CA25F",
    "baseline_only": "#e66101",
    "essential": "#2F6F4E",
}

rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]


def parse_args():
    parser = argparse.ArgumentParser(description="Build redesigned Figure 5 candidate prioritization panels.")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml")
    parser.add_argument("--upstream-root", default="outputs/Figure3a")
    parser.add_argument("--output-root", default=ROOT_DEFAULT)
    parser.add_argument("--protocol", default=PROTOCOL)
    parser.add_argument("--species", default=SPECIES)
    parser.add_argument("--model", default=MODEL)
    parser.add_argument("--baseline-feature", default=BASELINE)
    parser.add_argument("--esm2-feature", default=ESM2)
    parser.add_argument("--seeds", nargs="*", type=int, default=None)
    parser.add_argument("--preferred-seed", type=int, default=PREFERRED_SEED)
    parser.add_argument("--default-top-k", type=int, default=500)
    return parser.parse_args()


def ensure_dirs(*paths):
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def write_tsv(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def save_pair(fig, pdf_path, png_path):
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def compute_umap(matrix, random_state):
    import umap.umap_ as umap

    scaled = StandardScaler().fit_transform(matrix)
    metric = "cosine" if matrix.shape[1] > 10 else "euclidean"
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=30,
        min_dist=0.1,
        metric=metric,
        random_state=random_state,
    )
    return reducer.fit_transform(scaled), metric


def archive_legacy_hidden_umap(output_root):
    archive_dir = Path(output_root) / "archive_old_hidden_umap"
    archive_dir.mkdir(parents=True, exist_ok=True)
    legacy_roots = [
        Path("results/Figure5/plots"),
        Path("results/Figure5/data"),
        Path("results/Figure5/summary"),
        Path("results/Figure5/tables"),
    ]
    copied = []
    for root in legacy_roots:
        if not root.exists():
            continue
        for pattern in [
            "Figure5_new_A_hidden_space_rescue_umap*",
            "Figure5a_hidden_umap_error_transition_fgraminearum_seed1029*",
        ]:
            for src in root.glob(pattern):
                dst = archive_dir / src.name
                shutil.copy2(src, dst)
                copied.append(dst)
    note_path = archive_dir / "README.md"
    lines = [
        "# Archived Figure5A hidden-space outputs",
        "",
        "- These files were copied from legacy hidden-space / FN_to_TP-rescue Figure5A outputs.",
        "- The active Figure5A manuscript output is now prediction-space UMAP under this results root.",
    ]
    if copied:
        lines.append("")
        lines.append("## Archived files")
        lines.extend([f"- `{path.name}`" for path in sorted(copied)])
    note_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return archive_dir, copied


def resolve_run_seed(upstream_root, protocol, model, baseline_feature, esm2_feature, preferred_seed, runtime_seeds):
    available = []
    for seed in sorted(set(int(seed) for seed in runtime_seeds)):
        base = Path(upstream_root) / protocol / model / baseline_feature / f"run_{seed}" / "predictions.tsv"
        esm2 = Path(upstream_root) / protocol / model / esm2_feature / f"run_{seed}" / "predictions.tsv"
        if base.exists() and esm2.exists():
            available.append(seed)
    if not available:
        raise FileNotFoundError("No paired baseline/ESM2 prediction runs were found.")
    chosen = preferred_seed if preferred_seed in available else available[0]
    return chosen, available


def load_prediction(upstream_root, protocol, model, feature_setting, seed):
    path = Path(upstream_root) / protocol / model / feature_setting / f"run_{seed}" / "predictions.tsv"
    df = read_tsv(path).copy()
    if "pred_score" not in df.columns:
        raise ValueError(f"Missing pred_score in {path}")
    df["gene_id"] = df["graph_gene_id"].astype(str)
    df["canonical_gene_id"] = df["canonical_gene_id"].astype(str)
    df["score"] = pd.to_numeric(df["pred_score"], errors="raise").astype(float)
    return df, path


def rank_table_for_scores(df, score_col, prefix):
    ranked = df[["gene_id", score_col]].copy()
    ranked = ranked.sort_values([score_col, "gene_id"], ascending=[False, True]).reset_index(drop=True)
    ranked[f"rank_{prefix}"] = np.arange(1, len(ranked) + 1, dtype=int)
    if len(ranked) > 1:
        ranked[f"percentile_{prefix}"] = 1.0 - ((ranked[f"rank_{prefix}"] - 1) / float(len(ranked) - 1))
    else:
        ranked[f"percentile_{prefix}"] = 1.0
    return ranked


def build_prediction_space_seed_table(args, seed):
    base_df, base_path = load_prediction(args.upstream_root, args.protocol, args.model, args.baseline_feature, seed)
    esm2_df, esm2_path = load_prediction(args.upstream_root, args.protocol, args.model, args.esm2_feature, seed)
    base_df = base_df.rename(
        columns={
            "score": "score_ORT_EXP_SUB",
            "pred_label": "pred_label_ORT_EXP_SUB",
        }
    )
    esm2_df = esm2_df.rename(
        columns={
            "score": "score_ORT_EXP_SUB_ESM2",
            "pred_label": "pred_label_ORT_EXP_SUB_ESM2",
        }
    )
    merged = base_df.merge(
        esm2_df[
            [
                "gene_id",
                "canonical_gene_id",
                "split",
                "label",
                "is_labeled",
                "in_graph",
                "score_ORT_EXP_SUB_ESM2",
                "pred_label_ORT_EXP_SUB_ESM2",
            ]
        ],
        on="gene_id",
        how="inner",
        validate="one_to_one",
    )
    merged["canonical_gene_id"] = merged["canonical_gene_id_x"].where(
        merged["canonical_gene_id_x"].astype(str).str.len() > 0,
        merged["canonical_gene_id_y"],
    )
    merged["split"] = merged["split_x"].where(
        merged["split_x"].astype(str).str.len() > 0,
        merged["split_y"],
    )
    merged["label"] = pd.to_numeric(merged["label_x"], errors="coerce").where(
        pd.to_numeric(merged["label_x"], errors="coerce").notna(),
        pd.to_numeric(merged["label_y"], errors="coerce"),
    )
    merged["is_labeled"] = merged["is_labeled_x"].fillna(merged["is_labeled_y"])
    merged["in_graph"] = merged["in_graph_x"].fillna(merged["in_graph_y"])
    merged["score_ORT_EXP_SUB"] = pd.to_numeric(merged["score_ORT_EXP_SUB"], errors="raise").astype(float)
    merged["score_ORT_EXP_SUB_ESM2"] = pd.to_numeric(merged["score_ORT_EXP_SUB_ESM2"], errors="raise").astype(float)

    base_rank = rank_table_for_scores(merged, "score_ORT_EXP_SUB", "ORT_EXP_SUB")
    esm2_rank = rank_table_for_scores(merged, "score_ORT_EXP_SUB_ESM2", "ORT_EXP_SUB_ESM2")
    merged = merged.merge(
        base_rank[["gene_id", "rank_ORT_EXP_SUB", "percentile_ORT_EXP_SUB"]],
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    merged = merged.merge(
        esm2_rank[["gene_id", "rank_ORT_EXP_SUB_ESM2", "percentile_ORT_EXP_SUB_ESM2"]],
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    merged["delta_score"] = merged["score_ORT_EXP_SUB_ESM2"] - merged["score_ORT_EXP_SUB"]
    merged["delta_rank"] = merged["rank_ORT_EXP_SUB"] - merged["rank_ORT_EXP_SUB_ESM2"]
    merged["true_label"] = np.where(merged["is_labeled"].fillna(False), merged["label"], np.nan)
    merged["seed"] = int(seed)
    merged["prediction_path_ORT_EXP_SUB"] = str(base_path)
    merged["prediction_path_ORT_EXP_SUB_ESM2"] = str(esm2_path)
    return merged[
        [
            "gene_id",
            "canonical_gene_id",
            "split",
            "label",
            "true_label",
            "is_labeled",
            "in_graph",
            "score_ORT_EXP_SUB",
            "score_ORT_EXP_SUB_ESM2",
            "rank_ORT_EXP_SUB",
            "rank_ORT_EXP_SUB_ESM2",
            "percentile_ORT_EXP_SUB",
            "percentile_ORT_EXP_SUB_ESM2",
            "delta_score",
            "delta_rank",
            "seed",
            "prediction_path_ORT_EXP_SUB",
            "prediction_path_ORT_EXP_SUB_ESM2",
            "pred_label_ORT_EXP_SUB",
            "pred_label_ORT_EXP_SUB_ESM2",
        ]
    ].copy()


def find_embedding_candidates(run_dir):
    patterns = [
        "hidden_embeddings.tsv",
        "embeddings.tsv",
        "latent.tsv",
        "node_embeddings.tsv",
        "hidden_representation.tsv",
        "hidden_embedding.tsv",
    ]
    found = []
    for pattern in patterns:
        found.extend(sorted(run_dir.glob(pattern)))
    return found


def load_checkpoint_hidden_embedding(args, feature_setting, seed):
    run_dir = Path(args.upstream_root) / args.protocol / args.model / feature_setting / f"run_{seed}"
    candidate_files = find_embedding_candidates(run_dir)
    checkpoint_path = run_dir / "best_model.pt"
    if not checkpoint_path.exists():
        raise FileNotFoundError(
            "Hidden embedding not available: missing checkpoint at {0}. Candidate files: {1}".format(
                checkpoint_path, ", ".join(str(path) for path in candidate_files) if candidate_files else "none"
            )
        )
    try:
        bundle = load_protocol_dataset(args.runtime_config, args.protocol, feature_setting)
        config = yaml.safe_load(Path(args.runtime_config).read_text(encoding="utf-8"))
        model_cfg = resolve_graph_model_config(config, normalize_model_name(args.model))
        checkpoint = torch.load(str(checkpoint_path), map_location="cpu")
        x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
        edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
        hidden, _logits = graphsage_penultimate_and_logits(x, edge_index, checkpoint["state_dict"], model_cfg["aggregator_type"])
        node_df = bundle["node_manifest"].copy().reset_index(drop=True)
        node_df["gene_id"] = node_df["graph_gene_id"].astype(str)
        hidden_df = node_df[["gene_id"]].copy()
        hidden_matrix = hidden.detach().cpu().numpy().astype(np.float32, copy=False)
        hidden_df["embedding_index"] = np.arange(hidden_matrix.shape[0], dtype=int)
    except Exception as exc:
        raise RuntimeError(
            "Hidden embedding not available for Figure5A. Tried checkpoint-derived penultimate embedding from {0}. "
            "Run dir: {1}. Candidate files: {2}. Original error: {3}".format(
                checkpoint_path,
                run_dir,
                ", ".join(str(path) for path in candidate_files) if candidate_files else "none",
                exc,
            )
        ) from exc
    return {
        "run_dir": run_dir,
        "checkpoint_path": checkpoint_path,
        "candidate_files": candidate_files,
        "embedding_df": hidden_df,
        "hidden_matrix": hidden_matrix,
        "embedding_dim": int(hidden_matrix.shape[1]),
        "source_kind": "checkpoint_penultimate_hidden",
        "source_label": str(checkpoint_path),
    }


def topk_definitions(n_genes):
    top5 = max(1, int(math.ceil(n_genes * 0.05)))
    return [("100", 100), ("200", 200), ("500", 500), ("5%", top5)]


def add_topk_membership(df):
    n = len(df)
    for key, k in topk_definitions(n):
        safe = min(k, n)
        suffix = "top5percent" if key == "5%" else f"top{key}"
        df[f"in_{suffix}_ORT_EXP_SUB"] = df["rank_ORT_EXP_SUB"] <= safe
        df[f"in_{suffix}_ORT_EXP_SUB_ESM2"] = df["rank_ORT_EXP_SUB_ESM2"] <= safe
    return df


def build_figure5a_outputs(args, output_root, seed):
    df = build_prediction_space_seed_table(args, seed)
    df = add_topk_membership(df)
    baseline_embedding_info = load_checkpoint_hidden_embedding(args, args.baseline_feature, seed)
    esm2_embedding_info = load_checkpoint_hidden_embedding(args, args.esm2_feature, seed)

    df = df.merge(
        baseline_embedding_info["embedding_df"].rename(columns={"embedding_index": "embedding_index_ORT_EXP_SUB"}),
        on="gene_id",
        how="inner",
        validate="one_to_one",
    )
    df = df.merge(
        esm2_embedding_info["embedding_df"].rename(columns={"embedding_index": "embedding_index_ORT_EXP_SUB_ESM2"}),
        on="gene_id",
        how="inner",
        validate="one_to_one",
    ).sort_values("embedding_index_ORT_EXP_SUB_ESM2").reset_index(drop=True)
    if df.empty:
        raise RuntimeError(
            "Hidden embedding merge produced zero aligned genes. Candidate files baseline: {0}; ESM2: {1}".format(
                ", ".join(str(path) for path in baseline_embedding_info["candidate_files"]) if baseline_embedding_info["candidate_files"] else "none",
                ", ".join(str(path) for path in esm2_embedding_info["candidate_files"]) if esm2_embedding_info["candidate_files"] else "none",
            )
        )
    baseline_hidden_matrix = baseline_embedding_info["hidden_matrix"][df["embedding_index_ORT_EXP_SUB"].to_numpy(dtype=int)]
    esm2_hidden_matrix = esm2_embedding_info["hidden_matrix"][df["embedding_index_ORT_EXP_SUB_ESM2"].to_numpy(dtype=int)]
    baseline_coords, baseline_umap_metric = compute_umap(baseline_hidden_matrix, random_state=seed)
    esm2_coords, esm2_umap_metric = compute_umap(esm2_hidden_matrix, random_state=seed)
    df["UMAP1_ORT_EXP_SUB"] = baseline_coords[:, 0]
    df["UMAP2_ORT_EXP_SUB"] = baseline_coords[:, 1]
    df["UMAP1_ORT_EXP_SUB_ESM2"] = esm2_coords[:, 0]
    df["UMAP2_ORT_EXP_SUB_ESM2"] = esm2_coords[:, 1]

    summary_rows = []
    for key, k in topk_definitions(len(df)):
        safe = min(k, len(df))
        base_mask = df["rank_ORT_EXP_SUB"] <= safe
        esm2_mask = df["rank_ORT_EXP_SUB_ESM2"] <= safe
        shared_mask = base_mask & esm2_mask
        esm2_specific_mask = (~base_mask) & esm2_mask
        lost_mask = base_mask & (~esm2_mask)
        esm2_specific_delta = pd.to_numeric(df.loc[esm2_specific_mask, "delta_score"], errors="coerce")
        denom = int(shared_mask.sum() + esm2_specific_mask.sum() + lost_mask.sum())
        summary_rows.append(
            {
                "K": key,
                "n_top_noESM2": int(base_mask.sum()),
                "n_top_ESM2": int(esm2_mask.sum()),
                "n_shared": int(shared_mask.sum()),
                "n_ESM2_specific": int(esm2_specific_mask.sum()),
                "n_lost_after_ESM2": int(lost_mask.sum()),
                "jaccard": float(shared_mask.sum() / denom) if denom else 0.0,
                "median_delta_score_ESM2_specific": float(esm2_specific_delta.median()) if esm2_specific_delta.notna().any() else np.nan,
                "mean_delta_score_ESM2_specific": float(esm2_specific_delta.mean()) if esm2_specific_delta.notna().any() else np.nan,
            }
        )
        suffix = "top5percent" if key == "5%" else f"top{key}"
        df[f"is_shared_{suffix}"] = shared_mask
        df[f"is_esm2_specific_{suffix}"] = esm2_specific_mask
        df[f"is_lost_{suffix}"] = lost_mask

    summary_df = pd.DataFrame(summary_rows)
    main_suffix = f"top{args.default_top_k}"
    df["class_A1"] = np.where(df[f"in_{main_suffix}_ORT_EXP_SUB"], "Predicted essential", "Other genes")
    df["class_A2"] = np.select(
        [
            df[f"is_esm2_specific_{main_suffix}"],
            df[f"is_shared_{main_suffix}"],
        ],
        [
            "ESM2-specific essential",
            "Shared essential",
        ],
        default="Other genes",
    )

    plot_data_path = Path(output_root) / "Figure5A_prediction_umap_plot_data.tsv"
    summary_path = Path(output_root) / "Figure5A_prediction_shift_summary.tsv"
    shared_path = Path(output_root) / "Figure5A_shared_essential_genes.tsv"
    esm2_specific_path = Path(output_root) / "Figure5A_ESM2_specific_essential_genes.tsv"
    essential_pdf_path = Path(output_root) / "Figure5A_prediction_umap_essential_vs_nonessential.pdf"
    essential_png_path = Path(output_root) / "Figure5A_prediction_umap_essential_vs_nonessential.png"
    specific_pdf_path = Path(output_root) / "Figure5A_prediction_umap_esm2_specific.pdf"
    specific_png_path = Path(output_root) / "Figure5A_prediction_umap_esm2_specific.png"
    two_panel_pdf_path = Path(output_root) / "Figure5A_prediction_space_umap_two_panel.pdf"
    two_panel_png_path = Path(output_root) / "Figure5A_prediction_space_umap_two_panel.png"
    audit_path = Path(output_root) / "Figure5A_prediction_space_audit.md"

    write_tsv(df, plot_data_path)
    write_tsv(summary_df, summary_path)

    esm2_specific = df[df[f"is_esm2_specific_{main_suffix}"]].copy()
    esm2_specific = esm2_specific.sort_values(["delta_score", "rank_ORT_EXP_SUB_ESM2", "gene_id"], ascending=[False, True, True]).reset_index(drop=True)
    shared = df[df[f"is_shared_{main_suffix}"]].copy()
    shared = shared.sort_values(["rank_ORT_EXP_SUB_ESM2", "rank_ORT_EXP_SUB", "gene_id"], ascending=[True, True, True]).reset_index(drop=True)

    for subset_df, mode in [(esm2_specific, "esm2_specific"), (shared, "shared")]:
        subset_df["in_top100"] = subset_df["is_esm2_specific_top100"] if mode == "esm2_specific" else subset_df["is_shared_top100"]
        subset_df["in_top200"] = subset_df["is_esm2_specific_top200"] if mode == "esm2_specific" else subset_df["is_shared_top200"]
        subset_df["in_top500"] = subset_df["is_esm2_specific_top500"] if mode == "esm2_specific" else subset_df["is_shared_top500"]
        subset_df["in_top5percent"] = subset_df["is_esm2_specific_top5percent"] if mode == "esm2_specific" else subset_df["is_shared_top5percent"]

    gene_cols = [
        "gene_id",
        "score_ORT_EXP_SUB",
        "score_ORT_EXP_SUB_ESM2",
        "rank_ORT_EXP_SUB",
        "rank_ORT_EXP_SUB_ESM2",
        "percentile_ORT_EXP_SUB",
        "percentile_ORT_EXP_SUB_ESM2",
        "delta_score",
        "delta_rank",
        "true_label",
        "in_top100",
        "in_top200",
        "in_top500",
        "in_top5percent",
    ]
    write_tsv(esm2_specific[gene_cols], esm2_specific_path)
    write_tsv(shared[gene_cols], shared_path)

    base_xlim = (float(df["UMAP1_ORT_EXP_SUB"].min()), float(df["UMAP1_ORT_EXP_SUB"].max()))
    base_ylim = (float(df["UMAP2_ORT_EXP_SUB"].min()), float(df["UMAP2_ORT_EXP_SUB"].max()))
    esm2_xlim = (float(df["UMAP1_ORT_EXP_SUB_ESM2"].min()), float(df["UMAP1_ORT_EXP_SUB_ESM2"].max()))
    esm2_ylim_lower = 4.0
    esm2_ylim_upper = float(df["UMAP2_ORT_EXP_SUB_ESM2"].max())
    esm2_display_ylim = (esm2_ylim_lower, esm2_ylim_upper)

    def style_axis(ax, title, xlim, ylim):
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("UMAP1", fontsize=10)
        ax.set_ylabel("UMAP2", fontsize=10)
        ax.set_facecolor("white")
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.7)
        ax.spines["bottom"].set_linewidth(0.7)
        ax.tick_params(width=0.7, labelsize=9)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

    def draw_essential_panel(ax):
        other = df[df["pred_label_ORT_EXP_SUB_ESM2"].astype(int) == 0].copy()
        essential = df[df["pred_label_ORT_EXP_SUB_ESM2"].astype(int) == 1].copy()
        ax.scatter(other["UMAP1_ORT_EXP_SUB_ESM2"], other["UMAP2_ORT_EXP_SUB_ESM2"], s=1.0, c="#4C78A8", alpha=0.65, linewidths=0)
        ax.scatter(essential["UMAP1_ORT_EXP_SUB_ESM2"], essential["UMAP2_ORT_EXP_SUB_ESM2"], s=1.8, c="#D81B60", alpha=0.90, linewidths=0)
        style_axis(ax, "Prediction UMAP of essential and non-essential genes\nGraphSAGE with ORT+EXP+SUB+ESM2", esm2_xlim, esm2_display_ylim)

    def draw_specific_panel(ax):
        bg = df.copy()
        shared_fg = shared.copy()
        rescued_fg = esm2_specific.copy()
        ax.scatter(bg["UMAP1_ORT_EXP_SUB_ESM2"], bg["UMAP2_ORT_EXP_SUB_ESM2"], s=1.0, c=COLORS["other"], alpha=0.32, linewidths=0)
        ax.scatter(shared_fg["UMAP1_ORT_EXP_SUB_ESM2"], shared_fg["UMAP2_ORT_EXP_SUB_ESM2"], s=1.8, c="#D81B60", alpha=0.90, linewidths=0)
        ax.scatter(rescued_fg["UMAP1_ORT_EXP_SUB_ESM2"], rescued_fg["UMAP2_ORT_EXP_SUB_ESM2"], s=2.2, c=COLORS["esm2_specific"], alpha=0.95, linewidths=0)
        style_axis(ax, "Prediction UMAP of essential and non-essential genes\nGraphSAGE with ORT+EXP+SUB+ESM2", esm2_xlim, esm2_display_ylim)

    def draw_two_panel_left(ax):
        other = df[df["pred_label_ORT_EXP_SUB"].astype(int) == 0].copy()
        essential = df[df["pred_label_ORT_EXP_SUB"].astype(int) == 1].copy()
        ax.scatter(other["UMAP1_ORT_EXP_SUB"], other["UMAP2_ORT_EXP_SUB"], s=1.0, c="#4C78A8", alpha=0.65, linewidths=0)
        ax.scatter(essential["UMAP1_ORT_EXP_SUB"], essential["UMAP2_ORT_EXP_SUB"], s=1.8, c="#D81B60", alpha=0.90, linewidths=0)
        style_axis(ax, "ORT+EXP+SUB", base_xlim, base_ylim)

    def draw_two_panel_right(ax):
        other = df[df["pred_label_ORT_EXP_SUB_ESM2"].astype(int) == 0].copy()
        essential = df[df["pred_label_ORT_EXP_SUB_ESM2"].astype(int) == 1].copy()
        rescued_fg = esm2_specific.copy()
        essential_regular = essential[~essential["gene_id"].isin(set(rescued_fg["gene_id"]))].copy()
        ax.scatter(other["UMAP1_ORT_EXP_SUB_ESM2"], other["UMAP2_ORT_EXP_SUB_ESM2"], s=1.0, c="#4C78A8", alpha=0.65, linewidths=0)
        ax.scatter(essential_regular["UMAP1_ORT_EXP_SUB_ESM2"], essential_regular["UMAP2_ORT_EXP_SUB_ESM2"], s=1.8, c="#D81B60", alpha=0.90, linewidths=0)
        ax.scatter(rescued_fg["UMAP1_ORT_EXP_SUB_ESM2"], rescued_fg["UMAP2_ORT_EXP_SUB_ESM2"], s=2.2, c=COLORS["esm2_specific"], alpha=0.95, linewidths=0)
        style_axis(ax, "ORT+EXP+SUB+ESM2", esm2_xlim, esm2_display_ylim)

    fig, ax = plt.subplots(1, 1, figsize=(8.8, 4.6), facecolor="white")
    draw_essential_panel(ax)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="", markersize=4.8, markerfacecolor="#4C78A8", markeredgecolor="none", label="Predicted non-essential"),
            Line2D([0], [0], marker="o", linestyle="", markersize=5.6, markerfacecolor="#D81B60", markeredgecolor="none", label="Predicted essential"),
        ],
        loc="best",
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    save_pair(fig, essential_pdf_path, essential_png_path)

    fig, ax = plt.subplots(1, 1, figsize=(8.8, 4.6), facecolor="white")
    draw_specific_panel(ax)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="", markersize=4.8, markerfacecolor=COLORS["shared"], markeredgecolor="none", label="Shared Top500 essential"),
            Line2D([0], [0], marker="o", linestyle="", markersize=5.4, markerfacecolor=COLORS["esm2_specific"], markeredgecolor="none", label="ESM2-specific Top500 essential"),
            Line2D([0], [0], marker="o", linestyle="", markersize=4.2, markerfacecolor=COLORS["other"], markeredgecolor="none", label="Other genes"),
        ],
        loc="best",
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    save_pair(fig, specific_pdf_path, specific_png_path)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6), facecolor="white")
    draw_two_panel_left(axes[0])
    draw_two_panel_right(axes[1])
    fig.suptitle("Prediction UMAP reveals ESM2-specific essential candidates", fontsize=12, y=1.02)
    axes[1].legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="", markersize=4.8, markerfacecolor="#4C78A8", markeredgecolor="none", label="Predicted non-essential"),
            Line2D([0], [0], marker="o", linestyle="", markersize=5.4, markerfacecolor="#D81B60", markeredgecolor="none", label="Predicted essential"),
            Line2D([0], [0], marker="o", linestyle="", markersize=5.8, markerfacecolor=COLORS["esm2_specific"], markeredgecolor="none", label="ESM2-specific essential"),
        ],
        loc="best",
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    save_pair(fig, two_panel_pdf_path, two_panel_png_path)

    top20 = esm2_specific.head(20)[["gene_id", "score_ORT_EXP_SUB", "score_ORT_EXP_SUB_ESM2", "delta_score", "delta_rank", "rank_ORT_EXP_SUB", "rank_ORT_EXP_SUB_ESM2"]]
    audit_lines = [
        "# Figure5A prediction UMAP audit",
        "",
        f"- Input file path ORT_EXP_SUB: `{df['prediction_path_ORT_EXP_SUB'].iloc[0]}`.",
        f"- Input file path ORT_EXP_SUB_ESM2: `{df['prediction_path_ORT_EXP_SUB_ESM2'].iloc[0]}`.",
        f"- Using seed/run: `run_{seed}`.",
        f"- UMAP input file: baseline `{baseline_embedding_info['source_label']}`; ESM2 `{esm2_embedding_info['source_label']}`.",
        f"- UMAP input dimensionality: baseline `{baseline_hidden_matrix.shape[0]} x {baseline_hidden_matrix.shape[1]}`; ESM2 `{esm2_hidden_matrix.shape[0]} x {esm2_hidden_matrix.shape[1]}`.",
        "- Score column name used only for coloring/highlight and Top-K definition: `pred_score`.",
        f"- Aligned gene count: `{len(df)}`.",
        "- UMAP feature columns: `checkpoint-derived GraphSAGE penultimate hidden embedding`.",
        f"- UMAP metric: baseline `{baseline_umap_metric}`; ESM2 `{esm2_umap_metric}`.",
        f"- K definition for main figure: `Top{args.default_top_k}` by descending prediction score within each model.",
        "- Additional summaries also include `Top100`, `Top200`, and `Top5%`.",
        "- This figure uses hidden embedding for UMAP coordinates and prediction score only for coloring/highlight.",
        "- Hidden embedding source: `checkpoint_penultimate_hidden` for both panels.",
        "- ESM2-specific essential candidates colored green (`#2CA25F`).",
        f"- Left ORT+EXP+SUB panel keeps original y-axis limits `{base_ylim}`.",
        f"- Right ORT+EXP+SUB+ESM2 panel display is restricted to `UMAP2 >= {esm2_ylim_lower}` with y-axis limits `{esm2_display_ylim}`.",
        "- Underlying UMAP coordinates are unchanged unless a filtered plot version is explicitly used.",
        "- No score/rank/delta fallback was used.",
        "",
        "## Output files",
        "",
        f"- `{essential_pdf_path}`",
        f"- `{essential_png_path}`",
        f"- `{specific_pdf_path}`",
        f"- `{specific_png_path}`",
        f"- `{two_panel_pdf_path}`",
        f"- `{two_panel_png_path}`",
        f"- `{summary_path}`",
        f"- `{esm2_specific_path}`",
        f"- `{shared_path}`",
        f"- `{plot_data_path}`",
        "",
        "## Top 20 ESM2-specific genes by delta_score",
        "",
        top20.to_markdown(index=False),
        "",
    ]
    audit_path.write_text("\n".join(audit_lines), encoding="utf-8")

    return {
        "seed_df": df,
        "summary_df": summary_df,
        "esm2_specific_df": esm2_specific[gene_cols].copy(),
        "shared_df": shared[gene_cols].copy(),
        "plot_data_path": plot_data_path,
        "summary_path": summary_path,
        "shared_path": shared_path,
        "esm2_specific_path": esm2_specific_path,
        "essential_pdf_path": essential_pdf_path,
        "essential_png_path": essential_png_path,
        "specific_pdf_path": specific_pdf_path,
        "specific_png_path": specific_png_path,
        "two_panel_pdf_path": two_panel_pdf_path,
        "two_panel_png_path": two_panel_png_path,
        "audit_path": audit_path,
        "embedding_dim": int(esm2_hidden_matrix.shape[1]),
        "embedding_source": esm2_embedding_info["source_label"],
        "umap_metric": esm2_umap_metric,
    }


def build_model_ranking(upstream_root, protocol, model, feature_setting, seeds):
    rows = []
    for seed in seeds:
        df, path = load_prediction(upstream_root, protocol, model, feature_setting, seed)
        ranked = rank_table_for_scores(df.rename(columns={"score": "score_seed"}), "score_seed", "seed")
        df = df.merge(ranked[["gene_id", "rank_seed", "percentile_seed"]], on="gene_id", how="left", validate="one_to_one")
        df["seed"] = int(seed)
        df["prediction_path"] = str(path)
        rows.append(df)
    per_seed = pd.concat(rows, ignore_index=True)
    ranking = (
        per_seed.groupby("gene_id", as_index=False)
        .agg(
            canonical_gene_id=("canonical_gene_id", "first"),
            split=("split", "first"),
            label=("label", "first"),
            is_labeled=("is_labeled", "first"),
            in_graph=("in_graph", "first"),
            final_score=("percentile_seed", "mean"),
            mean_raw_score=("score", "mean"),
            n_seeds=("seed", "nunique"),
        )
        .sort_values(["final_score", "gene_id"], ascending=[False, True])
        .reset_index(drop=True)
    )
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    return ranking, per_seed


def build_candidate_table(args, seeds):
    base, base_seed = build_model_ranking(args.upstream_root, args.protocol, args.model, args.baseline_feature, seeds)
    esm2, esm2_seed = build_model_ranking(args.upstream_root, args.protocol, args.model, args.esm2_feature, seeds)
    merged = base.merge(
        esm2[["gene_id", "final_score", "mean_raw_score", "rank", "n_seeds"]],
        on="gene_id",
        how="inner",
        validate="one_to_one",
        suffixes=("_baseline", "_ESM2"),
    )
    merged = merged.rename(
        columns={
            "final_score_baseline": "score_baseline",
            "final_score_ESM2": "score_ESM2",
            "mean_raw_score_baseline": "mean_raw_score_baseline",
            "mean_raw_score_ESM2": "mean_raw_score_ESM2",
            "rank_baseline": "rank_baseline",
            "rank_ESM2": "rank_ESM2",
        }
    )
    merged["delta_score"] = merged["score_ESM2"] - merged["score_baseline"]
    merged["percentile_gain"] = merged["delta_score"]
    merged["rank_gain"] = merged["rank_baseline"] - merged["rank_ESM2"]
    for k in TOP_KS:
        merged[f"top{k}_baseline"] = merged["rank_baseline"] <= k
        merged[f"top{k}_ESM2"] = merged["rank_ESM2"] <= k
        merged[f"candidate_class_top{k}"] = np.select(
            [
                merged[f"top{k}_baseline"] & merged[f"top{k}_ESM2"],
                ~merged[f"top{k}_baseline"] & merged[f"top{k}_ESM2"],
                merged[f"top{k}_baseline"] & ~merged[f"top{k}_ESM2"],
            ],
            ["Shared", "ESM2_unique", "Baseline_only"],
            default="Neither",
        )
    final_cols = [
        "gene_id",
        "canonical_gene_id",
        "split",
        "label",
        "is_labeled",
        "in_graph",
        "score_baseline",
        "score_ESM2",
        "mean_raw_score_baseline",
        "mean_raw_score_ESM2",
        "rank_baseline",
        "rank_ESM2",
        "delta_score",
        "percentile_gain",
        "rank_gain",
        "n_seeds_baseline",
        "n_seeds_ESM2",
    ]
    for k in TOP_KS:
        final_cols.extend([f"top{k}_baseline", f"top{k}_ESM2", f"candidate_class_top{k}"])
    merged = merged[final_cols].sort_values(["rank_ESM2", "gene_id"]).reset_index(drop=True)
    base_seed["model_feature"] = args.baseline_feature
    esm2_seed["model_feature"] = args.esm2_feature
    return merged, pd.concat([base_seed, esm2_seed], ignore_index=True)


def build_overlap_summary(candidate):
    rows = []
    for k in TOP_KS:
        shared = int(((candidate[f"top{k}_baseline"]) & (candidate[f"top{k}_ESM2"])).sum())
        esm2_unique = int(((~candidate[f"top{k}_baseline"]) & (candidate[f"top{k}_ESM2"])).sum())
        baseline_only = int(((candidate[f"top{k}_baseline"]) & (~candidate[f"top{k}_ESM2"])).sum())
        denom = shared + esm2_unique + baseline_only
        rows.append(
            {
                "top_k": k,
                "shared_candidates": shared,
                "ESM2_unique_candidates": esm2_unique,
                "baseline_only_candidates": baseline_only,
                "jaccard": shared / float(denom) if denom else 0.0,
            }
        )
    return pd.DataFrame(rows)


def build_rescued(candidate, top_k):
    rescued = candidate[
        (~candidate[f"top{top_k}_baseline"])
        & (candidate[f"top{top_k}_ESM2"])
        & ((candidate["rank_gain"] >= 500) | (candidate["percentile_gain"] >= 0.2))
    ].copy()
    rescued["rescued_definition"] = (
        f"not Top{top_k} baseline; Top{top_k} ESM2; rank_gain >= 500 or percentile_gain >= 0.2"
    )
    return rescued.sort_values(["rank_gain", "delta_score"], ascending=[False, False]).reset_index(drop=True)


def plot_panel_b(output_root):
    coords_path = Path("results/Figure5/data/Figure5c_input_vs_hidden_compare_fgraminearum_seed1029_coords.tsv")
    labels_path = Path("results/Figure5/data/Figure5a_hidden_umap_error_transition_fgraminearum_seed1029_coords.tsv")
    if not coords_path.exists() or not labels_path.exists():
        return None
    coords = read_tsv(coords_path)
    labels = read_tsv(labels_path)[["node_id", "transition"]].drop_duplicates()
    coords = coords.merge(labels, on="node_id", how="left")
    coords["is_rescued"] = coords["transition"].astype(str).eq("FN_to_TP_rescued")
    out_data = Path(output_root) / "Figure5_new_B_input_vs_hidden_umap_plot_data.tsv"
    write_tsv(coords, out_data)

    spaces = [space for space in ["bio_input", "bio_hidden"] if space in set(coords["space"])]
    if not spaces:
        spaces = sorted(coords["space"].dropna().astype(str).unique())
    fig, axes = plt.subplots(1, len(spaces), figsize=(5.1 * len(spaces), 4.4), squeeze=False, facecolor="white")
    for ax, space in zip(axes[0], spaces):
        sub = coords[coords["space"] == space].copy()
        bg = sub[~sub["is_rescued"]]
        rescued = sub[sub["is_rescued"]]
        ax.scatter(bg["umap1"], bg["umap2"], s=7, c="#C8C8C8", alpha=0.45, linewidths=0)
        ax.scatter(rescued["umap1"], rescued["umap2"], s=28, c=COLORS["esm2_specific"], edgecolors="white", linewidths=0.35)
        title = "Input representation" if space == "bio_input" else "Hidden representation"
        ax.set_title(title)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Rescued genes shift from input to hidden representation", y=1.02)
    pdf_path = Path(output_root) / "Figure5_new_B_input_vs_hidden_umap.pdf"
    png_path = Path(output_root) / "Figure5_new_B_input_vs_hidden_umap.png"
    save_pair(fig, pdf_path, png_path)
    return {
        "inputs": [coords_path, labels_path],
        "outputs": [pdf_path, png_path, out_data],
    }


def build_grouped_importance(args, output_root):
    global_path = Path("results/Figure5/tables/Figure5d_group_ablation_global_summary.tsv")
    schema_base = Path(args.upstream_root) / args.protocol / args.model / args.baseline_feature / "run_1029" / "feature_schema.tsv"
    schema_esm2 = Path(args.upstream_root) / args.protocol / args.model / args.esm2_feature / "run_1029" / "feature_schema.tsv"
    group_map = {
        "orthologs": "ORT",
        "expression": "EXP",
        "sublocalization": "SUB",
        "esm2": "ESM2",
        "degree": "Graph topology / degree",
    }
    rows = []
    for model_label, schema_path, has_attr in [
        ("GraphSAGE + ORT_EXP_SUB", schema_base, False),
        ("GraphSAGE + ORT_EXP_SUB_ESM2", schema_esm2, True),
    ]:
        schema = read_tsv(schema_path)
        for item in schema.itertuples(index=False):
            rows.append(
                {
                    "model": model_label,
                    "feature_group": group_map.get(str(item.feature_block), str(item.feature_block)),
                    "source_feature_block": str(item.feature_block),
                    "n_features": int(item.dimension),
                    "metric": "delta_auprc",
                    "grouped_importance": np.nan,
                    "importance_error": np.nan,
                    "importance_source": "",
                    "status": "not_available_baseline_group_ablation_not_found" if not has_attr else "pending",
                }
            )
    out = pd.DataFrame(rows)
    if global_path.exists():
        imp = read_tsv(global_path)
        imp = imp[imp["protocol"].astype(str).eq(str(args.protocol))].copy()
        for idx, row in out[out["model"].eq("GraphSAGE + ORT_EXP_SUB_ESM2")].iterrows():
            hit = imp[imp["feature_group"].astype(str).eq(str(row["feature_group"]))]
            if not hit.empty:
                out.loc[idx, "grouped_importance"] = float(hit.iloc[0]["delta_auprc_mean"])
                out.loc[idx, "importance_error"] = float(hit.iloc[0]["delta_auprc_std"])
                out.loc[idx, "importance_source"] = str(global_path)
                out.loc[idx, "status"] = "available_group_ablation_delta_auprc"
            else:
                out.loc[idx, "status"] = "not_available_group_not_in_existing_ablation"
    out_path = Path(output_root) / "Figure5_new_C_grouped_feature_importance.tsv"
    write_tsv(out, out_path)

    plot_df = out[out["model"].eq("GraphSAGE + ORT_EXP_SUB_ESM2") & out["grouped_importance"].notna()].copy()
    order = ["ORT", "EXP", "SUB", "Graph topology / degree", "ESM2"]
    plot_df["feature_group"] = pd.Categorical(plot_df["feature_group"], order, ordered=True)
    plot_df = plot_df.sort_values("feature_group")
    fig, ax = plt.subplots(figsize=(5.8, 4.2), facecolor="white")
    ax.bar(
        plot_df["feature_group"].astype(str),
        plot_df["grouped_importance"],
        yerr=plot_df["importance_error"],
        color=["#4C78A8", "#59A14F", "#F28E2B", "#8A8A8A", "#C83E4D"][: len(plot_df)],
        capsize=3,
    )
    ax.set_ylabel("Grouped importance (Delta AUPRC)")
    ax.set_title("ESM2 adds independent grouped feature signal")
    ax.tick_params(axis="x", rotation=30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    pdf_path = Path(output_root) / "Figure5_new_C_grouped_feature_importance.pdf"
    png_path = Path(output_root) / "Figure5_new_C_grouped_feature_importance.png"
    save_pair(fig, pdf_path, png_path)
    return {
        "inputs": [global_path, schema_base, schema_esm2],
        "outputs": [out_path, pdf_path, png_path],
    }


def plot_panel_d(overlap, output_root):
    plot_df = overlap[overlap["top_k"].isin([200, 500])].copy()
    fig, ax = plt.subplots(figsize=(5.4, 4.3), facecolor="white")
    x = np.arange(len(plot_df))
    bottom = np.zeros(len(plot_df))
    segments = [
        ("shared_candidates", "Shared", COLORS["shared"]),
        ("ESM2_unique_candidates", "ESM2 unique", COLORS["esm2_specific"]),
        ("baseline_only_candidates", "Baseline only", COLORS["baseline_only"]),
    ]
    for col, label, color in segments:
        vals = plot_df[col].to_numpy()
        ax.bar(x, vals, bottom=bottom, label=label, color=color)
        bottom += vals
    ax.set_xticks(x, [f"Top{int(k)}" for k in plot_df["top_k"]])
    ax.set_ylabel("Number of candidates")
    ax.set_title("Top-K candidate overlap")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    pdf_path = Path(output_root) / "Figure5_new_D_topk_candidate_overlap.pdf"
    png_path = Path(output_root) / "Figure5_new_D_topk_candidate_overlap.png"
    save_pair(fig, pdf_path, png_path)
    return [pdf_path, png_path]


def plot_panel_e(candidate, rescued, output_root):
    fig, ax = plt.subplots(figsize=(5.4, 5.0), facecolor="white")
    sample = candidate.copy()
    ax.scatter(
        sample["rank_baseline"],
        sample["rank_ESM2"],
        s=8,
        c=np.where(sample["label"].fillna(-1).astype(float).eq(1.0), COLORS["essential"], "#C7C7C7"),
        alpha=0.35,
        linewidths=0,
    )
    if not rescued.empty:
        ax.scatter(
            rescued["rank_baseline"],
            rescued["rank_ESM2"],
            s=34,
            c=COLORS["esm2_specific"],
            edgecolors="white",
            linewidths=0.35,
            label=f"ESM2-rescued Top500 (n={len(rescued)})",
        )
    lim = max(candidate["rank_baseline"].max(), candidate["rank_ESM2"].max())
    ax.plot([1, lim], [1, lim], color="black", lw=1, ls="--", alpha=0.65)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel("Baseline rank (lower is better)")
    ax.set_ylabel("ESM2 rank (lower is better)")
    ax.set_title("ESM2 rank-shift rescue")
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    pdf_path = Path(output_root) / "Figure5_new_E_rank_shift_rescue_plot.pdf"
    png_path = Path(output_root) / "Figure5_new_E_rank_shift_rescue_plot.png"
    save_pair(fig, pdf_path, png_path)
    return [pdf_path, png_path]


def feature_group_profiles(args, candidate):
    bundle = load_protocol_dataset(args.runtime_config, args.protocol, args.esm2_feature)
    schema = bundle["feature_schema"].copy()
    manifest = bundle["node_manifest"].copy()
    manifest["gene_id"] = manifest["graph_gene_id"].astype(str)
    matrix = np.asarray(bundle["feature_matrix"], dtype=float)
    profile = manifest[["gene_id"]].copy()
    name_map = {
        "orthologs": "orthology_support",
        "expression": "expression_signal",
        "sublocalization": "subcellular_signal",
        "degree": "graph_degree_feature",
    }
    for item in schema.itertuples(index=False):
        block = str(item.feature_block)
        if block == "esm2":
            continue
        start = int(item.start_col)
        end = int(item.end_col) + 1
        values = matrix[:, start:end]
        signal = values[:, 0] if values.shape[1] == 1 else np.nanmean(np.abs(values), axis=1)
        profile[name_map.get(block, f"{block}_signal")] = signal
    if "expression_signal" in profile.columns:
        profile["expression_percentile"] = profile["expression_signal"].rank(method="average", pct=True)
    edge_path = Path(args.upstream_root) / args.protocol / args.model / args.esm2_feature / "run_1029" / "edge_table.tsv"
    if edge_path.exists():
        edge = read_tsv(edge_path)
        degree = pd.concat([edge["A"].astype(str), edge["B"].astype(str)]).value_counts().rename_axis("gene_id")
        degree = degree.reset_index(name="ppi_degree")
        profile = profile.merge(degree, on="gene_id", how="left")
        profile["ppi_degree"] = profile["ppi_degree"].fillna(0).astype(float)
    return candidate.merge(profile, on="gene_id", how="left"), schema


def build_biological_profile(args, candidate, rescued, output_root):
    prof, _schema = feature_group_profiles(args, candidate)
    prof["profile_group_top500"] = prof["candidate_class_top500"].replace(
        {"Shared": "Shared_top500", "ESM2_unique": "ESM2_unique_top500", "Baseline_only": "Baseline_unique_top500"}
    )
    prof.loc[prof["candidate_class_top500"].eq("Neither"), "profile_group_top500"] = "Neither"
    prof["is_esm2_rescued_top500"] = prof["gene_id"].isin(set(rescued["gene_id"]))
    profile_path = Path(output_root) / "Figure5_new_F_biological_profile_gene_level.tsv"
    write_tsv(prof, profile_path)

    metrics = [
        ("n_genes", None),
        ("known_essential_fraction", "label"),
        ("mean_ppi_degree", "ppi_degree"),
        ("median_ppi_degree", "ppi_degree"),
        ("mean_orthology_support", "orthology_support"),
        ("mean_expression_percentile", "expression_percentile"),
        ("mean_subcellular_signal", "subcellular_signal"),
        ("mean_rank_gain", "rank_gain"),
        ("mean_delta_score", "delta_score"),
    ]
    groups = ["Shared_top500", "ESM2_unique_top500", "Baseline_unique_top500", "Neither"]
    rows = []
    for group in groups:
        sub = prof[prof["profile_group_top500"] == group]
        for metric, col in metrics:
            if metric == "n_genes":
                value = len(sub)
                status = "available"
            elif metric == "known_essential_fraction":
                labeled = sub[pd.to_numeric(sub[col], errors="coerce").notna()].copy()
                value = float(pd.to_numeric(labeled[col], errors="coerce").eq(1).mean()) if len(labeled) else np.nan
                status = "available" if len(labeled) else "missing_no_labeled_genes"
            elif col in sub.columns:
                vals = pd.to_numeric(sub[col], errors="coerce")
                value = float(vals.median()) if metric.startswith("median") else float(vals.mean())
                status = "available" if vals.notna().any() else f"missing_no_values_for_{col}"
            else:
                value = np.nan
                status = f"missing_column_{col}"
            rows.append({"comparison_group": group, "metric": metric, "value": value, "status": status})
    summary = pd.DataFrame(rows)
    summary_path = Path(output_root) / "Figure5_new_F_biological_profile_summary.tsv"
    write_tsv(summary, summary_path)

    enrichment = pd.DataFrame(
        [
            {
                "analysis": "GO_enrichment",
                "status": "not_completed",
                "reason": "No local F. graminearum GO annotation/background table was found in the inspected project paths.",
            },
            {
                "analysis": "KEGG_enrichment",
                "status": "not_completed",
                "reason": "No local F. graminearum KEGG annotation/background table was found in the inspected project paths.",
            },
        ]
    )
    enrichment_path = Path(output_root) / "Figure5_new_F_enrichment_results.tsv"
    write_tsv(enrichment, enrichment_path)

    heat = summary[summary["status"].eq("available")].pivot(index="metric", columns="comparison_group", values="value")
    heat = heat.drop(index=["n_genes"], errors="ignore")
    z = heat.copy()
    for idx in z.index:
        vals = z.loc[idx].astype(float)
        denom = vals.max() - vals.min()
        z.loc[idx] = 0.0 if denom == 0 or np.isnan(denom) else (vals - vals.min()) / denom
    fig, ax = plt.subplots(figsize=(6.8, 4.7), facecolor="white")
    im = ax.imshow(z.fillna(0).to_numpy(dtype=float), cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(np.arange(z.shape[1]), z.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(z.shape[0]), z.index)
    ax.set_title("Biological profile of ESM2-unique candidates")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            label = "NA" if pd.isna(val) else f"{val:.2g}"
            ax.text(j, i, label, ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, label="Row-normalized value")
    pdf_path = Path(output_root) / "Figure5_new_F_biological_profile_heatmap.pdf"
    png_path = Path(output_root) / "Figure5_new_F_biological_profile_heatmap.png"
    save_pair(fig, pdf_path, png_path)
    return {
        "inputs": [],
        "outputs": [profile_path, summary_path, enrichment_path, pdf_path, png_path],
        "summary_df": summary,
        "enrichment_df": enrichment,
    }


def write_rank_shift_summary(candidate, rescued, output_root, top_k):
    rows = [
        {"metric": "n_genome_wide_genes", "value": len(candidate)},
        {"metric": f"n_esm2_rescued_top{top_k}", "value": len(rescued)},
        {"metric": "max_rank_gain", "value": float(candidate["rank_gain"].max())},
        {"metric": "median_rank_gain", "value": float(candidate["rank_gain"].median())},
        {"metric": "n_positive_rank_gain", "value": int((candidate["rank_gain"] > 0).sum())},
        {"metric": "n_negative_rank_gain", "value": int((candidate["rank_gain"] < 0).sum())},
    ]
    out = pd.DataFrame(rows)
    path = Path(output_root) / "Figure5_new_E_rank_shift_summary.tsv"
    write_tsv(out, path)
    return path


def write_summary(args, seeds, figure5a_result, inputs_by_panel, outputs_by_panel, overlap, candidate, rescued, bio_summary, enrichment, output_root, representative_seed):
    top20 = figure5a_result["esm2_specific_df"].head(20)
    lines = [
        "# Figure5 new design and results summary",
        "",
        'Theme: "ESM2 expands the essentiality prediction manifold and identifies ESM2-specific essential gene candidates."',
        "",
        f"Protocol: `{args.protocol}`; model: `{args.model}`; seeds: {', '.join(map(str, seeds))}.",
        f"Representative Figure5A seed: `{representative_seed}`.",
        "",
        "## Panel inputs and outputs",
    ]
    for panel in ["A", "B", "C", "D", "E", "F"]:
        if panel not in outputs_by_panel:
            continue
        lines.append(f"### Panel 5{panel}")
        lines.append("Inputs:")
        for path in inputs_by_panel.get(panel, []):
            lines.append(f"- `{path}`")
        lines.append("Outputs:")
        for path in outputs_by_panel.get(panel, []):
            lines.append(f"- `{path}`")
        lines.append("")
    lines.extend(["## Figure5A Top-K summary", "", figure5a_result["summary_df"].to_markdown(index=False), ""])
    lines.extend(["## Figure5A top 20 ESM2-specific genes by delta_score", "", top20.to_markdown(index=False), ""])
    lines.extend(["## Top-K overlap", "", overlap.to_markdown(index=False), ""])
    lines.extend(
        [
            "## ESM2-rescued candidates",
            "",
            (
                f"Definition: not in Top{args.default_top_k} baseline, in Top{args.default_top_k} ESM2, "
                "and rank_gain >= 500 or percentile_gain >= 0.2."
            ),
            f"Count: {len(rescued)}",
            "",
        ]
    )
    bio = bio_summary[bio_summary["status"].eq("available")].copy()
    if not bio.empty:
        lines.extend(["## Biological profile", "", bio.to_markdown(index=False), ""])
    incomplete_rows = enrichment[~enrichment["status"].eq("completed")]
    lines.extend(["## Incomplete analyses", ""])
    for row in incomplete_rows.itertuples(index=False):
        lines.append(f"- {row.analysis}: {row.status}; {row.reason}")
    lines.extend(
        [
            "",
            "## Draft manuscript figure legend",
            "",
            (
                "Figure 5. ESM2 expands the essentiality prediction manifold and identifies ESM2-specific "
                "essential gene candidates. (A) GraphSAGE penultimate hidden-embedding UMAP from the "
                "ORT+EXP+SUB+ESM2 model highlights predicted essential genes and separates shared versus "
                "ESM2-specific Top500 essential predictions by final score-based ranking. (B) Matched input- "
                "and hidden-representation UMAPs show how "
                "rescued genes reorganize across representation levels. (C) Grouped feature ablation shows "
                "that ESM2 provides an independent contribution beyond ORT, EXP, and SUB feature groups. "
                "(D) Genome-wide Top-K overlap compares candidates prioritized by baseline GraphSAGE and "
                "GraphSAGE+ESM2 across seeds. (E) Rank-shift analysis identifies ESM2-rescued candidates that "
                "enter Top500 only after ESM2 and show large rank or percentile gains. (F) Biological "
                "profiling summarizes graph, orthology, expression, localization, and label-support "
                "characteristics of shared, ESM2-unique, baseline-unique, and non-top500 genes."
            ),
            "",
        ]
    )
    path = Path(output_root) / "Figure5_new_design_and_results_summary.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main():
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dirs(output_root)
    archive_dir, _copied = archive_legacy_hidden_umap(output_root)

    runtime_seeds = args.seeds or load_runtime_seed_list(args.runtime_config)
    representative_seed, available_seeds = resolve_run_seed(
        args.upstream_root,
        args.protocol,
        args.model,
        args.baseline_feature,
        args.esm2_feature,
        args.preferred_seed,
        runtime_seeds,
    )
    seeds = [seed for seed in runtime_seeds if seed in set(available_seeds)]

    figure5a_result = build_figure5a_outputs(args, output_root, representative_seed)

    inputs_by_panel = {
        "A": [
            figure5a_result["seed_df"]["prediction_path_ORT_EXP_SUB"].iloc[0],
            figure5a_result["seed_df"]["prediction_path_ORT_EXP_SUB_ESM2"].iloc[0],
        ]
    }
    outputs_by_panel = {
        "A": [
            str(figure5a_result["essential_pdf_path"]),
            str(figure5a_result["essential_png_path"]),
            str(figure5a_result["specific_pdf_path"]),
            str(figure5a_result["specific_png_path"]),
            str(figure5a_result["two_panel_pdf_path"]),
            str(figure5a_result["two_panel_png_path"]),
            str(figure5a_result["summary_path"]),
            str(figure5a_result["esm2_specific_path"]),
            str(figure5a_result["shared_path"]),
            str(figure5a_result["audit_path"]),
        ]
    }

    panel_b = plot_panel_b(output_root)
    if panel_b is not None:
        inputs_by_panel["B"] = [str(path) for path in panel_b["inputs"]]
        outputs_by_panel["B"] = [str(path) for path in panel_b["outputs"]]

    panel_c = build_grouped_importance(args, output_root)
    inputs_by_panel["C"] = [str(path) for path in panel_c["inputs"]]
    outputs_by_panel["C"] = [str(path) for path in panel_c["outputs"]]

    candidate, per_seed = build_candidate_table(args, seeds)
    candidate_path = output_root / "Figure5_new_candidate_rank_table.tsv"
    per_seed_path = output_root / "Figure5_new_D_per_seed_rank_percentiles.tsv"
    write_tsv(candidate, candidate_path)
    write_tsv(per_seed, per_seed_path)
    overlap = build_overlap_summary(candidate)
    overlap_path = output_root / "Figure5_new_D_topk_overlap_summary.tsv"
    write_tsv(overlap, overlap_path)
    d_plots = plot_panel_d(overlap, output_root)
    inputs_by_panel["D"] = [str(per_seed_path)]
    outputs_by_panel["D"] = [str(candidate_path), str(overlap_path)] + [str(path) for path in d_plots]

    rescued = build_rescued(candidate, args.default_top_k)
    rescued_path = output_root / "Figure5_new_E_esm2_rescued_candidates.tsv"
    write_tsv(rescued, rescued_path)
    rank_summary_path = write_rank_shift_summary(candidate, rescued, output_root, args.default_top_k)
    e_plots = plot_panel_e(candidate, rescued, output_root)
    inputs_by_panel["E"] = [str(candidate_path)]
    outputs_by_panel["E"] = [str(rescued_path), str(rank_summary_path)] + [str(path) for path in e_plots]

    panel_f = build_biological_profile(args, candidate, rescued, output_root)
    bio_summary = panel_f["summary_df"]
    enrichment = panel_f["enrichment_df"]
    inputs_by_panel["F"] = [str(candidate_path), str(rescued_path)]
    outputs_by_panel["F"] = [str(path) for path in panel_f["outputs"]]

    summary_path = write_summary(
        args,
        seeds,
        figure5a_result,
        inputs_by_panel,
        outputs_by_panel,
        overlap,
        candidate,
        rescued,
        bio_summary,
        enrichment,
        output_root,
        representative_seed,
    )

    print(f"Figure5A essential-vs-nonessential PDF path: {figure5a_result['essential_pdf_path']}")
    print(f"Figure5A ESM2-specific PDF path: {figure5a_result['specific_pdf_path']}")
    print(f"Figure5A two-panel PDF path: {figure5a_result['two_panel_pdf_path']}")
    print(f"Figure5A summary table: {figure5a_result['summary_path']}")
    print("Figure5A first 20 ESM2-specific genes:")
    print(figure5a_result["esm2_specific_df"].head(20).to_csv(sep="\t", index=False).rstrip())
    print("Figure5A audit preview (first 80 lines):")
    audit_lines = figure5a_result["audit_path"].read_text(encoding="utf-8").splitlines()
    print("\n".join(audit_lines[:80]))
    print(f"Legacy hidden-space files archived under: {archive_dir}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
