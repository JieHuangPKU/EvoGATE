import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, matthews_corrcoef
import warnings

try:
    from scipy.stats import mannwhitneyu
except ImportError:  # pragma: no cover - optional dependency fallback
    mannwhitneyu = None

from src.analysis.figure4_representation_common import (
    error_transition,
    graphsage_penultimate_and_logits,
    locate_run_dir,
    species_title,
    subset_indices,
)
from src.analysis.figure5_final_common import (
    FIGURE5_PROTOCOLS,
    FUSARIUM_PROTOCOL,
    STABLE_SEED_THRESHOLD,
    SUBSET,
    confidence_tier_from_count,
    focus_transition_columns,
    load_runtime_seed_list,
    protocol_output_slug,
    save_plot_pair,
    write_manifest,
    write_markdown,
)
from src.data.frozen_protocol_loader import load_protocol_dataset
from src.train.run_frozen_protocol_feature_combo_model import resolve_graph_model_config
from src.train.run_frozen_protocol_model import normalize_model_name


MODEL = "GraphSAGE"
FULL_FEATURE_SETTING = "ORT_EXP_SUB_ESM2"
BASELINE_FEATURE_SETTING = "ORT_EXP_SUB"
FEATURE_GROUP_ORDER = ["ORT", "EXP", "SUB", "ESM2"]
FEATURE_GROUP_LABELS = {
    "ORT": "ORT",
    "EXP": "EXP",
    "SUB": "SUB",
    "ESM2": "ESM2",
}
FEATURE_GROUP_COLORS = {
    "ORT": "#4C78A8",
    "EXP": "#59A14F",
    "SUB": "#F28E2B",
    "ESM2": "#E15759",
}
FEATURE_BLOCK_TO_GROUP = {
    "orthologs": "ORT",
    "expression": "EXP",
    "sublocalization": "SUB",
    "esm2": "ESM2",
}
METRIC_ORDER = ["delta_auprc", "delta_mcc"]
METRIC_LABELS = {
    "delta_auprc": "ΔAUPRC",
    "delta_mcc": "ΔMCC",
}
METRIC_COLORS = {
    "delta_auprc": "#4C78A8",
    "delta_mcc": "#E15759",
}
GENE_SET_ORDER = ["stable_rescued", "always_correct_TP", "persistent_FN", "corrected_FP"]
GENE_SET_LABELS = {
    "stable_rescued": "Stable rescued",
    "always_correct_TP": "Always-correct TP",
    "persistent_FN": "Persistent FN",
    "corrected_FP": "Corrected FP",
}
GENE_SET_DEFINITION_TEXT = {
    "stable_rescued": "essential genes with `stable_rescued_ge2 = True` (Figure5a consensus threshold >= 2 seeds)",
    "always_correct_TP": "essential genes with `TP_stable == n_seeds_observed`",
    "persistent_FN": "essential genes with `FN_persistent == n_seeds_observed`",
    "corrected_FP": "non-essential genes with `FP_to_TN_corrected == n_seeds_observed`",
}
MASKING_RULE = "zero_out_selected_standardized_columns"
MASKING_RULE_DESCRIPTION = (
    "The saved frozen-protocol bundle is already z-scored by the training split; masking sets the selected standardized columns to 0.0."
)
VALUE_DEFINITION = "probability_full - probability_drop_group"

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


def parse_args():
    parser = argparse.ArgumentParser(description="Build Figure5d feature-group attribution outputs")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocols", nargs="+", default=FIGURE5_PROTOCOLS)
    parser.add_argument("--subset", default=SUBSET, type=str)
    parser.add_argument("--model", default=MODEL, type=str)
    parser.add_argument("--feature-setting", default=FULL_FEATURE_SETTING, type=str)
    parser.add_argument("--baseline-feature-setting", default=BASELINE_FEATURE_SETTING, type=str)
    parser.add_argument("--stable-seed-threshold", default=STABLE_SEED_THRESHOLD, type=int)
    parser.add_argument("--data-dir", default="results/Figure5/data", type=str)
    parser.add_argument("--table-dir", default="results/Figure5/tables", type=str)
    parser.add_argument("--plot-dir", default="results/Figure5/plots", type=str)
    parser.add_argument("--summary-dir", default="results/Figure5/summary", type=str)
    return parser.parse_args()


def manifest_row(category, path):
    return {"category": category, "path": str(Path(path).resolve())}


def load_prediction_subset(upstream_root, protocol, feature_setting, seed, subset):
    path = Path(upstream_root) / protocol / MODEL / feature_setting / f"run_{seed}" / "predictions.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df[df["split"].astype(str) == subset].copy().reset_index(drop=True)
    df["graph_gene_id"] = df["graph_gene_id"].astype(str)
    df["label"] = pd.to_numeric(df["label"], errors="raise").astype(int)
    df["pred_label"] = pd.to_numeric(df["pred_label"], errors="raise").astype(int)
    df["pred_score"] = pd.to_numeric(df["pred_score"], errors="raise").astype(float)
    return df


def build_gene_consensus(upstream_root, protocols, seeds, subset, stable_seed_threshold):
    transition_defs = focus_transition_columns(stable_seed_threshold)
    per_seed_rows = []
    for protocol in protocols:
        species_slug = protocol_output_slug(protocol)
        for seed in seeds:
            baseline_df = load_prediction_subset(upstream_root, protocol, BASELINE_FEATURE_SETTING, seed, subset)
            full_df = load_prediction_subset(upstream_root, protocol, FULL_FEATURE_SETTING, seed, subset)
            merged = baseline_df.merge(
                full_df[["graph_gene_id", "pred_score", "pred_label"]],
                on="graph_gene_id",
                how="inner",
                validate="one_to_one",
                suffixes=("_baseline", "_full"),
            )
            merged["transition"] = [
                error_transition(gold, baseline_pred, full_pred)
                for gold, baseline_pred, full_pred in zip(
                    merged["label"].astype(int),
                    merged["pred_label_baseline"].astype(int),
                    merged["pred_label_full"].astype(int),
                )
            ]
            merged["node_id"] = merged["graph_gene_id"].astype(str)
            merged["species"] = species_slug
            merged["protocol"] = protocol
            merged["seed"] = int(seed)
            merged["subset"] = subset
            merged["delta_probability_full_minus_baseline"] = (
                merged["pred_score_full"].astype(float) - merged["pred_score_baseline"].astype(float)
            )
            per_seed_rows.append(
                merged[
                    [
                        "node_id",
                        "graph_gene_id",
                        "species",
                        "protocol",
                        "label",
                        "split",
                        "seed",
                        "subset",
                        "pred_label_baseline",
                        "pred_score_baseline",
                        "pred_label_full",
                        "pred_score_full",
                        "delta_probability_full_minus_baseline",
                        "transition",
                    ]
                ].copy()
            )

    per_seed_df = pd.concat(per_seed_rows, ignore_index=True)
    group_cols = ["node_id", "species", "protocol", "label", "split"]
    metric_df = (
        per_seed_df.groupby(group_cols, dropna=False)
        .agg(
            n_seeds_observed=("seed", "nunique"),
            mean_delta_probability_full_minus_baseline=("delta_probability_full_minus_baseline", "mean"),
            mean_baseline_probability=("pred_score_baseline", "mean"),
            mean_full_probability=("pred_score_full", "mean"),
            mean_baseline_pred_label=("pred_label_baseline", "mean"),
            mean_full_pred_label=("pred_label_full", "mean"),
        )
        .reset_index()
    )
    transition_counts = (
        per_seed_df.groupby(group_cols + ["transition"], dropna=False).size().unstack(fill_value=0).reset_index()
    )
    gene_level_df = metric_df.merge(transition_counts, on=group_cols, how="left", validate="one_to_one")
    for item in transition_defs:
        if item["raw"] not in gene_level_df.columns:
            gene_level_df[item["raw"]] = 0
        gene_level_df[item["count_col"]] = gene_level_df[item["raw"]].fillna(0).astype(int)
        gene_level_df[item["frequency_col"]] = gene_level_df[item["count_col"]] / gene_level_df["n_seeds_observed"].astype(int)
        gene_level_df[item["stable_col"]] = gene_level_df[item["count_col"]] >= stable_seed_threshold
    gene_level_df["rescue_count"] = gene_level_df["fn_to_tp_rescued_count"].astype(int)
    gene_level_df["rescue_frequency"] = gene_level_df["fn_to_tp_rescued_frequency"].astype(float)
    stable_rescue_col = next(item["stable_col"] for item in transition_defs if item["slug"] == "fn_to_tp_rescued")
    gene_level_df["stable_rescued_ge2"] = gene_level_df[stable_rescue_col]
    gene_level_df["confidence_tier"] = gene_level_df["rescue_count"].astype(int).map(
        lambda value: confidence_tier_from_count(value, int(gene_level_df["n_seeds_observed"].max()), stable_seed_threshold)
    )

    def assign_gene_set(row):
        if int(row["label"]) == 1 and bool(row["stable_rescued_ge2"]):
            return "stable_rescued"
        if int(row["label"]) == 1 and int(row.get("TP_stable", 0)) == int(row["n_seeds_observed"]):
            return "always_correct_TP"
        if int(row["label"]) == 1 and int(row.get("FN_persistent", 0)) == int(row["n_seeds_observed"]):
            return "persistent_FN"
        if int(row["label"]) == 0 and int(row.get("FP_to_TN_corrected", 0)) == int(row["n_seeds_observed"]):
            return "corrected_FP"
        return "other"

    gene_level_df["gene_set"] = gene_level_df.apply(assign_gene_set, axis=1)
    gene_level_df["gene_set_definition"] = gene_level_df["gene_set"].map(GENE_SET_DEFINITION_TEXT).fillna("")
    return per_seed_df, gene_level_df


def build_feature_group_manifest(runtime_config_path, protocols):
    rows = []
    for protocol in protocols:
        bundle = load_protocol_dataset(runtime_config_path, protocol, FULL_FEATURE_SETTING)
        species_slug = protocol_output_slug(protocol, bundle["species"])
        schema = bundle["feature_schema"].copy()
        grouped_blocks = set(FEATURE_BLOCK_TO_GROUP)
        for feature_block, feature_group in FEATURE_BLOCK_TO_GROUP.items():
            block_df = schema[schema["feature_block"].astype(str) == feature_block].copy()
            if block_df.empty:
                raise ValueError(f"Required feature block '{feature_block}' was not found for protocol '{protocol}'")
            rows.append(
                {
                    "protocol": protocol,
                    "species": species_slug,
                    "feature_group": feature_group,
                    "source_feature_block": feature_block,
                    "start_col": int(block_df["start_col"].min()),
                    "end_col": int(block_df["end_col"].max()),
                    "n_features": int(block_df["dimension"].sum()),
                    "masking_rule_used": MASKING_RULE,
                    "attribution_status": "attributed",
                }
            )
        for remaining in schema[~schema["feature_block"].astype(str).isin(grouped_blocks)].itertuples(index=False):
            rows.append(
                {
                    "protocol": protocol,
                    "species": species_slug,
                    "feature_group": f"{str(remaining.feature_block).upper()}_UNATTRIBUTED",
                    "source_feature_block": str(remaining.feature_block),
                    "start_col": int(remaining.start_col),
                    "end_col": int(remaining.end_col),
                    "n_features": int(remaining.dimension),
                    "masking_rule_used": "held_constant_not_masked",
                    "attribution_status": "held_constant",
                }
            )
    manifest_df = pd.DataFrame(rows)
    protocol_order = {protocol: idx for idx, protocol in enumerate(protocols)}
    manifest_df["protocol_rank"] = manifest_df["protocol"].map(protocol_order)
    manifest_df = manifest_df.sort_values(
        ["protocol_rank", "start_col", "end_col", "feature_group"], kind="stable"
    ).drop(columns="protocol_rank")
    return manifest_df.reset_index(drop=True)


def feature_group_columns(feature_schema):
    group_to_columns = {}
    for feature_block, feature_group in FEATURE_BLOCK_TO_GROUP.items():
        block_df = feature_schema[feature_schema["feature_block"].astype(str) == feature_block].copy()
        if block_df.empty:
            raise ValueError(f"Feature block '{feature_block}' missing from feature_schema")
        columns = []
        for row in block_df.itertuples(index=False):
            columns.extend(range(int(row.start_col), int(row.end_col) + 1))
        group_to_columns[feature_group] = columns
    return group_to_columns


def predict_probabilities(bundle, checkpoint_path, model_cfg, masked_columns=None):
    aggregator_type = str(model_cfg["aggregator_type"]).strip().lower()
    if aggregator_type != "mean":
        raise ValueError(
            "Figure5d currently supports the frozen-protocol mean GraphSAGE inference path only; "
            f"got aggregator_type='{aggregator_type}'"
        )
    checkpoint = torch.load(str(checkpoint_path), map_location="cpu")
    state_dict = checkpoint["state_dict"]
    x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
    if masked_columns:
        x[:, masked_columns] = 0.0
    edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
    with torch.no_grad():
        _, logits = graphsage_penultimate_and_logits(x, edge_index, state_dict, aggregator_type)
        probabilities = torch.sigmoid(logits.view(-1)).cpu().numpy().astype(np.float32, copy=False)
    return probabilities


def summarize_metrics(labels, probabilities):
    labels = np.asarray(labels, dtype=int)
    probabilities = np.asarray(probabilities, dtype=float)
    pred_labels = (probabilities >= 0.5).astype(int)
    return {
        "auprc": float(average_precision_score(labels, probabilities)),
        "mcc": float(matthews_corrcoef(labels, pred_labels)),
        "pred_labels": pred_labels,
    }


def benjamini_hochberg(values):
    values = np.asarray(values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    mask = np.isfinite(values)
    if not mask.any():
        return adjusted
    valid = values[mask]
    order = np.argsort(valid)
    ranked = valid[order]
    n = float(len(ranked))
    scaled = ranked * n / np.arange(1.0, n + 1.0)
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    scaled = np.clip(scaled, 0.0, 1.0)
    adjusted_valid = np.empty_like(valid)
    adjusted_valid[order] = scaled
    adjusted[mask] = adjusted_valid
    return adjusted


def build_global_importance_plot(plot_df, protocols, plot_pdf, plot_png):
    protocol_order = [protocol for protocol in protocols if protocol in set(plot_df["protocol"])]
    fig, axes = plt.subplots(1, len(protocol_order), figsize=(5.2 * len(protocol_order), 4.2), facecolor="white", sharey=True)
    if len(protocol_order) == 1:
        axes = [axes]
    bar_width = 0.34
    metric_offsets = {
        "delta_auprc": -bar_width / 2.0,
        "delta_mcc": bar_width / 2.0,
    }
    for ax, protocol in zip(axes, protocol_order):
        species_df = plot_df[plot_df["protocol"] == protocol].copy()
        x_positions = np.arange(len(FEATURE_GROUP_ORDER), dtype=float)
        for metric_name in METRIC_ORDER:
            metric_df = (
                species_df[species_df["metric"] == metric_name]
                .set_index("feature_group")
                .reindex(FEATURE_GROUP_ORDER)
                .reset_index()
            )
            ax.bar(
                x_positions + metric_offsets[metric_name],
                metric_df["plot_value"].astype(float),
                width=bar_width,
                color=METRIC_COLORS[metric_name],
                edgecolor="white",
                linewidth=0.6,
                yerr=metric_df["plot_error"].astype(float),
                capsize=2.5,
                label=METRIC_LABELS[metric_name],
            )
        ax.axhline(0.0, color="#6B6B6B", linewidth=0.8)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([FEATURE_GROUP_LABELS[name] for name in FEATURE_GROUP_ORDER])
        ax.set_title(species_title(protocol), fontsize=11)
        ax.set_xlabel("Masked feature group")
        ax.set_facecolor("white")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Global performance drop")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, frameon=False, loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 1.03))
    fig.suptitle("Figure5d-1. Feature-group global importance", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_plot_pair(fig, plot_pdf, plot_png)


def build_rescued_vs_other_plot(plot_df, protocols, plot_pdf, plot_png):
    protocol_order = [protocol for protocol in protocols if protocol in set(plot_df["protocol"])]
    gene_sets = [name for name in GENE_SET_ORDER if name in set(plot_df["gene_set"])]
    fig, axes = plt.subplots(1, len(protocol_order), figsize=(6.2 * len(protocol_order), 4.6), facecolor="white", sharey=True)
    if len(protocol_order) == 1:
        axes = [axes]
    offsets = np.linspace(-0.3, 0.3, num=len(FEATURE_GROUP_ORDER))
    width = 0.14
    for ax, protocol in zip(axes, protocol_order):
        species_df = plot_df[plot_df["protocol"] == protocol].copy()
        centers = np.arange(len(gene_sets), dtype=float)
        legend_handles = []
        for feature_index, feature_group in enumerate(FEATURE_GROUP_ORDER):
            data = []
            for gene_set in gene_sets:
                values = species_df[
                    (species_df["gene_set"] == gene_set) & (species_df["feature_group"] == feature_group)
                ]["plot_value"].astype(float).to_numpy()
                data.append(values if len(values) else np.array([np.nan], dtype=float))
            artist = ax.boxplot(
                data,
                positions=centers + offsets[feature_index],
                widths=width,
                patch_artist=True,
                showfliers=False,
                manage_ticks=False,
                medianprops={"color": "#1F1F1F", "linewidth": 1.0},
                whiskerprops={"color": "#4A4A4A", "linewidth": 0.8},
                capprops={"color": "#4A4A4A", "linewidth": 0.8},
            )
            for patch in artist["boxes"]:
                patch.set_facecolor(FEATURE_GROUP_COLORS[feature_group])
                patch.set_edgecolor("white")
                patch.set_linewidth(0.6)
                patch.set_alpha(0.95)
            legend_handles.append(artist["boxes"][0])
        ax.axhline(0.0, color="#6B6B6B", linewidth=0.8)
        ax.set_xticks(centers)
        ax.set_xticklabels([GENE_SET_LABELS[name] for name in gene_sets], rotation=15, ha="right")
        ax.set_title(species_title(protocol), fontsize=11)
        ax.set_xlabel("Gene subset")
        ax.set_facecolor("white")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(legend_handles, [FEATURE_GROUP_LABELS[name] for name in FEATURE_GROUP_ORDER], frameon=False, loc="best")
    axes[0].set_ylabel("Mean per-gene contribution\n(full probability - masked probability)")
    fig.suptitle("Figure5d-2. Feature-group contribution by gene subset", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_plot_pair(fig, plot_pdf, plot_png)


def build_gene_set_heatmap_plot(plot_df, protocols, plot_pdf, plot_png):
    protocol_order = [protocol for protocol in protocols if protocol in set(plot_df["protocol"])]
    gene_sets = [name for name in GENE_SET_ORDER if name in set(plot_df["gene_set"])]
    matrices = []
    for protocol in protocol_order:
        pivot_df = (
            plot_df[plot_df["protocol"] == protocol]
            .pivot(index="gene_set", columns="feature_group", values="plot_value")
            .reindex(index=gene_sets, columns=FEATURE_GROUP_ORDER)
        )
        matrices.append(pivot_df.to_numpy(dtype=float))
    max_abs = float(np.nanmax(np.abs(np.concatenate(matrices, axis=None)))) if matrices else 1.0
    if not np.isfinite(max_abs) or max_abs <= 0.0:
        max_abs = 1.0
    fig, axes = plt.subplots(1, len(protocol_order), figsize=(4.2 * len(protocol_order), 4.6), facecolor="white", sharey=True)
    if len(protocol_order) == 1:
        axes = [axes]
    im = None
    for ax, protocol in zip(axes, protocol_order):
        pivot_df = (
            plot_df[plot_df["protocol"] == protocol]
            .pivot(index="gene_set", columns="feature_group", values="plot_value")
            .reindex(index=gene_sets, columns=FEATURE_GROUP_ORDER)
        )
        matrix = pivot_df.to_numpy(dtype=float)
        im = ax.imshow(matrix, cmap="RdBu_r", vmin=-max_abs, vmax=max_abs, aspect="auto")
        ax.set_title(species_title(protocol), fontsize=11)
        ax.set_xticks(np.arange(len(FEATURE_GROUP_ORDER)))
        ax.set_xticklabels([FEATURE_GROUP_LABELS[name] for name in FEATURE_GROUP_ORDER])
        ax.set_yticks(np.arange(len(gene_sets)))
        ax.set_yticklabels([GENE_SET_LABELS[name] for name in gene_sets])
        ax.set_facecolor("white")
        for row_index, gene_set in enumerate(gene_sets):
            for col_index, feature_group in enumerate(FEATURE_GROUP_ORDER):
                value = matrix[row_index, col_index]
                if np.isnan(value):
                    continue
                ax.text(
                    col_index,
                    row_index,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if abs(value) > (max_abs * 0.55) else "#1F1F1F",
                )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    if im is not None:
        colorbar = fig.colorbar(im, ax=axes, fraction=0.046, pad=0.04)
        colorbar.set_label("Mean contribution")
    fig.suptitle("Figure5d-3. Feature-group composition heatmap", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_plot_pair(fig, plot_pdf, plot_png)


def main():
    args = parse_args()
    data_dir = Path(args.data_dir)
    table_dir = Path(args.table_dir)
    plot_dir = Path(args.plot_dir)
    summary_dir = Path(args.summary_dir)
    for path in [data_dir, table_dir, plot_dir, summary_dir]:
        path.mkdir(parents=True, exist_ok=True)

    seeds = load_runtime_seed_list(args.runtime_config)
    per_seed_consensus_df, gene_level_df = build_gene_consensus(
        upstream_root=args.upstream_root,
        protocols=args.protocols,
        seeds=seeds,
        subset=args.subset,
        stable_seed_threshold=args.stable_seed_threshold,
    )

    manifest_rows = []
    feature_manifest_df = build_feature_group_manifest(args.runtime_config, args.protocols)
    feature_manifest_path = data_dir / "Figure5d_feature_group_manifest.tsv"
    feature_manifest_df.to_csv(feature_manifest_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", feature_manifest_path))

    import yaml

    base_config = yaml.safe_load(Path(args.runtime_config).read_text(encoding="utf-8"))
    model_cfg = resolve_graph_model_config(base_config, normalize_model_name(args.model))

    global_metric_rows = []
    per_gene_rows = []
    validation_rows = []

    for protocol in args.protocols:
        species_slug = protocol_output_slug(protocol)
        for seed in seeds:
            bundle = load_protocol_dataset(args.runtime_config, protocol, args.feature_setting)
            checkpoint_path = locate_run_dir(args.upstream_root, protocol, args.model, args.feature_setting, seed) / "best_model.pt"
            selected_idx = subset_indices(bundle, args.subset)
            subset_manifest = bundle["node_manifest"].iloc[selected_idx].copy().reset_index(drop=True)
            subset_manifest["graph_gene_id"] = subset_manifest["graph_gene_id"].astype(str)
            subset_manifest["label"] = pd.to_numeric(subset_manifest["label"], errors="raise").astype(int)
            subset_manifest["split"] = subset_manifest["split"].astype(str)

            full_probabilities = predict_probabilities(bundle, checkpoint_path, model_cfg)
            full_subset_probabilities = full_probabilities[selected_idx]
            full_metrics = summarize_metrics(subset_manifest["label"].to_numpy(dtype=int), full_subset_probabilities)

            saved_full_df = load_prediction_subset(args.upstream_root, protocol, args.feature_setting, seed, args.subset)
            validation_df = subset_manifest.merge(
                saved_full_df[["graph_gene_id", "pred_score", "pred_label"]],
                on="graph_gene_id",
                how="left",
                validate="one_to_one",
            )
            if validation_df["pred_score"].isna().any():
                raise ValueError(f"Missing saved predictions after merge for protocol={protocol} seed={seed}")
            validation_rows.append(
                {
                    "protocol": protocol,
                    "species": species_slug,
                    "seed": int(seed),
                    "subset": args.subset,
                    "max_abs_probability_diff_vs_saved": float(
                        np.max(np.abs(full_subset_probabilities - validation_df["pred_score"].astype(float).to_numpy()))
                    ),
                    "mean_abs_probability_diff_vs_saved": float(
                        np.mean(np.abs(full_subset_probabilities - validation_df["pred_score"].astype(float).to_numpy()))
                    ),
                    "saved_pred_label_match_fraction": float(
                        np.mean(full_metrics["pred_labels"] == validation_df["pred_label"].astype(int).to_numpy())
                    ),
                }
            )

            group_columns = feature_group_columns(bundle["feature_schema"])
            masked_probabilities = {}
            for feature_group in FEATURE_GROUP_ORDER:
                masked = predict_probabilities(bundle, checkpoint_path, model_cfg, masked_columns=group_columns[feature_group])
                masked_probabilities[feature_group] = masked
                masked_subset = masked[selected_idx]
                masked_metrics = summarize_metrics(subset_manifest["label"].to_numpy(dtype=int), masked_subset)
                global_metric_rows.append(
                    {
                        "protocol": protocol,
                        "species": species_slug,
                        "seed": int(seed),
                        "subset": args.subset,
                        "feature_group": feature_group,
                        "full_auprc": full_metrics["auprc"],
                        "drop_auprc": masked_metrics["auprc"],
                        "delta_auprc": float(full_metrics["auprc"] - masked_metrics["auprc"]),
                        "full_mcc": full_metrics["mcc"],
                        "drop_mcc": masked_metrics["mcc"],
                        "delta_mcc": float(full_metrics["mcc"] - masked_metrics["mcc"]),
                        "n_test_genes": int(len(selected_idx)),
                        "n_positive_test_genes": int(subset_manifest["label"].sum()),
                    }
                )

            consensus_lookup = gene_level_df[
                (gene_level_df["protocol"] == protocol) & (gene_level_df["split"].astype(str) == args.subset)
            ][["node_id", "gene_set", "gene_set_definition"]].drop_duplicates(subset=["node_id"], keep="first")
            seed_df = subset_manifest.rename(columns={"graph_gene_id": "gene_id", "label": "true_label"}).copy()
            seed_df["protocol"] = protocol
            seed_df["species"] = species_slug
            seed_df["seed"] = int(seed)
            seed_df["predicted_label_full"] = full_metrics["pred_labels"].astype(int)
            seed_df["probability_full"] = full_subset_probabilities.astype(float)
            for feature_group in FEATURE_GROUP_ORDER:
                drop_subset = masked_probabilities[feature_group][selected_idx].astype(float)
                seed_df[f"probability_drop_{feature_group}"] = drop_subset
                seed_df[f"contribution_{feature_group}"] = seed_df["probability_full"].astype(float) - drop_subset
            seed_df = seed_df.merge(
                consensus_lookup.rename(columns={"node_id": "gene_id"}),
                on="gene_id",
                how="left",
                validate="one_to_one",
            )
            seed_df["gene_set"] = seed_df["gene_set"].fillna("other")
            seed_df["gene_set_definition"] = seed_df["gene_set_definition"].fillna("")
            per_gene_rows.append(seed_df)

    per_gene_df = pd.concat(per_gene_rows, ignore_index=True)
    per_gene_path = data_dir / "Figure5d_group_ablation_per_gene.tsv"
    per_gene_df.to_csv(per_gene_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", per_gene_path))

    validation_df = pd.DataFrame(validation_rows).sort_values(["protocol", "seed"]).reset_index(drop=True)
    validation_path = table_dir / "Figure5d_inference_reproducibility_check.tsv"
    validation_df.to_csv(validation_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", validation_path))

    global_metrics_df = pd.DataFrame(global_metric_rows).sort_values(["protocol", "seed", "feature_group"]).reset_index(drop=True)
    global_metrics_path = table_dir / "Figure5d_group_ablation_global_metrics.tsv"
    global_metrics_df.to_csv(global_metrics_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", global_metrics_path))

    global_summary_df = (
        global_metrics_df.groupby(["protocol", "species", "feature_group"], dropna=False)
        .agg(
            n_seeds=("seed", "nunique"),
            full_auprc_mean=("full_auprc", "mean"),
            full_auprc_std=("full_auprc", "std"),
            drop_auprc_mean=("drop_auprc", "mean"),
            drop_auprc_std=("drop_auprc", "std"),
            delta_auprc_mean=("delta_auprc", "mean"),
            delta_auprc_std=("delta_auprc", "std"),
            full_mcc_mean=("full_mcc", "mean"),
            full_mcc_std=("full_mcc", "std"),
            drop_mcc_mean=("drop_mcc", "mean"),
            drop_mcc_std=("drop_mcc", "std"),
            delta_mcc_mean=("delta_mcc", "mean"),
            delta_mcc_std=("delta_mcc", "std"),
        )
        .reset_index()
    )
    for column in [name for name in global_summary_df.columns if name.endswith("_std")]:
        global_summary_df[column] = global_summary_df[column].fillna(0.0)
    global_summary_path = table_dir / "Figure5d_group_ablation_global_summary.tsv"
    global_summary_df.to_csv(global_summary_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", global_summary_path))

    global_plot_rows = []
    for row in global_summary_df.itertuples(index=False):
        for metric_name in METRIC_ORDER:
            base_metric = metric_name.replace("delta_", "")
            global_plot_rows.append(
                {
                    "protocol": row.protocol,
                    "species": row.species,
                    "feature_group": row.feature_group,
                    "metric": metric_name,
                    "metric_label": METRIC_LABELS[metric_name],
                    "plot_value": float(getattr(row, f"{metric_name}_mean")),
                    "plot_error": float(getattr(row, f"{metric_name}_std")),
                    "full_metric_mean": float(getattr(row, f"full_{base_metric}_mean")),
                    "drop_metric_mean": float(getattr(row, f"drop_{base_metric}_mean")),
                    "n_seeds": int(row.n_seeds),
                    "value_transform": "none",
                    "plot_value_definition": "full_metric_mean - masked_metric_mean",
                }
            )
    global_plot_df = pd.DataFrame(global_plot_rows)
    global_plot_path = data_dir / "Figure5d_feature_group_global_importance_plot_data.tsv"
    global_plot_df.to_csv(global_plot_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", global_plot_path))

    per_gene_mean_df = (
        per_gene_df.groupby(["gene_id", "protocol", "species", "true_label", "gene_set", "gene_set_definition"], dropna=False)
        .agg(
            n_seed_observations=("seed", "nunique"),
            mean_probability_full=("probability_full", "mean"),
            mean_probability_drop_ORT=("probability_drop_ORT", "mean"),
            mean_probability_drop_EXP=("probability_drop_EXP", "mean"),
            mean_probability_drop_SUB=("probability_drop_SUB", "mean"),
            mean_probability_drop_ESM2=("probability_drop_ESM2", "mean"),
            mean_contribution_ORT=("contribution_ORT", "mean"),
            mean_contribution_EXP=("contribution_EXP", "mean"),
            mean_contribution_SUB=("contribution_SUB", "mean"),
            mean_contribution_ESM2=("contribution_ESM2", "mean"),
        )
        .reset_index()
    )
    per_gene_mean_df = per_gene_mean_df[per_gene_mean_df["gene_set"].isin(GENE_SET_ORDER)].copy()

    rescued_plot_rows = []
    for row in per_gene_mean_df.itertuples(index=False):
        for feature_group in FEATURE_GROUP_ORDER:
            rescued_plot_rows.append(
                {
                    "gene_id": row.gene_id,
                    "protocol": row.protocol,
                    "species": row.species,
                    "true_label": int(row.true_label),
                    "gene_set": row.gene_set,
                    "gene_set_definition": row.gene_set_definition,
                    "feature_group": feature_group,
                    "plot_value": float(getattr(row, f"mean_contribution_{feature_group}")),
                    "raw_mean_contribution": float(getattr(row, f"mean_contribution_{feature_group}")),
                    "n_seed_observations": int(row.n_seed_observations),
                    "seed_aggregation": "mean_across_seeds",
                    "value_transform": "none",
                    "plot_value_definition": VALUE_DEFINITION,
                }
            )
    rescued_plot_df = pd.DataFrame(rescued_plot_rows)
    rescued_plot_path = data_dir / "Figure5d_feature_group_rescued_vs_other_plot_data.tsv"
    rescued_plot_df.to_csv(rescued_plot_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", rescued_plot_path))

    by_gene_set_df = (
        rescued_plot_df.groupby(["protocol", "species", "gene_set", "gene_set_definition", "feature_group"], dropna=False)
        .agg(
            n_genes=("gene_id", "nunique"),
            mean_contribution=("plot_value", "mean"),
            median_contribution=("plot_value", "median"),
            std_contribution=("plot_value", "std"),
            min_contribution=("plot_value", "min"),
            max_contribution=("plot_value", "max"),
        )
        .reset_index()
    )
    by_gene_set_df["std_contribution"] = by_gene_set_df["std_contribution"].fillna(0.0)
    by_gene_set_path = table_dir / "Figure5d_group_ablation_by_gene_set.tsv"
    by_gene_set_df.to_csv(by_gene_set_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", by_gene_set_path))

    stats_rows = []
    for protocol in args.protocols:
        species_slug = protocol_output_slug(protocol)
        protocol_df = rescued_plot_df[rescued_plot_df["protocol"] == protocol].copy()
        available_gene_sets = [name for name in GENE_SET_ORDER if name in set(protocol_df["gene_set"])]
        if "stable_rescued" not in available_gene_sets:
            continue
        for feature_group in FEATURE_GROUP_ORDER:
            rescue_values = protocol_df[
                (protocol_df["gene_set"] == "stable_rescued") & (protocol_df["feature_group"] == feature_group)
            ]["plot_value"].astype(float).to_numpy()
            for comparator in [name for name in available_gene_sets if name != "stable_rescued"]:
                comparator_values = protocol_df[
                    (protocol_df["gene_set"] == comparator) & (protocol_df["feature_group"] == feature_group)
                ]["plot_value"].astype(float).to_numpy()
                p_value = float("nan")
                statistic = float("nan")
                if len(rescue_values) > 0 and len(comparator_values) > 0 and mannwhitneyu is not None:
                    statistic, p_value = mannwhitneyu(rescue_values, comparator_values, alternative="two-sided")
                stats_rows.append(
                    {
                        "protocol": protocol,
                        "species": species_slug,
                        "feature_group": feature_group,
                        "reference_gene_set": "stable_rescued",
                        "comparison_gene_set": comparator,
                        "reference_n_genes": int(len(rescue_values)),
                        "comparison_n_genes": int(len(comparator_values)),
                        "reference_mean": float(np.mean(rescue_values)) if len(rescue_values) else float("nan"),
                        "comparison_mean": float(np.mean(comparator_values)) if len(comparator_values) else float("nan"),
                        "reference_median": float(np.median(rescue_values)) if len(rescue_values) else float("nan"),
                        "comparison_median": float(np.median(comparator_values)) if len(comparator_values) else float("nan"),
                        "mean_difference_reference_minus_comparison": (
                            float(np.mean(rescue_values) - np.mean(comparator_values))
                            if len(rescue_values) and len(comparator_values)
                            else float("nan")
                        ),
                        "median_difference_reference_minus_comparison": (
                            float(np.median(rescue_values) - np.median(comparator_values))
                            if len(rescue_values) and len(comparator_values)
                            else float("nan")
                        ),
                        "mannwhitney_u": statistic,
                        "mannwhitney_pvalue": p_value,
                        "test_available": str(mannwhitneyu is not None).lower(),
                    }
                )
    stats_df = pd.DataFrame(stats_rows)
    if not stats_df.empty:
        stats_df["mannwhitney_fdr_bh"] = benjamini_hochberg(stats_df["mannwhitney_pvalue"].to_numpy(dtype=float))
    else:
        stats_df["mannwhitney_fdr_bh"] = []
    stats_path = table_dir / "Figure5d_group_ablation_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", stats_path))

    heatmap_df = by_gene_set_df[
        ["protocol", "species", "gene_set", "gene_set_definition", "feature_group", "n_genes", "mean_contribution"]
    ].copy()
    heatmap_df = heatmap_df.rename(columns={"mean_contribution": "plot_value"})
    heatmap_df["value_transform"] = "none"
    heatmap_df["plot_value_definition"] = "mean_across_seeds(probability_full - probability_drop_group)"
    heatmap_plot_path = data_dir / "Figure5d_feature_group_gene_set_heatmap_plot_data.tsv"
    heatmap_df.to_csv(heatmap_plot_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", heatmap_plot_path))

    heatmap_table_path = table_dir / "Figure5d_group_ablation_gene_set_heatmap_values.tsv"
    heatmap_df.to_csv(heatmap_table_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", heatmap_table_path))

    global_pdf = plot_dir / "Figure5d_feature_group_global_importance.pdf"
    global_png = plot_dir / "Figure5d_feature_group_global_importance.png"
    build_global_importance_plot(global_plot_df, args.protocols, global_pdf, global_png)
    manifest_rows.append(manifest_row("plot", global_pdf))
    manifest_rows.append(manifest_row("plot", global_png))

    rescued_pdf = plot_dir / "Figure5d_feature_group_rescued_vs_other.pdf"
    rescued_png = plot_dir / "Figure5d_feature_group_rescued_vs_other.png"
    build_rescued_vs_other_plot(rescued_plot_df, args.protocols, rescued_pdf, rescued_png)
    manifest_rows.append(manifest_row("plot", rescued_pdf))
    manifest_rows.append(manifest_row("plot", rescued_png))

    heatmap_pdf = plot_dir / "Figure5d_feature_group_gene_set_heatmap.pdf"
    heatmap_png = plot_dir / "Figure5d_feature_group_gene_set_heatmap.png"
    build_gene_set_heatmap_plot(heatmap_df, args.protocols, heatmap_pdf, heatmap_png)
    manifest_rows.append(manifest_row("plot", heatmap_pdf))
    manifest_rows.append(manifest_row("plot", heatmap_png))

    fusarium_global = global_summary_df[global_summary_df["protocol"] == FUSARIUM_PROTOCOL].copy()
    fusarium_best_auprc = (
        fusarium_global.sort_values(["delta_auprc_mean", "delta_mcc_mean"], ascending=[False, False]).iloc[0]
        if not fusarium_global.empty
        else None
    )
    rescued_vs_tp = by_gene_set_df[
        (by_gene_set_df["protocol"] == FUSARIUM_PROTOCOL)
        & (by_gene_set_df["feature_group"] == "ESM2")
        & (by_gene_set_df["gene_set"].isin(["stable_rescued", "always_correct_TP"]))
    ].copy()
    rescued_vs_tp_map = dict(zip(rescued_vs_tp["gene_set"], rescued_vs_tp["mean_contribution"]))
    heatmap_fusarium = heatmap_df[heatmap_df["protocol"] == FUSARIUM_PROTOCOL].copy()
    top_by_gene_set = {}
    for gene_set, subset_df in heatmap_fusarium.groupby("gene_set", dropna=False):
        subset_df = subset_df.sort_values(["plot_value", "feature_group"], ascending=[False, True]).reset_index(drop=True)
        if not subset_df.empty:
            top_by_gene_set[gene_set] = (subset_df.loc[0, "feature_group"], float(subset_df.loc[0, "plot_value"]))

    summary_lines = [
        "# Figure5d Feature-Group Attribution",
        "",
        "## Method",
        "",
        "- Primary attribution method: fixed-model leave-one-group-out inference ablation on the saved `ORT_EXP_SUB_ESM2` GraphSAGE checkpoint.",
        f"- Masking rule: `{MASKING_RULE}`.",
        f"- Masking rule detail: {MASKING_RULE_DESCRIPTION}",
        "- Attribution score definition: `contribution(group) = probability_full - probability_drop_group`.",
        "- Stored prediction quantity: probabilities from `sigmoid(logit)`.",
        "- No model retraining was performed for the ablated groups.",
        "- The additional graph-degree scalar present in the frozen bundle was detected and held constant; Figure5d attributes only the four requested groups `ORT / EXP / SUB / ESM2`.",
        "",
        "## Gene subsets used",
        "",
        *[f"- `{name}`: {GENE_SET_DEFINITION_TEXT[name]}." for name in GENE_SET_ORDER if name in set(gene_level_df["gene_set"])],
        "",
        "## Main findings",
        "",
    ]
    if fusarium_best_auprc is not None:
        summary_lines.append(
            "- Fusarium global importance: the largest mean `ΔAUPRC` drop was observed for "
            f"`{fusarium_best_auprc['feature_group']}` (`{fusarium_best_auprc['delta_auprc_mean']:.4f}`), "
            f"with mean `ΔMCC = {fusarium_best_auprc['delta_mcc_mean']:.4f}`."
        )
    if {"stable_rescued", "always_correct_TP"}.issubset(rescued_vs_tp_map):
        direction = "higher" if rescued_vs_tp_map["stable_rescued"] > rescued_vs_tp_map["always_correct_TP"] else "not higher"
        summary_lines.append(
            "- Fusarium ESM2 dependence: stable rescued genes had "
            f"`{direction}` mean ESM2 contribution than always-correct TPs "
            f"(`stable_rescued = {rescued_vs_tp_map['stable_rescued']:.4f}`, "
            f"`always_correct_TP = {rescued_vs_tp_map['always_correct_TP']:.4f}`)."
        )
    for gene_set in GENE_SET_ORDER:
        if gene_set not in top_by_gene_set:
            continue
        top_group, top_value = top_by_gene_set[gene_set]
        summary_lines.append(
            f"- Fusarium `{gene_set}` heatmap pattern: highest mean contribution came from `{top_group}` (`{top_value:.4f}`)."
        )
    summary_lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- Figure5d is group-level attribution only; it does not identify individual raw dimensions, residues, or causal sequence motifs.",
            "- Zeroing standardized columns is an inference-time perturbation on the frozen normalized feature matrix, not a retrained drop-group model.",
            "- Contribution values can be negative when masking a group increases predicted probability or improves a metric for a specific gene/run.",
            "",
            "## Outputs",
            "",
            f"- Feature-group manifest: `{feature_manifest_path}`.",
            f"- Per-gene attribution table: `{per_gene_path}`.",
            f"- Global summary table: `{global_summary_path}`.",
            f"- Gene-set summary table: `{by_gene_set_path}`.",
            f"- Statistical comparison table: `{stats_path}`.",
            f"- Global importance plot: `{global_pdf}` / `{global_png}`.",
            f"- Rescued-vs-other plot: `{rescued_pdf}` / `{rescued_png}`.",
            f"- Gene-set heatmap: `{heatmap_pdf}` / `{heatmap_png}`.",
        ]
    )
    summary_path = summary_dir / "Figure5d_feature_group_attribution.md"
    write_markdown(summary_path, summary_lines)
    manifest_rows.append(manifest_row("summary", summary_path))

    manifest_path = table_dir / "Figure5d_output_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)


if __name__ == "__main__":
    main()
