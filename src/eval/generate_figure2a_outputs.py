from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize

from src.eval.aggregate_frozen_protocol_runs import aggregate_runs, collect_metrics, dataframe_to_markdown
from src.eval.publication_summary import build_publication_summary, publication_markdown


PROTOCOL_TO_REGIME = {
    "fgraminearum_oldlabel": "old440",
    "fgraminearum_newlabel": "new_label",
}
REGIME_ORDER = ["old440", "new_label"]
DEFAULT_MODELS = ["GraphSAGE", "GCN", "GAT", "GIN"]
DEFAULT_FEATURE_SETTINGS = ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB"]
FEATURE_LABELS = {
    "ORT": "ortholog",
    "EXP": "expression",
    "SUB": "subloc",
    "ORT_EXP": "ortholog+expression",
    "ORT_SUB": "ortholog+subloc",
    "EXP_SUB": "expression+subloc",
    "ORT_EXP_SUB": "ortholog+expression+subloc",
}
PUBLIC_METRIC_COLUMNS = {
    "test_auroc": "AUROC",
    "test_auprc": "AUPRC",
    "test_mcc": "MCC",
    "test_f1": "F1",
}
BOOL_TRUE = {"true", "t", "1", "yes", "y"}
OLD_NEW_PALETTE = {"old440": "#768B93", "new_label": "#BC240F"}
MODEL_PALETTE = {
    "GraphSAGE": "#855C75",
    "GCN": "#C97C5D",
    "GAT": "#ADC3C9",
    "GIN": "#7D9D90",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Figure2a summaries and plots.")
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    parser.add_argument("--prefix", default="Figure2a", type=str)
    parser.add_argument("--target-name", default="fgraminearum", type=str)
    parser.add_argument("--models", default=",".join(DEFAULT_MODELS), type=str)
    parser.add_argument("--feature-settings", default=",".join(DEFAULT_FEATURE_SETTINGS), type=str)
    parser.add_argument("--seeds", default="", type=str)
    return parser.parse_args()


def split_csv(raw_value: str) -> list[str]:
    return [item.strip() for item in raw_value.split(",") if item.strip()]


def configure_plotting() -> None:
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["font.size"] = 9
    plt.rcParams["axes.titlesize"] = 10
    plt.rcParams["axes.labelsize"] = 9
    plt.rcParams["xtick.labelsize"] = 8
    plt.rcParams["ytick.labelsize"] = 8
    plt.rcParams["legend.fontsize"] = 8
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1.0
    plt.rcParams["xtick.color"] = "black"
    plt.rcParams["ytick.color"] = "black"
    plt.rcParams["axes.grid"] = False
    sns.set_theme(style="white")


def save_pdf(fig: plt.Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    if output_path.stat().st_size <= 1024:
        raise RuntimeError(f"PDF output is unexpectedly small: {output_path}")


def apply_phase2b_axis_style(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.grid(False)
    ax.xaxis.grid(False, which="both")
    ax.yaxis.grid(False, which="both")
    ax.minorticks_off()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["left"].set_linewidth(1.0)
    ax.tick_params(axis="both", colors="black", width=0.9, length=3.5)


def normalize_label_regime(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["label_regime"] = out["protocol"].map(PROTOCOL_TO_REGIME).fillna(out.get("label_regime", ""))
    out["label_regime"] = pd.Categorical(out["label_regime"], categories=REGIME_ORDER, ordered=True)
    return out


def filter_expected_grid(per_run: pd.DataFrame, models: list[str], feature_settings: list[str], seeds: list[str]) -> pd.DataFrame:
    df = per_run.copy()
    df = df[df["protocol"].isin(PROTOCOL_TO_REGIME)].copy()
    df = df[df["model"].isin(models)].copy()
    df = df[df["feature_setting"].isin(feature_settings)].copy()
    if seeds:
        df = df[df["seed"].astype(str).isin(seeds)].copy()
    return df.reset_index(drop=True)


def add_public_metric_aliases(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for source, public_name in PUBLIC_METRIC_COLUMNS.items():
        mean_column = f"{source}_mean"
        std_column = f"{source}_std"
        if mean_column in out.columns:
            out[f"{public_name}_mean"] = out[mean_column]
        if std_column in out.columns:
            out[f"{public_name}_std"] = out[std_column]
    return out


def validate_grid(per_run: pd.DataFrame, models: list[str], feature_settings: list[str], seeds: list[str]) -> None:
    expected = len(PROTOCOL_TO_REGIME) * len(models) * len(feature_settings) * len(seeds)
    if len(per_run) != expected:
        observed = per_run.groupby(["protocol", "model", "feature_setting"], dropna=False).size().reset_index(name="n_runs")
        raise RuntimeError(
            f"Figure2a expected {expected} per-run rows "
            f"({len(models)} models x {len(feature_settings)} feature settings x {len(seeds)} seeds x 2 regimes), "
            f"found {len(per_run)}.\nObserved grid:\n{observed.to_string(index=False)}"
        )


def expected_aggregated(aggregated: pd.DataFrame, models: list[str], feature_settings: list[str]) -> pd.DataFrame:
    df = aggregated.copy()
    df["model"] = pd.Categorical(df["model"], categories=models, ordered=True)
    df["feature_setting"] = pd.Categorical(df["feature_setting"], categories=feature_settings, ordered=True)
    return df.sort_values(["label_regime", "model", "feature_setting"], kind="stable").reset_index(drop=True)


def build_comparison_table(aggregated: pd.DataFrame) -> pd.DataFrame:
    metrics = ["test_auroc_mean", "test_auprc_mean", "test_mcc_mean", "test_f1_mean"]
    pivot = aggregated.pivot_table(
        index=["model", "feature_setting"],
        columns="label_regime",
        values=metrics,
        aggfunc="first",
        observed=False,
    )
    pivot.columns = [f"{regime}_{metric}" for metric, regime in pivot.columns]
    out = pivot.reset_index()
    for metric in metrics:
        out[f"delta_{metric}_new_minus_old"] = out[f"new_label_{metric}"] - out[f"old440_{metric}"]
    return out.sort_values(["model", "feature_setting"], kind="stable").reset_index(drop=True)


def metric_matrix(aggregated: pd.DataFrame, label_regime: str, metric: str, models: list[str], feature_settings: list[str]) -> pd.DataFrame:
    df = aggregated[aggregated["label_regime"].astype(str) == label_regime].copy()
    matrix = df.pivot(index="model", columns="feature_setting", values=metric)
    return matrix.reindex(index=models, columns=feature_settings)


def plot_heatmap(
    matrix: pd.DataFrame,
    feature_settings: list[str],
    title: str,
    cbar_label: str,
    output_path: Path,
    cmap: str = "viridis",
    center: float | None = None,
    norm: Normalize | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 3.9))
    sns.heatmap(
        matrix,
        cmap=cmap,
        center=center,
        norm=norm,
        linewidths=0.5,
        linecolor="#d9d9d9",
        annot=True,
        fmt=".3f",
        cbar_kws={"label": cbar_label},
        ax=ax,
    )
    ax.set_xlabel("Feature combination")
    ax.set_ylabel("Model")
    ax.set_title(title)
    ax.set_xticklabels([FEATURE_LABELS.get(item, item) for item in feature_settings], rotation=35, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    save_pdf(fig, output_path)


def draw_errorbar_bars(ax: plt.Axes, errors: np.ndarray) -> None:
    for patch, err in zip(ax.patches, errors):
        if pd.isna(err):
            continue
        x = patch.get_x() + patch.get_width() / 2.0
        y = patch.get_height()
        ax.errorbar(x, y, yerr=float(err), color="black", capsize=2.5, linewidth=0.9, fmt="none", zorder=10)


def plot_regime_barplot(aggregated: pd.DataFrame, label_regime: str, models: list[str], feature_settings: list[str], output_path: Path) -> None:
    plot_df = aggregated[aggregated["label_regime"].astype(str) == label_regime].copy()
    plot_df["feature_label"] = plot_df["feature_setting"].map(FEATURE_LABELS)
    plot_df["feature_label"] = pd.Categorical(
        plot_df["feature_label"],
        categories=[FEATURE_LABELS[item] for item in feature_settings],
        ordered=True,
    )
    plot_df["model"] = pd.Categorical(plot_df["model"], categories=models, ordered=True)
    plot_df = plot_df.sort_values(["feature_label", "model"], kind="stable")
    fig, ax = plt.subplots(figsize=(10.0, 4.5))
    sns.barplot(
        data=plot_df,
        x="feature_label",
        y="test_auprc_mean",
        hue="model",
        palette=MODEL_PALETTE,
        errorbar=None,
        ax=ax,
    )
    draw_errorbar_bars(ax, plot_df["test_auprc_std"].to_numpy(dtype=float))
    ax.set_xlabel("Feature combination")
    ax.set_ylabel("AUPRC mean +/- seed SD")
    ax.set_title(f"Figure2a {label_regime}: 4 GNN x 7 feature combinations")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="Model", frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.20))
    apply_phase2b_axis_style(ax)
    save_pdf(fig, output_path)


def plot_combined_barplot(aggregated: pd.DataFrame, models: list[str], feature_settings: list[str], output_path: Path) -> None:
    plot_df = aggregated.copy()
    plot_df["feature_label"] = plot_df["feature_setting"].map(FEATURE_LABELS)
    plot_df["facet"] = plot_df["label_regime"].astype(str)
    plot_df["feature_label"] = pd.Categorical(
        plot_df["feature_label"],
        categories=[FEATURE_LABELS[item] for item in feature_settings],
        ordered=True,
    )
    plot_df["model"] = pd.Categorical(plot_df["model"], categories=models, ordered=True)
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.0), sharex=True, sharey=True)
    for ax, regime in zip(axes, REGIME_ORDER):
        sub = plot_df[plot_df["facet"] == regime].sort_values(["feature_label", "model"], kind="stable")
        sns.barplot(data=sub, x="feature_label", y="test_auprc_mean", hue="model", palette=MODEL_PALETTE, errorbar=None, ax=ax)
        draw_errorbar_bars(ax, sub["test_auprc_std"].to_numpy(dtype=float))
        ax.set_title(regime)
        ax.set_xlabel("")
        ax.set_ylabel("AUPRC mean +/- seed SD")
        if regime == REGIME_ORDER[0]:
            ax.legend(title="Model", frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.28))
        else:
            ax.get_legend().remove()
        apply_phase2b_axis_style(ax)
    axes[-1].set_xlabel("Feature combination")
    axes[-1].tick_params(axis="x", rotation=35)
    save_pdf(fig, output_path)


def build_best_combo_table(aggregated: pd.DataFrame, models: list[str]) -> pd.DataFrame:
    rows = []
    for model in models:
        for regime in REGIME_ORDER:
            subset = aggregated[(aggregated["model"].astype(str) == model) & (aggregated["label_regime"].astype(str) == regime)].copy()
            if subset.empty:
                raise RuntimeError(f"No aggregated rows for {model} / {regime}")
            best = subset.sort_values(["test_auprc_mean", "test_auroc_mean", "test_mcc_mean"], ascending=[False, False, False]).iloc[0]
            rows.append(
                {
                    "model": model,
                    "label_regime": regime,
                    "best_feature_setting": str(best["feature_setting"]),
                    "best_feature_label": FEATURE_LABELS.get(str(best["feature_setting"]), str(best["feature_setting"])),
                    "AUPRC_mean": float(best["test_auprc_mean"]),
                    "AUPRC_std": float(best["test_auprc_std"]),
                    "AUROC_mean": float(best["test_auroc_mean"]),
                    "AUROC_std": float(best["test_auroc_std"]),
                    "MCC_mean": float(best["test_mcc_mean"]),
                    "MCC_std": float(best["test_mcc_std"]),
                    "F1_mean": float(best["test_f1_mean"]),
                    "F1_std": float(best["test_f1_std"]),
                    "n_runs": int(best["n_runs"]),
                    "seed_list": str(best["seed_list"]),
                }
            )
    return pd.DataFrame(rows)


def plot_best_combo_barplot(best_df: pd.DataFrame, models: list[str], output_path: Path) -> None:
    plot_df = best_df.copy()
    plot_df["model"] = pd.Categorical(plot_df["model"], categories=models, ordered=True)
    plot_df["label_regime"] = pd.Categorical(plot_df["label_regime"], categories=REGIME_ORDER, ordered=True)
    plot_df = plot_df.sort_values(["model", "label_regime"], kind="stable")
    fig, ax = plt.subplots(figsize=(7.8, 4.0))
    sns.barplot(data=plot_df, x="model", y="AUPRC_mean", hue="label_regime", palette=OLD_NEW_PALETTE, errorbar=None, ax=ax)
    draw_errorbar_bars(ax, plot_df["AUPRC_std"].to_numpy(dtype=float))
    ax.set_xlabel("Model")
    ax.set_ylabel("Best-combo AUPRC mean +/- seed SD")
    ax.set_title("Figure2a best feature combination per model: old440 vs new label")
    ax.legend(title="Label regime", frameon=False)
    ymin = max(0.0, float(plot_df["AUPRC_mean"].min()) - float(plot_df["AUPRC_std"].max()) - 0.02)
    ymax = min(1.0, float(plot_df["AUPRC_mean"].max()) + float(plot_df["AUPRC_std"].max()) + 0.05)
    ax.set_ylim(ymin, ymax)
    apply_phase2b_axis_style(ax)
    save_pdf(fig, output_path)


def parse_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(BOOL_TRUE)


def pr_curve_points(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    order = np.argsort(-scores, kind="mergesort")
    labels_sorted = labels[order]
    positives = int(np.sum(labels_sorted == 1))
    if positives == 0:
        raise ValueError("PR requires at least one positive label.")
    tp = np.cumsum(labels_sorted == 1)
    fp = np.cumsum(labels_sorted == 0)
    recall = tp / positives
    precision = tp / np.maximum(tp + fp, 1)
    recall_curve = np.concatenate(([0.0], recall))
    precision_curve = np.concatenate(([1.0], precision))
    auprc = float(np.sum(np.diff(recall_curve) * precision_curve[1:]))
    return recall_curve, precision_curve, auprc


def load_test_predictions(predictions_path: Path) -> pd.DataFrame:
    df = pd.read_csv(predictions_path, sep="\t")
    df["is_labeled_flag"] = parse_bool_series(df["is_labeled"])
    df["label"] = pd.to_numeric(df["label"], errors="coerce")
    df["pred_score"] = pd.to_numeric(df["pred_score"], errors="coerce")
    filtered = df[(df["is_labeled_flag"]) & (df["split"] == "test")].copy()
    return filtered.dropna(subset=["label", "pred_score"]).copy()


def build_pr_curve_data(output_root: Path, aggregated: pd.DataFrame, models: list[str]) -> pd.DataFrame:
    rows = []
    for model in models:
        new_rows = aggregated[(aggregated["model"].astype(str) == model) & (aggregated["label_regime"].astype(str) == "new_label")].copy()
        best_new = new_rows.sort_values(["test_auprc_mean", "test_auroc_mean", "test_mcc_mean"], ascending=[False, False, False]).iloc[0]
        feature_setting = str(best_new["feature_setting"])
        for regime in REGIME_ORDER:
            protocol = {value: key for key, value in PROTOCOL_TO_REGIME.items()}[regime]
            predictions = []
            run_paths = sorted((output_root / protocol / model / feature_setting).glob("run_*/predictions.tsv"))
            if not run_paths:
                raise FileNotFoundError(f"No predictions.tsv files found for {protocol}/{model}/{feature_setting}")
            for path in run_paths:
                predictions.append(load_test_predictions(path))
            combined = pd.concat(predictions, ignore_index=True)
            labels = combined["label"].to_numpy(dtype=int)
            scores = combined["pred_score"].to_numpy(dtype=float)
            recalls, precisions, auprc = pr_curve_points(labels, scores)
            for idx, (recall, precision) in enumerate(zip(recalls, precisions)):
                rows.append(
                    {
                        "model": model,
                        "label_regime": regime,
                        "feature_setting": feature_setting,
                        "feature_label": FEATURE_LABELS.get(feature_setting, feature_setting),
                        "selection_strategy": "per_model_best_new_label_combo",
                        "curve_index": int(idx),
                        "recall": float(recall),
                        "precision": float(precision),
                        "pooled_test_auprc": float(auprc),
                        "n_test_predictions": int(len(combined)),
                        "n_source_runs": int(len(run_paths)),
                    }
                )
    return pd.DataFrame(rows)


def plot_pr_curves(pr_df: pd.DataFrame, models: list[str], output_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6), sharex=True, sharey=True)
    axes_flat = axes.ravel()
    colors = OLD_NEW_PALETTE
    for ax, model in zip(axes_flat, models):
        sub_model = pr_df[pr_df["model"] == model]
        feature_label = str(sub_model["feature_label"].iloc[0])
        for regime in REGIME_ORDER:
            sub = sub_model[sub_model["label_regime"] == regime]
            auprc = float(sub["pooled_test_auprc"].iloc[0])
            ax.plot(sub["recall"], sub["precision"], color=colors[regime], linewidth=1.8, label=f"{regime} AUPRC={auprc:.3f}")
        ax.set_title(f"{model} | {feature_label}")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.legend(frameon=False, loc="upper right")
        apply_phase2b_axis_style(ax)
    save_pdf(fig, output_path)


def main() -> None:
    args = parse_args()
    configure_plotting()
    output_root = Path(args.output_root).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    summary_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix
    models = split_csv(args.models)
    feature_settings = split_csv(args.feature_settings)
    seeds = split_csv(args.seeds)

    per_run = collect_metrics(output_root)
    per_run = filter_expected_grid(per_run, models, feature_settings, seeds)
    if seeds:
        validate_grid(per_run, models, feature_settings, seeds)
    per_run = normalize_label_regime(per_run)
    per_run["target"] = args.target_name
    per_run["feature_label"] = per_run["feature_setting"].map(FEATURE_LABELS)
    per_run = per_run.sort_values(["label_regime", "model", "feature_setting", "seed"], kind="stable").reset_index(drop=True)

    aggregated = aggregate_runs(per_run)
    aggregated = normalize_label_regime(aggregated)
    aggregated["target"] = args.target_name
    aggregated["feature_label"] = aggregated["feature_setting"].map(FEATURE_LABELS)
    aggregated = add_public_metric_aliases(expected_aggregated(aggregated, models, feature_settings))
    publication = build_publication_summary(aggregated)
    comparison = build_comparison_table(aggregated)
    best_combo = build_best_combo_table(aggregated, models)
    pr_curve_data = build_pr_curve_data(output_root, aggregated, models)

    per_run.to_csv(summary_dir / f"{prefix}_per_run_metrics.tsv", sep="\t", index=False)
    aggregated.to_csv(summary_dir / f"{prefix}_aggregated_metrics.tsv", sep="\t", index=False)
    publication.to_csv(summary_dir / f"{prefix}_final_summary.tsv", sep="\t", index=False)
    comparison.to_csv(summary_dir / f"{prefix}_old440_vs_newlabel_comparison.tsv", sep="\t", index=False)
    best_combo.to_csv(summary_dir / f"{prefix}_best_combo_selection.tsv", sep="\t", index=False)
    pr_curve_data.to_csv(summary_dir / f"{prefix}_best_combo_pr_curve_data.tsv", sep="\t", index=False)

    (summary_dir / f"{prefix}_per_run_metrics.md").write_text(
        dataframe_to_markdown(per_run, f"{prefix} Per-Run Metrics", "Figure2a per-seed metrics for 4 GNN x 7 feature combinations x old440/new_label."),
        encoding="utf-8",
    )
    (summary_dir / f"{prefix}_aggregated_metrics.md").write_text(
        dataframe_to_markdown(
            aggregated,
            f"{prefix} Aggregated Metrics",
            "Mean and seed-standard-deviation metrics grouped by model, label regime, and feature combination.",
        ),
        encoding="utf-8",
    )
    (summary_dir / f"{prefix}_final_summary.md").write_text(
        publication_markdown(publication, f"{prefix} Final Summary", "Publication-facing Figure2a summary."),
        encoding="utf-8",
    )

    old_auprc = metric_matrix(aggregated, "old440", "test_auprc_mean", models, feature_settings)
    new_auprc = metric_matrix(aggregated, "new_label", "test_auprc_mean", models, feature_settings)
    old_mcc = metric_matrix(aggregated, "old440", "test_mcc_mean", models, feature_settings)
    new_mcc = metric_matrix(aggregated, "new_label", "test_mcc_mean", models, feature_settings)
    delta_auprc = new_auprc - old_auprc
    delta_mcc = new_mcc - old_mcc
    delta_auprc_values = delta_auprc.to_numpy(dtype=float)
    delta_mcc_values = delta_mcc.to_numpy(dtype=float)
    delta_auprc_norm = Normalize(vmin=float(np.nanmin(delta_auprc_values)), vmax=float(np.nanmax(delta_auprc_values)))
    delta_mcc_norm = Normalize(vmin=float(np.nanmin(delta_mcc_values)), vmax=float(np.nanmax(delta_mcc_values)))
    plot_heatmap(old_auprc, feature_settings, "Figure2a old440 AUPRC mean", "AUPRC mean", summary_dir / f"{prefix}_heatmap_old440_auprc.pdf")
    plot_heatmap(new_auprc, feature_settings, "Figure2a new label AUPRC mean", "AUPRC mean", summary_dir / f"{prefix}_heatmap_newlabel_auprc.pdf")
    plot_heatmap(old_mcc, feature_settings, "Figure2a old440 MCC mean", "MCC mean", summary_dir / f"{prefix}_heatmap_old440_mcc.pdf")
    plot_heatmap(new_mcc, feature_settings, "Figure2a new label MCC mean", "MCC mean", summary_dir / f"{prefix}_heatmap_newlabel_mcc.pdf")
    plot_heatmap(delta_auprc, feature_settings, "Figure2a delta AUPRC mean (new - old)", "Delta AUPRC", summary_dir / f"{prefix}_delta_heatmap_auprc.pdf", cmap="coolwarm", norm=delta_auprc_norm)
    plot_heatmap(delta_mcc, feature_settings, "Figure2a delta MCC mean (new - old)", "Delta MCC", summary_dir / f"{prefix}_delta_heatmap_mcc.pdf", cmap="coolwarm", norm=delta_mcc_norm)
    plot_regime_barplot(aggregated, "old440", models, feature_settings, summary_dir / f"{prefix}_barplot_old440_auprc.pdf")
    plot_regime_barplot(aggregated, "new_label", models, feature_settings, summary_dir / f"{prefix}_barplot_newlabel_auprc.pdf")
    plot_combined_barplot(aggregated, models, feature_settings, summary_dir / f"{prefix}_barplot_old440_newlabel_auprc.pdf")
    plot_best_combo_barplot(best_combo, models, summary_dir / f"{prefix}_best_feature_old_vs_new_barplot.pdf")
    plot_pr_curves(pr_curve_data, models, summary_dir / f"{prefix}_best_combo_pr_curve.pdf")


if __name__ == "__main__":
    main()
