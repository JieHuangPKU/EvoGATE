import argparse
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


PROTOCOLS = ["fgraminearum_oldlabel", "fgraminearum_newlabel"]
MODELS = ["GraphSAGE", "GCN", "GAT", "GIN"]
FEATURE_COMBOS = ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB"]
METRIC_COLUMNS = ["test_auroc_mean", "test_auprc_mean", "test_mcc_mean", "test_f1_mean", "test_accuracy_mean"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze Fusarium old-vs-new label feature-combo benchmark outputs")
    parser.add_argument("--summary-dir", required=True, type=str)
    parser.add_argument("--analysis-dir", required=True, type=str)
    return parser.parse_args()


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
    sns.set_theme(style="whitegrid")


def save_pdf(fig, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def write_markdown_table(df: pd.DataFrame, title: str, intro: str, output_path: Path) -> None:
    lines = [f"# {title}", "", intro, ""]
    if df.empty:
        lines.extend(["No rows available.", ""])
    else:
        lines.append(df.to_markdown(index=False))
        lines.append("")
    output_path.write_text("\n".join(lines), encoding="utf-8")


def load_inputs(summary_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    aggregated = pd.read_csv(summary_dir / "aggregated_metrics.tsv", sep="\t")
    per_run = pd.read_csv(summary_dir / "per_run_metrics.tsv", sep="\t")
    aggregated = aggregated[aggregated["protocol"].isin(PROTOCOLS)].copy()
    aggregated = aggregated[aggregated["model"].isin(MODELS)].copy()
    aggregated = aggregated[aggregated["feature_setting"].isin(FEATURE_COMBOS)].copy()
    per_run = per_run[per_run["protocol"].isin(PROTOCOLS)].copy()
    per_run = per_run[per_run["model"].isin(MODELS)].copy()
    per_run = per_run[per_run["feature_setting"].isin(FEATURE_COMBOS)].copy()
    aggregated["protocol"] = pd.Categorical(aggregated["protocol"], categories=PROTOCOLS, ordered=True)
    aggregated["model"] = pd.Categorical(aggregated["model"], categories=MODELS, ordered=True)
    aggregated["feature_setting"] = pd.Categorical(aggregated["feature_setting"], categories=FEATURE_COMBOS, ordered=True)
    aggregated = aggregated.sort_values(["protocol", "model", "feature_setting"]).reset_index(drop=True)
    return aggregated, per_run


def build_old_vs_new_comparison(aggregated: pd.DataFrame) -> pd.DataFrame:
    pivot = aggregated.pivot_table(
        index=["model", "feature_setting"],
        columns="protocol",
        values=["test_auprc_mean", "test_auroc_mean", "test_mcc_mean", "n_runs"],
        aggfunc="first",
    )
    pivot.columns = [f"{metric}_{protocol}" for metric, protocol in pivot.columns]
    out = pivot.reset_index()
    out["delta_test_auprc_mean_new_minus_old"] = out["test_auprc_mean_fgraminearum_newlabel"] - out["test_auprc_mean_fgraminearum_oldlabel"]
    out["delta_test_auroc_mean_new_minus_old"] = out["test_auroc_mean_fgraminearum_newlabel"] - out["test_auroc_mean_fgraminearum_oldlabel"]
    out["delta_test_mcc_mean_new_minus_old"] = out["test_mcc_mean_fgraminearum_newlabel"] - out["test_mcc_mean_fgraminearum_oldlabel"]
    return out.sort_values(["delta_test_auprc_mean_new_minus_old", "delta_test_auroc_mean_new_minus_old"], ascending=[False, False]).reset_index(drop=True)


def build_best_by_model_and_feature(aggregated: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protocol in PROTOCOLS:
        protocol_df = aggregated[aggregated["protocol"] == protocol].copy()
        best_model = protocol_df.sort_values(["test_auprc_mean", "test_auroc_mean", "test_mcc_mean"], ascending=[False, False, False]).iloc[0]
        best_feature = (
            protocol_df.groupby("feature_setting", as_index=False)
            .agg(
                test_auprc_mean=("test_auprc_mean", "mean"),
                test_auroc_mean=("test_auroc_mean", "mean"),
                test_mcc_mean=("test_mcc_mean", "mean"),
            )
            .sort_values(["test_auprc_mean", "test_auroc_mean", "test_mcc_mean"], ascending=[False, False, False])
            .iloc[0]
        )
        rows.append(
            {
                "summary_type": "best_model_setting",
                "protocol": protocol,
                "model": str(best_model["model"]),
                "feature_setting": str(best_model["feature_setting"]),
                "test_auprc_mean": float(best_model["test_auprc_mean"]),
                "test_auroc_mean": float(best_model["test_auroc_mean"]),
                "test_mcc_mean": float(best_model["test_mcc_mean"]),
            }
        )
        rows.append(
            {
                "summary_type": "best_feature_overall",
                "protocol": protocol,
                "model": "",
                "feature_setting": str(best_feature["feature_setting"]),
                "test_auprc_mean": float(best_feature["test_auprc_mean"]),
                "test_auroc_mean": float(best_feature["test_auroc_mean"]),
                "test_mcc_mean": float(best_feature["test_mcc_mean"]),
            }
        )
    return pd.DataFrame(rows)


def build_publication_summary(aggregated: pd.DataFrame, old_vs_new: pd.DataFrame, best_table: pd.DataFrame) -> str:
    best_old = best_table[(best_table["summary_type"] == "best_model_setting") & (best_table["protocol"] == "fgraminearum_oldlabel")].iloc[0]
    best_new = best_table[(best_table["summary_type"] == "best_model_setting") & (best_table["protocol"] == "fgraminearum_newlabel")].iloc[0]
    best_feature_old = best_table[(best_table["summary_type"] == "best_feature_overall") & (best_table["protocol"] == "fgraminearum_oldlabel")].iloc[0]
    best_feature_new = best_table[(best_table["summary_type"] == "best_feature_overall") & (best_table["protocol"] == "fgraminearum_newlabel")].iloc[0]
    mean_delta_auprc = float(old_vs_new["delta_test_auprc_mean_new_minus_old"].mean())
    mean_delta_auroc = float(old_vs_new["delta_test_auroc_mean_new_minus_old"].mean())
    positive_new = int((old_vs_new["delta_test_auprc_mean_new_minus_old"] > 0).sum())
    total_rows = int(len(old_vs_new))
    lines = [
        "# Fusarium Old-vs-New Label Benchmark Summary",
        "",
        "## Scope",
        "",
        "- Compared `fgraminearum_oldlabel` versus `fgraminearum_newlabel` under the current processed-data frozen protocol pipeline.",
        "- Benchmarked 4 graph models: `GraphSAGE`, `GCN`, `GAT`, `GIN`.",
        "- Benchmarked 7 feature combinations: `ORT`, `EXP`, `SUB`, `ORT_EXP`, `ORT_SUB`, `EXP_SUB`, `ORT_EXP_SUB`.",
        "",
        "## Main Findings",
        "",
        "- Best setting for oldlabel: model `{}` with feature combination `{}` (test AUPRC mean {:.4f}, AUROC mean {:.4f}, MCC mean {:.4f}).".format(
            best_old["model"], best_old["feature_setting"], best_old["test_auprc_mean"], best_old["test_auroc_mean"], best_old["test_mcc_mean"]
        ),
        "- Best setting for newlabel: model `{}` with feature combination `{}` (test AUPRC mean {:.4f}, AUROC mean {:.4f}, MCC mean {:.4f}).".format(
            best_new["model"], best_new["feature_setting"], best_new["test_auprc_mean"], best_new["test_auroc_mean"], best_new["test_mcc_mean"]
        ),
        "- Best feature combination for oldlabel after averaging across models: `{}`.".format(best_feature_old["feature_setting"]),
        "- Best feature combination for newlabel after averaging across models: `{}`.".format(best_feature_new["feature_setting"]),
        "- Newlabel minus oldlabel mean delta across all model-feature settings: AUPRC {:.4f}, AUROC {:.4f}.".format(mean_delta_auprc, mean_delta_auroc),
        "- Newlabel is better than oldlabel on {}/{} model-feature configurations when ranked by mean test AUPRC.".format(positive_new, total_rows),
        "",
        "## Recommendation",
        "",
        "- Use the newlabel regime as the default Fusarium benchmark if the objective is current processed-data mainline evaluation.",
        "- Retain the oldlabel branch as an explicit historical comparison rather than merging the two regimes.",
        "- Prefer the best-scoring model-feature pair above as the publication-ready Fusarium setting.",
        "",
        "## Interpretation",
        "",
        "- Whether topology plus certain feature groups dominate performance should be judged from the per-regime heatmaps and the feature-combination ranking figure.",
        "- The benchmark remains fully explicit about label regime, model family, and feature combination; no regime is treated implicitly.",
    ]
    return "\n".join(lines) + "\n"


def plot_old_vs_new_overview(old_vs_new: pd.DataFrame, output_path: Path) -> None:
    plot_df = old_vs_new.copy()
    plot_df["label"] = plot_df["model"].astype(str) + " | " + plot_df["feature_setting"].astype(str)
    plot_df = plot_df.sort_values("delta_test_auprc_mean_new_minus_old", ascending=True)
    fig, ax = plt.subplots(figsize=(8.5, 6.8))
    ax.barh(plot_df["label"], plot_df["delta_test_auprc_mean_new_minus_old"], color="#1f77b4")
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("Newlabel - Oldlabel test AUPRC mean")
    ax.set_ylabel("")
    ax.set_title("Fusarium label-regime effect across models and feature combinations")
    save_pdf(fig, output_path)


def plot_protocol_heatmap(aggregated: pd.DataFrame, protocol: str, output_path: Path) -> None:
    protocol_df = aggregated[aggregated["protocol"] == protocol].copy()
    matrix = protocol_df.pivot(index="model", columns="feature_setting", values="test_auprc_mean").reindex(index=MODELS, columns=FEATURE_COMBOS)
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    sns.heatmap(matrix, cmap="viridis", linewidths=0.5, linecolor="#d9d9d9", annot=True, fmt=".3f", cbar_kws={"label": "test AUPRC mean"}, ax=ax)
    ax.set_xlabel("Feature combination")
    ax.set_ylabel("Model")
    ax.set_title(f"{protocol}: model-feature heatmap")
    save_pdf(fig, output_path)


def plot_feature_combo_rankings(aggregated: pd.DataFrame, output_path: Path) -> None:
    ranking = (
        aggregated.groupby(["protocol", "feature_setting"], as_index=False)
        .agg(
            test_auprc_mean=("test_auprc_mean", "mean"),
            test_auroc_mean=("test_auroc_mean", "mean"),
        )
        .sort_values(["protocol", "test_auprc_mean", "test_auroc_mean"], ascending=[True, False, False])
    )
    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    sns.barplot(data=ranking, x="feature_setting", y="test_auprc_mean", hue="protocol", palette=["#7f7f7f", "#4c78a8"], ax=ax)
    ax.set_xlabel("Feature combination")
    ax.set_ylabel("Mean test AUPRC across models")
    ax.set_title("Feature-combination ranking by label regime")
    ax.tick_params(axis="x", rotation=30)
    save_pdf(fig, output_path)


def main() -> None:
    args = parse_args()
    configure_plotting()
    summary_dir = Path(args.summary_dir).resolve()
    analysis_dir = Path(args.analysis_dir).resolve()
    analysis_dir.mkdir(parents=True, exist_ok=True)

    aggregated, per_run = load_inputs(summary_dir)
    old_vs_new = build_old_vs_new_comparison(aggregated)
    best_table = build_best_by_model_and_feature(aggregated)

    aggregated.to_csv(analysis_dir / "model_feature_regime_summary.tsv", sep="\t", index=False)
    per_run.to_csv(analysis_dir / "per_run_metrics.tsv", sep="\t", index=False)
    old_vs_new.to_csv(analysis_dir / "old_vs_new_comparison.tsv", sep="\t", index=False)
    best_table.to_csv(analysis_dir / "best_by_model_and_feature.tsv", sep="\t", index=False)

    write_markdown_table(
        old_vs_new,
        "Old vs New Comparison",
        "Each row compares `fgraminearum_newlabel` minus `fgraminearum_oldlabel` for one model-feature configuration.",
        analysis_dir / "old_vs_new_comparison.md",
    )
    write_markdown_table(
        best_table,
        "Best by Model and Feature",
        "Best-setting summaries per label regime under the frozen processed-data pipeline.",
        analysis_dir / "best_by_model_and_feature.md",
    )
    (analysis_dir / "publication_summary.md").write_text(build_publication_summary(aggregated, old_vs_new, best_table), encoding="utf-8")

    plot_old_vs_new_overview(old_vs_new, analysis_dir / "fig_old_vs_new_overview.pdf")
    plot_protocol_heatmap(aggregated, "fgraminearum_oldlabel", analysis_dir / "fig_model_feature_heatmap_oldlabel.pdf")
    plot_protocol_heatmap(aggregated, "fgraminearum_newlabel", analysis_dir / "fig_model_feature_heatmap_newlabel.pdf")
    plot_feature_combo_rankings(aggregated, analysis_dir / "fig_feature_combo_rankings.pdf")


if __name__ == "__main__":
    main()
