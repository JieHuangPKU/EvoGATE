"""
Build a complete F. graminearum benchmark summary from aggregated benchmark outputs.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


MODELS = ["GCN", "GAT", "GraphSAGE", "GIN"]
FEATURE_COMBOS = ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB"]
FEATURE_LABELS = {
    "ORT": "ortholog",
    "EXP": "expression",
    "SUB": "subloc",
    "ORT_EXP": "ortholog+expression",
    "ORT_SUB": "ortholog+subloc",
    "EXP_SUB": "expression+subloc",
    "ORT_EXP_SUB": "ortholog+expression+subloc",
}
METRIC_MAP = {
    "auroc_mean": "AUROC_mean",
    "auroc_std": "AUROC_std",
    "auprc_mean": "AUPRC_mean",
    "auprc_std": "AUPRC_std",
    "f1_mean": "F1_mean",
    "f1_std": "F1_std",
    "mcc_mean": "MCC_mean",
    "mcc_std": "MCC_std",
    "acc_mean": "ACC_mean",
    "acc_std": "ACC_std",
}
DELTA_SPECS = [
    ("ortholog_full_delta", "ortholog", "ORT_EXP_SUB", "EXP_SUB"),
    ("expression_full_delta", "expression", "ORT_EXP_SUB", "ORT_SUB"),
    ("subloc_full_delta", "subloc", "ORT_EXP_SUB", "ORT_EXP"),
    ("ortholog_to_ortholog_expression", "expression_on_ortholog_base", "ORT_EXP", "ORT"),
    ("ortholog_to_ortholog_subloc", "subloc_on_ortholog_base", "ORT_SUB", "ORT"),
    ("expression_to_expression_subloc", "subloc_on_expression_base", "EXP_SUB", "EXP"),
]


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize F. graminearum feature-combination benchmark results")
    parser.add_argument(
        "--input-root",
        default="outputs/fgraminearum_feature_combo_benchmark",
        help="Root containing {model}/{feature_combo}/thr_300/aggregated_metrics.tsv",
    )
    parser.add_argument(
        "--summary-root",
        default="outputs/fgraminearum_feature_combo_benchmark/summary",
        help="Directory for summary tables and figures",
    )
    return parser.parse_args()


def configure_matplotlib():
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 9
    plt.rcParams["axes.titlesize"] = 10
    plt.rcParams["axes.labelsize"] = 9
    plt.rcParams["xtick.labelsize"] = 8
    plt.rcParams["ytick.labelsize"] = 8
    plt.rcParams["legend.fontsize"] = 8
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    sns.set_theme(style="white")


def load_aggregated(input_root):
    input_root = Path(input_root)
    rows = []
    missing = []
    for model in MODELS:
        for combo in FEATURE_COMBOS:
            path = input_root / model / combo / "thr_300" / "aggregated_metrics.tsv"
            if not path.exists():
                missing.append(str(path))
                continue
            row = pd.read_csv(path, sep="\t").iloc[0].to_dict()
            row["source_path"] = str(path)
            rows.append(row)
    if missing:
        raise FileNotFoundError("Missing aggregated_metrics.tsv files:\n{}".format("\n".join(missing)))
    df = pd.DataFrame(rows)
    expected = len(MODELS) * len(FEATURE_COMBOS)
    if len(df) != expected:
        raise RuntimeError("Expected {} aggregated rows, found {}".format(expected, len(df)))
    df = df.rename(columns=METRIC_MAP)
    df["feature_label"] = df["feature_combo"].map(FEATURE_LABELS)
    df["model"] = pd.Categorical(df["model"], categories=MODELS, ordered=True)
    df["feature_combo"] = pd.Categorical(df["feature_combo"], categories=FEATURE_COMBOS, ordered=True)
    df = df.sort_values(["model", "feature_combo"]).reset_index(drop=True)
    return df


def build_model_feature_matrix(df):
    cols = [
        "model",
        "feature_combo",
        "feature_label",
        "AUROC_mean",
        "AUROC_std",
        "AUPRC_mean",
        "AUPRC_std",
        "F1_mean",
        "F1_std",
        "MCC_mean",
        "MCC_std",
        "ACC_mean",
        "ACC_std",
        "runs_completed",
    ]
    return df[cols].copy()


def build_best_tables(df):
    best_per_model = (
        df.sort_values(["model", "AUPRC_mean", "AUROC_mean", "MCC_mean"], ascending=[True, False, False, False])
        .groupby("model", as_index=False)
        .first()
    )
    model_mean = (
        df.groupby("model", as_index=False)
        .agg(
            AUPRC_mean_overall=("AUPRC_mean", "mean"),
            AUROC_mean_overall=("AUROC_mean", "mean"),
            MCC_mean_overall=("MCC_mean", "mean"),
            F1_mean_overall=("F1_mean", "mean"),
        )
        .sort_values(["AUPRC_mean_overall", "AUROC_mean_overall", "MCC_mean_overall"], ascending=[False, False, False])
    )
    best_overall = df.sort_values(["AUPRC_mean", "AUROC_mean", "MCC_mean"], ascending=[False, False, False]).iloc[0]
    out_rows = []
    for _, row in best_per_model.iterrows():
        out_rows.append(
            {
                "summary_type": "best_feature_combo_per_model",
                "model": str(row["model"]),
                "feature_combo": str(row["feature_combo"]),
                "feature_label": str(row["feature_label"]),
                "AUPRC_mean": float(row["AUPRC_mean"]),
                "AUROC_mean": float(row["AUROC_mean"]),
                "MCC_mean": float(row["MCC_mean"]),
                "F1_mean": float(row["F1_mean"]),
            }
        )
    out_rows.append(
        {
            "summary_type": "best_model_overall_by_mean_AUPRC",
            "model": str(model_mean.iloc[0]["model"]),
            "feature_combo": "",
            "feature_label": "",
            "AUPRC_mean": float(model_mean.iloc[0]["AUPRC_mean_overall"]),
            "AUROC_mean": float(model_mean.iloc[0]["AUROC_mean_overall"]),
            "MCC_mean": float(model_mean.iloc[0]["MCC_mean_overall"]),
            "F1_mean": float(model_mean.iloc[0]["F1_mean_overall"]),
        }
    )
    out_rows.append(
        {
            "summary_type": "best_feature_combo_overall",
            "model": str(best_overall["model"]),
            "feature_combo": str(best_overall["feature_combo"]),
            "feature_label": str(best_overall["feature_label"]),
            "AUPRC_mean": float(best_overall["AUPRC_mean"]),
            "AUROC_mean": float(best_overall["AUROC_mean"]),
            "MCC_mean": float(best_overall["MCC_mean"]),
            "F1_mean": float(best_overall["F1_mean"]),
        }
    )
    return pd.DataFrame(out_rows), best_per_model, model_mean, best_overall


def compute_feature_importance(df):
    index_df = df.set_index(["model", "feature_combo"])
    rows = []
    for model in MODELS:
        for delta_name, feature_name, numer_combo, denom_combo in DELTA_SPECS:
            num = index_df.loc[(model, numer_combo)]
            den = index_df.loc[(model, denom_combo)]
            row = {
                "model": model,
                "comparison": delta_name,
                "feature": feature_name,
                "numerator_combo": numer_combo,
                "denominator_combo": denom_combo,
            }
            for metric in ["AUPRC_mean", "AUROC_mean", "F1_mean", "MCC_mean"]:
                row[metric + "_delta"] = float(num[metric] - den[metric])
            rows.append(row)
    feature_importance = pd.DataFrame(rows)

    full_ablation = feature_importance[feature_importance["comparison"].isin(["ortholog_full_delta", "expression_full_delta", "subloc_full_delta"])].copy()
    ranking = (
        full_ablation.groupby("feature", as_index=False)
        .agg(
            AUPRC_mean_delta=("AUPRC_mean_delta", "mean"),
            AUROC_mean_delta=("AUROC_mean_delta", "mean"),
            MCC_mean_delta=("MCC_mean_delta", "mean"),
            F1_mean_delta=("F1_mean_delta", "mean"),
            positive_AUPRC_models=("AUPRC_mean_delta", lambda s: int((s > 0).sum())),
            positive_AUROC_models=("AUROC_mean_delta", lambda s: int((s > 0).sum())),
        )
        .sort_values(["AUPRC_mean_delta", "AUROC_mean_delta", "MCC_mean_delta"], ascending=[False, False, False])
        .reset_index(drop=True)
    )
    ranking["importance_rank"] = range(1, len(ranking) + 1)
    return feature_importance, ranking


def compute_stability(df):
    model_stability = (
        df.groupby("model", as_index=False)
        .agg(
            mean_AUPRC_std=("AUPRC_std", "mean"),
            mean_AUROC_std=("AUROC_std", "mean"),
            mean_MCC_std=("MCC_std", "mean"),
            mean_F1_std=("F1_std", "mean"),
        )
    )
    model_stability["stability_score"] = model_stability[
        ["mean_AUPRC_std", "mean_AUROC_std", "mean_MCC_std", "mean_F1_std"]
    ].mean(axis=1)
    model_stability = model_stability.sort_values(["stability_score", "mean_AUPRC_std"], ascending=[True, True]).reset_index(drop=True)

    combo_stability = (
        df.groupby("feature_combo", as_index=False)
        .agg(
            mean_AUPRC_std=("AUPRC_std", "mean"),
            mean_AUROC_std=("AUROC_std", "mean"),
            mean_MCC_std=("MCC_std", "mean"),
            mean_F1_std=("F1_std", "mean"),
        )
    )
    combo_stability["instability_score"] = combo_stability[
        ["mean_AUPRC_std", "mean_AUROC_std", "mean_MCC_std", "mean_F1_std"]
    ].mean(axis=1)
    combo_stability = combo_stability.sort_values(["instability_score", "mean_AUPRC_std"], ascending=[False, False]).reset_index(drop=True)
    return model_stability, combo_stability


def save_pdf(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="pdf", bbox_inches=None)
    plt.close(fig)
    size = path.stat().st_size
    if size <= 10 * 1024:
        raise RuntimeError("PDF output is unexpectedly small ({} bytes): {}".format(size, path))


def make_heatmap(df, output_path):
    pivot = (
        df.pivot(index="model", columns="feature_combo", values="AUPRC_mean")
        .reindex(index=MODELS, columns=FEATURE_COMBOS)
    )
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    sns.heatmap(
        pivot,
        cmap="viridis",
        linewidths=0.5,
        linecolor="#d9d9d9",
        square=True,
        cbar_kws={"label": "AUPRC mean"},
        annot=True,
        fmt=".3f",
        ax=ax,
    )
    ax.set_xlabel("Feature combination")
    ax.set_ylabel("Model")
    ax.set_title("Fusarium benchmark: AUPRC across models and feature combinations")
    ax.set_xticklabels([FEATURE_LABELS[c] for c in FEATURE_COMBOS], rotation=35, ha="right")
    ax.set_yticklabels(MODELS, rotation=0)
    fig.subplots_adjust(left=0.14, right=0.96, bottom=0.29, top=0.88)
    save_pdf(fig, output_path)


def make_feature_barplot(feature_importance, output_path):
    plot_df = feature_importance[
        feature_importance["comparison"].isin(["ortholog_full_delta", "expression_full_delta", "subloc_full_delta"])
    ].copy()
    name_map = {
        "ortholog_full_delta": "Δortholog",
        "expression_full_delta": "Δexpression",
        "subloc_full_delta": "Δsubloc",
    }
    plot_df["comparison_label"] = plot_df["comparison"].map(name_map)
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    sns.barplot(
        data=plot_df,
        x="model",
        y="AUPRC_mean_delta",
        hue="comparison_label",
        palette=["#1f77b4", "#ff7f0e", "#2ca02c"],
        ax=ax,
    )
    ax.set_xlabel("Model")
    ax.set_ylabel("AUPRC delta")
    ax.set_title("Feature contribution to AUPRC by model")
    ax.grid(False)
    sns.despine(ax=ax)
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.18, top=0.87)
    save_pdf(fig, output_path)


def make_model_best_barplot(best_per_model, output_path):
    plot_df = best_per_model.sort_values(["AUPRC_mean", "AUROC_mean"], ascending=[True, True]).copy()
    plot_df["label"] = plot_df["model"].astype(str) + " | " + plot_df["feature_label"].astype(str)
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.barh(plot_df["label"], plot_df["AUPRC_mean"], color="#4c78a8")
    ax.set_xlabel("AUPRC mean")
    ax.set_ylabel("")
    ax.set_title("Best feature combination per model")
    ax.grid(False)
    sns.despine(ax=ax)
    fig.subplots_adjust(left=0.29, right=0.97, bottom=0.19, top=0.86)
    save_pdf(fig, output_path)


def build_interpretation(best_per_model, model_mean, best_overall, ranking, feature_importance, model_stability, combo_stability):
    ortholog_consistent = bool((feature_importance[feature_importance["comparison"] == "ortholog_full_delta"]["AUPRC_mean_delta"] > 0).all())
    top_feature = str(ranking.iloc[0]["feature"])
    top_feature_auprc = float(ranking.iloc[0]["AUPRC_mean_delta"])
    ortholog_row = ranking[ranking["feature"] == "ortholog"].iloc[0]
    expression_row = ranking[ranking["feature"] == "expression"].iloc[0]
    subloc_row = ranking[ranking["feature"] == "subloc"].iloc[0]
    multi_best_count = int(best_per_model["feature_combo"].astype(str).str.contains("_").sum())
    most_stable_model = str(model_stability.iloc[0]["model"])
    most_unstable_combo = str(combo_stability.iloc[0]["feature_combo"])

    lines = [
        "# Fusarium Feature Combination Benchmark Interpretation",
        "",
        "The benchmark comprises 4 graph models, 7 feature combinations, and 3 runs per configuration under a fixed `STRING threshold = 300` and `include_degree = false`.",
        "",
        "## Main Findings",
        "",
        "1. Which feature contributes most?",
        "The largest average gain is attributable to `{}`. Across models, its full-ablation contrast yields a mean AUPRC gain of {:.4f}, exceeding the corresponding gains for expression ({:.4f}) and subloc ({:.4f}).".format(
            top_feature,
            top_feature_auprc,
            float(expression_row["AUPRC_mean_delta"]),
            float(subloc_row["AUPRC_mean_delta"]),
        ),
        "2. Does ortholog dominate?",
        "{}. Ortholog produces a mean AUPRC gain of {:.4f} and is positive in {}/4 models in the full-ablation comparison.".format(
            "Yes" if ortholog_consistent and top_feature == "ortholog" else "Not uniformly",
            float(ortholog_row["AUPRC_mean_delta"]),
            int(ortholog_row["positive_AUPRC_models"]),
        ),
        "3. Does expression improve performance?",
        "{}. Expression increases AUPRC by a mean of {:.4f}; the effect is positive in {}/4 models.".format(
            "Yes" if float(expression_row["AUPRC_mean_delta"]) > 0 else "Not reliably",
            float(expression_row["AUPRC_mean_delta"]),
            int(expression_row["positive_AUPRC_models"]),
        ),
        "4. Is subloc useful?",
        "{}. Subloc shows a mean AUPRC effect of {:.4f}; this indicates {} contribution rather than a universally beneficial signal.".format(
            "Yes, but weakly" if float(subloc_row["AUPRC_mean_delta"]) > 0 else "No",
            float(subloc_row["AUPRC_mean_delta"]),
            "a modest positive" if float(subloc_row["AUPRC_mean_delta"]) > 0 else "a neutral-to-negative",
        ),
        "5. Is multi-feature always better?",
        "No. Multi-feature combinations are frequently competitive, but they are not universally optimal. Only {}/4 models achieve their best AUPRC with a multi-feature combination, so feature addition should be treated as model-dependent rather than automatically beneficial.".format(
            multi_best_count
        ),
        "6. Which model is best for Fusarium?",
        "By mean AUPRC across all feature combinations, the strongest model is `{}` (mean AUPRC {:.4f}). The single best configuration is `{}` with `{}` (AUPRC {:.4f}, AUROC {:.4f}, MCC {:.4f}).".format(
            str(model_mean.iloc[0]["model"]),
            float(model_mean.iloc[0]["AUPRC_mean_overall"]),
            str(best_overall["model"]),
            str(best_overall["feature_label"]),
            float(best_overall["AUPRC_mean"]),
            float(best_overall["AUROC_mean"]),
            float(best_overall["MCC_mean"]),
        ),
        "",
        "## Stability",
        "",
        "The most stable model by mean metric standard deviation is `{}`. The most unstable feature combination is `{}`.".format(
            most_stable_model,
            most_unstable_combo,
        ),
        "",
        "## Evidence Summary",
        "",
        "- Best feature per model was selected by `AUPRC_mean`.",
        "- Feature contributions were quantified with explicit ablation deltas relative to the full `ortholog+expression+subloc` configuration and with pairwise gains on simpler bases.",
        "- All claims above are derived from the aggregated 3-run summaries under `outputs/fgraminearum_feature_combo_benchmark/*/*/thr_300/aggregated_metrics.tsv`.",
    ]
    return "\n".join(lines) + "\n"


def build_figure_legend():
    lines = [
        "# Figure Legends",
        "",
        "## Figure 1. Heatmap of model-by-feature-combination AUPRC.",
        "Cells report mean AUPRC across three independent runs for each graph model and feature combination in *F. graminearum*. Rows correspond to models and columns correspond to the seven tested feature combinations. Higher values indicate improved precision-recall performance.",
        "",
        "## Figure 2. Feature contribution to AUPRC by model.",
        "Bars show ablation-defined AUPRC deltas for ortholog, expression, and subloc features within each model. Positive values indicate that adding the focal feature improves the full three-feature configuration relative to the matched two-feature baseline.",
        "",
        "## Figure 3. Best feature combination per model.",
        "Horizontal bars compare the highest-AUPRC feature combination achieved by each graph model. Bars are ordered by AUPRC mean to facilitate direct model comparison under the best-performing feature setting for each model.",
    ]
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    configure_matplotlib()

    summary_root = Path(args.summary_root)
    figures_root = summary_root / "figures"
    summary_root.mkdir(parents=True, exist_ok=True)
    figures_root.mkdir(parents=True, exist_ok=True)

    df = load_aggregated(args.input_root)
    model_feature_matrix = build_model_feature_matrix(df)
    best_config, best_per_model, model_mean, best_overall = build_best_tables(df)
    feature_importance, feature_ranking = compute_feature_importance(df)
    model_stability, combo_stability = compute_stability(df)

    model_feature_matrix.to_csv(summary_root / "model_feature_matrix.tsv", sep="\t", index=False)
    best_config.to_csv(summary_root / "best_config.tsv", sep="\t", index=False)
    feature_importance.to_csv(summary_root / "feature_importance.tsv", sep="\t", index=False)
    feature_ranking.to_csv(summary_root / "feature_importance_ranking.tsv", sep="\t", index=False)
    model_stability.to_csv(summary_root / "model_stability.tsv", sep="\t", index=False)
    combo_stability.to_csv(summary_root / "feature_combo_stability.tsv", sep="\t", index=False)

    interpretation = build_interpretation(
        best_per_model,
        model_mean,
        best_overall,
        feature_ranking,
        feature_importance,
        model_stability,
        combo_stability,
    )
    (summary_root / "interpretation.md").write_text(interpretation, encoding="utf-8")
    (summary_root / "figure_legend.md").write_text(build_figure_legend(), encoding="utf-8")

    make_heatmap(df, figures_root / "heatmap_model_feature_auprc.pdf")
    make_feature_barplot(feature_importance, figures_root / "feature_contribution_barplot.pdf")
    make_model_best_barplot(best_per_model, figures_root / "model_best_comparison.pdf")


if __name__ == "__main__":
    main()
