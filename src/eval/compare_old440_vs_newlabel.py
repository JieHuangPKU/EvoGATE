"""Compare the completed old440 benchmark against the new-label companion benchmark."""

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


def parse_args():
    parser = argparse.ArgumentParser(description="Compare old440 and new-label Fusarium feature-combo benchmarks")
    parser.add_argument("--old-root", default="outputs/fgraminearum_feature_combo_benchmark")
    parser.add_argument("--new-root", default="outputs/fgraminearum_feature_combo_newlabel_benchmark")
    parser.add_argument("--results-root", default="results/phase2b_new_label")
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


def _load_grid(root, prefix):
    rows = []
    missing = []
    root = Path(root)
    for model in MODELS:
        for combo in FEATURE_COMBOS:
            path = root / model / combo / "thr_300" / "aggregated_metrics.tsv"
            if not path.exists():
                missing.append(str(path))
                continue
            row = pd.read_csv(path, sep="\t").iloc[0].to_dict()
            rows.append(row)
    if missing:
        raise FileNotFoundError("Missing aggregated benchmark files:\n{}".format("\n".join(missing)))
    df = pd.DataFrame(rows)
    rename_map = {
        "auroc_mean": f"{prefix}_auroc_mean",
        "auprc_mean": f"{prefix}_auprc_mean",
        "mcc_mean": f"{prefix}_mcc_mean",
        "f1_mean": f"{prefix}_f1_mean",
        "auroc_std": f"{prefix}_auroc_std",
        "auprc_std": f"{prefix}_auprc_std",
        "mcc_std": f"{prefix}_mcc_std",
        "f1_std": f"{prefix}_f1_std",
    }
    df = df.rename(columns=rename_map)
    keep = ["model", "feature_combo"] + list(rename_map.values())
    return df[keep].copy()


def _save_pdf(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="pdf", bbox_inches=None)
    plt.close(fig)
    if path.stat().st_size <= 10 * 1024:
        raise RuntimeError(f"PDF output is unexpectedly small: {path}")


def _make_heatmap(df, value_col, title, output_path, cmap="viridis"):
    pivot = df.pivot(index="model", columns="feature_combo", values=value_col).reindex(index=MODELS, columns=FEATURE_COMBOS)
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    sns.heatmap(
        pivot,
        cmap=cmap,
        linewidths=0.5,
        linecolor="#d9d9d9",
        square=True,
        cbar_kws={"label": value_col},
        annot=True,
        fmt=".3f",
        ax=ax,
    )
    ax.set_xlabel("Feature combination")
    ax.set_ylabel("Model")
    ax.set_title(title)
    ax.set_xticklabels([FEATURE_LABELS[c] for c in FEATURE_COMBOS], rotation=35, ha="right")
    ax.set_yticklabels(MODELS, rotation=0)
    fig.subplots_adjust(left=0.14, right=0.96, bottom=0.29, top=0.88)
    _save_pdf(fig, output_path)


def _make_best_barplot(best_old, best_new, output_path):
    old_df = best_old[["model", "feature_combo", "old_auprc_mean"]].copy()
    old_df["regime"] = "old440"
    old_df["auprc_mean"] = old_df["old_auprc_mean"]
    new_df = best_new[["model", "feature_combo", "new_auprc_mean"]].copy()
    new_df["regime"] = "new_label"
    new_df["auprc_mean"] = new_df["new_auprc_mean"]
    plot_df = pd.concat([old_df[["model", "feature_combo", "regime", "auprc_mean"]], new_df[["model", "feature_combo", "regime", "auprc_mean"]]], ignore_index=True)
    plot_df["label"] = plot_df["model"] + " | " + plot_df["feature_combo"].map(FEATURE_LABELS)
    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    sns.barplot(data=plot_df, x="auprc_mean", y="model", hue="regime", palette=["#4c78a8", "#f58518"], ax=ax)
    ax.set_xlabel("Best-combo AUPRC mean")
    ax.set_ylabel("Model")
    ax.set_title("Best feature combination per model: old440 vs new label")
    ax.grid(False)
    sns.despine(ax=ax)
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.18, top=0.87)
    _save_pdf(fig, output_path)


def _compute_newlabel_feature_effects(new_df):
    rows = []
    indexed = new_df.set_index(["model", "feature_combo"])
    for feature_name, full_combo, baseline_combo in [
        ("ortholog", "ORT_EXP_SUB", "EXP_SUB"),
        ("expression", "ORT_EXP_SUB", "ORT_SUB"),
        ("subloc", "ORT_EXP_SUB", "ORT_EXP"),
    ]:
        deltas = []
        for model in MODELS:
            num = indexed.loc[(model, full_combo)]
            den = indexed.loc[(model, baseline_combo)]
            deltas.append(
                {
                    "feature_group": feature_name,
                    "model": model,
                    "delta_auprc": float(num["new_auprc_mean"] - den["new_auprc_mean"]),
                    "delta_auroc": float(num["new_auroc_mean"] - den["new_auroc_mean"]),
                    "delta_mcc": float(num["new_mcc_mean"] - den["new_mcc_mean"]),
                }
            )
        rows.extend(deltas)
    detail_df = pd.DataFrame(rows)
    summary_df = detail_df.groupby("feature_group", as_index=False).agg(
        mean_delta_auprc=("delta_auprc", "mean"),
        mean_delta_auroc=("delta_auroc", "mean"),
        mean_delta_mcc=("delta_mcc", "mean"),
        positive_auprc_models=("delta_auprc", lambda s: int((s > 0).sum())),
    ).sort_values(["mean_delta_auprc", "mean_delta_auroc", "mean_delta_mcc"], ascending=[False, False, False]).reset_index(drop=True)
    return detail_df, summary_df


def main():
    args = parse_args()
    configure_matplotlib()

    results_root = Path(args.results_root)
    figures_root = results_root / "figures"
    results_root.mkdir(parents=True, exist_ok=True)
    figures_root.mkdir(parents=True, exist_ok=True)

    old_df = _load_grid(args.old_root, "old")
    new_df = _load_grid(args.new_root, "new")
    merged = old_df.merge(new_df, on=["model", "feature_combo"], how="inner")
    expected = len(MODELS) * len(FEATURE_COMBOS)
    if len(merged) != expected:
        raise RuntimeError(f"Expected {expected} comparable rows, found {len(merged)}")

    merged["delta_auroc"] = merged["new_auroc_mean"] - merged["old_auroc_mean"]
    merged["delta_auprc"] = merged["new_auprc_mean"] - merged["old_auprc_mean"]
    merged["delta_mcc"] = merged["new_mcc_mean"] - merged["old_mcc_mean"]
    merged["delta_f1"] = merged["new_f1_mean"] - merged["old_f1_mean"]
    merged = merged.sort_values(["model", "feature_combo"]).reset_index(drop=True)
    merged.to_csv(results_root / "old440_vs_newlabel_comparison.tsv", sep="\t", index=False)

    new_final_summary = new_df.rename(columns={
        "new_auroc_mean": "auroc_mean",
        "new_auprc_mean": "auprc_mean",
        "new_mcc_mean": "mcc_mean",
        "new_f1_mean": "f1_mean",
        "new_auroc_std": "auroc_std",
        "new_auprc_std": "auprc_std",
        "new_mcc_std": "mcc_std",
        "new_f1_std": "f1_std",
    })
    new_final_summary.to_csv(results_root / "final_summary.tsv", sep="\t", index=False)

    best_old = merged.sort_values(["model", "old_auprc_mean", "old_auroc_mean"], ascending=[True, False, False]).groupby("model", as_index=False).first()
    best_new = merged.sort_values(["model", "new_auprc_mean", "new_auroc_mean"], ascending=[True, False, False]).groupby("model", as_index=False).first()
    new_model_mean = merged.groupby("model", as_index=False).agg(new_auprc_mean=("new_auprc_mean", "mean")).sort_values("new_auprc_mean", ascending=False)
    new_best = merged.sort_values(["new_auprc_mean", "new_auroc_mean", "new_mcc_mean"], ascending=[False, False, False]).iloc[0]
    old_best = merged.sort_values(["old_auprc_mean", "old_auroc_mean", "old_mcc_mean"], ascending=[False, False, False]).iloc[0]
    feature_detail, feature_summary = _compute_newlabel_feature_effects(merged)
    feature_detail.to_csv(results_root / "newlabel_feature_main_effects_by_model.tsv", sep="\t", index=False)
    feature_summary.to_csv(results_root / "newlabel_feature_main_effects.tsv", sep="\t", index=False)

    _make_heatmap(merged.assign(auprc_mean=merged["old_auprc_mean"]), "old_auprc_mean", "old440 benchmark AUPRC", figures_root / "heatmap_old440_auprc.pdf")
    _make_heatmap(merged.assign(auprc_mean=merged["new_auprc_mean"]), "new_auprc_mean", "new-label benchmark AUPRC", figures_root / "heatmap_newlabel_auprc.pdf")
    _make_heatmap(merged.assign(delta_auprc=merged["delta_auprc"]), "delta_auprc", "delta AUPRC (new - old)", figures_root / "delta_heatmap_auprc.pdf", cmap="coolwarm")
    _make_heatmap(merged.assign(delta_mcc=merged["delta_mcc"]), "delta_mcc", "delta MCC (new - old)", figures_root / "delta_heatmap_mcc.pdf", cmap="coolwarm")
    _make_best_barplot(best_old, best_new, figures_root / "best_feature_old_vs_new_barplot.pdf")

    overall_delta_auprc = float(merged["delta_auprc"].mean())
    overall_delta_mcc = float(merged["delta_mcc"].mean())
    top_feature = feature_summary.iloc[0]
    lines = [
        "# Phase 2B New-Label vs old440 Summary",
        "",
        "## Key Answers",
        "",
        "1. new label 是否优于 old440",
        "{}。Across all 28 matched configurations, the mean delta is AUPRC={:.4f} and MCC={:.4f}.".format("是" if overall_delta_auprc > 0 else "否", overall_delta_auprc, overall_delta_mcc),
        "2. 哪个 model 在 new label 下最强",
        "`{}` is the strongest model under the new label regime by mean AUPRC across all feature combinations ({:.4f}).".format(str(new_model_mean.iloc[0]["model"]), float(new_model_mean.iloc[0]["new_auprc_mean"])),
        "3. 哪个 feature combo 在 new label 下最强",
        "The best single configuration under the new label regime is `{}` with `{}` (AUPRC={:.4f}, AUROC={:.4f}, MCC={:.4f}).".format(str(new_best["model"]), FEATURE_LABELS[str(new_best["feature_combo"])], float(new_best["new_auprc_mean"]), float(new_best["new_auroc_mean"]), float(new_best["new_mcc_mean"])),
        "4. ortholog / expression / subloc 哪类特征在新 label 下增益最大",
        "`{}` shows the largest full-ablation gain under the new label regime (mean delta AUPRC={:.4f}, mean delta MCC={:.4f}).".format(str(top_feature["feature_group"]), float(top_feature["mean_delta_auprc"]), float(top_feature["mean_delta_mcc"])),
        "5. 如果 ortholog 或 expression 增益明显，是否说明后续优化这些特征值得做",
        "{}".format("值得优先做 ortholog/expression 方向优化。" if str(top_feature["feature_group"]) in ["ortholog", "expression"] and float(top_feature["mean_delta_auprc"]) > 0 else "当前结果不支持把 ortholog/expression 作为唯一优先优化方向。"),
        "",
        "## Additional Evidence",
        "",
        "- old440 best configuration: `{}` + `{}` with AUPRC={:.4f}.".format(str(old_best["model"]), FEATURE_LABELS[str(old_best["feature_combo"])], float(old_best["old_auprc_mean"])),
        "- new-label best configuration: `{}` + `{}` with AUPRC={:.4f}.".format(str(new_best["model"]), FEATURE_LABELS[str(new_best["feature_combo"])], float(new_best["new_auprc_mean"])),
        "- Comparison table: `results/phase2b_new_label/old440_vs_newlabel_comparison.tsv`.",
    ]
    (results_root / "phase2b_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
