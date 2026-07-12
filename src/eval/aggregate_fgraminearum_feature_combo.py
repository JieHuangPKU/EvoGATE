"""
Aggregate and analyze fixed-config F. graminearum feature-combination benchmarks.
"""

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SPECIES = "fgraminearum"
MODELS = ["GAT", "GCN", "GIN", "GraphSAGE"]
FEATURE_COMBOS = ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB"]
METRIC_COLUMNS = ["auroc", "auprc", "mcc", "f1", "accuracy"]
HEATMAP_METRICS = {
    "auroc_mean": ("AUROC", "fgraminearum_auroc_heatmap.png", "YlGn"),
    "auprc_mean": ("AUPRC", "fgraminearum_auprc_heatmap.png", "YlOrBr"),
    "mcc_mean": ("MCC", "fgraminearum_mcc_heatmap.png", "PuBu"),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate F. graminearum feature-combination benchmark outputs")
    parser.add_argument("--run-root", type=str)
    parser.add_argument("--aggregated-output", type=str)
    parser.add_argument("--output-root", type=str)
    parser.add_argument("--final-summary-output", type=str)
    parser.add_argument("--results-root", type=str)
    return parser.parse_args()


def _contains_feature(combo_name, feature_name):
    tokens = set(str(combo_name).split("_"))
    return feature_name in tokens


def _as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ["true", "1", "yes"]


def _load_run_metrics(run_dir):
    metrics_path = os.path.join(run_dir, "metrics.tsv")
    feature_path = os.path.join(run_dir, "feature_summary.tsv")
    if not os.path.exists(metrics_path):
        raise FileNotFoundError("Missing metrics.tsv: {}".format(metrics_path))
    if not os.path.exists(feature_path):
        raise FileNotFoundError("Missing feature_summary.tsv: {}".format(feature_path))
    metric_row = pd.read_csv(metrics_path, sep="\t").iloc[0].to_dict()
    feature_row = pd.read_csv(feature_path, sep="\t").iloc[0].to_dict()
    return metric_row, feature_row


def _aggregate_one(run_root, aggregated_output):
    run_root = os.path.abspath(run_root)
    aggregated_output = os.path.abspath(aggregated_output)
    run_dirs = []
    for name in sorted(os.listdir(run_root)):
        path = os.path.join(run_root, name)
        if os.path.isdir(path) and name.startswith("run_"):
            run_dirs.append(path)
    if not run_dirs:
        raise RuntimeError("No run directories found under {}".format(run_root))

    metric_rows = []
    reference = None
    for run_dir in run_dirs:
        metric_row, feature_row = _load_run_metrics(run_dir)
        metric_row["run_id"] = os.path.basename(run_dir)
        metric_rows.append(metric_row)
        if reference is None:
            reference = {
                "species": str(feature_row["species"]),
                "model": str(feature_row["model"]),
                "feature_combo": str(feature_row["feature_combo"]),
                "string_thr": int(feature_row["string_thr"]),
                "include_degree": _as_bool(feature_row["include_degree"]),
                "orthologs_enabled": _as_bool(feature_row.get("orthologs_enabled", _contains_feature(feature_row["feature_combo"], "ORT"))),
                "expression_enabled": _as_bool(feature_row.get("expression_enabled", _contains_feature(feature_row["feature_combo"], "EXP"))),
                "sublocalization_enabled": _as_bool(feature_row.get("sublocalization_enabled", _contains_feature(feature_row["feature_combo"], "SUB"))),
            }

    metric_df = pd.DataFrame(metric_rows)
    out = dict(reference)
    out["runs_completed"] = int(len(metric_df))
    for metric in METRIC_COLUMNS:
        metric_prefix = "acc" if metric == "accuracy" else metric
        out[metric_prefix + "_mean"] = float(metric_df[metric].mean())
        out[metric_prefix + "_std"] = float(metric_df[metric].std(ddof=0))
    Path(os.path.dirname(aggregated_output)).mkdir(parents=True, exist_ok=True)
    pd.DataFrame([out]).to_csv(aggregated_output, sep="\t", index=False)


def _collect_aggregated_rows(output_root):
    rows = []
    missing = []
    for model in MODELS:
        for feature_combo in FEATURE_COMBOS:
            agg_path = os.path.join(output_root, model, feature_combo, "thr_300", "aggregated_metrics.tsv")
            if not os.path.exists(agg_path):
                missing.append(agg_path)
                continue
            row = pd.read_csv(agg_path, sep="\t").iloc[0].to_dict()
            row["include_degree"] = _as_bool(row.get("include_degree", False))
            row["orthologs_enabled"] = _as_bool(row.get("orthologs_enabled", _contains_feature(row["feature_combo"], "ORT")))
            row["expression_enabled"] = _as_bool(row.get("expression_enabled", _contains_feature(row["feature_combo"], "EXP")))
            row["sublocalization_enabled"] = _as_bool(row.get("sublocalization_enabled", _contains_feature(row["feature_combo"], "SUB")))
            rows.append(row)
    summary = pd.DataFrame(rows)
    if summary.empty:
        summary = pd.DataFrame(
            columns=[
                "species",
                "model",
                "feature_combo",
                "string_thr",
                "include_degree",
                "orthologs_enabled",
                "expression_enabled",
                "sublocalization_enabled",
                "runs_completed",
                "auroc_mean",
                "auroc_std",
                "auprc_mean",
                "auprc_std",
                "mcc_mean",
                "mcc_std",
                "f1_mean",
                "f1_std",
                "acc_mean",
                "acc_std",
            ]
        )
    else:
        summary["model_rank"] = summary["model"].map({name: idx for idx, name in enumerate(MODELS)})
        summary["combo_rank"] = summary["feature_combo"].map({name: idx for idx, name in enumerate(FEATURE_COMBOS)})
        summary = summary.sort_values(["model_rank", "combo_rank"]).drop(columns=["model_rank", "combo_rank"]).reset_index(drop=True)
    return summary, missing


def _plot_heatmap(summary_df, metric_col, title, output_path, cmap):
    matrix = summary_df.pivot(index="feature_combo", columns="model", values=metric_col).reindex(index=FEATURE_COMBOS, columns=MODELS)
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    image = ax.imshow(matrix.values, aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(MODELS)))
    ax.set_xticklabels(MODELS)
    ax.set_yticks(range(len(FEATURE_COMBOS)))
    ax.set_yticklabels(FEATURE_COMBOS)
    ax.set_title(title)
    numeric_values = matrix.values[~pd.isna(matrix.values)]
    mean_value = float(numeric_values.mean()) if len(numeric_values) else 0.0
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            if pd.isna(value):
                text = "NA"
                text_color = "black"
            else:
                text = "{:.3f}".format(float(value))
                text_color = "white" if float(value) < mean_value else "black"
            ax.text(j, i, text, ha="center", va="center", color=text_color, fontsize=9)
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.ax.set_ylabel(metric_col, rotation=270, labelpad=16)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _compute_feature_main_effects(summary_df):
    rows = []
    for feature_name in ["ORT", "EXP", "SUB"]:
        present_mask = summary_df["feature_combo"].astype(str).map(lambda combo: _contains_feature(combo, feature_name))
        present_df = summary_df[present_mask].copy()
        absent_df = summary_df[~present_mask].copy()
        row = {
            "feature_group": feature_name,
            "n_with_feature": int(len(present_df)),
            "n_without_feature": int(len(absent_df)),
        }
        for metric_prefix in ["auroc", "auprc", "mcc"]:
            with_mean = float(present_df[metric_prefix + "_mean"].mean())
            without_mean = float(absent_df[metric_prefix + "_mean"].mean())
            row[metric_prefix + "_with_feature_mean"] = with_mean
            row[metric_prefix + "_without_feature_mean"] = without_mean
            row[metric_prefix + "_average_gain"] = with_mean - without_mean
        positive_model_count = 0
        per_model_gains = []
        for model in MODELS:
            model_df = summary_df[summary_df["model"] == model].copy()
            model_present = model_df[model_df["feature_combo"].astype(str).map(lambda combo: _contains_feature(combo, feature_name))]
            model_absent = model_df[~model_df["feature_combo"].astype(str).map(lambda combo: _contains_feature(combo, feature_name))]
            model_gain = float(model_present["auroc_mean"].mean()) - float(model_absent["auroc_mean"].mean())
            per_model_gains.append(model_gain)
            if model_gain > 0:
                positive_model_count += 1
        row["auroc_positive_models"] = int(positive_model_count)
        row["auroc_positive_model_fraction"] = float(positive_model_count / len(MODELS))
        row["auroc_model_gain_std"] = float(pd.Series(per_model_gains).std(ddof=0))
        rows.append(row)
    effect_df = pd.DataFrame(rows).sort_values(["auroc_average_gain", "auprc_average_gain", "feature_group"], ascending=[False, False, True]).reset_index(drop=True)
    return effect_df


def _build_report(summary_df, effect_df):
    if summary_df.empty:
        return "\n".join(
            [
                "# F. graminearum Feature Combo Benchmark Summary",
                "",
                "No aggregated benchmark rows were found, so the benchmark conclusions are not available.",
            ]
        )

    best_row = summary_df.sort_values(["auroc_mean", "auprc_mean", "mcc_mean", "auroc_std"], ascending=[False, False, False, True]).iloc[0]
    best_model_summary = (
        summary_df.groupby("model", as_index=False)
        .agg(
            best_auroc=("auroc_mean", "max"),
            mean_auroc=("auroc_mean", "mean"),
            mean_auprc=("auprc_mean", "mean"),
            mean_mcc=("mcc_mean", "mean"),
        )
        .sort_values(["best_auroc", "mean_auroc", "mean_auprc", "mean_mcc"], ascending=[False, False, False, False])
        .iloc[0]
    )
    top_effect = effect_df.sort_values(["auroc_average_gain", "auprc_average_gain", "mcc_average_gain"], ascending=[False, False, False]).iloc[0]
    exp_row = effect_df[effect_df["feature_group"] == "EXP"].iloc[0]
    ort_row = effect_df[effect_df["feature_group"] == "ORT"].iloc[0]

    exp_stable = bool(exp_row["auroc_average_gain"] > 0 and exp_row["auprc_average_gain"] > 0 and exp_row["auroc_positive_models"] >= 3)
    ort_dominant = bool(top_effect["feature_group"] == "ORT" and ort_row["auroc_average_gain"] > 0 and ort_row["auprc_average_gain"] > 0)

    optimize_ortholog = "值得继续优化" if ort_row["auroc_average_gain"] > 0 or ort_row["auprc_average_gain"] > 0 else "暂不建议优先优化"
    optimize_expression = "值得继续优化" if exp_stable else "目前不建议优先作为主优化方向"

    lines = [
        "# F. graminearum Feature Combo Benchmark Summary",
        "",
        "## Key Answers",
        "",
        "- Fusarium 上整体最优图模型是 `{}`，其最佳组合表现出现在 `{}` 上（AUROC={:.4f}, AUPRC={:.4f}, MCC={:.4f}）。".format(
            str(best_model_summary["model"]),
            str(best_row["feature_combo"]),
            float(best_row["auroc_mean"]),
            float(best_row["auprc_mean"]),
            float(best_row["mcc_mean"]),
        ),
        "- Fusarium 上最优 feature combo 是 `{}`，对应模型为 `{}`。".format(str(best_row["feature_combo"]), str(best_row["model"])),
        "- 组合级主效应分析显示贡献最大的特征类型是 `{}`（AUROC 平均增益={:.4f}, AUPRC 平均增益={:.4f}, MCC 平均增益={:.4f}）。".format(
            str(top_effect["feature_group"]),
            float(top_effect["auroc_average_gain"]),
            float(top_effect["auprc_average_gain"]),
            float(top_effect["mcc_average_gain"]),
        ),
        "- expression 是否提供稳定增益: {}。依据是 EXP 的 AUROC/AUPRC 平均增益分别为 {:.4f}/{:.4f}，且在 {}/{} 个模型上带来 AUROC 正增益。".format(
            "是" if exp_stable else "否",
            float(exp_row["auroc_average_gain"]),
            float(exp_row["auprc_average_gain"]),
            int(exp_row["auroc_positive_models"]),
            len(MODELS),
        ),
        "- ortholog 是否是主导信号: {}。依据是 ORT 的 AUROC/AUPRC/MCC 平均增益为 {:.4f}/{:.4f}/{:.4f}。".format(
            "是" if ort_dominant else "否",
            float(ort_row["auroc_average_gain"]),
            float(ort_row["auprc_average_gain"]),
            float(ort_row["mcc_average_gain"]),
        ),
        "- 是否值得进一步优化 ortholog 或 expression: ortholog {}, expression {}。".format(optimize_ortholog, optimize_expression),
        "",
        "## Notes",
        "",
        "- 结论基于固定配置 benchmark: species=`fgraminearum`, string threshold=`300`, include_degree=`false`, 4 个图模型 x 7 个 feature combos x 3 runs。",
        "- feature contribution analysis 为组合级主效应估计，不是 node-level explainer。",
    ]
    return "\n".join(lines) + "\n"


def _analyze_all(output_root, final_summary_output, results_root):
    output_root = os.path.abspath(output_root)
    final_summary_output = os.path.abspath(final_summary_output)
    results_root = Path(results_root).resolve()
    results_root.mkdir(parents=True, exist_ok=True)

    summary_df, missing = _collect_aggregated_rows(output_root)
    summary_df.to_csv(final_summary_output, sep="\t", index=False)
    summary_df.to_csv(results_root / "model_feature_combo_summary.tsv", sep="\t", index=False)

    if summary_df.empty:
        effect_df = pd.DataFrame(
            columns=[
                "feature_group",
                "n_with_feature",
                "n_without_feature",
                "auroc_with_feature_mean",
                "auroc_without_feature_mean",
                "auroc_average_gain",
                "auprc_with_feature_mean",
                "auprc_without_feature_mean",
                "auprc_average_gain",
                "mcc_with_feature_mean",
                "mcc_without_feature_mean",
                "mcc_average_gain",
                "auroc_positive_models",
                "auroc_positive_model_fraction",
                "auroc_model_gain_std",
            ]
        )
    else:
        effect_df = _compute_feature_main_effects(summary_df)
        for metric_col, meta in HEATMAP_METRICS.items():
            title, filename, cmap = meta
            _plot_heatmap(summary_df, metric_col, "F. graminearum {} by model x feature combo".format(title), results_root / filename, cmap)

    effect_df.to_csv(results_root / "feature_main_effects.tsv", sep="\t", index=False)

    report_text = _build_report(summary_df, effect_df)
    if missing:
        report_text += "\n## Missing Aggregates\n\n"
        for path in missing:
            report_text += "- `{}`\n".format(path)
    (results_root / "fgraminearum_feature_combo_summary.md").write_text(report_text, encoding="utf-8")


def main():
    args = parse_args()
    if args.run_root and args.aggregated_output:
        _aggregate_one(args.run_root, args.aggregated_output)
        return
    if args.output_root and args.final_summary_output and args.results_root:
        _analyze_all(args.output_root, args.final_summary_output, args.results_root)
        return
    raise SystemExit(
        "Provide either --run-root/--aggregated-output or --output-root/--final-summary-output/--results-root"
    )


if __name__ == "__main__":
    main()
