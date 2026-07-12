import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.analysis.figure5_representation_common import (
    LABEL_COLORS,
    load_hidden_case,
    pair_cases,
    save_png,
    save_pdf,
    separation_metrics,
    species_title,
    write_json,
    write_markdown,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build Figure5B single-species hidden-space quantitative summary")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--model", default="GraphSAGE", type=str)
    parser.add_argument("--seed", default=1029, type=int)
    parser.add_argument("--subset", default="test", type=str)
    parser.add_argument("--plot-dir", default="results/Figure5/plots", type=str)
    parser.add_argument("--data-dir", default="results/Figure5/data", type=str)
    parser.add_argument("--table-dir", default="results/Figure5/tables", type=str)
    parser.add_argument("--summary-dir", default="results/Figure5/summary", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    plot_dir = Path(args.plot_dir).resolve()
    data_dir = Path(args.data_dir).resolve()
    table_dir = Path(args.table_dir).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    for path in [plot_dir, data_dir, table_dir, summary_dir]:
        path.mkdir(parents=True, exist_ok=True)

    baseline_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB", args.model, args.seed, args.subset)
    esm2_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB_ESM2", args.model, args.seed, args.subset)
    paired = pair_cases(baseline_case, esm2_case)
    species = esm2_case["species"]
    stem = f"Figure5b_hidden_quant_summary_{species}_seed{args.seed}"

    raw_rows = []
    for case in [baseline_case, esm2_case]:
        metrics = separation_metrics(case["hidden_matrix"], case["labels"])
        for metric_name, metric_value in metrics.items():
            raw_rows.append(
                {
                    "species": case["species"],
                    "protocol": case["protocol"],
                    "feature_setting": case["feature_setting"],
                    "summary_type": "separation",
                    "item": metric_name,
                    "value": metric_value,
                }
            )

    labels = paired["label"].astype(int)
    baseline_pred = paired["baseline_pred_label"].astype(int)
    esm2_pred = paired["esm2_pred_label"].astype(int)
    raw_rows.extend(
        [
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "total_TP", "value": int(((labels == 1) & (esm2_pred == 1)).sum())},
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "total_TN", "value": int(((labels == 0) & (esm2_pred == 0)).sum())},
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "FN_to_TP_rescued", "value": int((paired["transition"] == "FN_to_TP_rescued").sum())},
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "FP_to_TN_corrected", "value": int((paired["transition"] == "FP_to_TN_corrected").sum())},
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "persistent_FN", "value": int((paired["transition"] == "FN_persistent").sum())},
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "persistent_FP", "value": int((paired["transition"] == "FP_persistent").sum())},
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB", "summary_type": "error_baseline", "item": "total_TP", "value": int(((labels == 1) & (baseline_pred == 1)).sum())},
            {"species": species, "protocol": args.protocol, "feature_setting": "ORT_EXP_SUB", "summary_type": "error_baseline", "item": "total_TN", "value": int(((labels == 0) & (baseline_pred == 0)).sum())},
        ]
    )

    raw_df = pd.DataFrame(raw_rows)
    raw_path = data_dir / f"{stem}.tsv"
    stats_path = table_dir / f"{stem}_stats.tsv"
    pdf_path = plot_dir / f"{stem}.pdf"
    png_path = plot_dir / f"{stem}.png"
    summary_path = summary_dir / f"{stem}.md"
    json_path = table_dir / f"{stem}.json"
    raw_df.to_csv(raw_path, sep="\t", index=False)

    stats_df = raw_df.copy()
    stats_df["n"] = 1
    stats_df["mean"] = stats_df["value"].astype(float)
    stats_df["std"] = 0.0
    stats_df["variance"] = 0.0
    stats_df.to_csv(stats_path, sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), facecolor="white")
    separation_df = raw_df[raw_df["summary_type"] == "separation"].copy()
    metric_order = ["centroid_distance", "silhouette_score", "davies_bouldin_index"]
    x = np.arange(len(metric_order))
    width = 0.36
    for idx, feature_setting in enumerate(["ORT_EXP_SUB", "ORT_EXP_SUB_ESM2"]):
        subset = separation_df[separation_df["feature_setting"] == feature_setting].set_index("item").reindex(metric_order)
        axes[0].bar(x + (idx - 0.5) * width, subset["value"], width=width, label=feature_setting, color="#9AA5B1" if feature_setting == "ORT_EXP_SUB" else "#D62728")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(metric_order, rotation=20, ha="right")
    axes[0].set_title("{0}\nHidden separation summary".format(species_title(args.protocol)), fontsize=10)
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_facecolor("white")

    error_df = raw_df[raw_df["summary_type"] == "error"].copy()
    error_order = ["total_TP", "total_TN", "FN_to_TP_rescued", "FP_to_TN_corrected", "persistent_FN", "persistent_FP"]
    colors = [LABEL_COLORS[1], LABEL_COLORS[0], "#D62728", "#9467BD", "#FF9896", "#C5B0D5"]
    axes[1].bar(np.arange(len(error_order)), error_df.set_index("item").reindex(error_order)["value"], color=colors)
    axes[1].set_xticks(np.arange(len(error_order)))
    axes[1].set_xticklabels(error_order, rotation=25, ha="right")
    axes[1].set_title("{0}\nError rescue summary".format(species_title(args.protocol)), fontsize=10)
    axes[1].set_facecolor("white")
    for ax in axes:
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
    fig.tight_layout()
    save_png(fig, png_path)
    save_pdf(fig, pdf_path)

    rescued = int(error_df[error_df["item"] == "FN_to_TP_rescued"]["value"].iloc[0])
    write_markdown(
        summary_path,
        [
            "# Figure5b Hidden-Space Quantitative Summary",
            "",
            f"- Species: `{species}`.",
            f"- Plot PDF: `{pdf_path}`.",
            f"- Plot PNG: `{png_path}`.",
            f"- Raw plot data: `{raw_path}`.",
            f"- Statistics table: `{stats_path}`.",
            f"- Rescued false negatives: `{rescued}`.",
            "- `n=1`, so `mean=value`, `std=0`, and `variance=0` for this mainline seed-specific table.",
        ],
    )

    write_json(
        json_path,
        {
            "panel": "Figure5b",
            "protocol": args.protocol,
            "species": species,
            "model": args.model,
            "seed": args.seed,
            "subset": args.subset,
            "raw_tsv": str(raw_path),
            "stats_tsv": str(stats_path),
            "plot_pdf": str(pdf_path),
            "plot_png": str(png_path),
            "summary_md": str(summary_path),
        },
    )


if __name__ == "__main__":
    main()
