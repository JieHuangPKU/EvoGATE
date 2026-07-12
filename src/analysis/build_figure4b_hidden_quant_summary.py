import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.analysis.figure4_representation_common import (
    LABEL_COLORS,
    load_hidden_case,
    pair_cases,
    save_pdf,
    separation_metrics,
    species_title,
    write_json,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build Figure4B single-species hidden quantitative summary")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--model", default="GraphSAGE", type=str)
    parser.add_argument("--seed", default=1029, type=int)
    parser.add_argument("--subset", default="test", type=str)
    parser.add_argument("--output-dir", required=True, type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    baseline_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB", args.model, args.seed, args.subset)
    esm2_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB_ESM2", args.model, args.seed, args.subset)
    paired = pair_cases(baseline_case, esm2_case)

    separation_rows = []
    for case in [baseline_case, esm2_case]:
        metrics = separation_metrics(case["hidden_matrix"], case["labels"])
        separation_rows.extend(
            [
                {"species": case["species"], "feature_setting": case["feature_setting"], "summary_type": "separation", "item": "centroid_distance", "value": metrics["centroid_distance"]},
                {"species": case["species"], "feature_setting": case["feature_setting"], "summary_type": "separation", "item": "silhouette_score", "value": metrics["silhouette_score"]},
                {"species": case["species"], "feature_setting": case["feature_setting"], "summary_type": "separation", "item": "davies_bouldin_index", "value": metrics["davies_bouldin_index"]},
            ]
        )

    labels = paired["label"].astype(int)
    baseline_pred = paired["baseline_pred_label"].astype(int)
    esm2_pred = paired["esm2_pred_label"].astype(int)
    error_rows = [
        {"species": esm2_case["species"], "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "total_TP", "value": int(((labels == 1) & (esm2_pred == 1)).sum())},
        {"species": esm2_case["species"], "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "total_TN", "value": int(((labels == 0) & (esm2_pred == 0)).sum())},
        {"species": esm2_case["species"], "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "FN_to_TP_rescued", "value": int((paired["transition"] == "FN_to_TP_rescued").sum())},
        {"species": esm2_case["species"], "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "FP_to_TN_corrected", "value": int((paired["transition"] == "FP_to_TN_corrected").sum())},
        {"species": esm2_case["species"], "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "persistent_FN", "value": int((paired["transition"] == "FN_persistent").sum())},
        {"species": esm2_case["species"], "feature_setting": "ORT_EXP_SUB_ESM2", "summary_type": "error", "item": "persistent_FP", "value": int((paired["transition"] == "FP_persistent").sum())},
        {"species": baseline_case["species"], "feature_setting": "ORT_EXP_SUB", "summary_type": "error_baseline", "item": "total_TP", "value": int(((labels == 1) & (baseline_pred == 1)).sum())},
        {"species": baseline_case["species"], "feature_setting": "ORT_EXP_SUB", "summary_type": "error_baseline", "item": "total_TN", "value": int(((labels == 0) & (baseline_pred == 0)).sum())},
    ]
    summary_df = pd.DataFrame(separation_rows + error_rows)
    tsv_path = output_dir / "Figure4B_{0}_hidden_quant_summary_seed{1}.tsv".format(esm2_case["species"], args.seed)
    pdf_path = output_dir / "Figure4B_{0}_hidden_quant_summary_seed{1}.pdf".format(esm2_case["species"], args.seed)
    summary_df.to_csv(tsv_path, sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), facecolor="white")
    separation_df = summary_df[summary_df["summary_type"] == "separation"].copy()
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

    error_df = summary_df[summary_df["summary_type"] == "error"].copy()
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
    save_pdf(fig, pdf_path)

    write_json(
        output_dir / "Figure4B_{0}_hidden_quant_summary_seed{1}.json".format(esm2_case["species"], args.seed),
        {
            "protocol": args.protocol,
            "species": esm2_case["species"],
            "model": args.model,
            "seed": args.seed,
            "subset": args.subset,
            "baseline_checkpoint_path": baseline_case["checkpoint_path"],
            "esm2_checkpoint_path": esm2_case["checkpoint_path"],
            "split_manifest_path": esm2_case["split_manifest_path"],
            "label_manifest_path": esm2_case["label_manifest_path"],
            "summary_tsv": str(tsv_path),
            "pdf_path": str(pdf_path),
        },
    )


if __name__ == "__main__":
    main()
