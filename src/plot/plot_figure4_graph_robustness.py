import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


COLORS = {
    "STRING_300": "#1f77b4",
    "STRING_700": "#2ca02c",
    "eFG": "#d95f02",
    "eFG_HIGH": "#6a3d9a",
    "eFG_HIGH_MEDIUM": "#ff7f0e",
    "eFG_ALL": "#8c564b",
    "AUPRC": "#1f77b4",
    "MCC": "#d95f02",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Figure 4 Graph Robustness panels and export redraw TSVs")
    parser.add_argument("--summary-dir", required=True, type=str)
    parser.add_argument("--output-dir", required=True, type=str)
    return parser.parse_args()


def style_axes(ax) -> None:
    ax.set_facecolor("white")
    for spine_name in ["left", "bottom"]:
        ax.spines[spine_name].set_color("black")
        ax.spines[spine_name].set_linewidth(1.0)
    ax.tick_params(axis="both", colors="black", width=0.8)


def save_pdf_png(fig, output_stub: Path) -> None:
    fig.savefig(output_stub.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stub.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    summary_dir = Path(args.summary_dir)
    output_dir = Path(args.output_dir)
    data_dir = output_dir / "data"
    supplementary_dir = output_dir / "supplementary"
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    supplementary_dir.mkdir(parents=True, exist_ok=True)

    threshold_per_run_df = pd.read_csv(summary_dir / "Figure4_threshold_per_run_metrics.tsv", sep="\t")
    threshold_df = pd.read_csv(summary_dir / "Figure4_threshold_aggregated_metrics.tsv", sep="\t")
    density_df = pd.read_csv(summary_dir / "Figure4_network_density_summary.tsv", sep="\t")
    source_per_run_df = pd.read_csv(summary_dir / "Figure4_source_comparison_per_run_metrics.tsv", sep="\t")
    source_df = pd.read_csv(summary_dir / "Figure4_source_comparison_metrics.tsv", sep="\t")
    overlap_df = pd.read_csv(summary_dir / "Figure4_edge_overlap_summary.tsv", sep="\t").iloc[0]
    supplementary_per_run_df = pd.read_csv(summary_dir / "Figure4_supplementary_source_comparison_per_run_metrics.tsv", sep="\t")
    supplementary_summary_df = pd.read_csv(summary_dir / "Figure4_supplementary_source_comparison_metrics.tsv", sep="\t")

    feature_setting = str(threshold_df["feature_setting"].iloc[0])
    model = str(threshold_df["model"].iloc[0])
    model_feature = f"{model} + {feature_setting}"

    # Figure4A data exports
    figure4a_seed_df = threshold_per_run_df[threshold_per_run_df["graph_family"].astype(str).eq("STRING")].copy()
    figure4a_seed_df = figure4a_seed_df[
        [
            "protocol",
            "species",
            "model",
            "feature_setting",
            "model_feature_setting",
            "graph_condition",
            "string_threshold",
            "seed",
            "test_auprc",
            "test_mcc",
            "val_auprc",
            "val_mcc",
            "test_auroc",
            "test_accuracy",
            "metrics_path",
        ]
    ].sort_values(["string_threshold", "seed"]).reset_index(drop=True)
    figure4a_seed_df.to_csv(data_dir / "Figure4A_threshold_performance_line.tsv", sep="\t", index=False)
    figure4a_summary_df = threshold_df.copy().sort_values("string_threshold").reset_index(drop=True)
    figure4a_summary_df.to_csv(data_dir / "Figure4A_threshold_performance_line_summary.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.2), facecolor="white")
    x = figure4a_summary_df["string_threshold"].astype(int)
    ax.errorbar(x, figure4a_summary_df["test_auprc_mean"], yerr=figure4a_summary_df["test_auprc_sd"], marker="o", linewidth=2, capsize=3, color=COLORS["AUPRC"], label="Test AUPRC")
    ax.errorbar(x, figure4a_summary_df["test_mcc_mean"], yerr=figure4a_summary_df["test_mcc_sd"], marker="s", linewidth=2, capsize=3, color=COLORS["MCC"], label="Test MCC")
    ax.set_xlabel("STRING threshold")
    ax.set_ylabel("Performance")
    ax.set_xticks(x)
    ax.set_title(f"Figure4A. {model_feature} performance vs STRING threshold")
    style_axes(ax)
    ax.legend(frameon=False, loc="best")
    save_pdf_png(fig, output_dir / "Figure4A_threshold_performance_line")

    # Figure4B data exports
    figure4b_df = density_df.copy().sort_values("string_threshold").reset_index(drop=True)
    figure4b_df.to_csv(data_dir / "Figure4B_network_density_plot.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.0), sharex=True, facecolor="white")
    axes[0].plot(figure4b_df["string_threshold"], figure4b_df["edge_count"], marker="o", color="#1f77b4", linewidth=2)
    axes[0].set_ylabel("Edge count")
    axes[1].plot(figure4b_df["string_threshold"], figure4b_df["graph_node_count"], marker="o", color="#2ca02c", linewidth=2)
    axes[1].set_ylabel("Node count")
    axes[2].plot(figure4b_df["string_threshold"], figure4b_df["isolated_node_count"], marker="o", color="#d95f02", linewidth=2)
    axes[2].set_ylabel("Isolated nodes")
    axes[2].set_xlabel("STRING threshold")
    axes[0].set_title(f"Figure4B. {model_feature} graph density vs STRING threshold")
    for ax in axes:
        style_axes(ax)
    save_pdf_png(fig, output_dir / "Figure4B_network_density_plot")

    # Figure4C main data exports
    figure4c_seed_df = source_per_run_df[
        [
            "protocol",
            "species",
            "model",
            "feature_setting",
            "model_feature_setting",
            "graph_condition",
            "seed",
            "test_auprc",
            "test_mcc",
            "val_auprc",
            "val_mcc",
            "metrics_path",
        ]
    ].sort_values(["graph_condition", "seed"]).reset_index(drop=True)
    figure4c_seed_df.to_csv(data_dir / "Figure4C_source_comparison_barplot.tsv", sep="\t", index=False)
    figure4c_summary_df = source_df.copy().reset_index(drop=True)
    figure4c_summary_df.to_csv(data_dir / "Figure4C_source_comparison_barplot_summary.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.4), sharey=False, facecolor="white")
    conditions = figure4c_summary_df["graph_condition"].tolist()
    xpos = range(len(conditions))
    axes[0].bar(
        list(xpos),
        figure4c_summary_df["test_auprc_mean"],
        yerr=figure4c_summary_df["test_auprc_sd"],
        capsize=4,
        color=[COLORS.get(condition, "#4c4c4c") for condition in conditions],
    )
    axes[0].set_xticks(list(xpos), conditions, rotation=20)
    axes[0].set_ylabel("Test AUPRC")
    axes[0].set_title("AUPRC")
    axes[1].bar(
        list(xpos),
        figure4c_summary_df["test_mcc_mean"],
        yerr=figure4c_summary_df["test_mcc_sd"],
        capsize=4,
        color=[COLORS.get(condition, "#4c4c4c") for condition in conditions],
    )
    axes[1].set_xticks(list(xpos), conditions, rotation=20)
    axes[1].set_ylabel("Test MCC")
    axes[1].set_title("MCC")
    fig.suptitle(f"Figure4C. {model_feature} source comparison", y=0.98)
    for ax in axes:
        style_axes(ax)
    fig.subplots_adjust(top=0.8, wspace=0.35)
    save_pdf_png(fig, output_dir / "Figure4C_source_comparison_barplot")

    # Figure4D data exports
    figure4d_rows = [
        {"plot_section": "total_edge_counts", "condition": "STRING_300", "category": "total", "value": int(overlap_df["string_300_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "total_edge_counts", "condition": "STRING_700", "category": "total", "value": int(overlap_df["string_700_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "total_edge_counts", "condition": "eFG", "category": "total", "value": int(overlap_df["efg_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "pairwise_overlap", "condition": "STRING_300 vs eFG", "category": "shared", "value": int(overlap_df["string_300_intersect_efg_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "pairwise_overlap", "condition": "STRING_300 vs eFG", "category": "STRING_only", "value": int(overlap_df["string_300_only_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "pairwise_overlap", "condition": "STRING_300 vs eFG", "category": "eFG_only", "value": int(overlap_df["efg_only_vs_string_300_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "pairwise_overlap", "condition": "STRING_700 vs eFG", "category": "shared", "value": int(overlap_df["string_700_intersect_efg_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "pairwise_overlap", "condition": "STRING_700 vs eFG", "category": "STRING_only", "value": int(overlap_df["string_700_only_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
        {"plot_section": "pairwise_overlap", "condition": "STRING_700 vs eFG", "category": "eFG_only", "value": int(overlap_df["efg_only_vs_string_700_edge_count"]), "model": model, "feature_setting": feature_setting, "model_feature_setting": model_feature},
    ]
    figure4d_df = pd.DataFrame(figure4d_rows)
    figure4d_df.to_csv(data_dir / "Figure4D_edge_overlap_plot.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4), facecolor="white")
    totals = pd.Series({"STRING_300": overlap_df["string_300_edge_count"], "STRING_700": overlap_df["string_700_edge_count"], "eFG": overlap_df["efg_edge_count"]})
    axes[0].bar(totals.index.tolist(), totals.values.tolist(), color=[COLORS["STRING_300"], COLORS["STRING_700"], COLORS["eFG"]])
    axes[0].set_ylabel("Edge count")
    axes[0].set_title("Total edge counts")
    overlap_parts = pd.DataFrame(
        {
            "condition": ["STRING_300 vs eFG", "STRING_700 vs eFG"],
            "shared": [overlap_df["string_300_intersect_efg_edge_count"], overlap_df["string_700_intersect_efg_edge_count"]],
            "string_only": [overlap_df["string_300_only_edge_count"], overlap_df["string_700_only_edge_count"]],
            "efg_only": [overlap_df["efg_only_vs_string_300_edge_count"], overlap_df["efg_only_vs_string_700_edge_count"]],
        }
    )
    y = range(len(overlap_parts))
    axes[1].barh(list(y), overlap_parts["shared"], color="#6a3d9a", label="Shared")
    axes[1].barh(list(y), overlap_parts["string_only"], left=overlap_parts["shared"], color="#1f77b4", label="STRING-only")
    axes[1].barh(list(y), overlap_parts["efg_only"], left=overlap_parts["shared"] + overlap_parts["string_only"], color="#d95f02", label="eFG-only")
    axes[1].set_yticks(list(y), overlap_parts["condition"].tolist())
    axes[1].set_xlabel("Edge count")
    axes[1].set_title("Pairwise overlap decomposition")
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle(f"Figure4D. {model_feature} STRING/eFG edge overlap decomposition", y=0.98)
    for ax in axes:
        style_axes(ax)
    fig.subplots_adjust(top=0.82, wspace=0.4)
    save_pdf_png(fig, output_dir / "Figure4D_edge_overlap_plot")

    # Supplementary source comparison
    supplementary_seed_export = supplementary_per_run_df[
        [
            "protocol",
            "species",
            "model",
            "feature_setting",
            "model_feature_setting",
            "graph_condition",
            "seed",
            "test_auprc",
            "test_mcc",
            "val_auprc",
            "val_mcc",
            "metrics_path",
        ]
    ].sort_values(["graph_condition", "seed"]).reset_index(drop=True)
    supplementary_seed_export.to_csv(supplementary_dir / "Figure4S1_eFG_confidence_source_comparison_barplot.tsv", sep="\t", index=False)
    supplementary_summary_df.to_csv(supplementary_dir / "Figure4S1_eFG_confidence_source_comparison_barplot_summary.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4), sharey=False, facecolor="white")
    supp_conditions = supplementary_summary_df["graph_condition"].tolist()
    supp_x = range(len(supp_conditions))
    supp_colors = [COLORS.get(condition, "#4c4c4c") for condition in supp_conditions]
    axes[0].bar(list(supp_x), supplementary_summary_df["test_auprc_mean"], yerr=supplementary_summary_df["test_auprc_sd"], capsize=4, color=supp_colors)
    axes[0].set_xticks(list(supp_x), supp_conditions, rotation=20)
    axes[0].set_ylabel("Test AUPRC")
    axes[0].set_title("AUPRC")
    axes[1].bar(list(supp_x), supplementary_summary_df["test_mcc_mean"], yerr=supplementary_summary_df["test_mcc_sd"], capsize=4, color=supp_colors)
    axes[1].set_xticks(list(supp_x), supp_conditions, rotation=20)
    axes[1].set_ylabel("Test MCC")
    axes[1].set_title("MCC")
    fig.suptitle(f"Figure4S1. {model_feature} eFG confidence-layered source comparison", y=0.98)
    for ax in axes:
        style_axes(ax)
    fig.subplots_adjust(top=0.8, wspace=0.35)
    save_pdf_png(fig, supplementary_dir / "Figure4S1_eFG_confidence_source_comparison_barplot")


if __name__ == "__main__":
    main()
