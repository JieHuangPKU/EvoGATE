import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.analysis.figure4_representation_common import TRANSITION_ORDER, fine_scatter, save_pdf, species_title


def parse_args():
    parser = argparse.ArgumentParser(description="Assemble dual-species Figure4 panels from single-species artifacts")
    parser.add_argument("--panel", required=True, choices=["a", "b", "c"])
    parser.add_argument("--left-species", required=True, type=str)
    parser.add_argument("--right-species", required=True, type=str)
    parser.add_argument("--left-input", required=True, type=str)
    parser.add_argument("--right-input", required=True, type=str)
    parser.add_argument("--output", required=True, type=str)
    return parser.parse_args()


def assemble_panel_a(left_species, right_species, left_input, right_input, output):
    left = pd.read_csv(left_input, sep="\t")
    right = pd.read_csv(right_input, sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), facecolor="white")
    fine_scatter(axes[0], left[["umap1", "umap2"]].to_numpy(), transitions=left["transition"].astype(str).to_numpy(), title=species_title(left_species), by="transition", legend=True)
    fine_scatter(axes[1], right[["umap1", "umap2"]].to_numpy(), transitions=right["transition"].astype(str).to_numpy(), title=species_title(right_species), by="transition", legend=True)
    fig.suptitle("Figure4A. GraphSAGE hidden UMAP colored by baseline to +ESM2 error transition", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_pdf(fig, output)


def assemble_panel_b(left_species, right_species, left_input, right_input, output):
    left = pd.read_csv(left_input, sep="\t")
    right = pd.read_csv(right_input, sep="\t")
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.2), facecolor="white")
    for row_idx, (species_name, df) in enumerate([(left_species, left), (right_species, right)]):
        sep_df = df[df["summary_type"] == "separation"].copy()
        err_df = df[df["summary_type"] == "error"].copy()
        metric_order = ["centroid_distance", "silhouette_score", "davies_bouldin_index"]
        x = np.arange(len(metric_order))
        width = 0.36
        for idx, feature_setting in enumerate(["ORT_EXP_SUB", "ORT_EXP_SUB_ESM2"]):
            subset = sep_df[sep_df["feature_setting"] == feature_setting].set_index("item").reindex(metric_order)
            axes[row_idx, 0].bar(x + (idx - 0.5) * width, subset["value"], width=width, label=feature_setting, color="#9AA5B1" if feature_setting == "ORT_EXP_SUB" else "#D62728")
        axes[row_idx, 0].set_xticks(x)
        axes[row_idx, 0].set_xticklabels(metric_order, rotation=20, ha="right")
        axes[row_idx, 0].set_title("{0} separation".format(species_title(species_name)), fontsize=10)
        axes[row_idx, 0].legend(frameon=False, fontsize=8)

        error_order = ["total_TP", "total_TN", "FN_to_TP_rescued", "FP_to_TN_corrected", "persistent_FN", "persistent_FP"]
        colors = ["#E45756", "#4C78A8", "#D62728", "#9467BD", "#FF9896", "#C5B0D5"]
        axes[row_idx, 1].bar(np.arange(len(error_order)), err_df.set_index("item").reindex(error_order)["value"], color=colors)
        axes[row_idx, 1].set_xticks(np.arange(len(error_order)))
        axes[row_idx, 1].set_xticklabels(error_order, rotation=25, ha="right")
        axes[row_idx, 1].set_title("{0} rescue summary".format(species_title(species_name)), fontsize=10)
    for ax in axes.ravel():
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.set_facecolor("white")
    fig.suptitle("Figure4B. Hidden-space quantitative summary and error rescue counts", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_pdf(fig, output)


def assemble_panel_c(left_species, right_species, left_input, right_input, output):
    left = pd.read_csv(left_input, sep="\t")
    right = pd.read_csv(right_input, sep="\t")
    fig = plt.figure(figsize=(14, 7), facecolor="white")
    outer = fig.add_gridspec(1, 2, wspace=0.18)
    for species_name, df, gs in [(left_species, left, outer[0]), (right_species, right, outer[1])]:
        inner = gs.subgridspec(2, 2, wspace=0.15, hspace=0.18)
        for ax, object_name in zip(
            [fig.add_subplot(inner[0, 0]), fig.add_subplot(inner[0, 1]), fig.add_subplot(inner[1, 0]), fig.add_subplot(inner[1, 1])],
            ["Input ORT_EXP_SUB", "Input ORT_EXP_SUB_ESM2", "Hidden ORT_EXP_SUB", "Hidden ORT_EXP_SUB_ESM2"],
        ):
            subset = df[df["feature_object"] == object_name].copy()
            fine_scatter(ax, subset[["umap1", "umap2"]].to_numpy(), labels=subset["label"].astype(int).to_numpy(), title=object_name, by="label", legend=True)
        title_ax = fig.add_subplot(gs)
        title_ax.axis("off")
        title_ax.set_title(species_title(species_name), fontsize=12, pad=16)
    fig.suptitle("Figure4C. Input-level vs hidden-level UMAP comparison", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_pdf(fig, output)


def main():
    args = parse_args()
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.panel == "a":
        assemble_panel_a(args.left_species, args.right_species, args.left_input, args.right_input, output)
    elif args.panel == "b":
        assemble_panel_b(args.left_species, args.right_species, args.left_input, args.right_input, output)
    else:
        assemble_panel_c(args.left_species, args.right_species, args.left_input, args.right_input, output)


if __name__ == "__main__":
    main()
