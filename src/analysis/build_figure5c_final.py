import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.analysis.figure5_final_common import (
    CONFIDENCE_TIER_ORDER,
    FIGURE5_PROTOCOLS,
    SUBSET,
    UMAP_PARAMS,
    compute_neighbor_fraction,
    compute_umap,
    fine_scatter,
    load_paired_cases,
    protocol_output_slug,
    save_plot_pair,
    separation_metrics,
    species_title,
    write_manifest,
    write_markdown,
)


SPACE_ORDER = ["bio_input", "bio_esm2_input", "bio_hidden", "bio_esm2_hidden"]
SPACE_LABELS = {
    "bio_input": "Bio input",
    "bio_esm2_input": "Bio+ESM2 input",
    "bio_hidden": "Bio hidden",
    "bio_esm2_hidden": "Bio+ESM2 hidden",
}
TIER_COLORS = {
    "high_confidence_rescued": "#D62728",
    "stable_rescued": "#F28E2B",
    "seed_specific_rescued": "#4C78A8",
    "not_rescued": "#9AA5B1",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Build final Figure5c local neighborhood mechanism outputs")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocols", nargs="+", default=FIGURE5_PROTOCOLS)
    parser.add_argument("--subset", default=SUBSET, type=str)
    parser.add_argument("--selection-table", required=True, type=str)
    parser.add_argument("--consensus-gene-table", required=True, type=str)
    parser.add_argument("--data-dir", default="results/Figure5/data", type=str)
    parser.add_argument("--table-dir", default="results/Figure5/tables", type=str)
    parser.add_argument("--plot-dir", default="results/Figure5/plots", type=str)
    parser.add_argument("--summary-dir", default="results/Figure5/summary", type=str)
    parser.add_argument("--k-neighbors", default=15, type=int)
    return parser.parse_args()


def manifest_row(category, path):
    return {"category": category, "path": str(Path(path).resolve())}


def build_species_compare_plot(coords_df, protocol, plot_pdf, plot_png):
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 7.8), facecolor="white")
    for ax, space_name in zip(axes.ravel(), SPACE_ORDER):
        subset = coords_df[coords_df["space"] == space_name].copy()
        fine_scatter(
            ax,
            subset[["umap1", "umap2"]].to_numpy(),
            labels=subset["label"].astype(int).to_numpy(),
            title=SPACE_LABELS[space_name],
            by="label",
            legend=True,
        )
    fig.suptitle(f"Figure5c input vs hidden comparison ({species_title(protocol)})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_plot_pair(fig, plot_pdf, plot_png)


def build_dual_compare_plot(species_coords, protocol_order, plot_pdf, plot_png):
    fig = plt.figure(figsize=(14.2, 7.2), facecolor="white")
    outer = fig.add_gridspec(1, 2, wspace=0.18)
    for protocol, grid in zip(protocol_order, outer):
        species_slug = protocol_output_slug(protocol)
        coords_df = species_coords[species_slug]
        inner = grid.subgridspec(2, 2, wspace=0.16, hspace=0.18)
        for ax, space_name in zip(
            [fig.add_subplot(inner[0, 0]), fig.add_subplot(inner[0, 1]), fig.add_subplot(inner[1, 0]), fig.add_subplot(inner[1, 1])],
            SPACE_ORDER,
        ):
            subset = coords_df[coords_df["space"] == space_name].copy()
            fine_scatter(
                ax,
                subset[["umap1", "umap2"]].to_numpy(),
                labels=subset["label"].astype(int).to_numpy(),
                title=SPACE_LABELS[space_name],
                by="label",
                legend=True,
            )
        title_ax = fig.add_subplot(grid)
        title_ax.axis("off")
        title_ax.set_title(species_title(protocol), fontsize=12, pad=18)
    fig.suptitle("Figure5c input vs hidden comparison (dual species)", fontsize=12)
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

    representative_seed = int(pd.read_csv(args.selection_table, sep="\t").query("is_representative_seed")["seed"].iloc[0])
    consensus_df = pd.read_csv(args.consensus_gene_table, sep="\t")
    consensus_df["label"] = consensus_df["label"].astype(int)
    consensus_df = consensus_df[consensus_df["label"] == 1].copy()
    umap_params = dict(UMAP_PARAMS)
    umap_params["random_state"] = representative_seed

    species_coords = {}
    metric_frames = []
    analysis_frames = []
    manifest_rows = []

    for protocol in args.protocols:
        baseline_case, esm2_case, _ = load_paired_cases(args.runtime_config, args.upstream_root, protocol, representative_seed, args.subset)
        species_slug = protocol_output_slug(protocol, baseline_case["species"])
        spaces = {
            "bio_input": baseline_case["input_matrix"],
            "bio_esm2_input": esm2_case["input_matrix"],
            "bio_hidden": baseline_case["hidden_matrix"],
            "bio_esm2_hidden": esm2_case["hidden_matrix"],
        }
        labels = baseline_case["labels"]
        node_ids = baseline_case["node_ids"]

        coords_rows = []
        metric_rows = []
        for space_name, matrix in spaces.items():
            coords = compute_umap(matrix, umap_params)
            metrics = separation_metrics(matrix, labels)
            for idx, node_id in enumerate(node_ids):
                coords_rows.append(
                    {
                        "node_id": node_id,
                        "species": species_slug,
                        "protocol": protocol,
                        "seed": representative_seed,
                        "space": space_name,
                        "label": int(labels[idx]),
                        "umap1": float(coords[idx, 0]),
                        "umap2": float(coords[idx, 1]),
                    }
                )
            for metric_name, metric_value in metrics.items():
                metric_rows.append(
                    {
                        "protocol": protocol,
                        "species": species_slug,
                        "seed": representative_seed,
                        "space": space_name,
                        "metric": metric_name,
                        "value": float(metric_value),
                        "n": 1,
                        "mean": float(metric_value),
                        "std": 0.0,
                        "variance": 0.0,
                    }
                )

        coords_df = pd.DataFrame(coords_rows)
        metrics_df = pd.DataFrame(metric_rows)
        species_coords[species_slug] = coords_df.copy()
        metric_frames.append(metrics_df.copy())

        coords_path = data_dir / f"Figure5c_input_vs_hidden_compare_{species_slug}_seed{representative_seed}_coords.tsv"
        metrics_path = table_dir / f"Figure5c_input_vs_hidden_compare_{species_slug}_seed{representative_seed}_metrics.tsv"
        compare_pdf = plot_dir / f"Figure5c_input_vs_hidden_compare_{species_slug}_seed{representative_seed}.pdf"
        compare_png = plot_dir / f"Figure5c_input_vs_hidden_compare_{species_slug}_seed{representative_seed}.png"
        summary_path = summary_dir / f"Figure5c_input_vs_hidden_compare_{species_slug}_seed{representative_seed}.md"
        coords_df.to_csv(coords_path, sep="\t", index=False)
        metrics_df.to_csv(metrics_path, sep="\t", index=False)
        build_species_compare_plot(coords_df, protocol, compare_pdf, compare_png)
        write_markdown(
            summary_path,
            [
                f"# Figure5c Input-Vs-Hidden Comparison ({species_slug})",
                "",
                f"- Representative seed: `{representative_seed}`.",
                f"- Coordinates: `{coords_path}`.",
                f"- Comparison metrics: `{metrics_path}`.",
                f"- Plot PDF: `{compare_pdf}`.",
                f"- Plot PNG: `{compare_png}`.",
            ],
        )
        for output_path in [coords_path, metrics_path, compare_pdf, compare_png, summary_path]:
            category = output_path.parent.name.rstrip("s")
            manifest_rows.append(manifest_row(category, output_path))

        essential_neighbor = {
            space_name: compute_neighbor_fraction(matrix, labels, args.k_neighbors)
            for space_name, matrix in spaces.items()
        }
        node_df = pd.DataFrame(
            {
                "node_id": node_ids,
                "label": labels.astype(int),
                "essential_neighbor_fraction_bio_input": essential_neighbor["bio_input"],
                "essential_neighbor_fraction_bio_esm2_input": essential_neighbor["bio_esm2_input"],
                "essential_neighbor_fraction_bio_hidden": essential_neighbor["bio_hidden"],
                "essential_neighbor_fraction_bio_esm2_hidden": essential_neighbor["bio_esm2_hidden"],
            }
        )
        species_consensus_df = consensus_df[consensus_df["protocol"] == protocol].copy()
        analysis_df = species_consensus_df.merge(node_df, on=["node_id", "label"], how="inner", validate="one_to_one")
        analysis_df["protocol"] = protocol
        analysis_df["species"] = species_slug
        analysis_df["seed"] = representative_seed
        analysis_df["k_neighbors"] = args.k_neighbors
        analysis_frames.append(analysis_df)

    combined_coords_df = pd.concat([species_coords["fgraminearum"], species_coords["scerevisiae"]], ignore_index=True)
    combined_metrics_df = pd.concat(metric_frames, ignore_index=True)

    coords_path = data_dir / "Figure5c_input_vs_hidden_compare_coords.tsv"
    metrics_path = table_dir / "Figure5c_input_vs_hidden_compare_metrics.tsv"
    all_species_metrics_path = table_dir / "Figure5c_input_vs_hidden_compare_all_species_metrics.tsv"
    dual_coords_path = data_dir / f"Figure5c_input_vs_hidden_compare_dual_species_seed{representative_seed}_coords.tsv"
    dual_metrics_path = table_dir / f"Figure5c_input_vs_hidden_compare_dual_species_seed{representative_seed}_metrics.tsv"
    dual_compare_pdf = plot_dir / f"Figure5c_input_vs_hidden_compare_dual_species_seed{representative_seed}.pdf"
    dual_compare_png = plot_dir / f"Figure5c_input_vs_hidden_compare_dual_species_seed{representative_seed}.png"
    dual_summary_path = summary_dir / f"Figure5c_input_vs_hidden_compare_dual_species_seed{representative_seed}.md"
    combined_coords_df.to_csv(coords_path, sep="\t", index=False)
    combined_coords_df.to_csv(dual_coords_path, sep="\t", index=False)
    combined_metrics_df.to_csv(metrics_path, sep="\t", index=False)
    combined_metrics_df.to_csv(all_species_metrics_path, sep="\t", index=False)
    combined_metrics_df.to_csv(dual_metrics_path, sep="\t", index=False)
    build_dual_compare_plot(species_coords, args.protocols, dual_compare_pdf, dual_compare_png)
    write_markdown(
        dual_summary_path,
        [
            "# Figure5c Input-Vs-Hidden Comparison (dual_species)",
            "",
            f"- Representative seed: `{representative_seed}`.",
            f"- Combined coordinates: `{dual_coords_path}`.",
            f"- Combined comparison metrics: `{dual_metrics_path}`.",
            f"- Plot PDF: `{dual_compare_pdf}`.",
            f"- Plot PNG: `{dual_compare_png}`.",
        ],
    )
    for output_path in [
        coords_path,
        metrics_path,
        all_species_metrics_path,
        dual_coords_path,
        dual_metrics_path,
        dual_compare_pdf,
        dual_compare_png,
        dual_summary_path,
    ]:
        category = output_path.parent.name.rstrip("s")
        manifest_rows.append(manifest_row(category, output_path))

    analysis_df = pd.concat(analysis_frames, ignore_index=True)
    analysis_path = data_dir / "Figure5c_local_neighborhood_analysis.tsv"
    analysis_df.to_csv(analysis_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", analysis_path))

    summary_rows = []
    for (protocol, species_slug, tier_name), subset in analysis_df.groupby(["protocol", "species", "confidence_tier"], dropna=False):
        for space_name in SPACE_ORDER:
            column = f"essential_neighbor_fraction_{space_name}"
            summary_rows.append(
                {
                    "protocol": protocol,
                    "species": species_slug,
                    "gene_set": tier_name,
                    "space": space_name,
                    "n": int(subset.shape[0]),
                    "mean": float(subset[column].mean()),
                    "std": float(subset[column].std(ddof=0)),
                    "variance": float(subset[column].var(ddof=0)),
                }
            )
    neighborhood_summary_df = pd.DataFrame(summary_rows)
    neighborhood_summary_path = table_dir / "Figure5c_local_neighborhood_summary.tsv"
    neighborhood_summary_df.to_csv(neighborhood_summary_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", neighborhood_summary_path))

    local_pdf = plot_dir / "Figure5c_local_neighborhood_summary.pdf"
    local_png = plot_dir / "Figure5c_local_neighborhood_summary.png"
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 4.8), facecolor="white")
    for ax, protocol in zip(axes, args.protocols):
        species_slug = protocol_output_slug(protocol)
        plot_df = neighborhood_summary_df[neighborhood_summary_df["protocol"] == protocol].copy()
        tiers_present = [tier for tier in CONFIDENCE_TIER_ORDER if tier in set(plot_df["gene_set"])]
        x = np.arange(len(SPACE_ORDER))
        width = 0.18
        for idx, tier_name in enumerate(tiers_present):
            tier_df = plot_df[plot_df["gene_set"] == tier_name].set_index("space").reindex(SPACE_ORDER).fillna(0.0).reset_index()
            ax.bar(
                x + (idx - (len(tiers_present) - 1) / 2) * width,
                tier_df["mean"],
                width=width,
                yerr=tier_df["std"],
                capsize=3,
                label=tier_name,
                color=TIER_COLORS.get(tier_name, "#7F7F7F"),
            )
        ax.set_xticks(x)
        ax.set_xticklabels([SPACE_LABELS[space_name] for space_name in SPACE_ORDER], rotation=20, ha="right")
        ax.set_ylabel("Mean essential-neighbor fraction")
        ax.set_title(species_title(protocol))
        ax.legend(frameon=False, fontsize=8)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
    fig.suptitle("Figure5c local neighborhood summary", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(local_pdf, format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(local_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    manifest_rows.extend([manifest_row("plot", local_pdf), manifest_row("plot", local_png)])

    summary_path = summary_dir / "Figure5c_input_vs_hidden_compare.md"
    write_markdown(
        summary_path,
        [
            "# Figure5c Input-Vs-Hidden Mechanism Summary",
            "",
            f"- Representative seed: `{representative_seed}`.",
            f"- Local mechanism analysis uses `k={args.k_neighbors}` nearest-neighbor essential-label cohesion on essential genes only.",
            "- Figure5c compares input-space and hidden-space layouts for the representative seed, then summarizes local essential-neighbor structure by rescue tier.",
            f"- Combined coordinates: `{coords_path}`.",
            f"- Combined comparison metrics: `{metrics_path}`.",
            f"- Local neighborhood analysis: `{analysis_path}`.",
            f"- Local neighborhood summary: `{neighborhood_summary_path}`.",
            f"- Dual-species compare plot: `{dual_compare_pdf}` / `{dual_compare_png}`.",
            f"- Local neighborhood plot: `{local_pdf}` / `{local_png}`.",
        ],
    )
    manifest_rows.append(manifest_row("summary", summary_path))

    manifest_path = table_dir / "Figure5c_output_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)


if __name__ == "__main__":
    main()
