import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.analysis.figure5_final_common import (
    CONFIDENCE_TIER_LABELS,
    CONFIDENCE_TIER_ORDER,
    FIGURE5_PROTOCOLS,
    STABLE_SEED_THRESHOLD,
    SUBSET,
    UMAP_PARAMS,
    compute_umap,
    confidence_tier_from_count,
    fine_scatter,
    focus_transition_columns,
    load_paired_cases,
    load_runtime_seed_list,
    protocol_output_slug,
    save_plot_pair,
    species_title,
    write_manifest,
    write_markdown,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build final Figure5a outputs")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocols", nargs="+", default=FIGURE5_PROTOCOLS)
    parser.add_argument("--subset", default=SUBSET, type=str)
    parser.add_argument("--selection-table", required=True, type=str)
    parser.add_argument("--stable-seed-threshold", default=STABLE_SEED_THRESHOLD, type=int)
    parser.add_argument("--data-dir", default="results/Figure5/data", type=str)
    parser.add_argument("--table-dir", default="results/Figure5/tables", type=str)
    parser.add_argument("--plot-dir", default="results/Figure5/plots", type=str)
    parser.add_argument("--summary-dir", default="results/Figure5/summary", type=str)
    return parser.parse_args()


def manifest_row(category, path):
    return {"category": category, "path": str(Path(path).resolve())}


def main():
    args = parse_args()
    data_dir = Path(args.data_dir)
    table_dir = Path(args.table_dir)
    plot_dir = Path(args.plot_dir)
    summary_dir = Path(args.summary_dir)
    for path in [data_dir, table_dir, plot_dir, summary_dir]:
        path.mkdir(parents=True, exist_ok=True)

    selection_df = pd.read_csv(args.selection_table, sep="\t")
    representative_seed = int(selection_df.loc[selection_df["is_representative_seed"], "seed"].iloc[0])
    seeds = load_runtime_seed_list(args.runtime_config)
    transition_defs = focus_transition_columns(args.stable_seed_threshold)
    transition_display_map = {item["raw"]: item["display"] for item in transition_defs}
    umap_params = dict(UMAP_PARAMS)
    umap_params["random_state"] = representative_seed

    per_seed_rows = []
    representative_frames = {}
    representative_counts_tables = []
    manifest_rows = []

    for protocol in args.protocols:
        for seed in seeds:
            baseline_case, esm2_case, paired = load_paired_cases(args.runtime_config, args.upstream_root, protocol, seed, args.subset)
            species_slug = protocol_output_slug(protocol, esm2_case["species"])
            per_seed = paired.copy().rename(columns={"graph_gene_id": "node_id"})
            per_seed["protocol"] = protocol
            per_seed["species"] = species_slug
            per_seed["seed"] = seed
            per_seed["subset"] = args.subset
            per_seed["delta_probability"] = per_seed["esm2_pred_score"].astype(float) - per_seed["baseline_pred_score"].astype(float)
            per_seed["transition_display"] = per_seed["transition"].map(transition_display_map).fillna(per_seed["transition"])
            per_seed_rows.append(per_seed)

            if seed != representative_seed:
                continue

            rep_coords = compute_umap(esm2_case["hidden_matrix"], umap_params)
            rep_df = per_seed.copy()
            rep_df["umap1"] = rep_coords[:, 0]
            rep_df["umap2"] = rep_coords[:, 1]
            representative_frames[species_slug] = rep_df.copy()

            stem = f"Figure5a_hidden_umap_error_transition_{species_slug}_seed{representative_seed}"
            coords_path = data_dir / f"{stem}_coords.tsv"
            transition_path = data_dir / f"{stem}_transition_labels.tsv"
            counts_path = table_dir / f"{stem}_transition_counts.tsv"
            pdf_path = plot_dir / f"{stem}.pdf"
            png_path = plot_dir / f"{stem}.png"
            summary_path = summary_dir / f"{stem}.md"

            rep_df.to_csv(coords_path, sep="\t", index=False)
            rep_df[
                [
                    "node_id",
                    "species",
                    "protocol",
                    "label",
                    "split",
                    "seed",
                    "baseline_pred_label",
                    "esm2_pred_label",
                    "baseline_pred_score",
                    "esm2_pred_score",
                    "transition",
                    "transition_display",
                ]
            ].to_csv(transition_path, sep="\t", index=False)
            rep_counts = (
                rep_df.groupby(["species", "protocol", "transition", "transition_display"], dropna=False)
                .size()
                .reset_index(name="gene_count")
            )
            rep_counts["seed"] = representative_seed
            rep_counts.to_csv(counts_path, sep="\t", index=False)
            representative_counts_tables.append(rep_counts)

            fig, ax = plt.subplots(figsize=(5.8, 5.0), facecolor="white")
            fine_scatter(
                ax,
                rep_df[["umap1", "umap2"]].to_numpy(),
                transitions=rep_df["transition"].astype(str).to_numpy(),
                title="{0}\nRepresentative hidden UMAP by error transition".format(species_title(protocol)),
                by="transition",
                legend=True,
            )
            save_plot_pair(fig, pdf_path, png_path)
            write_markdown(
                summary_path,
                [
                    f"# Figure5a Hidden UMAP Error Transition ({species_slug})",
                    "",
                    f"- Representative seed selected by the cross-species median-deviation rule: `{representative_seed}`.",
                    f"- UMAP random state for this representative view: `{representative_seed}`.",
                    "- This panel is a representative visualization; manuscript-facing stability summaries are computed across all five seeds.",
                    f"- Coordinates: `{coords_path}`.",
                    f"- Transition labels: `{transition_path}`.",
                    f"- Transition counts: `{counts_path}`.",
                    f"- Plot PDF: `{pdf_path}`.",
                    f"- Plot PNG: `{png_path}`.",
                ],
            )
            for output_path in [coords_path, transition_path, counts_path, pdf_path, png_path, summary_path]:
                category = output_path.parent.name.rstrip("s")
                manifest_rows.append(manifest_row(category, output_path))

    per_seed_df = pd.concat(per_seed_rows, ignore_index=True)
    per_seed_path = data_dir / "Figure5a_rescue_consensus_per_seed.tsv"
    per_seed_df.to_csv(per_seed_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", per_seed_path))

    group_cols = ["node_id", "species", "protocol", "label", "split"]
    metric_df = (
        per_seed_df.groupby(group_cols, dropna=False)
        .agg(
            n_seeds_observed=("seed", "nunique"),
            mean_delta_probability=("delta_probability", "mean"),
            std_delta_probability=("delta_probability", "std"),
            mean_baseline_probability=("baseline_pred_score", "mean"),
            mean_esm2_probability=("esm2_pred_score", "mean"),
            mean_baseline_pred_label=("baseline_pred_label", "mean"),
            mean_esm2_pred_label=("esm2_pred_label", "mean"),
        )
        .reset_index()
    )
    metric_df["std_delta_probability"] = metric_df["std_delta_probability"].fillna(0.0)

    transition_counts = (
        per_seed_df.groupby(group_cols + ["transition"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    gene_level_df = metric_df.merge(transition_counts, on=group_cols, how="left", validate="one_to_one")
    for item in transition_defs:
        if item["raw"] not in gene_level_df.columns:
            gene_level_df[item["raw"]] = 0
        gene_level_df[item["count_col"]] = gene_level_df[item["raw"]].fillna(0).astype(int)
        gene_level_df[item["frequency_col"]] = gene_level_df[item["count_col"]] / len(seeds)
        gene_level_df[item["stable_col"]] = gene_level_df[item["count_col"]] >= args.stable_seed_threshold
    gene_level_df["rescue_count"] = gene_level_df["fn_to_tp_rescued_count"]
    gene_level_df["rescue_frequency"] = gene_level_df["fn_to_tp_rescued_frequency"]
    stable_rescue_col = next(item["stable_col"] for item in transition_defs if item["slug"] == "fn_to_tp_rescued")
    gene_level_df["stable_rescued_ge2"] = gene_level_df[stable_rescue_col]
    gene_level_df["confidence_tier"] = gene_level_df["rescue_count"].astype(int).map(
        lambda count: confidence_tier_from_count(count, len(seeds), args.stable_seed_threshold)
    )
    gene_level_path = data_dir / "Figure5a_rescue_consensus_gene_level.tsv"
    gene_level_df.to_csv(gene_level_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", gene_level_path))

    summary_rows = []
    stable_bar_rows = []
    for protocol in args.protocols:
        species_slug = protocol_output_slug(protocol)
        species_gene_df = gene_level_df[gene_level_df["protocol"] == protocol].copy()
        essential_df = species_gene_df[species_gene_df["label"].astype(int) == 1].copy()
        row = {
            "protocol": protocol,
            "species": species_slug,
            "n_seeds": len(seeds),
            "stable_seed_threshold": args.stable_seed_threshold,
            "representative_seed": representative_seed,
            "total_test_genes": int(species_gene_df.shape[0]),
            "total_essential_test_genes": int(essential_df.shape[0]),
            "ever_rescued_genes": int((essential_df["rescue_count"] > 0).sum()),
            f"stable_rescued_genes_ge{args.stable_seed_threshold}": int((essential_df["rescue_count"] >= args.stable_seed_threshold).sum()),
            "high_confidence_rescued_genes_ge4": int((essential_df["confidence_tier"] == "high_confidence_rescued").sum()),
            f"stable_rescued_tier_genes_ge{args.stable_seed_threshold}_lt4": int((essential_df["confidence_tier"] == "stable_rescued").sum()),
            "seed_specific_rescued_genes_1seed": int((essential_df["confidence_tier"] == "seed_specific_rescued").sum()),
            "not_rescued_essential_genes": int((essential_df["confidence_tier"] == "not_rescued").sum()),
        }
        for item in transition_defs:
            stable_value = int(species_gene_df[item["stable_col"]].sum())
            row[f"stable_{item['display']}_genes_ge{args.stable_seed_threshold}"] = stable_value
            row[f"ever_{item['display']}_genes"] = int((species_gene_df[item["count_col"]] > 0).sum())
            stable_bar_rows.append(
                {
                    "protocol": protocol,
                    "species": species_slug,
                    "transition_display": item["display"],
                    "gene_count": stable_value,
                }
            )
        summary_rows.append(row)
    summary_table = pd.DataFrame(summary_rows).sort_values(["species"]).reset_index(drop=True)
    summary_table_path = table_dir / "Figure5a_rescue_consensus_summary.tsv"
    summary_table.to_csv(summary_table_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", summary_table_path))

    essential_gene_level_df = gene_level_df[gene_level_df["label"].astype(int) == 1].copy()
    tier_counts = (
        essential_gene_level_df.groupby(["protocol", "species", "confidence_tier"], dropna=False)
        .size()
        .reset_index(name="gene_count")
    )
    tier_counts["tier_label"] = tier_counts["confidence_tier"].map(CONFIDENCE_TIER_LABELS)
    tier_counts["stable_seed_threshold"] = args.stable_seed_threshold
    tier_counts["n_seeds"] = len(seeds)
    tier_counts["tier_rank"] = tier_counts["confidence_tier"].map({name: idx for idx, name in enumerate(CONFIDENCE_TIER_ORDER)})
    tier_counts = tier_counts.sort_values(["species", "tier_rank"]).drop(columns="tier_rank").reset_index(drop=True)
    tier_counts_path = table_dir / "Figure5a_rescue_confidence_tiers.tsv"
    tier_counts.to_csv(tier_counts_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", tier_counts_path))

    stability_rows = []
    for protocol in args.protocols:
        species_slug = protocol_output_slug(protocol)
        species_gene_df = gene_level_df[gene_level_df["protocol"] == protocol].copy()
        for item in transition_defs:
            counts = species_gene_df[item["count_col"]].astype(int)
            for occurrence_count in range(len(seeds) + 1):
                gene_count = int((counts == occurrence_count).sum())
                stability_rows.append(
                    {
                        "protocol": protocol,
                        "species": species_slug,
                        "transition": item["raw"],
                        "transition_display": item["display"],
                        "occurrence_count": occurrence_count,
                        "occurrence_frequency": occurrence_count / len(seeds),
                        "gene_count": gene_count,
                        "gene_fraction": float(gene_count / max(len(species_gene_df), 1)),
                        "stable_seed_threshold": args.stable_seed_threshold,
                        "stable_retained": occurrence_count >= args.stable_seed_threshold,
                    }
                )
    stability_dist_df = pd.DataFrame(stability_rows)
    stability_dist_path = table_dir / "Figure5a_rescue_transition_stability_distribution.tsv"
    stability_dist_df.to_csv(stability_dist_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", stability_dist_path))

    rescue_frequency_data = essential_gene_level_df[
        [
            "node_id",
            "species",
            "protocol",
            "label",
            "split",
            "rescue_count",
            "rescue_frequency",
            "stable_rescued_ge2",
            "confidence_tier",
            "mean_delta_probability",
            "std_delta_probability",
        ]
    ].copy()
    rescue_frequency_data_path = data_dir / "Figure5a_rescue_frequency_distribution_data.tsv"
    rescue_frequency_data.to_csv(rescue_frequency_data_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("data", rescue_frequency_data_path))

    rescue_frequency_table = (
        rescue_frequency_data.groupby(["protocol", "species", "rescue_count"], dropna=False)
        .size()
        .reset_index(name="gene_count")
    )
    rescue_frequency_table["rescue_frequency"] = rescue_frequency_table["rescue_count"] / len(seeds)
    rescue_frequency_table["gene_fraction"] = rescue_frequency_table["gene_count"] / rescue_frequency_table.groupby(
        ["protocol", "species"]
    )["gene_count"].transform("sum")
    rescue_frequency_table["stable_retained"] = rescue_frequency_table["rescue_count"] >= args.stable_seed_threshold
    rescue_frequency_table_path = table_dir / "Figure5a_rescue_frequency_distribution.tsv"
    rescue_frequency_table.to_csv(rescue_frequency_table_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", rescue_frequency_table_path))

    all_transition_counts_path = table_dir / "Figure5a_transition_counts_all_species.tsv"
    pd.concat(representative_counts_tables, ignore_index=True).to_csv(all_transition_counts_path, sep="\t", index=False)
    manifest_rows.append(manifest_row("table", all_transition_counts_path))

    dual_species_df = pd.concat(
        [representative_frames["fgraminearum"], representative_frames["scerevisiae"]],
        ignore_index=True,
    )
    dual_stem = f"Figure5a_hidden_umap_error_transition_dual_species_seed{representative_seed}"
    dual_coords_path = data_dir / f"{dual_stem}_coords.tsv"
    dual_transition_path = data_dir / f"{dual_stem}_transition_labels.tsv"
    dual_counts_path = table_dir / f"{dual_stem}_transition_counts.tsv"
    dual_pdf_path = plot_dir / f"{dual_stem}.pdf"
    dual_png_path = plot_dir / f"{dual_stem}.png"
    dual_summary_path = summary_dir / f"{dual_stem}.md"
    dual_species_df.to_csv(dual_coords_path, sep="\t", index=False)
    dual_species_df[
        [
            "node_id",
            "species",
            "protocol",
            "label",
            "split",
            "seed",
            "baseline_pred_label",
            "esm2_pred_label",
            "baseline_pred_score",
            "esm2_pred_score",
            "transition",
            "transition_display",
        ]
    ].to_csv(dual_transition_path, sep="\t", index=False)
    dual_counts = (
        dual_species_df.groupby(["species", "protocol", "transition", "transition_display"], dropna=False)
        .size()
        .reset_index(name="gene_count")
    )
    dual_counts["seed"] = representative_seed
    dual_counts.to_csv(dual_counts_path, sep="\t", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.9), facecolor="white")
    for ax, protocol in zip(axes, args.protocols):
        species_slug = protocol_output_slug(protocol)
        species_df = representative_frames[species_slug]
        fine_scatter(
            ax,
            species_df[["umap1", "umap2"]].to_numpy(),
            transitions=species_df["transition"].astype(str).to_numpy(),
            title=species_title(protocol),
            by="transition",
            legend=True,
        )
    fig.suptitle("Figure5a dual-species representative hidden UMAPs", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_plot_pair(fig, dual_pdf_path, dual_png_path)
    write_markdown(
        dual_summary_path,
        [
            "# Figure5a Hidden UMAP Error Transition (dual_species)",
            "",
            f"- Representative seed selected by the cross-species median-deviation rule: `{representative_seed}`.",
            f"- Combined representative coordinate table: `{dual_coords_path}`.",
            f"- Combined transition label table: `{dual_transition_path}`.",
            f"- Combined transition counts table: `{dual_counts_path}`.",
            f"- Plot PDF: `{dual_pdf_path}`.",
            f"- Plot PNG: `{dual_png_path}`.",
        ],
    )
    for output_path in [dual_coords_path, dual_transition_path, dual_counts_path, dual_pdf_path, dual_png_path, dual_summary_path]:
        category = output_path.parent.name.rstrip("s")
        manifest_rows.append(manifest_row(category, output_path))

    representative_visual_summary_path = summary_dir / "Figure5a_representative_visualization_summary.md"
    write_markdown(
        representative_visual_summary_path,
        [
            "# Figure5a Representative Visualization Summary",
            "",
            f"- Representative seed selected by the cross-species median-deviation rule: `{representative_seed}`.",
            f"- Representative species plots use seed-tagged filenames under `{plot_dir}` and `{summary_dir}`.",
            f"- Fusarium summary: `{summary_dir / f'Figure5a_hidden_umap_error_transition_fgraminearum_seed{representative_seed}.md'}`.",
            f"- Saccharomyces summary: `{summary_dir / f'Figure5a_hidden_umap_error_transition_scerevisiae_seed{representative_seed}.md'}`.",
            f"- Dual-species summary: `{dual_summary_path}`.",
        ],
    )
    manifest_rows.append(manifest_row("summary", representative_visual_summary_path))

    bar_pdf = plot_dir / "Figure5a_rescue_consensus_barplot.pdf"
    bar_png = plot_dir / "Figure5a_rescue_consensus_barplot.png"
    bar_df = pd.DataFrame(stable_bar_rows)
    transition_order = [item["display"] for item in transition_defs]
    species_order = [protocol_output_slug(protocol) for protocol in args.protocols]
    bar_df["transition_rank"] = bar_df["transition_display"].map({name: idx for idx, name in enumerate(transition_order)})
    bar_df["species_rank"] = bar_df["species"].map({name: idx for idx, name in enumerate(species_order)})
    bar_df = bar_df.sort_values(["transition_rank", "species_rank"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.2, 4.6), facecolor="white")
    x = np.arange(len(transition_order))
    width = 0.32
    colors = {"fgraminearum": "#D62728", "scerevisiae": "#4C78A8"}
    for idx, species_slug in enumerate(species_order):
        subset = bar_df[bar_df["species"] == species_slug].set_index("transition_display").reindex(transition_order).fillna(0.0).reset_index()
        ax.bar(
            x + (idx - (len(species_order) - 1) / 2) * width,
            subset["gene_count"],
            width=width,
            color=colors.get(species_slug, "#9AA5B1"),
            label=species_slug,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(transition_order, rotation=20, ha="right")
    ax.set_ylabel("Genes retained as stable")
    ax.set_title(f"Figure5a stability-retained transition counts (>= {args.stable_seed_threshold} of {len(seeds)} seeds)")
    ax.legend(frameon=False)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    save_plot_pair(fig, bar_pdf, bar_png)
    manifest_rows.extend([manifest_row("plot", bar_pdf), manifest_row("plot", bar_png)])

    freq_pdf = plot_dir / "Figure5a_rescue_frequency_distribution.pdf"
    freq_png = plot_dir / "Figure5a_rescue_frequency_distribution.png"
    fig, ax = plt.subplots(figsize=(7.4, 4.6), facecolor="white")
    rescue_count_order = list(range(len(seeds) + 1))
    width = 0.32
    for idx, species_slug in enumerate(species_order):
        subset = rescue_frequency_table[rescue_frequency_table["species"] == species_slug].set_index("rescue_count").reindex(rescue_count_order).fillna(0.0).reset_index()
        ax.bar(
            np.arange(len(rescue_count_order)) + (idx - (len(species_order) - 1) / 2) * width,
            subset["gene_count"],
            width=width,
            color=colors.get(species_slug, "#9AA5B1"),
            label=species_slug,
        )
    ax.axvline(args.stable_seed_threshold - 0.5, linestyle="--", linewidth=1.2, color="#5D5D5D")
    ax.set_xticks(np.arange(len(rescue_count_order)))
    ax.set_xticklabels([str(value) for value in rescue_count_order])
    ax.set_xlabel("FN_to_TP rescue count across five seeds")
    ax.set_ylabel("Essential genes")
    ax.set_title(f"Figure5a rescue frequency distribution (stable >= {args.stable_seed_threshold} seeds)")
    ax.legend(frameon=False)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    save_plot_pair(fig, freq_pdf, freq_png)
    manifest_rows.extend([manifest_row("plot", freq_pdf), manifest_row("plot", freq_png)])

    consensus_summary_path = summary_dir / "Figure5a_rescue_consensus_summary.md"
    write_markdown(
        consensus_summary_path,
        [
            "# Figure5a Rescue Consensus Summary",
            "",
            f"- Seeds considered: `{', '.join(str(seed) for seed in seeds)}`.",
            f"- Representative seed used for visualization: `{representative_seed}`.",
            f"- Manuscript-facing stability threshold: `>= {args.stable_seed_threshold}` seeds.",
            "- Rescued gene definition: essential test genes with `baseline_pred_label = 0` and `esm2_pred_label = 1` under matched seed pairing.",
            "- Confidence tiers: `high_confidence_rescued = 4-5 seeds`, `stable_rescued = 2-3 seeds`, `seed_specific_rescued = 1 seed`, `not_rescued = 0 seeds`.",
            f"- Per-seed transition table: `{per_seed_path}`.",
            f"- Gene-level consensus table: `{gene_level_path}`.",
            f"- Rescue consensus summary table: `{summary_table_path}`.",
            f"- Rescue confidence tiers table: `{tier_counts_path}`.",
            f"- Transition stability distribution table: `{stability_dist_path}`.",
            f"- Stability-retained consensus barplot: `{bar_pdf}` / `{bar_png}`.",
        ],
    )
    manifest_rows.append(manifest_row("summary", consensus_summary_path))

    stable_counts_for_summary = summary_table[["species", f"stable_rescued_genes_ge{args.stable_seed_threshold}"]].copy()
    stable_counts_for_summary = stable_counts_for_summary.rename(columns={f"stable_rescued_genes_ge{args.stable_seed_threshold}": "stable_rescued_gene_count"})
    frequency_summary_path = summary_dir / "Figure5a_rescue_frequency_distribution.md"
    frequency_lines = [
        "# Figure5a Rescue Frequency Distribution",
        "",
        f"- Frequency plot support data: `{rescue_frequency_data_path}`.",
        f"- Aggregated rescue frequency table: `{rescue_frequency_table_path}`.",
        f"- Plot PDF: `{freq_pdf}`.",
        f"- Plot PNG: `{freq_png}`.",
        f"- Stability threshold used for manuscript-facing retention: `>= {args.stable_seed_threshold}` seeds.",
        "",
        "## Stable rescued genes by species",
        "",
    ]
    for row in stable_counts_for_summary.itertuples(index=False):
        frequency_lines.append(f"- `{row.species}` stable rescued genes (>= {args.stable_seed_threshold} seeds): `{int(row.stable_rescued_gene_count)}`.")
    write_markdown(frequency_summary_path, frequency_lines)
    manifest_rows.append(manifest_row("summary", frequency_summary_path))

    manifest_path = table_dir / "Figure5a_output_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)


if __name__ == "__main__":
    main()
