import argparse
from pathlib import Path

import pandas as pd


def write_markdown(path, lines):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Write Figure5 build summary markdown")
    parser.add_argument("--results-root", default="results/Figure5", type=str)
    parser.add_argument("--output", default="results/Figure5/summary/Figure5_build_summary.md", type=str)
    return parser.parse_args()


def load_manifest(path):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=["category", "path"])
    return pd.read_csv(path, sep="\t")


def has_pattern(manifest_df, pattern):
    if manifest_df.empty:
        return False
    return manifest_df["path"].astype(str).str.contains(pattern, regex=False).any()


def panel_status(panel_name, manifest_df, required_patterns):
    missing = [pattern for pattern in required_patterns if not has_pattern(manifest_df, pattern)]
    category_counts = manifest_df.groupby("category", dropna=False).size().to_dict()
    status = "complete" if not missing else "incomplete"
    return {
        "panel": panel_name,
        "plots": int(category_counts.get("plot", 0)),
        "data": int(category_counts.get("data", 0)),
        "tables": int(category_counts.get("table", 0)),
        "summary": int(category_counts.get("summary", 0)),
        "status": status,
        "missing_patterns": "; ".join(missing),
    }


def main():
    args = parse_args()
    root = Path(args.results_root).resolve()
    tables_dir = root / "tables"
    summary_dir = root / "summary"
    plot_dir = root / "plots"
    data_dir = root / "data"

    representative_table = pd.read_csv(tables_dir / "Figure5_representative_seed_selection.tsv", sep="\t")
    representative_seed = int(representative_table.loc[representative_table["is_representative_seed"], "seed"].iloc[0])

    manifest_a = load_manifest(tables_dir / "Figure5a_output_manifest.tsv")
    manifest_b = load_manifest(tables_dir / "Figure5b_output_manifest.tsv")
    manifest_c = load_manifest(tables_dir / "Figure5c_output_manifest.tsv")
    manifest_d = load_manifest(tables_dir / "Figure5d_output_manifest.tsv")

    contract_rows = [
        panel_status(
            "Figure5a",
            manifest_a,
            [
                "Figure5a_hidden_umap_error_transition_dual_species_seed",
                "Figure5a_hidden_umap_error_transition_fgraminearum_seed",
                "Figure5a_hidden_umap_error_transition_scerevisiae_seed",
                "Figure5a_rescue_consensus_per_seed.tsv",
                "Figure5a_rescue_consensus_gene_level.tsv",
                "Figure5a_rescue_consensus_summary.tsv",
                "Figure5a_rescue_confidence_tiers.tsv",
                "Figure5a_rescue_transition_stability_distribution.tsv",
                "Figure5a_rescue_frequency_distribution.tsv",
                "Figure5a_representative_visualization_summary.md",
                "Figure5a_rescue_consensus_summary.md",
                "Figure5a_rescue_frequency_distribution.md",
            ],
        ),
        panel_status(
            "Figure5b",
            manifest_b,
            [
                "Figure5b_hidden_quant_summary_dual_species_seed",
                "Figure5b_hidden_quant_summary_fgraminearum_seed",
                "Figure5b_hidden_quant_summary_scerevisiae_seed",
                "Figure5b_hidden_quant_summary_raw.tsv",
                "Figure5b_hidden_quant_summary_per_seed.tsv",
                "Figure5b_hidden_quant_summary_stats.tsv",
                "Figure5b_hidden_quant_summary_scaling_metadata.tsv",
                "Figure5b_hidden_quant_summary.md",
            ],
        ),
        panel_status(
            "Figure5c",
            manifest_c,
            [
                "Figure5c_input_vs_hidden_compare_dual_species_seed",
                "Figure5c_input_vs_hidden_compare_fgraminearum_seed",
                "Figure5c_input_vs_hidden_compare_scerevisiae_seed",
                "Figure5c_input_vs_hidden_compare_coords.tsv",
                "Figure5c_local_neighborhood_analysis.tsv",
                "Figure5c_input_vs_hidden_compare_metrics.tsv",
                "Figure5c_local_neighborhood_summary.tsv",
                "Figure5c_input_vs_hidden_compare.md",
            ],
        ),
        panel_status(
            "Figure5d",
            manifest_d,
            [
                "Figure5d_feature_group_manifest.tsv",
                "Figure5d_group_ablation_per_gene.tsv",
                "Figure5d_feature_group_global_importance_plot_data.tsv",
                "Figure5d_feature_group_rescued_vs_other_plot_data.tsv",
                "Figure5d_feature_group_gene_set_heatmap_plot_data.tsv",
                "Figure5d_group_ablation_global_summary.tsv",
                "Figure5d_group_ablation_global_metrics.tsv",
                "Figure5d_group_ablation_by_gene_set.tsv",
                "Figure5d_group_ablation_stats.tsv",
                "Figure5d_group_ablation_gene_set_heatmap_values.tsv",
                "Figure5d_feature_group_global_importance.pdf",
                "Figure5d_feature_group_rescued_vs_other.pdf",
                "Figure5d_feature_group_gene_set_heatmap.pdf",
                "Figure5d_feature_group_attribution.md",
            ],
        ),
    ]
    contract_df = pd.DataFrame(contract_rows)
    contract_path = tables_dir / "Figure5_contract_check.tsv"
    contract_df.to_csv(contract_path, sep="\t", index=False)

    write_markdown(
        Path(args.output),
        [
            "# Figure5 Build Summary",
            "",
            f"- Representative seed selected for seed-tagged exports: `{representative_seed}`.",
            "- The active Figure5 DAG now contains only the manuscript-facing mechanism / representation outputs under `results/Figure5/`.",
            "- Legacy export-style representation plots are intentionally excluded from the active Figure5 workflow.",
            "",
            "## Output roots",
            "",
            f"- Plots: `{plot_dir}`.",
            f"- Data: `{data_dir}`.",
            f"- Tables: `{tables_dir}`.",
            f"- Summary: `{summary_dir}`.",
            "",
            "## Contract check",
            "",
            *[
                f"- {row['panel']}: plots={row['plots']}, data={row['data']}, tables={row['tables']}, summary={row['summary']}, status=`{row['status']}`."
                for row in contract_rows
            ],
            "",
            "## Key tables",
            "",
            f"- Representative seed selection: `{tables_dir / 'Figure5_representative_seed_selection.tsv'}` ({len(representative_table)} rows).",
            f"- Figure5a transition counts: `{tables_dir / 'Figure5a_transition_counts_all_species.tsv'}`.",
            f"- Figure5b hidden summary stats: `{tables_dir / 'Figure5b_hidden_quant_summary_all_species_stats.tsv'}`.",
            f"- Figure5c input-vs-hidden metrics: `{tables_dir / 'Figure5c_input_vs_hidden_compare_all_species_metrics.tsv'}`.",
            f"- Figure5d global attribution summary: `{tables_dir / 'Figure5d_group_ablation_global_summary.tsv'}`.",
            f"- Figure5d per-gene attribution table: `{data_dir / 'Figure5d_group_ablation_per_gene.tsv'}`.",
            f"- Contract check table: `{contract_path}`.",
        ],
    )

    write_markdown(
        summary_dir / "Figure5_rerun_notes.md",
        [
            "# Figure5 Rerun Notes",
            "",
            "1. Run `scripts/run_Figure5_fusarium_graphsage_bio_esm2_representation_mechanism.sh -j 48` to rebuild the full Figure5 module.",
            "2. Run `scripts/run_Figure5a_fusarium_graphsage_bio_esm2_hidden_umap_error_transition.sh -j 48`, `scripts/run_Figure5b_fusarium_graphsage_bio_esm2_hidden_quant_summary.sh -j 48`, `scripts/run_Figure5c_fusarium_graphsage_bio_esm2_input_vs_hidden_compare.sh -j 48`, or `scripts/run_Figure5d_fusarium_graphsage_bio_esm2_feature_group_attribution.sh -j 48` for panel-specific rebuilds.",
            "3. The legacy export-style Figure5 representation plots are no longer part of the active Figure5 DAG.",
        ],
    )


if __name__ == "__main__":
    main()
