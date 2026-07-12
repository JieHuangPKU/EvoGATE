import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.analysis.figure5_final_common import (
    FIGURE5_PROTOCOLS,
    SUBSET,
    load_paired_cases,
    load_runtime_seed_list,
    protocol_output_slug,
    separation_metrics,
    species_title,
    write_manifest,
    write_markdown,
)


METRIC_ORDER = ["centroid_distance", "silhouette_score", "davies_bouldin_index"]
SPACE_ORDER = ["bio_input", "bio_esm2_input", "bio_hidden", "bio_esm2_hidden"]
SPACE_COLORS = {
    "bio_input": "#9AA5B1",
    "bio_esm2_input": "#F28E2B",
    "bio_hidden": "#4C78A8",
    "bio_esm2_hidden": "#D62728",
}
SPECIES_DISPLAY = {
    "fgraminearum": "Fusarium graminearum",
    "scerevisiae": "Saccharomyces cerevisiae",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Build final Figure5b hidden-space quantitative summary")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocols", nargs="+", default=FIGURE5_PROTOCOLS)
    parser.add_argument("--subset", default=SUBSET, type=str)
    parser.add_argument("--selection-table", required=True, type=str)
    parser.add_argument("--data-dir", default="results/Figure5/data", type=str)
    parser.add_argument("--table-dir", default="results/Figure5/tables", type=str)
    parser.add_argument("--plot-dir", default="results/Figure5/plots", type=str)
    parser.add_argument("--summary-dir", default="results/Figure5/summary", type=str)
    return parser.parse_args()


def manifest_row(category, path):
    return {"category": category, "path": str(Path(path).resolve())}


def determine_scaling_metadata(raw_df, stats_df):
    rows = []
    for metric in METRIC_ORDER:
        metric_values = raw_df.loc[raw_df["metric"] == metric, "value"].astype(float).abs()
        nonzero_values = metric_values[metric_values > 1e-12]
        max_abs_value = float(metric_values.max()) if not metric_values.empty else 0.0
        min_nonzero_value = float(nonzero_values.min()) if not nonzero_values.empty else 0.0
        dynamic_range = float(max_abs_value / min_nonzero_value) if min_nonzero_value > 0 else 0.0
        max_std = float(stats_df.loc[stats_df["metric"] == metric, "std"].fillna(0.0).max())
        apply_scale = bool(max_abs_value >= 80.0 or (dynamic_range >= 40.0 and max_abs_value >= 20.0 and max_std >= 5.0))
        rows.append(
            {
                "metric": metric,
                "scale_factor": 0.1 if apply_scale else 1.0,
                "scaling_applied": apply_scale,
                "max_abs_value": max_abs_value,
                "min_nonzero_abs_value": min_nonzero_value,
                "dynamic_range": dynamic_range,
                "max_std": max_std,
                "scaling_reason": (
                    "display values scaled by 0.1 because the raw metric range and variance exceed the readability threshold"
                    if apply_scale
                    else "no scaling applied because the metric remains readable in raw units"
                ),
            }
        )
    return pd.DataFrame(rows)


def metric_title(metric, scaling_lookup):
    scale_factor = float(scaling_lookup.loc[metric, "scale_factor"])
    if abs(scale_factor - 0.1) < 1e-12:
        return f"{metric} (display x0.1)"
    return metric


def scaled_values(df, scaling_lookup):
    scaled_df = df.copy()
    scaled_df["scale_factor"] = scaled_df["metric"].map(scaling_lookup["scale_factor"])
    if "value" in scaled_df.columns:
        scaled_df["display_value"] = scaled_df["value"].astype(float) * scaled_df["scale_factor"].astype(float)
    if "mean" in scaled_df.columns:
        scaled_df["display_mean"] = scaled_df["mean"].astype(float) * scaled_df["scale_factor"].astype(float)
    if "std" in scaled_df.columns:
        scaled_df["display_std"] = scaled_df["std"].astype(float) * scaled_df["scale_factor"].astype(float)
    return scaled_df


def display_species_title(species_slug):
    return SPECIES_DISPLAY.get(species_slug, str(species_slug))


def build_species_plot(plot_df, species_slug, plot_pdf, plot_png, title, use_error_bars):
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.2), facecolor="white")
    x = np.arange(len(SPACE_ORDER))
    scaling_lookup = plot_df[["metric", "scale_factor"]].drop_duplicates().set_index("metric")
    for ax, metric in zip(axes, METRIC_ORDER):
        subset = plot_df[plot_df["metric"] == metric].set_index("space").reindex(SPACE_ORDER).reset_index()
        if use_error_bars:
            ax.bar(
                x,
                subset["display_mean"],
                yerr=subset["display_std"],
                color=[SPACE_COLORS[space] for space in SPACE_ORDER],
                capsize=4,
            )
        else:
            ax.bar(
                x,
                subset["display_value"],
                color=[SPACE_COLORS[space] for space in SPACE_ORDER],
            )
        ax.set_xticks(x)
        ax.set_xticklabels(SPACE_ORDER, rotation=20, ha="right")
        ax.set_title(metric_title(metric, scaling_lookup))
        if metric == METRIC_ORDER[0]:
            ax.set_ylabel("Representative value" if not use_error_bars else "Mean +/- std across seeds")
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
    fig.suptitle(title.format(species=display_species_title(species_slug)), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(plot_pdf, format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(plot_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_dual_plot(plot_df, plot_pdf, plot_png, title, use_error_bars):
    species_order = ["fgraminearum", "scerevisiae"]
    fig, axes = plt.subplots(2, 3, figsize=(13.6, 7.2), facecolor="white")
    x = np.arange(len(SPACE_ORDER))
    scaling_lookup = plot_df[["metric", "scale_factor"]].drop_duplicates().set_index("metric")
    for row_idx, species_slug in enumerate(species_order):
        species_df = plot_df[plot_df["species"] == species_slug].copy()
        for col_idx, metric in enumerate(METRIC_ORDER):
            ax = axes[row_idx, col_idx]
            subset = species_df[species_df["metric"] == metric].set_index("space").reindex(SPACE_ORDER).reset_index()
            if use_error_bars:
                ax.bar(
                    x,
                    subset["display_mean"],
                    yerr=subset["display_std"],
                    color=[SPACE_COLORS[space] for space in SPACE_ORDER],
                    capsize=4,
                )
            else:
                ax.bar(
                    x,
                    subset["display_value"],
                    color=[SPACE_COLORS[space] for space in SPACE_ORDER],
                )
            ax.set_xticks(x)
            ax.set_xticklabels(SPACE_ORDER, rotation=20, ha="right")
            ax.set_title(f"{display_species_title(species_slug)}\n{metric_title(metric, scaling_lookup)}", fontsize=10)
            if col_idx == 0:
                ax.set_ylabel("Representative value" if not use_error_bars else "Mean +/- std across seeds")
            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(plot_pdf, format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(plot_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    data_dir = Path(args.data_dir)
    table_dir = Path(args.table_dir)
    plot_dir = Path(args.plot_dir)
    summary_dir = Path(args.summary_dir)
    for path in [data_dir, table_dir, plot_dir, summary_dir]:
        path.mkdir(parents=True, exist_ok=True)

    representative_seed = int(pd.read_csv(args.selection_table, sep="\t").query("is_representative_seed")["seed"].iloc[0])
    seeds = load_runtime_seed_list(args.runtime_config)

    raw_rows = []
    per_seed_rows = []
    for protocol in args.protocols:
        for seed in seeds:
            baseline_case, esm2_case, _ = load_paired_cases(args.runtime_config, args.upstream_root, protocol, seed, args.subset)
            species_slug = protocol_output_slug(protocol, baseline_case["species"])
            spaces = {
                "bio_input": baseline_case["input_matrix"],
                "bio_esm2_input": esm2_case["input_matrix"],
                "bio_hidden": baseline_case["hidden_matrix"],
                "bio_esm2_hidden": esm2_case["hidden_matrix"],
            }
            row = {"seed": seed, "protocol": protocol, "species": species_slug}
            for space_name, matrix in spaces.items():
                metrics = separation_metrics(matrix, baseline_case["labels"])
                for metric_name, metric_value in metrics.items():
                    raw_rows.append(
                        {
                            "seed": seed,
                            "protocol": protocol,
                            "species": species_slug,
                            "space": space_name,
                            "metric": metric_name,
                            "value": float(metric_value),
                        }
                    )
                    row[f"{space_name}__{metric_name}"] = float(metric_value)
            per_seed_rows.append(row)

    raw_df = pd.DataFrame(raw_rows)
    per_seed_df = pd.DataFrame(per_seed_rows).sort_values(["species", "seed"]).reset_index(drop=True)
    stats_df = (
        raw_df.groupby(["protocol", "species", "space", "metric"], dropna=False)["value"]
        .agg(["count", "mean", "std", "var"])
        .reset_index()
        .rename(columns={"count": "n", "var": "variance"})
    )
    stats_df["std"] = stats_df["std"].fillna(0.0)
    stats_df["variance"] = stats_df["variance"].fillna(0.0)
    scaling_df = determine_scaling_metadata(raw_df, stats_df)
    scaling_lookup = scaling_df.set_index("metric")

    raw_path = data_dir / "Figure5b_hidden_quant_summary_raw.tsv"
    per_seed_path = table_dir / "Figure5b_hidden_quant_summary_per_seed.tsv"
    stats_path = table_dir / "Figure5b_hidden_quant_summary_stats.tsv"
    all_species_stats_path = table_dir / "Figure5b_hidden_quant_summary_all_species_stats.tsv"
    scaling_path = table_dir / "Figure5b_hidden_quant_summary_scaling_metadata.tsv"
    raw_df.to_csv(raw_path, sep="\t", index=False)
    per_seed_df.to_csv(per_seed_path, sep="\t", index=False)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    stats_df.to_csv(all_species_stats_path, sep="\t", index=False)
    scaling_df.to_csv(scaling_path, sep="\t", index=False)

    manifest_rows = [
        manifest_row("data", raw_path),
        manifest_row("table", per_seed_path),
        manifest_row("table", stats_path),
        manifest_row("table", all_species_stats_path),
        manifest_row("table", scaling_path),
    ]

    representative_df = raw_df[raw_df["seed"] == representative_seed].copy()
    representative_df = scaled_values(representative_df, scaling_lookup)
    aggregated_plot_df = scaled_values(stats_df, scaling_lookup)

    representative_summary_paths = []
    representative_metric_tables = []
    for protocol in args.protocols:
        species_slug = protocol_output_slug(protocol)
        species_rep_df = representative_df[representative_df["protocol"] == protocol].copy()
        species_data_path = data_dir / f"Figure5b_hidden_quant_summary_{species_slug}_seed{representative_seed}.tsv"
        species_table_path = table_dir / f"Figure5b_hidden_quant_summary_{species_slug}_seed{representative_seed}_stats.tsv"
        species_plot_pdf = plot_dir / f"Figure5b_hidden_quant_summary_{species_slug}_seed{representative_seed}.pdf"
        species_plot_png = plot_dir / f"Figure5b_hidden_quant_summary_{species_slug}_seed{representative_seed}.png"
        species_summary_path = summary_dir / f"Figure5b_hidden_quant_summary_{species_slug}_seed{representative_seed}.md"
        species_rep_df.to_csv(species_data_path, sep="\t", index=False)

        species_metric_df = species_rep_df[["protocol", "species", "seed", "space", "metric", "value"]].copy()
        species_metric_df["n"] = 1
        species_metric_df["mean"] = species_metric_df["value"]
        species_metric_df["std"] = 0.0
        species_metric_df["variance"] = 0.0
        species_metric_df["scale_factor"] = species_metric_df["metric"].map(scaling_lookup["scale_factor"])
        species_metric_df["display_value"] = species_metric_df["value"] * species_metric_df["scale_factor"]
        species_metric_df.to_csv(species_table_path, sep="\t", index=False)
        representative_metric_tables.append(species_metric_df)

        build_species_plot(
            species_rep_df,
            species_slug,
            species_plot_pdf,
            species_plot_png,
            "Figure5b representative hidden-space separability ({species})",
            use_error_bars=False,
        )
        scaling_note = "No Figure5b metric required 0.1 display scaling." if not scaling_df["scaling_applied"].any() else "At least one Figure5b metric is displayed at x0.1; see the scaling metadata table."
        write_markdown(
            species_summary_path,
            [
                f"# Figure5b Hidden Quantitative Summary ({species_slug})",
                "",
                f"- Representative seed: `{representative_seed}`.",
                "- This file reports representative-seed values for the left-side Figure5b panel family.",
                f"- Quant summary data: `{species_data_path}`.",
                f"- Quant summary table: `{species_table_path}`.",
                f"- Plot PDF: `{species_plot_pdf}`.",
                f"- Plot PNG: `{species_plot_png}`.",
                f"- Scaling metadata: `{scaling_path}`.",
                f"- Scaling note: {scaling_note}",
            ],
        )
        representative_summary_paths.append(species_summary_path)
        for output_path in [species_data_path, species_table_path, species_plot_pdf, species_plot_png, species_summary_path]:
            category = output_path.parent.name.rstrip("s")
            manifest_rows.append(manifest_row(category, output_path))

    dual_data_path = data_dir / f"Figure5b_hidden_quant_summary_dual_species_seed{representative_seed}.tsv"
    dual_table_path = table_dir / f"Figure5b_hidden_quant_summary_dual_species_seed{representative_seed}_stats.tsv"
    dual_plot_pdf = plot_dir / f"Figure5b_hidden_quant_summary_dual_species_seed{representative_seed}.pdf"
    dual_plot_png = plot_dir / f"Figure5b_hidden_quant_summary_dual_species_seed{representative_seed}.png"
    dual_summary_path = summary_dir / f"Figure5b_hidden_quant_summary_dual_species_seed{representative_seed}.md"
    representative_df.to_csv(dual_data_path, sep="\t", index=False)
    dual_metric_df = pd.concat(representative_metric_tables, ignore_index=True)
    dual_metric_df.to_csv(dual_table_path, sep="\t", index=False)
    build_dual_plot(
        representative_df,
        dual_plot_pdf,
        dual_plot_png,
        "Figure5b representative hidden-space separability",
        use_error_bars=False,
    )
    write_markdown(
        dual_summary_path,
        [
            "# Figure5b Hidden Quantitative Summary (dual_species)",
            "",
            f"- Representative seed: `{representative_seed}`.",
            f"- Combined representative data: `{dual_data_path}`.",
            f"- Combined representative stats: `{dual_table_path}`.",
            f"- Plot PDF: `{dual_plot_pdf}`.",
            f"- Plot PNG: `{dual_plot_png}`.",
            f"- Scaling metadata: `{scaling_path}`.",
        ],
    )
    for output_path in [dual_data_path, dual_table_path, dual_plot_pdf, dual_plot_png, dual_summary_path]:
        category = output_path.parent.name.rstrip("s")
        manifest_rows.append(manifest_row(category, output_path))

    plot_pdf = plot_dir / "Figure5b_hidden_quant_summary.pdf"
    plot_png = plot_dir / "Figure5b_hidden_quant_summary.png"
    build_dual_plot(
        aggregated_plot_df,
        plot_pdf,
        plot_png,
        "Figure5b hidden-space separability across five seeds",
        use_error_bars=True,
    )
    manifest_rows.extend([manifest_row("plot", plot_pdf), manifest_row("plot", plot_png)])

    scaling_note = "No metric required 0.1 display scaling; all Figure5b plots remain in raw units." if not scaling_df["scaling_applied"].any() else "At least one Figure5b metric is displayed at x0.1 for readability; the raw values remain unchanged in all stored tables."
    summary_path = summary_dir / "Figure5b_hidden_quant_summary.md"
    write_markdown(
        summary_path,
        [
            "# Figure5b Hidden Quantitative Summary",
            "",
            f"- Seeds considered: `{', '.join(str(seed) for seed in seeds)}`.",
            f"- Representative seed used for seed-tagged panel exports: `{representative_seed}`.",
            "- Spaces compared: `bio_input`, `bio_esm2_input`, `bio_hidden`, `bio_esm2_hidden`.",
            "- Metrics reported with variance-aware statistics: `centroid_distance`, `silhouette_score`, `davies_bouldin_index`.",
            f"- Raw data: `{raw_path}`.",
            f"- Per-seed table: `{per_seed_path}`.",
            f"- Aggregated stats with mean/std/variance: `{stats_path}`.",
            f"- Scaling metadata: `{scaling_path}`.",
            f"- Aggregated manuscript-facing plot: `{plot_pdf}` / `{plot_png}`.",
            f"- Scaling note: {scaling_note}",
        ],
    )
    manifest_rows.append(manifest_row("summary", summary_path))

    manifest_path = table_dir / "Figure5b_output_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)


if __name__ == "__main__":
    main()
