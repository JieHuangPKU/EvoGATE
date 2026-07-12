#!/usr/bin/env python3
"""Plot a three-panel figure from precomputed label-transfer Sankey inputs.

This script reads existing files under results/label_transfer_sankey/ and does
display-only filtering/relabeling:
  - Positive_High -> Positive_Final
  - Positive (n=7) is omitted from the displayed Sankey and examples table

It does not rebuild any upstream data and does not modify input TSV files.
"""

from __future__ import annotations

import os
from pathlib import Path as FilePath

import matplotlib as mpl

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd


REPO_ROOT = FilePath(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results" / "label_transfer_sankey"

GENE_LEVEL_TSV = RESULTS_DIR / "sankey_gene_level_long.tsv"
EDGE_TSV = RESULTS_DIR / "sankey_aggregated_edges.tsv"
STAGE_COUNT_TSV = RESULTS_DIR / "sankey_stage_counts.tsv"
EXAMPLE_TSV = RESULTS_DIR / "representative_mapping_examples.tsv"
README_MD = RESULTS_DIR / "README.md"

PDF_OUT = RESULTS_DIR / "label_transfer_sankey_figure.pdf"
SVG_OUT = RESULTS_DIR / "label_transfer_sankey_figure.svg"
PNG_OUT = RESULTS_DIR / "label_transfer_sankey_figure.png"
LEGEND_OUT = RESULTS_DIR / "figure_legend.txt"

SOURCE_ORDER = ["scer_only", "spom_only", "both"]
SOURCE_LABELS = {
    "scer_only": "scer_only",
    "spom_only": "spom_only",
    "both": "both",
}
SOURCE_COLORS = {
    "scer_only": "#3C78A8",
    "spom_only": "#5B9E6D",
    "both": "#7A68A6",
}
FINAL_COLORS = {
    "Positive_Final": "#9E3D31",
    "Excluded": "#9A9A9A",
}
FRAME_COLOR = "#3A3A3A"
LIGHT_FILL = "#F4F4F4"
TEXT_COLOR = "#222222"


def require_inputs() -> None:
    missing = [path for path in [GENE_LEVEL_TSV, EDGE_TSV, STAGE_COUNT_TSV, EXAMPLE_TSV, README_MD] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required plotting inputs:\n- " + "\n- ".join(map(str, missing)))


def configure_matplotlib() -> None:
    FilePath(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_level = pd.read_csv(GENE_LEVEL_TSV, sep="\t", dtype=str).fillna("")
    edges = pd.read_csv(EDGE_TSV, sep="\t", dtype=str).fillna("")
    stage_counts = pd.read_csv(STAGE_COUNT_TSV, sep="\t", dtype=str).fillna("")
    examples = pd.read_csv(EXAMPLE_TSV, sep="\t", dtype=str).fillna("")
    return gene_level, edges, stage_counts, examples


def prepare_display_data(gene_level: pd.DataFrame, examples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    display = gene_level.loc[gene_level["bridge_resolved"].eq("yes")].copy()
    omitted_positive_n = int(display["final_label_class"].eq("Positive").sum())
    if omitted_positive_n != 7:
        raise ValueError(f"Expected display-only omitted Positive class size of 7, observed {omitted_positive_n}")
    display = display.loc[~display["final_label_class"].eq("Positive")].copy()
    display["display_final_class"] = display["final_label_class"].replace({"Positive_High": "Positive_Final"})

    example_display = examples.loc[~examples["final_label_class"].eq("Positive")].copy()
    example_display["final_label_class"] = example_display["final_label_class"].replace(
        {"Positive_High": "Positive_Final"}
    )
    example_display["yeast_species"] = example_display["yeast_species"].replace(
        {"scer": "S. cerevisiae", "spom": "S. pombe", "both": "Both"}
    )
    example_display["bridge_method"] = example_display["bridge_method"].str.replace("_", " ", regex=False)
    example_display["copy_status"] = example_display["copy_status"].str.replace("_", " ", regex=False)
    return display, example_display


def summarize_display(display: pd.DataFrame) -> dict[str, int]:
    return {
        "records_total": int(len(display)),
        "orthogroups_total": int(display["orthogroup_id"].nunique()),
        "genes_total": int(display["canonical_fusarium_gene_id"].nunique()),
        "positive_final_total": int(display["display_final_class"].eq("Positive_Final").sum()),
        "excluded_total": int(display["display_final_class"].eq("Excluded").sum()),
    }


def draw_flowchart(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Panel A  Workflow", loc="left", fontweight="bold", pad=6)

    boxes = [
        (0.05, 0.38, 0.2, 0.24, "Yeast essential\nsupport"),
        (0.29, 0.38, 0.2, 0.24, "OrthoFinder\northogroups"),
        (0.53, 0.38, 0.2, 0.24, "Canonical\nFusarium gene mapping"),
        (0.77, 0.38, 0.18, 0.24, "Final retained\npositives"),
    ]
    for x, y, w, h, text in boxes:
        ax.add_patch(Rectangle((x, y), w, h, facecolor="white", edgecolor=FRAME_COLOR, linewidth=1.2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=TEXT_COLOR, fontsize=9)

    for start, end in [(0.25, 0.29), (0.49, 0.53), (0.73, 0.77)]:
        arrow = FancyArrowPatch(
            (start, 0.50),
            (end, 0.50),
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.1,
            color=FRAME_COLOR,
        )
        ax.add_patch(arrow)

    ax.text(
        0.05,
        0.18,
        "Display uses precomputed canonical bridge outputs and final materialized labels only.",
        ha="left",
        va="center",
        fontsize=8,
        color="#555555",
    )


def _allocate_vertical_blocks(counts: list[int], y0: float = 0.08, y1: float = 0.92, gap: float = 0.03) -> list[tuple[float, float]]:
    total = sum(counts)
    usable = y1 - y0 - gap * (len(counts) - 1)
    heights = [usable * (count / total) if total else 0 for count in counts]
    blocks = []
    cursor = y1
    for height in heights:
        top = cursor
        bottom = cursor - height
        blocks.append((bottom, top))
        cursor = bottom - gap
    return blocks


def _flow_patch(x0: float, x1: float, y0a: float, y0b: float, y1a: float, y1b: float, color: str, alpha: float = 0.88) -> PathPatch:
    dx = (x1 - x0) * 0.45
    vertices = [
        (x0, y0b),
        (x0 + dx, y0b),
        (x1 - dx, y1b),
        (x1, y1b),
        (x1, y1a),
        (x1 - dx, y1a),
        (x0 + dx, y0a),
        (x0, y0a),
        (x0, y0b),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    return PathPatch(MplPath(vertices, codes), facecolor=color, edgecolor="none", alpha=alpha)


def draw_collapsed_alluvial(ax: plt.Axes, display: pd.DataFrame) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Panel B  Aggregated Transfer Alluvial", loc="left", fontweight="bold", pad=6)

    source_counts = (
        display.groupby("yeast_support_source")
        .size()
        .reindex(SOURCE_ORDER, fill_value=0)
        .astype(int)
        .to_dict()
    )
    final_counts = (
        display.groupby("display_final_class")
        .size()
        .reindex(["Positive_Final", "Excluded"], fill_value=0)
        .astype(int)
        .to_dict()
    )
    stats = summarize_display(display)

    x_source, x_og, x_gene, x_final = 0.08, 0.38, 0.65, 0.90
    node_w = 0.04

    source_blocks = _allocate_vertical_blocks([source_counts[key] for key in SOURCE_ORDER], gap=0.04)
    final_blocks = _allocate_vertical_blocks([final_counts["Positive_Final"], final_counts["Excluded"]], gap=0.05)

    og_block = _allocate_vertical_blocks([stats["records_total"]], gap=0.0)[0]
    gene_block = og_block

    source_to_final = (
        display.groupby(["yeast_support_source", "display_final_class"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=SOURCE_ORDER, columns=["Positive_Final", "Excluded"], fill_value=0)
        .astype(int)
    )

    source_running = {source: source_blocks[i][0] for i, source in enumerate(SOURCE_ORDER)}
    og_running = og_block[0]
    source_second_segments = {}

    # Stage 1: source -> orthogroup pool
    for idx, source in enumerate(SOURCE_ORDER):
        count = source_counts[source]
        bottom, top = source_blocks[idx]
        target_bottom = og_running
        target_top = og_running + (top - bottom)
        source_second_segments[source] = (target_bottom, target_top)
        ax.add_patch(_flow_patch(x_source + node_w, x_og, bottom, top, target_bottom, target_top, SOURCE_COLORS[source]))
        og_running = target_top

    # Stage 2: orthogroup pool -> mapped gene pool
    og_running = og_block[0]
    gene_running = gene_block[0]
    final_split_from_gene = {}
    for source in SOURCE_ORDER:
        source_bottom, source_top = source_second_segments[source]
        ax.add_patch(
            _flow_patch(x_og + node_w, x_gene, source_bottom, source_top, gene_running, gene_running + (source_top - source_bottom), SOURCE_COLORS[source], alpha=0.75)
        )
        final_split_from_gene[source] = (gene_running, gene_running + (source_top - source_bottom))
        gene_running += source_top - source_bottom

    # Stage 3: mapped gene pool -> final classes
    gene_offsets = {source: final_split_from_gene[source][0] for source in SOURCE_ORDER}
    final_offsets = {
        "Positive_Final": final_blocks[0][0],
        "Excluded": final_blocks[1][0],
    }
    for source in SOURCE_ORDER:
        for final_class in ["Positive_Final", "Excluded"]:
            count = int(source_to_final.loc[source, final_class])
            if count == 0:
                continue
            g0 = gene_offsets[source]
            g1 = g0 + count / stats["records_total"] * (gene_block[1] - gene_block[0])
            f0 = final_offsets[final_class]
            f1 = f0 + count / stats["records_total"] * (gene_block[1] - gene_block[0])
            flow_color = SOURCE_COLORS[source] if final_class == "Excluded" else FINAL_COLORS["Positive_Final"]
            flow_alpha = 0.45 if final_class == "Excluded" else 0.88
            ax.add_patch(_flow_patch(x_gene + node_w, x_final, g0, g1, f0, f1, flow_color, alpha=flow_alpha))
            gene_offsets[source] = g1
            final_offsets[final_class] = f1

    # Draw nodes
    for idx, source in enumerate(SOURCE_ORDER):
        bottom, top = source_blocks[idx]
        ax.add_patch(Rectangle((x_source, bottom), node_w, top - bottom, facecolor=SOURCE_COLORS[source], edgecolor=FRAME_COLOR, linewidth=0.8))
        ax.text(x_source - 0.01, (bottom + top) / 2, f"{SOURCE_LABELS[source]}\n(n={source_counts[source]})", ha="right", va="center", fontsize=8, color=TEXT_COLOR)

    for x, block, label in [
        (x_og, og_block, f"Orthogroups with\nessential support\n{stats['orthogroups_total']} unique"),
        (x_gene, gene_block, f"Mapped Fusarium\ngenes\n{stats['genes_total']} unique"),
    ]:
        bottom, top = block
        ax.add_patch(Rectangle((x, bottom), node_w, top - bottom, facecolor=LIGHT_FILL, edgecolor=FRAME_COLOR, linewidth=0.9))
        ax.text(x + node_w / 2, (bottom + top) / 2, label, ha="center", va="center", fontsize=8.5, color=TEXT_COLOR)

    final_labels = ["Positive_Final", "Excluded"]
    for idx, final_class in enumerate(final_labels):
        bottom, top = final_blocks[idx]
        ax.add_patch(Rectangle((x_final, bottom), node_w, top - bottom, facecolor=FINAL_COLORS[final_class], edgecolor=FRAME_COLOR, linewidth=0.8))
        ax.text(x_final + node_w + 0.01, (bottom + top) / 2, f"{final_class}\n(n={final_counts[final_class]})", ha="left", va="center", fontsize=8, color=TEXT_COLOR)

    headers = [
        (x_source + node_w / 2, "Yeast support\nsource"),
        (x_og + node_w / 2, "Orthogroups"),
        (x_gene + node_w / 2, "Canonical genes"),
        (x_final + node_w / 2, "Final class"),
    ]
    for x, header in headers:
        ax.text(x, 0.98, header, ha="center", va="top", fontsize=8.5, color=TEXT_COLOR, fontweight="bold")

    ax.text(
        0.02,
        0.015,
        "Display filter: Positive_High renamed to Positive_Final; minor Positive class (n=7) omitted.",
        ha="left",
        va="bottom",
        fontsize=7.8,
        color="#666666",
    )


def draw_examples_table(ax: plt.Axes, examples: pd.DataFrame) -> None:
    ax.axis("off")
    ax.set_title("Panel C  Representative Real Mapping Examples", loc="left", fontweight="bold", pad=6)

    preferred = pd.concat(
        [
            examples.loc[examples["final_label_class"].eq("Positive_Final")].head(6),
            examples.loc[examples["final_label_class"].eq("Excluded")].head(4),
        ],
        ignore_index=True,
    )
    if len(preferred) < 8:
        preferred = examples.head(8).copy()

    preferred = preferred.rename(
        columns={
            "yeast_gene": "Yeast gene",
            "yeast_species": "Species",
            "orthogroup_id": "Orthogroup",
            "canonical_fusarium_gene_id": "FGRAMPH1 gene ID",
            "occupancy_fraction": "Occupancy",
            "copy_status": "Copy status",
            "bridge_method": "Bridge method",
            "final_label_class": "Final class",
        }
    )
    preferred["Occupancy"] = preferred["Occupancy"].map(lambda value: f"{float(value):.2f}" if value != "" else "")
    preferred["Species"] = preferred["Species"].replace({"Both": "Both", "S. cerevisiae": "S. cerevisiae", "S. pombe": "S. pombe"})

    table_df = preferred[
        [
            "Yeast gene",
            "Species",
            "Orthogroup",
            "FGRAMPH1 gene ID",
            "Occupancy",
            "Copy status",
            "Bridge method",
            "Final class",
        ]
    ]

    col_widths = [0.11, 0.10, 0.11, 0.22, 0.08, 0.14, 0.16, 0.08]
    table = ax.table(
        cellText=table_df.values.tolist(),
        colLabels=table_df.columns.tolist(),
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.3)
    table.scale(1.0, 1.18)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#C8C8C8")
        cell.set_linewidth(0.5)
        if row == 0:
            cell.set_facecolor("#EDEDED")
            cell.set_text_props(weight="bold", color=TEXT_COLOR)
        else:
            cell.set_facecolor("white")
            final_class = table_df.iloc[row - 1]["Final class"]
            if col == len(table_df.columns) - 1:
                color = FINAL_COLORS.get(final_class, "#FFFFFF")
                cell.set_facecolor(color)
                cell.set_text_props(color="white" if final_class == "Positive_Final" else TEXT_COLOR, weight="bold")


def build_figure(display: pd.DataFrame, examples: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(15.5, 8.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.95], height_ratios=[1.0, 1.2])

    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[:, 1])
    ax_c = fig.add_subplot(grid[1, 0])

    draw_flowchart(ax_a)
    draw_collapsed_alluvial(ax_b, display)
    draw_examples_table(ax_c, examples)

    fig.suptitle(
        "Transfer of yeast essentiality evidence to Fusarium canonical gene IDs and final retained positives",
        x=0.03,
        y=0.995,
        ha="left",
        fontsize=13,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    fig.text(
        0.03,
        0.012,
        "A minor positive class (n=7) was omitted from the Sankey display because of its very small size; "
        "the final positive endpoint shown here corresponds to the high-confidence retained positive set.",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#555555",
    )
    return fig


def write_legend() -> None:
    legend_text = (
        "Panel A outlines the evidence-transfer workflow from yeast essential support through OrthoFinder "
        "orthogroups and canonical Fusarium gene mapping to final retained positives.\n\n"
        "Panel B shows an aggregated alluvial display derived from the real project Sankey inputs under "
        "results/label_transfer_sankey/. The displayed final positive endpoint, labeled Positive_Final, "
        "corresponds to the high-confidence retained positive set. A minor positive class (n=7) was omitted "
        "from the Sankey display because of its very small size.\n\n"
        "Panel C lists representative real mapping examples linking yeast support, orthogroup assignment, "
        "canonical FGRAMPH1 gene IDs, bridge method, and displayed final class.\n\n"
        "All counts are derived from real project data and canonical bridge outputs."
    )
    LEGEND_OUT.write_text(legend_text, encoding="utf-8")


def main() -> None:
    require_inputs()
    configure_matplotlib()
    gene_level, _edges, _stage_counts, examples = load_inputs()
    display, example_display = prepare_display_data(gene_level, examples)
    figure = build_figure(display, example_display)
    figure.savefig(PDF_OUT, bbox_inches="tight")
    figure.savefig(SVG_OUT, bbox_inches="tight")
    figure.savefig(PNG_OUT, dpi=300, bbox_inches="tight")
    plt.close(figure)
    write_legend()
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {SVG_OUT}")
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {LEGEND_OUT}")


if __name__ == "__main__":
    main()
