#!/usr/bin/env python3
"""Build the final hybrid A/B/C publication figure for Fusarium label transfer.

Stack used:
  - matplotlib: master canvas, Panel A, final vector export
  - floweaver: Sankey topology definition / validation for Panel B
  - plottable: styled representative example table for Panel C

The script reads existing real data products from results/label_transfer_sankey/
and performs display-only relabeling:
  - Positive_High -> Positive_Final
  - omit the minor Positive class (n=7) from the Sankey display
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import pandas as pd
from palettable.cartocolors.qualitative import Safe_10
from floweaver import Bundle, Dataset, Partition, ProcessGroup, SankeyDefinition, Waypoint, weave
from plottable import ColumnDefinition, Table


REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results" / "label_transfer_sankey"

GENE_LEVEL_TSV = RESULTS_DIR / "sankey_gene_level_long.tsv"
EDGE_TSV = RESULTS_DIR / "sankey_aggregated_edges.tsv"
STAGE_COUNT_TSV = RESULTS_DIR / "sankey_stage_counts.tsv"
EXAMPLE_TSV = RESULTS_DIR / "representative_mapping_examples.tsv"
README_MD = RESULTS_DIR / "README.md"

PH1_TRANSFER_TSV = REPO_ROOT / "data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv"
ORTHOGROUPS_TSV = Path(
    "/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/"
    "run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.tsv"
)
PROTEIN_BRIDGE_TSV = REPO_ROOT / "data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv"

PDF_OUT = RESULTS_DIR / "label_transfer_sankey_figure.pdf"
SVG_OUT = RESULTS_DIR / "label_transfer_sankey_figure.svg"
PNG_OUT = RESULTS_DIR / "label_transfer_sankey_figure.png"
LEGEND_OUT = RESULTS_DIR / "figure_legend.txt"
PALETTE_OUT = RESULTS_DIR / "palette_mapping.tsv"

SAFE10 = Safe_10.hex_colors
PALETTE = {
    "scer_only": SAFE10[0],
    "spom_only": SAFE10[1],
    "both": SAFE10[2],
    "Positive_Final": SAFE10[3],
    "Excluded": "#8E8E8E",
    "panelA_source_fill": "#EAF4FB",
    "panelA_bridge_fill": "#F8F1DA",
    "panelA_target_fill": "#E7F4E8",
    "panelA_summary_fill": "#F6F7F9",
}

SOURCE_ORDER = ["scer_only", "spom_only", "both"]
FINAL_ORDER = ["Positive_Final", "Excluded"]

TEXT = "#222222"
SUBTEXT = "#6B6B6B"
EDGE = "#495057"
MID_FILL = "#F4F5F7"
HEADER_BG = "#F4F5F7"


@dataclass
class DisplaySummary:
    total_support_universe: int
    total_orthogroups: int
    total_canonical_genes: int
    resolved_display_universe: int
    tracked_orthogroups: int
    tracked_genes: int
    positive_final: int
    excluded: int
    omitted_positive: int


def configure_matplotlib() -> None:
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def require_inputs() -> None:
    required = [
        GENE_LEVEL_TSV,
        EDGE_TSV,
        STAGE_COUNT_TSV,
        EXAMPLE_TSV,
        README_MD,
        PH1_TRANSFER_TSV,
        ORTHOGROUPS_TSV,
        PROTEIN_BRIDGE_TSV,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n- " + "\n- ".join(missing))


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_level = read_tsv(GENE_LEVEL_TSV)
    examples = read_tsv(EXAMPLE_TSV)
    return gene_level, examples


def prepare_display_tables(gene_level: pd.DataFrame, examples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, DisplaySummary]:
    display = gene_level.loc[gene_level["bridge_resolved"].eq("yes")].copy()
    omitted_positive = int(display["final_label_class"].eq("Positive").sum())
    if omitted_positive != 7:
        raise ValueError(f"Expected omitted Positive class size 7, observed {omitted_positive}")
    display = display.loc[~display["final_label_class"].eq("Positive")].copy()
    display["display_final_class"] = display["final_label_class"].replace({"Positive_High": "Positive_Final"})

    examples_display = examples.loc[~examples["final_label_class"].eq("Positive")].copy()
    examples_display["final_label_class"] = examples_display["final_label_class"].replace(
        {"Positive_High": "Positive_Final"}
    )

    ph1_transfer = read_tsv(PH1_TRANSFER_TSV)
    orthogroups = read_tsv(ORTHOGROUPS_TSV)
    bridge = read_tsv(PROTEIN_BRIDGE_TSV)
    resolved_bridge = bridge.loc[bridge["bridge_status"].eq("resolved") & bridge["resolved_canonical_gene_id"].ne("")]

    summary = DisplaySummary(
        total_support_universe=int(len(ph1_transfer)),
        total_orthogroups=int(orthogroups["Orthogroup"].nunique()),
        total_canonical_genes=int(resolved_bridge["resolved_canonical_gene_id"].nunique()),
        resolved_display_universe=int(len(gene_level.loc[gene_level["bridge_resolved"].eq("yes")])),
        tracked_orthogroups=int(display["orthogroup_id"].nunique()),
        tracked_genes=int(display["canonical_fusarium_gene_id"].nunique()),
        positive_final=int(display["display_final_class"].eq("Positive_Final").sum()),
        excluded=int(display["display_final_class"].eq("Excluded").sum()),
        omitted_positive=omitted_positive,
    )
    return display, examples_display, summary


def build_floweaver_sankey(display: pd.DataFrame):
    """Use floweaver to define and validate the Sankey topology for Panel B."""
    source_to_final = (
        display.groupby(["yeast_support_source", "display_final_class"], as_index=False)
        .size()
        .rename(columns={"yeast_support_source": "source", "display_final_class": "target", "size": "value"})
    )
    source_to_final["og_stage"] = "Yeast-supported OGs"
    source_to_final["gene_stage"] = "Mapped genes"

    dataset = Dataset(source_to_final[["source", "target", "og_stage", "gene_stage", "value"]])
    nodes = {
        "source_stage": ProcessGroup(
            SOURCE_ORDER,
            partition=Partition.Simple("source", SOURCE_ORDER),
            title="Yeast essential support",
        ),
        "og_stage": Waypoint(
            partition=Partition.Simple("og_stage", ["Yeast-supported OGs"]),
            title="OrthoFinder orthogroups",
        ),
        "gene_stage": Waypoint(
            partition=Partition.Simple("gene_stage", ["Mapped genes"]),
            title="Canonical Fusarium genes",
        ),
        "final_stage": ProcessGroup(
            FINAL_ORDER,
            partition=Partition.Simple("target", FINAL_ORDER),
            title="Final class",
        ),
    }
    bundles = [Bundle("source_stage", "final_stage", waypoints=["og_stage", "gene_stage"])]
    ordering = [["source_stage"], ["og_stage"], ["gene_stage"], ["final_stage"]]
    sankey_data = weave(SankeyDefinition(nodes, bundles, ordering), dataset, measures="value")
    return sankey_data, source_to_final


def _flow_patch(x0: float, x1: float, y0a: float, y0b: float, y1a: float, y1b: float, color: str, alpha: float = 0.95) -> PathPatch:
    dx = (x1 - x0) * 0.44
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


def allocate_segments(counts: list[int], y0: float, total_height: float) -> list[tuple[float, float]]:
    total = sum(counts)
    cursor = y0
    blocks = []
    for count in counts:
        h = total_height * count / total if total else 0.0
        blocks.append((cursor, cursor + h))
        cursor += h
    return blocks


def draw_panel_a(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Panel A  Workflow", loc="left", fontweight="bold", fontsize=12, pad=8)

    badges = [
        (0.08, 0.61, 0.22, 0.20, PALETTE["panelA_source_fill"], "Yeast essential genes", "source support space"),
        (0.39, 0.61, 0.22, 0.20, PALETTE["panelA_bridge_fill"], "OrthoFinder orthogroups", "cross-species bridge"),
        (0.70, 0.61, 0.22, 0.20, PALETTE["panelA_target_fill"], "Fusarium PH-1 gene IDs", "canonical target space"),
    ]
    for x, y, w, h, fill, title, subtitle in badges:
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.016,rounding_size=0.035",
            linewidth=1.2,
            edgecolor=EDGE,
            facecolor=fill,
        )
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h * 0.62, title, ha="center", va="center", fontsize=9.1, fontweight="bold", color=TEXT)
        ax.text(x + w / 2, y + h * 0.30, subtitle, ha="center", va="center", fontsize=8.0, color=SUBTEXT)

    connector_y = 0.71
    for x0, x1 in [(0.30, 0.38), (0.61, 0.69)]:
        ax.annotate(
            "",
            xy=(x1 - 0.012, connector_y),
            xytext=(x0 + 0.012, connector_y),
            arrowprops=dict(arrowstyle="-|>", lw=1.2, color="#7A7A7A", shrinkA=4, shrinkB=4, mutation_scale=11),
        )

    summary_box = FancyBboxPatch(
        (0.14, 0.17),
        0.72,
        0.22,
        boxstyle="round,pad=0.018,rounding_size=0.032",
        linewidth=1.0,
        edgecolor="#C7CCD1",
        facecolor=PALETTE["panelA_summary_fill"],
    )
    ax.add_patch(summary_box)
    ax.text(
        0.50,
        0.33,
        "Yeast essential support  \u2192  orthogroup filtering  \u2192  canonical mapping  \u2192  final retained positives",
        ha="center",
        va="center",
        fontsize=9.1,
        fontweight="bold",
        color=TEXT,
    )
    ax.text(
        0.50,
        0.24,
        "Real project evidence is transferred through OrthoFinder-defined orthogroups and bridge-resolved to canonical Fusarium PH-1 genes.",
        ha="center",
        va="center",
        fontsize=7.9,
        color=SUBTEXT,
    )


def draw_panel_b(ax: plt.Axes, display: pd.DataFrame, summary: DisplaySummary, source_to_final: pd.DataFrame) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Panel B  Tracked Transfer Subset", loc="left", fontweight="bold", fontsize=12, pad=8)

    x_positions = {"source": 0.10, "og": 0.40, "gene": 0.67, "final": 0.90}
    node_w = 0.055
    y0 = 0.16
    total_h = 0.63

    headers = [
        ("source", "Yeast essential\nsupport", f"Total n = {summary.total_support_universe:,}"),
        ("og", "OrthoFinder\northogroups", f"Total n = {summary.total_orthogroups:,}"),
        ("gene", "Canonical Fusarium\ngenes", f"Total n = {summary.total_canonical_genes:,}"),
        ("final", "Final class", f"Resolved n = {summary.resolved_display_universe:,}"),
    ]
    for key, title, subtitle in headers:
        x = x_positions[key] + node_w / 2
        ax.text(x, 0.94, title, ha="center", va="top", fontsize=10.2, fontweight="bold", color=TEXT)
        ax.text(x, 0.872, subtitle, ha="center", va="top", fontsize=8.2, color=SUBTEXT)

    source_counts = (
        display.groupby("yeast_support_source")
        .size()
        .reindex(SOURCE_ORDER, fill_value=0)
        .astype(int)
        .to_dict()
    )
    source_segments = allocate_segments([source_counts[src] for src in SOURCE_ORDER], y0, total_h)
    og_segments = allocate_segments([source_counts[src] for src in SOURCE_ORDER], y0, total_h)
    gene_segments = allocate_segments([source_counts[src] for src in SOURCE_ORDER], y0, total_h)
    final_segments = allocate_segments([summary.positive_final, summary.excluded], y0, total_h)

    # source nodes
    for idx, source in enumerate(SOURCE_ORDER):
        yb, yt = source_segments[idx]
        x = x_positions["source"]
        ax.add_patch(Rectangle((x, yb), node_w, yt - yb, facecolor=PALETTE[source], edgecolor=EDGE, linewidth=0.75))
        ax.text(x - 0.018, (yb + yt) / 2, f"{source}\n{source_counts[source]:,}", ha="right", va="center", fontsize=8.3, color=TEXT)

    # middle nodes
    for key, label, value in [
        ("og", "Yeast-supported OGs", summary.tracked_orthogroups),
        ("gene", "Mapped genes", summary.tracked_genes),
    ]:
        x = x_positions[key]
        ax.add_patch(Rectangle((x, y0), node_w, total_h, facecolor=MID_FILL, edgecolor=EDGE, linewidth=0.85))
        ax.text(x + node_w / 2, y0 + total_h / 2, f"{label}\n{value:,}", ha="center", va="center", fontsize=8.8, color=TEXT)

    # final nodes
    final_label_map = {"Positive_Final": summary.positive_final, "Excluded": summary.excluded}
    for idx, final_class in enumerate(FINAL_ORDER):
        yb, yt = final_segments[idx]
        x = x_positions["final"]
        ax.add_patch(Rectangle((x, yb), node_w, yt - yb, facecolor=PALETTE[final_class], edgecolor=EDGE, linewidth=0.85))
        ax.text(x + node_w + 0.018, (yb + yt) / 2, f"{final_class}\n{final_label_map[final_class]:,}", ha="left", va="center", fontsize=8.8, color=TEXT, fontweight="bold")

    # source-colored flows to gene pool
    x1 = x_positions["source"] + node_w
    x2 = x_positions["og"]
    x3 = x_positions["og"] + node_w
    x4 = x_positions["gene"]
    x5 = x_positions["gene"] + node_w
    x6 = x_positions["final"]
    for idx, source in enumerate(SOURCE_ORDER):
        ax.add_patch(_flow_patch(x1, x2, source_segments[idx][0], source_segments[idx][1], og_segments[idx][0], og_segments[idx][1], PALETTE[source], alpha=0.97))
        ax.add_patch(_flow_patch(x3, x4, og_segments[idx][0], og_segments[idx][1], gene_segments[idx][0], gene_segments[idx][1], PALETTE[source], alpha=0.93))

    # final-segment flows colored by endpoint
    source_to_final_map = (
        source_to_final.set_index(["source", "target"])["value"].astype(int).to_dict()
    )
    gene_offsets = {src: gene_segments[i][0] for i, src in enumerate(SOURCE_ORDER)}
    final_offsets = {"Positive_Final": final_segments[0][0], "Excluded": final_segments[1][0]}
    final_total_map = {"Positive_Final": summary.positive_final, "Excluded": summary.excluded}
    final_height_map = {
        "Positive_Final": final_segments[0][1] - final_segments[0][0],
        "Excluded": final_segments[1][1] - final_segments[1][0],
    }
    for idx, source in enumerate(SOURCE_ORDER):
        src_total = source_counts[source]
        src_h = gene_segments[idx][1] - gene_segments[idx][0]
        for final_class in FINAL_ORDER:
            count = source_to_final_map.get((source, final_class), 0)
            if count == 0:
                continue
            h0 = src_h * count / src_total if src_total else 0.0
            y_start0 = gene_offsets[source]
            y_start1 = y_start0 + h0
            h1 = final_height_map[final_class] * count / final_total_map[final_class] if final_total_map[final_class] else 0.0
            y_end0 = final_offsets[final_class]
            y_end1 = y_end0 + h1
            ax.add_patch(_flow_patch(x5, x6, y_start0, y_start1, y_end0, y_end1, PALETTE[final_class], alpha=0.94))
            gene_offsets[source] = y_start1
            final_offsets[final_class] = y_end1

    ax.text(
        0.10,
        0.06,
        "A minor Positive class (n=7) is omitted from display; only the tracked transfer subset is shown in the Sankey.",
        ha="left",
        va="center",
        fontsize=7.8,
        color=SUBTEXT,
    )


def build_panel_c_table(examples: pd.DataFrame, gene_level: pd.DataFrame) -> pd.DataFrame:
    examples = examples.copy()
    examples["final_label_class"] = examples["final_label_class"].replace({"Positive_High": "Positive_Final"})
    examples = examples.loc[~examples["final_label_class"].eq("Positive")].copy()

    lookup = (
        gene_level.drop_duplicates(subset=["canonical_fusarium_gene_id"])
        .set_index("canonical_fusarium_gene_id")
    )
    examples["ph1_copy_count"] = examples["canonical_fusarium_gene_id"].map(lookup["ph1_copy_count_from_membership"])
    examples["genome_occupancy"] = examples["canonical_fusarium_gene_id"].map(lookup["occupancy_count"])
    examples["phi_lethal_evidence"] = examples["canonical_fusarium_gene_id"].map(
        lambda g: "Yes" if str(lookup.loc[g, "is_lethal_supported_positive"]).lower() == "yes" else "No"
        if g in lookup.index
        else "No"
    )
    examples["yeast_species"] = examples["yeast_species"].replace({"scer": "S. cerevisiae", "spom": "S. pombe", "both": "Both"})
    examples["bridge_method"] = examples["bridge_method"].str.replace("_", " ", regex=False)
    examples["copy_status"] = examples["copy_status"].str.replace("_", " ", regex=False)
    examples["occupancy_fraction"] = examples["occupancy_fraction"].map(lambda x: f"{float(x):.2f}" if x else "")

    curated = pd.concat(
        [
            examples.loc[examples["final_label_class"].eq("Positive_Final")].head(6),
            examples.loc[examples["final_label_class"].eq("Excluded")].head(4),
        ],
        ignore_index=True,
    ).drop_duplicates(subset=["orthogroup_id", "canonical_fusarium_gene_id"])

    curated = curated.rename(
        columns={
            "yeast_gene": "Yeast essential gene",
            "yeast_species": "Yeast species",
            "orthogroup_id": "OrthoFinder orthogroup",
            "canonical_fusarium_gene_id": "Fusarium (PH-1) gene ID",
            "ph1_copy_count": "Copy number in PH-1",
            "occupancy_fraction": "Genome occupancy (18 spp.)",
            "phi_lethal_evidence": "PHI-base lethal evidence",
            "final_label_class": "Final label",
        }
    )
    return curated[
        [
            "Yeast essential gene",
            "Yeast species",
            "OrthoFinder orthogroup",
            "Fusarium (PH-1) gene ID",
            "Copy number in PH-1",
            "Genome occupancy (18 spp.)",
            "PHI-base lethal evidence",
            "Final label",
        ]
    ]


def draw_panel_c(ax: plt.Axes, table_df: pd.DataFrame) -> None:
    ax.axis("off")
    ax.set_title("Panel C  Representative Real Mapping Examples", loc="left", fontweight="bold", fontsize=12, pad=8)

    column_definitions = [
        ColumnDefinition("Yeast essential gene", title="Yeast gene", width=1.10, group="Yeast", textprops={"fontsize": 8.4}),
        ColumnDefinition("Yeast species", title="Species", width=0.90, group="Yeast", textprops={"fontsize": 8.4}),
        ColumnDefinition("OrthoFinder orthogroup", title="Orthogroup", width=1.00, group="Orthogroup", textprops={"fontsize": 8.4}),
        ColumnDefinition("Fusarium (PH-1) gene ID", title="FGRAMPH1 gene ID", width=1.65, group="Fusarium", textprops={"fontsize": 8.4}),
        ColumnDefinition("Copy number in PH-1", title="PH-1 copy", width=0.82, group="Fusarium", textprops={"fontsize": 8.4}),
        ColumnDefinition("Genome occupancy (18 spp.)", title="Occupancy", width=0.92, group="Fusarium", textprops={"fontsize": 8.4}),
        ColumnDefinition("PHI-base lethal evidence", title="PHI lethal", width=0.88, group="Evidence", textprops={"fontsize": 8.4}),
        ColumnDefinition(
            "Final label",
            title="Final label",
            width=0.95,
            group="Outcome",
            textprops={"fontsize": 8.4, "fontweight": "bold"},
            cmap=lambda val: PALETTE["Positive_Final"] if val == "Positive_Final" else PALETTE["Excluded"],
            text_cmap=lambda val: "#FFFFFF" if val == "Positive_Final" else "#111111",
        ),
    ]

    Table(
        table_df,
        ax=ax,
        index_col="Yeast essential gene",
        column_definitions=column_definitions,
        row_dividers=True,
        odd_row_color="#FBFBFC",
        even_row_color="#FFFFFF",
        col_label_cell_kw={"facecolor": HEADER_BG, "edgecolor": "#D0D4D9", "linewidth": 0.6},
        textprops={"fontsize": 8.4, "color": TEXT},
        row_divider_kw={"linewidth": 0.45, "color": "#D9DDE2"},
        column_border_kw={"linewidth": 0.55, "color": "#D9DDE2"},
    )


def write_palette_mapping() -> None:
    rows = [
        {"class": "scer_only", "hex_color": PALETTE["scer_only"], "source_palette": "palettable.cartocolors.qualitative.Safe_10[0]"},
        {"class": "spom_only", "hex_color": PALETTE["spom_only"], "source_palette": "palettable.cartocolors.qualitative.Safe_10[1]"},
        {"class": "both", "hex_color": PALETTE["both"], "source_palette": "palettable.cartocolors.qualitative.Safe_10[2]"},
        {"class": "Positive_Final", "hex_color": PALETTE["Positive_Final"], "source_palette": "palettable.cartocolors.qualitative.Safe_10[3]"},
        {"class": "Excluded", "hex_color": PALETTE["Excluded"], "source_palette": "manual_neutral_gray"},
    ]
    pd.DataFrame(rows).to_csv(PALETTE_OUT, sep="\t", index=False)


def write_legend(summary: DisplaySummary) -> None:
    legend = (
        "Panel A illustrates the conceptual workflow linking yeast essential support, OrthoFinder orthogroups, "
        "canonical Fusarium PH-1 gene mapping, and the final retained positives.\n\n"
        "Panel B shows the tracked transfer subset as a compact Sankey with lightweight totals above each stage. "
        f"Stage annotations summarize the yeast-support screen universe (n={summary.total_support_universe:,}), the "
        f"orthogroup universe (n={summary.total_orthogroups:,}), the bridge-resolved canonical gene universe "
        f"(n={summary.total_canonical_genes:,}), and the resolved transfer subset entering final display "
        f"(n={summary.resolved_display_universe:,}). Sankey ribbons show only the tracked transfer subset.\n\n"
        "Panel C lists representative real mapping examples derived from the project evidence tables.\n\n"
        "Positive_Final corresponds to the retained high-confidence positive transfer set. "
        "A minor positive class (n=7) was omitted from the Sankey display because of its very small size. "
        "All counts are derived from real project data and canonical bridge outputs."
    )
    LEGEND_OUT.write_text(legend, encoding="utf-8")


def build_figure(display: pd.DataFrame, gene_level: pd.DataFrame, examples_table: pd.DataFrame, summary: DisplaySummary, source_to_final: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(16.4, 10.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.55], height_ratios=[0.88, 1.12])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b, display, summary, source_to_final)
    draw_panel_c(ax_c, examples_table)

    return fig


def main() -> None:
    require_inputs()
    configure_matplotlib()
    gene_level, examples = load_data()
    display, examples_display, summary = prepare_display_tables(gene_level, examples)
    _sankey_data, source_to_final = build_floweaver_sankey(display)
    examples_table = build_panel_c_table(examples_display, gene_level)
    write_palette_mapping()
    write_legend(summary)

    fig = build_figure(display, gene_level, examples_table, summary, source_to_final)
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(SVG_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {SVG_OUT}")
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {LEGEND_OUT}")
    print(f"Wrote {PALETTE_OUT}")


if __name__ == "__main__":
    main()
