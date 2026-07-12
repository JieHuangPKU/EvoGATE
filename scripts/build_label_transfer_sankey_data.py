#!/usr/bin/env python3
"""Build real Sankey/alluvial input tables for Fusarium yeast-transfer label flow."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
BRIDGE_DIR = REPO_ROOT / "data/processed/essential_gene/fgraminearum/bridge"
NEWLABEL_DIR = REPO_ROOT / "data/processed/essential_gene/fgraminearum/newlabel"
DERIVED_DIR = REPO_ROOT / "data/derived_labels"
INTERIM_DIR = REPO_ROOT / "data/interim/protocol_refactor/fgraminearum_label_materialization"
OUTPUT_DIR = REPO_ROOT / "results/label_transfer_sankey"

ORTHO_ROOT = Path(
    "/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/"
    "run_20260405T213342_139369/Results_Apr05/Orthogroups"
)
ORTHOGROUPS_TSV = ORTHO_ROOT / "Orthogroups.tsv"
ORTHOGROUPS_GENECOUNT_TSV = ORTHO_ROOT / "Orthogroups.GeneCount.tsv"
ORTHOGROUPS_SINGLECOPY_TXT = ORTHO_ROOT / "Orthogroups_SingleCopyOrthologues.txt"

PROTEOME_MANIFEST_TSV = DERIVED_DIR / "proteome_manifest.tsv"
YEAST_TRANSFER_TSV = DERIVED_DIR / "ph1_yeast_essential_ortholog_labels.tsv"
YEAST_TRANSFER_SUMMARY_MD = DERIVED_DIR / "ph1_yeast_essential_ortholog_labels.summary.md"
SCER_ESSENTIAL_TSV = DERIVED_DIR / "yeast_essential/Scerevisiae.essential_genes.tsv"
SPOM_ESSENTIAL_TSV = DERIVED_DIR / "yeast_essential/Spombe.essential_genes.tsv"

PROTEIN_BRIDGE_TSV = BRIDGE_DIR / "protein_to_canonical_bridge.tsv"
HIGH_CONFIDENCE_TSV = BRIDGE_DIR / "high_confidence_yeast_transfer_candidates.tsv"
BRIDGE_SUMMARY_TSV = BRIDGE_DIR / "bridge_summary.tsv"
BRIDGE_MANIFEST_TSV = BRIDGE_DIR / "bridge_source_manifest.tsv"
UNRESOLVED_HIGH_TSV = BRIDGE_DIR / "unresolved_high_confidence_ids.tsv"

LABELS_TSV = NEWLABEL_DIR / "labels.tsv"
POSITIVE_TSV = NEWLABEL_DIR / "positive_genes.tsv"
LABEL_AUDIT_TSV = NEWLABEL_DIR / "label_construction_audit.tsv"
NEWLABEL_SOURCE_MANIFEST_TSV = NEWLABEL_DIR / "source_manifest.tsv"
NEWLABEL_SUMMARY_TSV = NEWLABEL_DIR / "summary.tsv"
LETHAL_POSITIVE_TSV = INTERIM_DIR / "lethal_positive_gene_list.tsv"

PH1_ACCESSION = "GCF_000240135.3"
VALID_SUPPORT_TYPES = {"scer_only", "spom_only", "both"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Output directory for Sankey-ready TSVs and README.",
    )
    return parser.parse_args()


def log(message: str) -> None:
    print(f"[build_label_transfer_sankey_data] {message}")


def require_files(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input files:\n- " + "\n- ".join(missing))


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def read_singlecopy_set(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def split_items(cell: str) -> list[str]:
    text = str(cell).strip()
    if not text:
        return []
    return [item.strip() for item in text.split(",") if item.strip()]


def first_nonempty(*values: str) -> str:
    for value in values:
        if str(value).strip():
            return str(value).strip()
    return ""


def classify_copy_status(row: pd.Series) -> str:
    is_exact_single_copy = str(row["is_exact_single_copy_core_18"]).strip() == "1"
    is_strict_core = str(row["is_strict_core_18"]).strip() == "1"
    occupancy = float(row["fg_occupancy_18"])
    single_copy_fraction = float(row["fg_single_copy_fraction_18"])
    max_copy = float(row["fg_max_copy_18"])
    if is_exact_single_copy:
        return "exact_single_copy_core"
    if is_strict_core and max_copy <= 2 and single_copy_fraction >= 0.90:
        return "near_single_copy_core"
    if is_strict_core:
        return "strict_core_multicopy"
    if occupancy >= 0.80 and single_copy_fraction >= 0.75:
        return "high_occupancy_mostly_single_copy"
    return "noncore_or_multicopy"


def load_orthogroup_context() -> tuple[pd.DataFrame, pd.DataFrame]:
    proteome_manifest = read_tsv(PROTEOME_MANIFEST_TSV)
    anchor_rows = proteome_manifest.loc[proteome_manifest["anchor_role"].eq("anchor")]
    if anchor_rows.empty:
        raise ValueError(f"No anchor row found in {PROTEOME_MANIFEST_TSV}")
    anchor_accession = anchor_rows["sample_id"].iloc[0]
    if anchor_accession != PH1_ACCESSION:
        raise ValueError(
            f"Unexpected PH-1 OrthoFinder accession: expected {PH1_ACCESSION}, observed {anchor_accession}"
        )

    orthogroups = read_tsv(ORTHOGROUPS_TSV)[["Orthogroup", PH1_ACCESSION, "Scerevisiae", "Spombe"]].rename(
        columns={
            "Orthogroup": "orthogroup_id",
            PH1_ACCESSION: "ph1_membership_cell",
            "Scerevisiae": "scer_membership_cell",
            "Spombe": "spom_membership_cell",
        }
    )
    gene_counts = read_tsv(ORTHOGROUPS_GENECOUNT_TSV).rename(columns={"Orthogroup": "orthogroup_id"})
    singlecopy_set = read_singlecopy_set(ORTHOGROUPS_SINGLECOPY_TXT)
    orthogroups["single_copy_orthofinder_listed"] = orthogroups["orthogroup_id"].isin(singlecopy_set).map(
        {True: "yes", False: "no"}
    )
    return orthogroups, gene_counts


def build_gene_level_table() -> tuple[pd.DataFrame, dict[str, object]]:
    yeast = read_tsv(YEAST_TRANSFER_TSV)
    bridge = read_tsv(PROTEIN_BRIDGE_TSV).rename(
        columns={
            "resolved_canonical_gene_id": "canonical_fusarium_gene_id",
        }
    )
    high_conf = read_tsv(HIGH_CONFIDENCE_TSV)
    positives = read_tsv(POSITIVE_TSV)
    audit = read_tsv(LABEL_AUDIT_TSV)
    lethal = read_tsv(LETHAL_POSITIVE_TSV)
    orthogroups, gene_counts = load_orthogroup_context()

    yeast = yeast.loc[yeast["has_any_yeast_essential"].eq("1")].copy()
    invalid_supports = sorted(set(yeast["yeast_essential_support_class"]) - VALID_SUPPORT_TYPES - {""})
    if invalid_supports:
        raise ValueError(f"Unexpected yeast support types: {invalid_supports}")

    bridge_cols = [
        "source_protein_id",
        "canonical_fusarium_gene_id",
        "bridge_status",
        "bridge_method",
        "mapping_confidence",
    ]
    merged = yeast.merge(bridge[bridge_cols], left_on="ph1_gene_id", right_on="source_protein_id", how="left")
    merged = merged.merge(
        orthogroups,
        on="orthogroup_id",
        how="left",
        validate="many_to_one",
    )
    merged = merged.merge(
        gene_counts[["orthogroup_id", "Total", PH1_ACCESSION]].rename(
            columns={"Total": "orthogroup_total_member_count", PH1_ACCESSION: "ph1_copy_count_from_membership"}
        ),
        on="orthogroup_id",
        how="left",
        validate="many_to_one",
    )

    merged["yeast_support_source"] = merged["yeast_essential_support_class"]
    merged["support_type"] = merged["yeast_essential_support_class"]
    merged["bridge_resolved"] = (
        merged["bridge_status"].eq("resolved") & merged["canonical_fusarium_gene_id"].astype(str).ne("")
    ).map({True: "yes", False: "no"})
    merged["passes_fusarium_filtering"] = merged["weak_positive_confidence"].ne("none").map(
        {True: "yes", False: "no"}
    )
    merged["occupancy_count"] = merged["fg_presence_count_18"]
    merged["occupancy_fraction"] = merged["fg_occupancy_18"]
    merged["copy_status"] = merged.apply(classify_copy_status, axis=1)
    merged["single_copy_or_near_single_copy"] = merged["copy_status"].isin(
        {"exact_single_copy_core", "near_single_copy_core", "high_occupancy_mostly_single_copy"}
    ).map({True: "yes", False: "no"})
    merged["orthogroup_member_count_total"] = merged["orthogroup_total_member_count"]
    merged["bridge_evidence_path"] = str(PROTEIN_BRIDGE_TSV.relative_to(REPO_ROOT))
    merged["supporting_xp_ids"] = merged["ph1_gene_id"]
    merged["supporting_scer_gene_ids"] = merged["scer_essential_gene_ids"]
    merged["supporting_spom_gene_ids"] = merged["spom_essential_gene_ids"]

    high_conf_set = set(high_conf["canonical_gene_id"])
    positive_set = set(positives["canonical_gene_id"])
    final_label_set = set(audit["canonical_gene_id"])
    lethal_set = set(lethal["canonical_gene_id"])
    audit_lookup = audit.set_index("canonical_gene_id").to_dict(orient="index")

    def final_class(row: pd.Series) -> str:
        canonical_gene = row["canonical_fusarium_gene_id"]
        if row["bridge_resolved"] != "yes":
            return "Excluded"
        if canonical_gene in positive_set and canonical_gene in high_conf_set:
            return "Positive_High"
        if canonical_gene in positive_set:
            return "Positive"
        return "Excluded"

    merged["in_high_confidence_set"] = merged["canonical_fusarium_gene_id"].isin(high_conf_set).map(
        {True: "yes", False: "no"}
    )
    merged["in_final_positive_set"] = merged["canonical_fusarium_gene_id"].isin(positive_set).map(
        {True: "yes", False: "no"}
    )
    merged["in_final_label_set"] = merged["canonical_fusarium_gene_id"].isin(final_label_set).map(
        {True: "yes", False: "no"}
    )
    merged["is_lethal_supported_positive"] = merged["canonical_fusarium_gene_id"].isin(lethal_set).map(
        {True: "yes", False: "no"}
    )
    merged["final_label_class"] = merged.apply(final_class, axis=1)
    merged["final_label_bucket"] = merged["canonical_fusarium_gene_id"].map(
        lambda gene: audit_lookup.get(gene, {}).get("construction_bucket", "") if gene else ""
    )

    group_cols = [
        "yeast_support_source",
        "support_type",
        "orthogroup_id",
        "canonical_fusarium_gene_id",
        "bridge_resolved",
        "bridge_method",
        "mapping_confidence",
        "bridge_evidence_path",
        "passes_fusarium_filtering",
        "weak_positive_confidence",
        "occupancy_count",
        "occupancy_fraction",
        "fg_mean_copy_18",
        "fg_median_copy_18",
        "fg_max_copy_18",
        "fg_single_copy_fraction_18",
        "copy_status",
        "single_copy_or_near_single_copy",
        "single_copy_orthofinder_listed",
        "orthogroup_member_count_total",
        "ph1_copy_count_from_membership",
        "in_final_positive_set",
        "in_high_confidence_set",
        "in_final_label_set",
        "is_lethal_supported_positive",
        "final_label_class",
        "final_label_bucket",
    ]
    aggregated = (
        merged.groupby(group_cols, dropna=False, as_index=False)
        .agg(
            ph1_gene_ids=("ph1_gene_id", lambda values: ";".join(sorted(set(map(str, values))))),
            scer_essential_gene_ids=("scer_essential_gene_ids", lambda values: ";".join(sorted({v for v in values if v}))),
            spom_essential_gene_ids=("spom_essential_gene_ids", lambda values: ";".join(sorted({v for v in values if v}))),
            source_row_count=("ph1_gene_id", "size"),
        )
        .sort_values(["orthogroup_id", "canonical_fusarium_gene_id", "yeast_support_source"], kind="stable")
        .reset_index(drop=True)
    )

    aggregated["canonical_fusarium_gene_id"] = aggregated["canonical_fusarium_gene_id"].fillna("")
    aggregated["yeast_gene"] = aggregated.apply(
        lambda row: first_nonempty(
            row["scer_essential_gene_ids"].split(";")[0] if row["scer_essential_gene_ids"] else "",
            row["spom_essential_gene_ids"].split(";")[0] if row["spom_essential_gene_ids"] else "",
        ),
        axis=1,
    )
    aggregated["yeast_species"] = aggregated["yeast_support_source"].map(
        {"scer_only": "scer", "spom_only": "spom", "both": "both"}
    )

    metadata = {
        "support_rows_total": int(len(yeast)),
        "gene_level_rows_total": int(len(aggregated)),
        "bridge_resolved_rows": int((aggregated["bridge_resolved"] == "yes").sum()),
        "bridge_unresolved_rows": int((aggregated["bridge_resolved"] == "no").sum()),
        "bridge_resolved_unique_genes": int(
            aggregated.loc[aggregated["bridge_resolved"].eq("yes"), "canonical_fusarium_gene_id"].nunique()
        ),
        "positive_high_rows": int((aggregated["final_label_class"] == "Positive_High").sum()),
        "positive_rows": int((aggregated["final_label_class"] == "Positive").sum()),
        "excluded_rows": int((aggregated["final_label_class"] == "Excluded").sum()),
    }
    return aggregated, metadata


def build_edge_tables(gene_level: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sankey_rows = gene_level.loc[gene_level["bridge_resolved"].eq("yes")].copy()
    edges = []

    def add_edges(source_layer: str, source_col: str, target_layer: str, target_col: str) -> None:
        grouped = (
            sankey_rows.groupby([source_col, target_col], as_index=False)
            .size()
            .rename(columns={source_col: "source_node", target_col: "target_node", "size": "count"})
        )
        grouped["source_layer"] = source_layer
        grouped["target_layer"] = target_layer
        edges.append(grouped[["source_layer", "source_node", "target_layer", "target_node", "count"]])

    add_edges("Layer1_YeastSupportSource", "yeast_support_source", "Layer2_Orthogroup", "orthogroup_id")
    add_edges("Layer2_Orthogroup", "orthogroup_id", "Layer3_MappedFusariumGene", "canonical_fusarium_gene_id")
    add_edges("Layer3_MappedFusariumGene", "canonical_fusarium_gene_id", "Layer4_FinalLabelClass", "final_label_class")

    edge_table = pd.concat(edges, ignore_index=True)
    stage_counts = pd.concat(
        [
            sankey_rows.groupby("yeast_support_source", as_index=False)
            .size()
            .rename(columns={"yeast_support_source": "node", "size": "count"})
            .assign(layer="Layer1_YeastSupportSource"),
            sankey_rows.groupby("orthogroup_id", as_index=False)
            .size()
            .rename(columns={"orthogroup_id": "node", "size": "count"})
            .assign(layer="Layer2_Orthogroup"),
            sankey_rows.groupby("canonical_fusarium_gene_id", as_index=False)
            .size()
            .rename(columns={"canonical_fusarium_gene_id": "node", "size": "count"})
            .assign(layer="Layer3_MappedFusariumGene"),
            sankey_rows.groupby("final_label_class", as_index=False)
            .size()
            .rename(columns={"final_label_class": "node", "size": "count"})
            .assign(layer="Layer4_FinalLabelClass"),
        ],
        ignore_index=True,
    )[["layer", "node", "count"]]

    return edge_table.sort_values(["source_layer", "source_node", "target_node"], kind="stable"), stage_counts.sort_values(
        ["layer", "node"], kind="stable"
    )


def build_examples(gene_level: pd.DataFrame) -> pd.DataFrame:
    resolved = gene_level.loc[gene_level["bridge_resolved"].eq("yes")].copy()
    resolved = resolved.sort_values(
        ["final_label_class", "yeast_support_source", "occupancy_fraction", "canonical_fusarium_gene_id"],
        ascending=[True, True, False, True],
        kind="stable",
    )

    selections = []
    for label_class, cap in [("Positive_High", 8), ("Positive", 4), ("Excluded", 8)]:
        subset = resolved.loc[resolved["final_label_class"].eq(label_class)].head(cap)
        selections.append(subset)

    examples = pd.concat(selections, ignore_index=True).drop_duplicates(
        subset=["orthogroup_id", "canonical_fusarium_gene_id"], keep="first"
    )
    examples = examples[
        [
            "yeast_gene",
            "yeast_species",
            "support_type",
            "orthogroup_id",
            "canonical_fusarium_gene_id",
            "occupancy_fraction",
            "copy_status",
            "bridge_method",
            "final_label_class",
        ]
    ].rename(columns={"support_type": "support_type"})
    return examples


def build_readme(
    output_dir: Path,
    gene_level: pd.DataFrame,
    metadata: dict[str, object],
    edges: pd.DataFrame,
    stage_counts: pd.DataFrame,
) -> str:
    high_positive = int((gene_level["final_label_class"] == "Positive_High").sum())
    positive = int((gene_level["final_label_class"] == "Positive").sum())
    excluded = int((gene_level["final_label_class"] == "Excluded").sum())
    unresolved = int((gene_level["bridge_resolved"] == "no").sum())
    resolved = int((gene_level["bridge_resolved"] == "yes").sum())
    final_positive_gene_count = int(gene_level.loc[gene_level["in_final_positive_set"].eq("yes"), "canonical_fusarium_gene_id"].nunique())
    high_conf_gene_count = int(gene_level.loc[gene_level["in_high_confidence_set"].eq("yes"), "canonical_fusarium_gene_id"].nunique())

    used_files = [
        str(YEAST_TRANSFER_TSV.relative_to(REPO_ROOT)),
        str(PROTEIN_BRIDGE_TSV.relative_to(REPO_ROOT)),
        str(HIGH_CONFIDENCE_TSV.relative_to(REPO_ROOT)),
        str(BRIDGE_SUMMARY_TSV.relative_to(REPO_ROOT)),
        str(BRIDGE_MANIFEST_TSV.relative_to(REPO_ROOT)),
        str(UNRESOLVED_HIGH_TSV.relative_to(REPO_ROOT)),
        str(LABELS_TSV.relative_to(REPO_ROOT)),
        str(POSITIVE_TSV.relative_to(REPO_ROOT)),
        str(LABEL_AUDIT_TSV.relative_to(REPO_ROOT)),
        str(NEWLABEL_SOURCE_MANIFEST_TSV.relative_to(REPO_ROOT)),
        str(NEWLABEL_SUMMARY_TSV.relative_to(REPO_ROOT)),
        str(LETHAL_POSITIVE_TSV.relative_to(REPO_ROOT)),
        str(PROTEOME_MANIFEST_TSV.relative_to(REPO_ROOT)),
        str(YEAST_TRANSFER_SUMMARY_MD.relative_to(REPO_ROOT)),
        str(SCER_ESSENTIAL_TSV.relative_to(REPO_ROOT)),
        str(SPOM_ESSENTIAL_TSV.relative_to(REPO_ROOT)),
        str(ORTHOGROUPS_TSV),
        str(ORTHOGROUPS_GENECOUNT_TSV),
        str(ORTHOGROUPS_SINGLECOPY_TXT),
    ]

    readme = f"""# Label Transfer Sankey Data

## Actual Inputs Used
- Canonical yeast-transfer source: `{YEAST_TRANSFER_TSV.relative_to(REPO_ROOT)}`
- Canonical bridge source: `{HIGH_CONFIDENCE_TSV.relative_to(REPO_ROOT)}` and `{PROTEIN_BRIDGE_TSV.relative_to(REPO_ROOT)}`
- Canonical final label source: `{LABEL_AUDIT_TSV.relative_to(REPO_ROOT)}`
- Canonical final positive subset: `{POSITIVE_TSV.relative_to(REPO_ROOT)}`
- Canonical final label audit source: `{LABEL_AUDIT_TSV.relative_to(REPO_ROOT)}`
- Lethal positive provenance source: `{LETHAL_POSITIVE_TSV.relative_to(REPO_ROOT)}`
- Orthogroup membership source: `{ORTHOGROUPS_TSV}`
- Orthogroup gene-count source: `{ORTHOGROUPS_GENECOUNT_TSV}`
- OrthoFinder single-copy list: `{ORTHOGROUPS_SINGLECOPY_TXT}`
- OrthoFinder source provenance summary used for path validation: `{YEAST_TRANSFER_SUMMARY_MD.relative_to(REPO_ROOT)}`

### Full File List
{chr(10).join(f"- `{path}`" for path in used_files)}

## Canonical Source Decisions
- Canonical final label source was assigned to `{LABEL_AUDIT_TSV.relative_to(REPO_ROOT)}` because it is the full materialized `newlabel` construction audit spanning all retained and excluded genes in the final label universe, while `{LABELS_TSV.relative_to(REPO_ROOT)}` is a narrower downstream label manifest.
- Canonical orthogroup membership source was assigned to `{ORTHOGROUPS_TSV}` because the repo-local summary `{YEAST_TRANSFER_SUMMARY_MD.relative_to(REPO_ROOT)}` records that exact OrthoFinder results directory as the upstream source for `ph1_yeast_essential_ortholog_labels.tsv`.
- Canonical XP to `fgraminearum::FGRAMPH1_*` mapping source was assigned to `{PROTEIN_BRIDGE_TSV.relative_to(REPO_ROOT)}` plus the protocolized positive subset `{HIGH_CONFIDENCE_TSV.relative_to(REPO_ROOT)}`.

## Layer Definitions And Counting Unit
- Counting unit for `sankey_gene_level_long.tsv`: one aggregated real evidence record per `(yeast_support_source, orthogroup_id, canonical_fusarium_gene_id)` after joining the upstream yeast-transfer table to the protocolized bridge. Multiple `XP_*` rows collapsing to the same orthogroup and canonical gene are aggregated, with `source_row_count` preserving the original row multiplicity.
- Layer 1 `Yeast support source`: taken directly from `yeast_essential_support_class` in `{YEAST_TRANSFER_TSV.relative_to(REPO_ROOT)}` and restricted to `scer_only`, `spom_only`, `both`.
- Layer 2 `Orthogroups with essential support`: `orthogroup_id` from `{YEAST_TRANSFER_TSV.relative_to(REPO_ROOT)}`; occupancy and copy metrics come from the same table (`fg_presence_count_18`, `fg_occupancy_18`, `fg_single_copy_fraction_18`, `is_strict_core_18`, `is_exact_single_copy_core_18`) and were cross-linked to `{ORTHOGROUPS_TSV}` / `{ORTHOGROUPS_GENECOUNT_TSV}` / `{ORTHOGROUPS_SINGLECOPY_TXT}`.
- Layer 3 `Mapped Fusarium genes`: `canonical_fusarium_gene_id` from the protocolized bridge join to `{PROTEIN_BRIDGE_TSV.relative_to(REPO_ROOT)}`; only bridge-resolved rows can enter the drawn Sankey edges because layer 3 is constrained to canonical `fgraminearum::FGRAMPH1_*` IDs.
- Layer 4 `Final label class`: assigned by joining canonical genes to the materialized newlabel outputs.

## Final Class Rules
- `Positive_High`: bridge-resolved canonical gene is present in `{POSITIVE_TSV.relative_to(REPO_ROOT)}` and also present in `{HIGH_CONFIDENCE_TSV.relative_to(REPO_ROOT)}`.
- `Positive`: bridge-resolved canonical gene is present in `{POSITIVE_TSV.relative_to(REPO_ROOT)}` but not present in the protocolized high-confidence transfer set. In practice these are yeast-supported mapped genes retained in final positives through non-high-confidence logic, typically lethal-backed positives.
- `Excluded`: everything else in the yeast-support pool. This includes bridge-resolved genes that ended up in final negatives or outside the final label table, plus unresolved bridge rows in the long table.

## Fusarium Filtering Definition
- `passes_fusarium_filtering = yes` means `weak_positive_confidence != none` in `{YEAST_TRANSFER_TSV.relative_to(REPO_ROOT)}`.
- `passes_fusarium_filtering = no` means `weak_positive_confidence == none`.
- This script does not re-infer the filter. It reuses the upstream materialized confidence field exactly.

## Copy-Status Definition
- `exact_single_copy_core`: `is_exact_single_copy_core_18 == 1`
- `near_single_copy_core`: strict core plus `fg_max_copy_18 <= 2` and `fg_single_copy_fraction_18 >= 0.90`
- `strict_core_multicopy`: strict core but not exact/near single-copy
- `high_occupancy_mostly_single_copy`: `fg_occupancy_18 >= 0.80` and `fg_single_copy_fraction_18 >= 0.75`
- `noncore_or_multicopy`: remaining rows

## Output Summary
- Upstream yeast-supported rows (`has_any_yeast_essential == 1`): {metadata["support_rows_total"]}
- Aggregated gene-level evidence rows: {metadata["gene_level_rows_total"]}
- Bridge-resolved aggregated rows: {resolved}
- Bridge-unresolved aggregated rows: {unresolved}
- `Positive_High` aggregated rows: {high_positive}
- `Positive` aggregated rows: {positive}
- `Excluded` aggregated rows: {excluded}
- Unique canonical genes in final positive set within this Sankey pool: {final_positive_gene_count}
- Unique canonical genes in protocolized high-confidence set within this Sankey pool: {high_conf_gene_count}
- Sankey edge rows written: {len(edges)}
- Stage-count rows written: {len(stage_counts)}

## Caveats
- The repo does not carry a repo-local copy of `Orthogroups.tsv`; the canonical membership table was read directly from the upstream OrthoFinder path documented in `{YEAST_TRANSFER_SUMMARY_MD.relative_to(REPO_ROOT)}`.
- Unresolved bridge rows are retained in `sankey_gene_level_long.tsv` with `bridge_resolved = no` and `final_label_class = Excluded`, but they are not drawn into `sankey_aggregated_edges.tsv` because layer 3 is restricted to canonical `fgraminearum::FGRAMPH1_*` gene IDs.
- `Orthogroups_SingleCopyOrthologues.txt` is used as supplemental context only. The primary occupancy and copy metrics still come from the materialized yeast-transfer table so the Sankey remains traceable to the exact project pipeline outputs.
"""
    return readme


def main() -> None:
    args = parse_args()
    require_files(
        [
            BRIDGE_DIR,
            NEWLABEL_DIR,
            YEAST_TRANSFER_TSV,
            PROTEOME_MANIFEST_TSV,
            PROTEIN_BRIDGE_TSV,
            HIGH_CONFIDENCE_TSV,
            BRIDGE_SUMMARY_TSV,
            BRIDGE_MANIFEST_TSV,
            UNRESOLVED_HIGH_TSV,
            LABELS_TSV,
            POSITIVE_TSV,
            LABEL_AUDIT_TSV,
            NEWLABEL_SOURCE_MANIFEST_TSV,
            NEWLABEL_SUMMARY_TSV,
            LETHAL_POSITIVE_TSV,
            YEAST_TRANSFER_SUMMARY_MD,
            SCER_ESSENTIAL_TSV,
            SPOM_ESSENTIAL_TSV,
            ORTHOGROUPS_TSV,
            ORTHOGROUPS_GENECOUNT_TSV,
            ORTHOGROUPS_SINGLECOPY_TXT,
        ]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log("Loading canonical bridge, orthogroup, and final label assets")
    gene_level, metadata = build_gene_level_table()
    log(f"Built {metadata['gene_level_rows_total']} aggregated gene-level records")

    log("Building Sankey edge and stage-count tables from bridge-resolved rows")
    edge_table, stage_counts = build_edge_tables(gene_level)

    log("Selecting representative real examples")
    examples = build_examples(gene_level)

    gene_level_path = args.output_dir / "sankey_gene_level_long.tsv"
    edge_path = args.output_dir / "sankey_aggregated_edges.tsv"
    stage_count_path = args.output_dir / "sankey_stage_counts.tsv"
    example_path = args.output_dir / "representative_mapping_examples.tsv"
    readme_path = args.output_dir / "README.md"

    gene_level.to_csv(gene_level_path, sep="\t", index=False)
    edge_table.to_csv(edge_path, sep="\t", index=False)
    stage_counts.to_csv(stage_count_path, sep="\t", index=False)
    examples.to_csv(example_path, sep="\t", index=False)
    readme_path.write_text(build_readme(args.output_dir, gene_level, metadata, edge_table, stage_counts), encoding="utf-8")

    log(f"Wrote {gene_level_path.relative_to(REPO_ROOT)}")
    log(f"Wrote {edge_path.relative_to(REPO_ROOT)}")
    log(f"Wrote {stage_count_path.relative_to(REPO_ROOT)}")
    log(f"Wrote {example_path.relative_to(REPO_ROOT)}")
    log(f"Wrote {readme_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
