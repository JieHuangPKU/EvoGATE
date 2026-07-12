import argparse
import hashlib
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the protocolized Fusarium XP-to-canonical bridge for newlabel")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    seq: list[str] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_ncbi_records(path: str | Path) -> pd.DataFrame:
    rows = []
    for header, sequence in read_fasta(path):
        xp_match = re.search(r"(XP_[0-9]+\.[0-9]+)", header)
        if xp_match is None:
            continue
        fgsg_match = re.search(r"(FGSG_[0-9]+)", header)
        rows.append(
            {
                "source_protein_id": xp_match.group(1),
                "source_header": header,
                "header_fgsg_id": fgsg_match.group(1) if fgsg_match else "",
                "sequence_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
                "sequence_length": len(sequence),
            }
        )
    return pd.DataFrame(rows)


def parse_legacy_proteome(path: str | Path) -> pd.DataFrame:
    rows = []
    for header, sequence in read_fasta(path):
        parts = header.split("|")
        legacy_acc = parts[1] if len(parts) > 1 else header.split()[0]
        gn_match = re.search(r"\bGN=([^ ]+)", header)
        rows.append(
            {
                "legacy_protein_accession": legacy_acc,
                "legacy_header": header,
                "legacy_gene_symbol": gn_match.group(1) if gn_match else "",
                "sequence_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
                "sequence_length": len(sequence),
            }
        )
    return pd.DataFrame(rows)


def build_fgsg_to_canonical_table(unified_map_path: str | Path) -> pd.DataFrame:
    rows = []
    master = read_tsv("data/interim/protocol_refactor/master_evidence_table.preliminary.tsv")
    fg_master = master[master["species"].eq("fgraminearum")].copy()
    fg_master = fg_master[
        fg_master["raw_gene_id"].astype(str).str.startswith("FGSG_")
        & fg_master["canonical_gene_id"].astype(str).str.startswith("fgraminearum::FGRAMPH1_")
    ][["raw_gene_id", "canonical_gene_id"]].drop_duplicates()
    for _, row in fg_master.iterrows():
        rows.append(
            {
                "source_id": row["raw_gene_id"],
                "source_type": "FGSG",
                "canonical_gene_id": row["canonical_gene_id"],
                "evidence_path": "data/interim/protocol_refactor/master_evidence_table.preliminary.tsv",
                "evidence_kind": "master_evidence_raw_gene_id",
            }
        )

    unified = read_tsv(unified_map_path)
    for _, row in unified.iterrows():
        canonical_gene = str(row.get("ph1_canonical_gene_id", "")).strip()
        if not canonical_gene.startswith("FGRAMPH1_"):
            continue
        canonical_gene = f"fgraminearum::{canonical_gene}"
        fgsg_ids = [item.strip() for item in str(row.get("fgsg_id", "")).split(";") if item.strip()]
        for source_id in fgsg_ids:
            if not source_id.startswith("FGSG_"):
                continue
            rows.append(
                {
                    "source_id": source_id,
                    "source_type": "FGSG",
                    "canonical_gene_id": canonical_gene,
                    "evidence_path": str(unified_map_path),
                    "evidence_kind": "unified_map_fgsg_id",
                }
            )

    mapping_specs = [
        ("data/processed/PPI/fgraminearum/string_id_mapping.tsv", "matched_aliases", "mapped_gene_id", "processed_ppi_alias_mapping"),
        ("data/processed/EXP/fgraminearum/exp_id_mapping.tsv", "raw_probe_id", "mapped_gene_id", "processed_expression_id_mapping"),
        ("data/processed/LC/fgraminearum/subloc_id_mapping.tsv", "raw_id", "mapped_gene_id", "processed_sublocalization_id_mapping"),
    ]
    for path, raw_col, mapped_col, kind in mapping_specs:
        frame = read_tsv(path)
        for _, row in frame.iterrows():
            raw_text = str(row[raw_col])
            source_ids = set(re.findall(r"FGSG_[0-9]+", raw_text))
            if raw_col == "raw_probe_id" and raw_text.startswith("FGSG_"):
                source_ids.add(raw_text)
            mapped_gene = str(row[mapped_col]).strip()
            if not mapped_gene.startswith("FGRAMPH1_"):
                continue
            canonical_gene = f"fgraminearum::{mapped_gene}"
            for source_id in source_ids:
                rows.append(
                    {
                        "source_id": source_id,
                        "source_type": "FGSG",
                        "canonical_gene_id": canonical_gene,
                        "evidence_path": path,
                        "evidence_kind": kind,
                    }
                )

    out = pd.DataFrame(rows).drop_duplicates()
    return out.sort_values(["source_id", "canonical_gene_id", "evidence_kind"], kind="stable").reset_index(drop=True)


def build_xp_to_canonical_table(path: str | Path) -> dict[str, list[str]]:
    unified = read_tsv(path)
    mapping: dict[str, set[str]] = defaultdict(set)
    for _, row in unified.iterrows():
        canonical_gene = str(row.get("ph1_canonical_gene_id", "")).strip()
        if not canonical_gene.startswith("FGRAMPH1_"):
            continue
        for protein_id in [item.strip() for item in str(row.get("ncbi_protein_ids", "")).split(";") if item.strip()]:
            mapping[protein_id].add(f"fgraminearum::{canonical_gene}")
    return {key: sorted(values) for key, values in mapping.items()}


def build_bridge(config: dict[str, Any]) -> dict[str, pd.DataFrame | str]:
    paths = config["paths"]
    bridge_dir = Path(paths["bridge_dir"])
    bridge_dir.mkdir(parents=True, exist_ok=True)

    ncbi_df = parse_ncbi_records(paths["anchor_ncbi_protein_fasta"])
    legacy_df = parse_legacy_proteome(paths["ph1_legacy_protein_fasta"])
    mapping_tab = read_tsv(paths["ph1_legacy_mapping_tab"])
    xp_to_canonical = build_xp_to_canonical_table(paths["ph1_unified_id_map"])
    mapping_tab.columns = ["query", "canonical_gene_id"]
    mapping_tab["legacy_protein_accession"] = mapping_tab["query"].astype(str).str.replace("229533.", "", regex=False)

    seq_to_legacy = legacy_df[["sequence_sha256", "legacy_protein_accession"]].drop_duplicates()
    seq_to_legacy_counts = seq_to_legacy.groupby("sequence_sha256")["legacy_protein_accession"].nunique().rename("legacy_seq_match_count")
    ncbi_bridge = ncbi_df.merge(seq_to_legacy, on="sequence_sha256", how="left")
    ncbi_bridge = ncbi_bridge.merge(seq_to_legacy_counts, on="sequence_sha256", how="left")
    ncbi_bridge["legacy_seq_match_count"] = pd.to_numeric(ncbi_bridge["legacy_seq_match_count"], errors="coerce").fillna(0).astype(int)
    ncbi_bridge = ncbi_bridge.merge(mapping_tab[["legacy_protein_accession", "canonical_gene_id"]], on="legacy_protein_accession", how="left")
    ncbi_bridge["canonical_gene_id"] = ncbi_bridge["canonical_gene_id"].map(
        lambda value: f"fgraminearum::{value}" if str(value).startswith("FGRAMPH1_") else ""
    )

    fgsg_to_canonical = build_fgsg_to_canonical_table(paths["ph1_unified_id_map"])
    fgsg_unique = (
        fgsg_to_canonical.groupby("source_id")["canonical_gene_id"]
        .agg(lambda values: sorted(set(value for value in values if str(value).startswith("fgraminearum::FGRAMPH1_"))))
        .reset_index(name="fgsg_candidates")
    )
    fgsg_unique["fgsg_candidate_count"] = fgsg_unique["fgsg_candidates"].map(len)
    fgsg_unique["fgsg_canonical_gene_id"] = fgsg_unique["fgsg_candidates"].map(lambda values: values[0] if len(values) == 1 else "")
    ncbi_bridge = ncbi_bridge.merge(fgsg_unique[["source_id", "fgsg_candidate_count", "fgsg_canonical_gene_id"]], left_on="header_fgsg_id", right_on="source_id", how="left")
    ncbi_bridge["fgsg_candidate_count"] = pd.to_numeric(ncbi_bridge["fgsg_candidate_count"], errors="coerce").fillna(0).astype(int)
    ncbi_bridge["fgsg_canonical_gene_id"] = ncbi_bridge["fgsg_canonical_gene_id"].fillna("")
    ncbi_bridge["unified_canonical_gene_id"] = ncbi_bridge["source_protein_id"].map(
        lambda value: xp_to_canonical[value][0] if value in xp_to_canonical and len(xp_to_canonical[value]) == 1 else ""
    )
    ncbi_bridge["unified_candidate_count"] = ncbi_bridge["source_protein_id"].map(
        lambda value: len(xp_to_canonical.get(value, []))
    )

    def resolve_row(row: pd.Series) -> tuple[str, str, str]:
        seq_canonical = str(row["canonical_gene_id"]).strip()
        fgsg_canonical = str(row["fgsg_canonical_gene_id"]).strip()
        unified_canonical = str(row["unified_canonical_gene_id"]).strip()
        seq_count = int(row["legacy_seq_match_count"])
        fgsg_count = int(row["fgsg_candidate_count"])
        unified_count = int(row["unified_candidate_count"])

        if unified_count > 1:
            return "", "ambiguous", "unified_map_multiple_candidate_targets"

        if seq_canonical and fgsg_canonical:
            if seq_canonical == fgsg_canonical:
                return seq_canonical, "resolved", "exact_sequence_plus_fgsg_support"
            if unified_canonical and unified_canonical == seq_canonical:
                return seq_canonical, "resolved", "exact_sequence_plus_unified_support"
            if unified_canonical and unified_canonical == fgsg_canonical:
                return unified_canonical, "resolved", "unified_plus_fgsg_support"
            return "", "ambiguous", "sequence_vs_fgsg_conflict"
        if unified_canonical and seq_canonical:
            if unified_canonical == seq_canonical:
                return seq_canonical, "resolved", "exact_sequence_plus_unified_support"
            return unified_canonical, "resolved", "unified_map_preferred_over_legacy_sequence_bridge"
        if unified_canonical and fgsg_canonical:
            if unified_canonical == fgsg_canonical:
                return unified_canonical, "resolved", "unified_plus_fgsg_support"
            return unified_canonical, "resolved", "unified_map_preferred_over_fgsg_header_bridge"
        if unified_canonical:
            return unified_canonical, "resolved", "ncbi_protein_id_via_unified_map"
        if seq_canonical:
            return seq_canonical, "resolved", "exact_sequence_to_legacy_protein"
        if fgsg_canonical:
            return fgsg_canonical, "resolved", "fgsg_header_repo_local_mapping"
        if seq_count > 1 or fgsg_count > 1:
            return "", "ambiguous", "multiple_candidate_targets"
        if row["header_fgsg_id"]:
            return "", "unresolved", "header_fgsg_not_mapped_to_canonical"
        if seq_count == 0:
            return "", "unresolved", "no_sequence_or_fgsg_bridge"
        return "", "unresolved", "unclassified"

    resolved = ncbi_bridge.apply(resolve_row, axis=1, result_type="expand")
    resolved.columns = ["resolved_canonical_gene_id", "bridge_status", "bridge_method"]
    bridge_df = pd.concat([ncbi_bridge, resolved], axis=1)
    bridge_df["mapping_confidence"] = bridge_df["bridge_method"].map(
        {
            "exact_sequence_plus_fgsg_support": "exact",
            "exact_sequence_to_legacy_protein": "exact",
            "fgsg_header_repo_local_mapping": "inferred",
            "sequence_vs_fgsg_conflict": "ambiguous",
            "multiple_candidate_targets": "ambiguous",
            "header_fgsg_not_mapped_to_canonical": "unresolved",
            "no_sequence_or_fgsg_bridge": "unresolved",
            "unclassified": "unresolved",
        }
    )
    bridge_df = bridge_df[
        [
            "source_protein_id",
            "source_header",
            "sequence_length",
            "sequence_sha256",
            "header_fgsg_id",
            "legacy_protein_accession",
            "legacy_seq_match_count",
            "fgsg_candidate_count",
            "unified_candidate_count",
            "resolved_canonical_gene_id",
            "bridge_status",
            "bridge_method",
            "mapping_confidence",
        ]
    ].sort_values("source_protein_id", kind="stable").reset_index(drop=True)

    yeast = read_tsv(paths["yeast_transfer_table"])
    yeast_bridge = yeast.merge(bridge_df, left_on="ph1_gene_id", right_on="source_protein_id", how="left")
    yeast_bridge["has_resolved_canonical"] = yeast_bridge["bridge_status"].eq("resolved")
    yeast_bridge["canonical_gene_id"] = yeast_bridge["resolved_canonical_gene_id"].fillna("")

    high_candidates = yeast_bridge[
        yeast_bridge["weak_positive_confidence"].eq("high") & yeast_bridge["bridge_status"].eq("resolved")
    ].copy()
    high_candidates = high_candidates.sort_values(["canonical_gene_id", "ph1_gene_id"], kind="stable")
    high_candidates["candidate_rank_within_gene"] = high_candidates.groupby("canonical_gene_id").cumcount() + 1
    high_protocolized = high_candidates.drop_duplicates("canonical_gene_id", keep="first").copy()
    high_protocolized["selection_rule"] = "keep_first_by_canonical_gene_after_sorted_xp"
    high_protocolized["construction_role"] = "high_confidence_yeast_transfer_positive_component"

    unresolved_high = yeast_bridge[
        yeast_bridge["weak_positive_confidence"].eq("high") & yeast_bridge["bridge_status"].ne("resolved")
    ].copy()
    unresolved_high["dropped_reason"] = unresolved_high["bridge_method"].fillna("missing_bridge_result")

    bridge_summary = pd.DataFrame(
        [
            {
                "ncbi_anchor_proteins_total": int(len(bridge_df)),
                "bridge_resolved_count": int((bridge_df["bridge_status"] == "resolved").sum()),
                "bridge_unresolved_count": int((bridge_df["bridge_status"] == "unresolved").sum()),
                "bridge_ambiguous_count": int((bridge_df["bridge_status"] == "ambiguous").sum()),
                "high_confidence_candidates_total": int((yeast_bridge["weak_positive_confidence"] == "high").sum()),
                "high_confidence_candidates_resolved_rows": int(len(high_candidates)),
                "high_confidence_candidates_protocolized_unique_genes": int(high_protocolized["canonical_gene_id"].nunique()),
                "high_confidence_candidates_unresolved_rows": int(len(unresolved_high)),
                "high_confidence_candidates_unified_map_resolved_rows": int(
                    (
                        yeast_bridge["weak_positive_confidence"].eq("high")
                        & yeast_bridge["bridge_method"].astype(str).str.contains("unified_map", regex=False)
                    ).sum()
                ),
                "yeast_none_candidates_resolved_unique_genes": int(
                    yeast_bridge.loc[
                        yeast_bridge["weak_positive_confidence"].eq("none") & yeast_bridge["bridge_status"].eq("resolved"),
                        "canonical_gene_id",
                    ].nunique()
                ),
                "source_anchor_ncbi_proteome": str(paths["anchor_ncbi_protein_fasta"]),
                "source_legacy_ph1_proteome": str(paths["ph1_legacy_protein_fasta"]),
                "source_legacy_mapping_tab": str(paths["ph1_legacy_mapping_tab"]),
            }
        ]
    )

    summary_md = "\n".join(
        [
            "# Fusarium Newlabel Protein-to-Canonical Bridge Summary",
            "",
            "## Purpose",
            "This bridge protocolizes the previously implicit XP-to-canonical mapping step required to convert the PH-1 yeast-transfer table into the final canonical Fusarium gene ID space.",
            "",
            "## Inputs",
            f"- repo-local proteome manifest: `{paths['proteome_manifest']}`",
            f"- PH-1 NCBI protein FASTA: `{paths['anchor_ncbi_protein_fasta']}`",
            f"- PH-1 legacy 229533 proteome: `{paths['ph1_legacy_protein_fasta']}`",
            f"- legacy accession-to-canonical table: `{paths['ph1_legacy_mapping_tab']}`",
            f"- unified XP/FGSG-to-PH1 map: `{paths['ph1_unified_id_map']}`",
            "- repo-local supplemental FGSG-to-canonical evidence: `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`, `data/processed/PPI/fgraminearum/string_id_mapping.tsv`, `data/processed/EXP/fgraminearum/exp_id_mapping.tsv`, `data/processed/LC/fgraminearum/subloc_id_mapping.tsv`",
            "",
            "## Bridge Logic",
            "The primary path is now the direct `XP_*` to `FGRAMPH1_*` mapping preserved in the shared unified ID map.",
            "That route is preferred when it provides a unique canonical PH-1 gene because it directly records the NCBI protein accession lineage used by the historical Fusarium orthology workflow.",
            "Sequence-based matching against the legacy `229533` proteome remains a second protocolized evidence path.",
            "When the NCBI header itself contains an `FGSG_*` locus tag, a repo-local FGSG-to-canonical bridge is also evaluated and recorded as supporting evidence.",
            "",
            "## Output Space",
            "The final canonical ID space is `fgraminearum::FGRAMPH1_*`.",
            "",
            "## Current Coverage",
            f"- total PH-1 NCBI proteins parsed: {int(len(bridge_df))}",
            f"- resolved protein-to-canonical mappings: {int((bridge_df['bridge_status'] == 'resolved').sum())}",
            f"- unresolved mappings: {int((bridge_df['bridge_status'] == 'unresolved').sum())}",
            f"- ambiguous mappings: {int((bridge_df['bridge_status'] == 'ambiguous').sum())}",
            f"- high-confidence yeast-transfer XP rows: {int((yeast_bridge['weak_positive_confidence'] == 'high').sum())}",
            f"- protocolized high-confidence transfer genes after canonical deduplication: {int(high_protocolized['canonical_gene_id'].nunique())}",
            f"- unresolved high-confidence transfer rows: {int(len(unresolved_high))}",
            f"- high-confidence transfer rows resolved via the unified XP/FGSG mapping path: {int((yeast_bridge['weak_positive_confidence'].eq('high') & yeast_bridge['bridge_method'].astype(str).str.contains('unified_map', regex=False)).sum())}",
            "",
            "## Caveats",
            "This bridge removes the direct dependence on `positive_set_P1.tsv` by making the XP-to-canonical mapping explicit.",
            "However, not every high-confidence XP protein currently resolves through the repo-local bridge evidence. Any unresolved XP rows are preserved in the audit outputs and are excluded from the protocolized rebuilt newlabel positive component until additional annotation evidence is protocolized.",
        ]
    ) + "\n"

    return {
        "bridge_df": bridge_df,
        "fgsg_to_canonical": fgsg_to_canonical,
        "high_protocolized": high_protocolized,
        "unresolved_high": unresolved_high,
        "summary_df": bridge_summary,
        "summary_md": summary_md,
    }


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    paths = config["paths"]
    bridge_dir = Path(paths["bridge_dir"])
    bridge_dir.mkdir(parents=True, exist_ok=True)

    bundle = build_bridge(config)
    bundle["bridge_df"].to_csv(bridge_dir / "protein_to_canonical_bridge.tsv", sep="\t", index=False)
    bundle["fgsg_to_canonical"].to_csv(bridge_dir / "source_to_canonical_mapping.tsv", sep="\t", index=False)
    bundle["high_protocolized"].to_csv(bridge_dir / "high_confidence_yeast_transfer_candidates.tsv", sep="\t", index=False)
    bundle["unresolved_high"].to_csv(bridge_dir / "unresolved_high_confidence_ids.tsv", sep="\t", index=False)
    bundle["summary_df"].to_csv(bridge_dir / "bridge_summary.tsv", sep="\t", index=False)
    (bridge_dir / "bridge_summary.md").write_text(bundle["summary_md"], encoding="utf-8")

    source_manifest = pd.DataFrame(
        [
            {
                "source_role": "proteome_manifest",
                "path": str(paths["proteome_manifest"]),
                "sha256": sha256_file(paths["proteome_manifest"]),
            },
            {
                "source_role": "anchor_ncbi_protein_fasta",
                "path": str(paths["anchor_ncbi_protein_fasta"]),
                "sha256": sha256_file(paths["anchor_ncbi_protein_fasta"]),
            },
            {
                "source_role": "legacy_ph1_protein_fasta",
                "path": str(paths["ph1_legacy_protein_fasta"]),
                "sha256": sha256_file(paths["ph1_legacy_protein_fasta"]),
            },
            {
                "source_role": "legacy_mapping_tab",
                "path": str(paths["ph1_legacy_mapping_tab"]),
                "sha256": sha256_file(paths["ph1_legacy_mapping_tab"]),
            },
            {
                "source_role": "ph1_unified_id_map",
                "path": str(paths["ph1_unified_id_map"]),
                "sha256": sha256_file(paths["ph1_unified_id_map"]),
            },
        ]
    )
    source_manifest.to_csv(bridge_dir / "bridge_source_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
