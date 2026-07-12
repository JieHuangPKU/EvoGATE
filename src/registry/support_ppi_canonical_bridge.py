from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


BRIDGE_COLUMNS = [
    "species",
    "raw_node_id",
    "raw_node_id_type",
    "epgat_processed_node_id",
    "epgat_processed_node_id_type",
    "current_embedding_gene_id",
    "current_embedding_gene_id_type",
    "proposed_canonical_gene_id",
    "bridge_rule",
    "confidence",
    "bridge_status",
    "needs_manual_review",
    "evidence_source",
    "notes",
]


def _master_labels(species: str) -> pd.DataFrame:
    df = pd.read_csv("data_registry/master_label_table.preliminary.tsv", sep="\t", dtype=str).fillna("")
    return df[df["species"] == species].copy()


def _human_bridge() -> pd.DataFrame:
    species = "human"
    labels = _master_labels(species)
    protein_to_canonical = dict(
        zip(labels["raw_protein_id"].astype(str), labels["canonical_gene_id"].astype(str))
    )
    ppi = pd.read_csv(
        "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/human/PPI/STRING/string.csv",
        usecols=["A", "B"],
        dtype=str,
    ).fillna("")
    raw_nodes = sorted(set(ppi["A"].astype(str)) | set(ppi["B"].astype(str)))
    rows = []
    for raw in raw_nodes:
        canonical = protein_to_canonical.get(raw, "")
        if canonical:
            rows.append(
                {
                    "species": species,
                    "raw_node_id": raw,
                    "raw_node_id_type": "uniprot_like_protein_accession",
                    "epgat_processed_node_id": raw,
                    "epgat_processed_node_id_type": "uniprot_like_protein_accession",
                    "current_embedding_gene_id": canonical.split("::", 1)[1],
                    "current_embedding_gene_id_type": "ENSG_identifier",
                    "proposed_canonical_gene_id": canonical,
                    "bridge_rule": "raw_protein_id_exact_match",
                    "confidence": "high",
                    "bridge_status": "exact",
                    "needs_manual_review": "false",
                    "evidence_source": "master_label_table.raw_protein_id",
                    "notes": "",
                }
            )
        else:
            rows.append(
                {
                    "species": species,
                    "raw_node_id": raw,
                    "raw_node_id_type": "uniprot_like_protein_accession",
                    "epgat_processed_node_id": raw,
                    "epgat_processed_node_id_type": "uniprot_like_protein_accession",
                    "current_embedding_gene_id": "",
                    "current_embedding_gene_id_type": "ENSG_identifier",
                    "proposed_canonical_gene_id": "",
                    "bridge_rule": "raw_protein_id_exact_match",
                    "confidence": "low",
                    "bridge_status": "missing",
                    "needs_manual_review": "true",
                    "evidence_source": "master_label_table.raw_protein_id",
                    "notes": "",
                }
            )
    return pd.DataFrame(rows, columns=BRIDGE_COLUMNS)


def _scerevisiae_bridge() -> pd.DataFrame:
    species = "scerevisiae"
    labels = _master_labels(species)
    protein_to_canonical = dict(
        zip(labels["raw_protein_id"].astype(str), labels["canonical_gene_id"].astype(str))
    )
    ppi = pd.read_csv(
        "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/yeast/PPI/STRING/string.csv",
        usecols=["A", "B", "protein1", "protein2"],
        dtype=str,
    ).fillna("")
    raw_nodes = sorted(set(ppi["A"].astype(str)) | set(ppi["B"].astype(str)))
    rows = []
    for raw in raw_nodes:
        canonical = protein_to_canonical.get(raw, "")
        if canonical:
            rows.append(
                {
                    "species": species,
                    "raw_node_id": raw,
                    "raw_node_id_type": "uniprot_like_protein_accession",
                    "epgat_processed_node_id": raw,
                    "epgat_processed_node_id_type": "uniprot_like_protein_accession",
                    "current_embedding_gene_id": canonical.split("::", 1)[1],
                    "current_embedding_gene_id_type": "S000_identifier",
                    "proposed_canonical_gene_id": canonical,
                    "bridge_rule": "string_A_B_to_raw_protein_id_exact_match",
                    "confidence": "high",
                    "bridge_status": "exact",
                    "needs_manual_review": "false",
                    "evidence_source": "master_label_table.raw_protein_id",
                    "notes": "STRING protein1/protein2 remain unresolved, but STRING A/B columns bridge directly to raw_protein_id.",
                }
            )
        else:
            rows.append(
                {
                    "species": species,
                    "raw_node_id": raw,
                    "raw_node_id_type": "uniprot_like_protein_accession",
                    "epgat_processed_node_id": raw,
                    "epgat_processed_node_id_type": "uniprot_like_protein_accession",
                    "current_embedding_gene_id": "",
                    "current_embedding_gene_id_type": "S000_identifier",
                    "proposed_canonical_gene_id": "",
                    "bridge_rule": "string_A_B_to_raw_protein_id_exact_match",
                    "confidence": "low",
                    "bridge_status": "missing",
                    "needs_manual_review": "true",
                    "evidence_source": "master_label_table.raw_protein_id",
                    "notes": "",
                }
            )
    return pd.DataFrame(rows, columns=BRIDGE_COLUMNS)


def _celegans_bridge() -> pd.DataFrame:
    species = "celegans"
    labels = _master_labels(species)
    wbgene_to_canonical = dict(
        zip(labels["raw_gene_id"].astype(str), labels["canonical_gene_id"].astype(str))
    )
    wb_cols = ["taxid", "species_name", "WBGene_id", "gene_symbol", "sequence_id"]
    wb = pd.read_csv(
        "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/celegans/PPI/STRING/wormbase.WS240.gene_ids.txt",
        sep="\t",
        comment="#",
        names=wb_cols,
        dtype=str,
    ).fillna("")
    alias_rows = []
    alias_to_wbgene: dict[tuple[str, str], str] = {}
    ambiguous_aliases: set[str] = set()

    def register_alias(alias: str, alias_type: str, wbgene_id: str) -> None:
        if not alias:
            return
        key = (alias_type, alias)
        existing = alias_to_wbgene.get(key)
        if existing is None:
            alias_to_wbgene[key] = wbgene_id
        elif existing != wbgene_id:
            ambiguous_aliases.add(f"{alias_type}::{alias}")
        alias_rows.append(
            {
                "alias": alias,
                "alias_type": alias_type,
                "WBGene_id": wbgene_id,
            }
        )

    for _, row in wb.iterrows():
        register_alias(row["gene_symbol"], "symbol", row["WBGene_id"])
        register_alias(row["sequence_id"], "sequence_id", row["WBGene_id"])

    alias_df = pd.DataFrame(alias_rows).drop_duplicates()
    Path("data_registry").mkdir(parents=True, exist_ok=True)
    # Add raw STRING-derived alias audit rows from the raw network file.
    raw_string = pd.read_csv(
        "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/celegans/PPI/STRING/6239.protein.links.detailed.v12.0.txt",
        sep=" ",
        usecols=["protein1", "protein2"],
        dtype=str,
    ).fillna("")
    raw_nodes = sorted(set(raw_string["protein1"].astype(str)) | set(raw_string["protein2"].astype(str)))

    def normalize_string_node(raw_id: str) -> tuple[str, str]:
        no_taxon = raw_id.split(".", 1)[1] if raw_id.startswith("6239.") else raw_id
        trimmed = re.sub(r"\.\d+$", "", no_taxon) if re.search(r"\.\d+$", no_taxon) else no_taxon
        return no_taxon, trimmed

    extra_alias_rows = []
    ppi_to_wbgene_rows = []
    rows = []
    for raw in raw_nodes:
        no_taxon, trimmed = normalize_string_node(raw)
        wbgene_id = ""
        canonical = ""
        rule = ""
        status = "missing"
        confidence = "low"
        review = "true"
        notes = ""

        checks = [
            ("string_node_raw", raw, "symbol->WBGene->canonical"),
            ("string_node_raw", raw, "sequence_id->WBGene->canonical"),
            ("string_node_no_taxon", no_taxon, "string_no_taxon->sequence_id->WBGene->canonical"),
            ("string_node_trimmed_isoform", trimmed, "string_trimmed_isoform->sequence_id->WBGene->canonical"),
            ("symbol", no_taxon, "symbol->WBGene->canonical"),
            ("sequence_id", no_taxon, "sequence_id->WBGene->canonical"),
            ("sequence_id", trimmed, "string_trimmed_isoform->sequence_id->WBGene->canonical"),
        ]

        # audit alias rows
        extra_alias_rows.extend(
            [
                {"alias": raw, "alias_type": "string_node_raw", "WBGene_id": ""},
                {"alias": no_taxon, "alias_type": "string_node_no_taxon", "WBGene_id": ""},
                {"alias": trimmed, "alias_type": "string_node_trimmed_isoform", "WBGene_id": ""},
            ]
        )

        candidate_wb = []
        for alias_type, alias_value, mapping_path in checks:
            key = (alias_type, alias_value)
            amb_key = f"{alias_type}::{alias_value}"
            if amb_key in ambiguous_aliases:
                status = "ambiguous"
                rule = mapping_path
                notes = "Alias maps to multiple WBGene IDs."
                break
            match = alias_to_wbgene.get(key, "")
            if match:
                candidate_wb = [match]
                rule = mapping_path
                break
            # allow symbol/sequence lookup when raw/no_taxon/trimmed equal those values
            if alias_type.startswith("string_node"):
                for fallback_type in ["symbol", "sequence_id"]:
                    amb_key = f"{fallback_type}::{alias_value}"
                    if amb_key in ambiguous_aliases:
                        status = "ambiguous"
                        rule = f"{mapping_path} via {fallback_type}"
                        notes = "Alias maps to multiple WBGene IDs."
                        break
                    match = alias_to_wbgene.get((fallback_type, alias_value), "")
                    if match:
                        candidate_wb = [match]
                        rule = mapping_path if fallback_type == "sequence_id" else f"{mapping_path} via symbol"
                        break
                if candidate_wb or status == "ambiguous":
                    break

        if status == "ambiguous":
            wbgene_id = ""
            canonical = ""
            confidence = "low"
            review = "true"
        elif candidate_wb:
            wbgene_id = candidate_wb[0]
            canonical = wbgene_to_canonical.get(wbgene_id, "")
            if canonical:
                status = "bridged" if rule.startswith("string_") else "alias"
                confidence = "high"
                review = "false"
            else:
                status = "conflict"
                confidence = "low"
                review = "true"
                notes = "WBGene resolved, but no canonical_gene_id was found in current registry."
        else:
            status = "missing"
            confidence = "low"
            review = "true"
            notes = "No exact alias match found in WormBase alias inventory."
            rule = "no_match"

        ppi_to_wbgene_rows.append(
            {
                "raw_ppi_node_id": raw,
                "wbgene_id": wbgene_id,
                "mapping_status": status,
                "mapping_path": rule,
                "needs_manual_review": review,
            }
        )
        rows.append(
            {
                "species": species,
                "raw_node_id": raw,
                "raw_node_id_type": "string_internal_protein_id",
                "epgat_processed_node_id": no_taxon,
                "epgat_processed_node_id_type": "string_node_no_taxon",
                "current_embedding_gene_id": wbgene_id,
                "current_embedding_gene_id_type": "WBGene_identifier",
                "proposed_canonical_gene_id": canonical,
                "bridge_rule": rule,
                "confidence": confidence,
                "bridge_status": status,
                "needs_manual_review": review,
                "evidence_source": "wormbase.WS240.gene_ids.txt + raw STRING node normalization + master_label_table.raw_gene_id",
                "notes": notes,
            }
        )

    alias_df = pd.concat([alias_df, pd.DataFrame(extra_alias_rows)], ignore_index=True).drop_duplicates()
    alias_df.to_csv("data_registry/wb_alias_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(ppi_to_wbgene_rows).to_csv(
        "outputs/support_graphs/celegans_ppi_to_wbgene.tsv", sep="\t", index=False
    )
    return pd.DataFrame(rows, columns=BRIDGE_COLUMNS)


def build_support_ppi_bridge_tables() -> dict[str, pd.DataFrame]:
    return {
        "scerevisiae": _scerevisiae_bridge(),
        "human": _human_bridge(),
        "celegans": _celegans_bridge(),
    }


def write_support_ppi_bridge_tables() -> None:
    registry_dir = Path("data_registry")
    registry_dir.mkdir(parents=True, exist_ok=True)
    for species, df in build_support_ppi_bridge_tables().items():
        df.to_csv(registry_dir / f"{species}_ppi_canonical_bridge.tsv", sep="\t", index=False)


if __name__ == "__main__":
    write_support_ppi_bridge_tables()
