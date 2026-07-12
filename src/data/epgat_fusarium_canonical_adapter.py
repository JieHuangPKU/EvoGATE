"""
Fusarium canonical adapter for Phase 1.5.

This adapter keeps the original-compatible legacy feature contract but
projects legacy Fusarium IDs into current ProGATE_v2 canonical_gene_id.
"""

import pandas as pd


def normalize_token(value):
    text = str(value).strip() if value is not None else ""
    if not text or text.lower() == "nan":
        return ""
    return text


def build_fusarium_canonical_lookup(registry_path):
    df = pd.read_csv(registry_path, sep="\t", dtype=str).fillna("")
    fg = df[df["species"] == "fgraminearum"].copy()
    rows = []
    alias_to_records = {}
    for _, row in fg.iterrows():
        canonical = row["canonical_gene_id"]
        aliases = [
            canonical,
            row.get("raw_gene_id", ""),
            row.get("raw_protein_id", ""),
            row.get("raw_transcript_id", ""),
        ]
        for alias in aliases:
            alias = normalize_token(alias)
            if not alias:
                continue
            alias_to_records.setdefault(alias, set()).add(canonical)
    resolved = {}
    for alias, canonicals in alias_to_records.items():
        if len(canonicals) == 1:
            resolved[alias] = list(canonicals)[0]
    return resolved


def map_legacy_ids_to_canonical(values, alias_lookup):
    records = []
    for value in values:
        token = normalize_token(value)
        if not token:
            records.append(
                {
                    "legacy_id": token,
                    "mapped_canonical_gene_id": "",
                    "mapping_status": "unresolved",
                    "mapping_rule": "empty_id",
                    "notes": "",
                }
            )
        elif token in alias_lookup:
            records.append(
                {
                    "legacy_id": token,
                    "mapped_canonical_gene_id": alias_lookup[token],
                    "mapping_status": "exact",
                    "mapping_rule": "exact_alias_lookup",
                    "notes": "",
                }
            )
        else:
            records.append(
                {
                    "legacy_id": token,
                    "mapped_canonical_gene_id": "",
                    "mapping_status": "unresolved",
                    "mapping_rule": "not_found_in_registry_alias_lookup",
                    "notes": "",
                }
            )
    return pd.DataFrame(records)
