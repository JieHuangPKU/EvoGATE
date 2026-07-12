"""
EPGAT original-compatible ID adapter.

This adapter makes legacy gene IDs explicit inside ProGATE_v2 Phase 1.
It is intentionally minimal and only supports the original gene-level
EPGAT contract, not the later extended canonical/multimodal layers.
"""

import pandas as pd


def normalize_legacy_gene_id(value):
    text = str(value).strip() if value is not None else ""
    if not text or text.lower() == "nan":
        return ""
    return text


def build_id_alignment_manifest(node_ids, label_ids, feature_tables):
    rows = []
    label_set = set([normalize_legacy_gene_id(v) for v in label_ids])
    feature_sets = {
        name: set([normalize_legacy_gene_id(v) for v in values])
        for name, values in feature_tables.items()
    }
    for node_id in sorted(set([normalize_legacy_gene_id(v) for v in node_ids])):
        rows.append(
            {
                "legacy_gene_id": node_id,
                "internal_gene_id": node_id,
                "has_label": node_id in label_set,
                "has_expression": node_id in feature_sets.get("expression", set()),
                "has_orthologs": node_id in feature_sets.get("orthologs", set()),
                "has_sublocalization": node_id in feature_sets.get("sublocalization", set()),
                "alignment_status": "aligned" if node_id else "invalid",
                "unaligned_reason": "" if node_id else "empty_gene_id",
            }
        )
    return pd.DataFrame(rows)
