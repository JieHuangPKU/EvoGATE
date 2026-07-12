"""
EPGAT original-compatible label adapter.

Reads either the legacy `ogee.csv` contract or explicit override
positive/negative label sets and returns the binary label table keyed by
legacy gene ID.
"""

import os

import pandas as pd

from src.data.epgat_legacy_id_adapter import normalize_legacy_gene_id


def _normalize_override_gene_id(value):
    gene_id = normalize_legacy_gene_id(value)
    if "::" in gene_id:
        gene_id = gene_id.split("::", 1)[1].strip()
    return normalize_legacy_gene_id(gene_id)


def load_legacy_labels(ogee_path):
    df = pd.read_csv(ogee_path, dtype={"Gene": str, "Label": int}).fillna("")
    if "Gene" not in df.columns or "Label" not in df.columns:
        raise ValueError("Legacy ogee.csv must contain Gene and Label columns")
    df["legacy_gene_id"] = df["Gene"].map(normalize_legacy_gene_id)
    df = df[df["legacy_gene_id"].astype(str).ne("")].copy()
    df = df.drop_duplicates(subset=["legacy_gene_id"], keep="first")
    df["label"] = df["Label"].astype(int)
    return df[["legacy_gene_id", "label"]].copy()

def load_override_labels(positive_path, negative_path):
    if not positive_path or not negative_path:
        raise ValueError("Both positive_path and negative_path are required for override labels")
    if not os.path.exists(positive_path):
        raise FileNotFoundError("Override positive label file not found: {}".format(positive_path))
    if not os.path.exists(negative_path):
        raise FileNotFoundError("Override negative label file not found: {}".format(negative_path))

    pos_df = pd.read_csv(positive_path, sep="\t", dtype=str).fillna("")
    neg_df = pd.read_csv(negative_path, sep="\t", dtype=str).fillna("")
    if "canonical_gene_id" not in pos_df.columns:
        raise ValueError("Override positive label file must contain canonical_gene_id column: {}".format(positive_path))
    if "canonical_gene_id" not in neg_df.columns:
        raise ValueError("Override negative label file must contain canonical_gene_id column: {}".format(negative_path))

    pos = pd.DataFrame({"legacy_gene_id": pos_df["canonical_gene_id"].map(_normalize_override_gene_id)})
    pos = pos[pos["legacy_gene_id"].astype(str).ne("")].drop_duplicates(subset=["legacy_gene_id"], keep="first").copy()
    pos["label"] = 1

    neg = pd.DataFrame({"legacy_gene_id": neg_df["canonical_gene_id"].map(_normalize_override_gene_id)})
    neg = neg[neg["legacy_gene_id"].astype(str).ne("")].drop_duplicates(subset=["legacy_gene_id"], keep="first").copy()
    neg["label"] = 0

    overlap = sorted(set(pos["legacy_gene_id"]).intersection(set(neg["legacy_gene_id"])))
    if overlap:
        raise ValueError("Override label sets overlap for {} genes; first examples: {}".format(len(overlap), ", ".join(overlap[:10])))

    labels = pd.concat([pos, neg], ignore_index=True)
    labels = labels.drop_duplicates(subset=["legacy_gene_id"], keep="first").reset_index(drop=True)
    return labels[["legacy_gene_id", "label"]].copy()

