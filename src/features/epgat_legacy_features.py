"""
EPGAT original-compatible feature assembly.

Important:
- This preserves the original feature block order observed in
  `EPGAT/utils/utils.py`:
  1. orthologs
  2. expression
  3. sublocalization
  4. degree
- This is Phase 1 only and intentionally excludes ESM/ESMC/PLM blocks.
"""

import os

import numpy as np
import pandas as pd


def load_feature_table(path, prefix):
    df = pd.read_csv(path, dtype=str).fillna("")
    if "Gene" not in df.columns:
        raise ValueError("Legacy feature table missing Gene column: {}".format(path))
    gene_ids = df["Gene"].astype(str)
    values = df.drop(columns=["Gene"]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.columns = ["{}_{}".format(prefix, idx) for idx in range(values.shape[1])]
    out = pd.DataFrame({"legacy_gene_id": gene_ids})
    for col in values.columns:
        out[col] = values[col].to_numpy()
    out = out.drop_duplicates(subset=["legacy_gene_id"], keep="first")
    return out


def append_degree_block(feature_df, edges_df):
    degree_map = {}
    for column in ["A", "B"]:
        counts = edges_df[column].astype(str).value_counts()
        for gene_id, count in counts.items():
            degree_map[gene_id] = degree_map.get(gene_id, 0) + int(count)
    feature_df["degree_0"] = feature_df["legacy_gene_id"].astype(str).map(degree_map).fillna(0).astype(float)
    return feature_df


def zscore_matrix(x):
    mean = x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    std[std < 1e-8] = 1.0
    return (x - mean) / std


def build_feature_schema(expression_df, ortholog_df, subloc_df, include_degree):
    rows = []
    start = 0
    blocks = [
        ("orthologs", ortholog_df),
        ("expression", expression_df),
        ("sublocalization", subloc_df),
    ]
    for block_name, df in blocks:
        cols = [c for c in df.columns if c != "legacy_gene_id"]
        if not cols:
            continue
        end = start + len(cols) - 1
        rows.append(
            {
                "feature_block": block_name,
                "start_col": start,
                "end_col": end,
                "dimension": len(cols),
                "data_source": block_name,
                "missing_strategy": "left_join_fillna_zero",
            }
        )
        start = end + 1
    if include_degree:
        rows.append(
            {
                "feature_block": "degree",
                "start_col": start,
                "end_col": start,
                "dimension": 1,
                "data_source": "ppi_degree_from_string_network",
                "missing_strategy": "computed_zero_if_absent",
            }
        )
    return pd.DataFrame(rows)
