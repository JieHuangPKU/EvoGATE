"""
Build a minimal original-compatible EPGAT legacy dataset inside ProGATE_v2.

This builder preserves an explicit legacy-compatible feature contract
driven by config flags:
- orthologs (optional)
- expression (optional)
- sublocalization (optional)
- degree (optional)

It outputs explicit manifests and schema files so the old implicit
assumptions are auditable in ProGATE_v2.
"""

import argparse
import os

import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import train_test_split

from src.data.epgat_legacy_id_adapter import build_id_alignment_manifest
from src.data.epgat_legacy_label_adapter import load_legacy_labels, load_override_labels
from src.features.epgat_legacy_features import (
    append_degree_block,
    build_feature_schema,
    load_feature_table,
    zscore_matrix,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build original-compatible EPGAT legacy dataset")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_legacy_ppi(ppi_path, string_thr, use_weights):
    edges = pd.read_csv(ppi_path)
    edge_weights = None
    if "combined_score" in edges.columns:
        edges = edges[edges["combined_score"].astype(float) > float(string_thr)].reset_index(drop=True)
        edge_weights = edges["combined_score"].astype(float).to_numpy() / 1000.0
    edge_df = edges[["A", "B"]].dropna().drop_duplicates().reset_index(drop=True)
    if edge_weights is not None and len(edge_weights) != len(edge_df):
        edge_weights = edges.loc[edge_df.index, "combined_score"].astype(float).to_numpy() / 1000.0
    if not use_weights:
        edge_weights = None
    return edge_df, edge_weights


def build_input_completeness(data_root, organism):
    files = {
        "labels": os.path.join(data_root, "EssentialGenes", "ogee.csv"),
        "ppi": os.path.join(data_root, "PPI", "STRING", "string.csv"),
        "expression": os.path.join(data_root, "Expression", "profile.csv"),
        "orthologs": os.path.join(data_root, "Orthologs", "orthologs.csv"),
        "sublocalization": os.path.join(data_root, "SubLocalizations", "subloc.csv"),
    }
    rows = []
    for name, path in files.items():
        exists = os.path.exists(path)
        row_count = ""
        column_count = ""
        usable = "false"
        notes = ""
        if exists:
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as handle:
                    row_count = max(sum(1 for _ in handle) - 1, 0)
                header = pd.read_csv(path, nrows=0, dtype=str)
                column_count = len(header.columns)
                usable = "true"
            except Exception as exc:
                notes = "read_error: {}".format(exc)
        else:
            notes = "missing_file"
        rows.append(
            {
                "species": organism,
                "input_block": name,
                "source_path": path,
                "exists": exists,
                "row_count": row_count,
                "column_count": column_count,
                "directly_usable": usable,
                "needs_species_specific_adapter": "false",
                "notes": notes,
            }
        )
    return pd.DataFrame(rows)


def build_dataset(config):
    organism = config["legacy"]["organism"]
    legacy_species = config["legacy"]["legacy_species_dir"]
    data_root = os.path.join(config["paths"]["legacy_epgat_root"], "data", "essential_genes", legacy_species)
    out_dir = os.path.join(config["paths"]["output_root"], config["run"]["name"])
    os.makedirs(out_dir, exist_ok=True)

    completeness = build_input_completeness(data_root, organism)
    completeness.to_csv(os.path.join(out_dir, "input_completeness.tsv"), sep="\t", index=False)
    completeness_lines = [
        "# Input Completeness",
        "",
    ]
    for _, row in completeness.iterrows():
        completeness_lines.append(
            "- {}: exists={}, row_count={}, column_count={}, directly_usable={}, notes={}".format(
                row["input_block"],
                row["exists"],
                row["row_count"],
                row["column_count"],
                row["directly_usable"],
                row["notes"],
            )
        )
    with open(os.path.join(out_dir, "input_completeness.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(completeness_lines))

    positive_override = str(config["legacy"].get("positive_set_path", "")).strip()
    negative_override = str(config["legacy"].get("negative_set_path", "")).strip()
    if positive_override or negative_override:
        if not (positive_override and negative_override):
            raise ValueError("Both positive_set_path and negative_set_path must be provided together")
        labels = load_override_labels(positive_override, negative_override)
    else:
        labels = load_legacy_labels(os.path.join(data_root, "EssentialGenes", "ogee.csv"))
    edge_df, edge_weights = load_legacy_ppi(
        os.path.join(data_root, "PPI", "STRING", "string.csv"),
        config["legacy"]["string_threshold"],
        bool(config["legacy"].get("use_weights", False)),
    )
    ppi_genes = sorted(set(edge_df["A"].astype(str)) | set(edge_df["B"].astype(str)))
    labels = labels[labels["legacy_gene_id"].isin(ppi_genes)].copy()
    genes = sorted(set(ppi_genes) | set(labels["legacy_gene_id"].astype(str)))

    orth_enabled = bool(config["legacy"].get("orthologs", False))
    expr_enabled = bool(config["legacy"].get("expression", False))
    sub_enabled = bool(config["legacy"].get("sublocalization", False))
    degree_enabled = bool(config["legacy"].get("include_degree", False))

    orth_df = load_feature_table(os.path.join(data_root, "Orthologs", "orthologs.csv"), "ortholog") if orth_enabled else pd.DataFrame({"legacy_gene_id": []})
    expr_df = load_feature_table(os.path.join(data_root, "Expression", "profile.csv"), "expression") if expr_enabled else pd.DataFrame({"legacy_gene_id": []})
    sub_df = load_feature_table(os.path.join(data_root, "SubLocalizations", "subloc.csv"), "subloc") if sub_enabled else pd.DataFrame({"legacy_gene_id": []})

    feature_df = pd.DataFrame({"legacy_gene_id": genes})
    if orth_enabled:
        feature_df = feature_df.merge(orth_df, on="legacy_gene_id", how="left")
    if expr_enabled:
        feature_df = feature_df.merge(expr_df, on="legacy_gene_id", how="left")
    if sub_enabled:
        feature_df = feature_df.merge(sub_df, on="legacy_gene_id", how="left")
    feature_df = feature_df.fillna(0.0)
    if degree_enabled:
        feature_df = append_degree_block(feature_df, edge_df)

    feature_cols = [c for c in feature_df.columns if c != "legacy_gene_id"]
    x = feature_df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
    x = zscore_matrix(x).astype(np.float32)
    mapping = dict(zip(genes, range(len(genes))))
    edge_index = np.vectorize(mapping.__getitem__)(edge_df[["A", "B"]].to_numpy(dtype=object)).astype(np.int64)

    label_map = dict(zip(labels["legacy_gene_id"], labels["label"]))
    label_manifest = pd.DataFrame({"legacy_gene_id": genes})
    label_manifest["label"] = label_manifest["legacy_gene_id"].map(label_map)
    label_manifest["is_labeled"] = label_manifest["label"].notna()

    labeled = label_manifest[label_manifest["is_labeled"]].copy()
    train_df, test_df = train_test_split(
        labeled,
        test_size=float(config["legacy"]["test_fraction"]),
        random_state=int(config["legacy"]["seed"]),
        stratify=labeled["label"].astype(int),
    )
    train_df, val_df = train_test_split(
        train_df,
        test_size=float(config["legacy"]["val_fraction"]) / (1.0 - float(config["legacy"]["test_fraction"])),
        random_state=int(config["legacy"]["seed"]),
        stratify=train_df["label"].astype(int),
    )
    split_map = {}
    for split_name, frame in [("train", train_df), ("val", val_df), ("test", test_df)]:
        for gene in frame["legacy_gene_id"].astype(str):
            split_map[gene] = split_name
    label_manifest["split"] = label_manifest["legacy_gene_id"].map(split_map).fillna("")

    id_manifest = build_id_alignment_manifest(
        node_ids=genes,
        label_ids=labels["legacy_gene_id"].astype(str).tolist(),
        feature_tables={
            "orthologs": orth_df["legacy_gene_id"].astype(str).tolist() if orth_enabled else [],
            "expression": expr_df["legacy_gene_id"].astype(str).tolist() if expr_enabled else [],
            "sublocalization": sub_df["legacy_gene_id"].astype(str).tolist() if sub_enabled else [],
        },
    )
    node_manifest = id_manifest.copy()
    node_manifest["in_ppi"] = node_manifest["legacy_gene_id"].isin(set(ppi_genes))
    node_manifest["legacy_species"] = organism

    schema = build_feature_schema(expr_df, orth_df, sub_df, include_degree=degree_enabled)

    np.save(os.path.join(out_dir, "feature_matrix.npy"), x)
    np.save(os.path.join(out_dir, "edge_index.npy"), edge_index)
    if edge_weights is not None:
        np.save(os.path.join(out_dir, "edge_weights.npy"), np.asarray(edge_weights, dtype=np.float32))

    node_manifest.to_csv(os.path.join(out_dir, "node_manifest.tsv"), sep="\t", index=False)
    label_manifest.to_csv(os.path.join(out_dir, "label_manifest.tsv"), sep="\t", index=False)
    edge_df.to_csv(os.path.join(out_dir, "edge_table.tsv"), sep="\t", index=False)
    schema.to_csv(os.path.join(out_dir, "feature_schema.tsv"), sep="\t", index=False)

    audit_rows = [
        {
            "species": organism,
            "row_count_after_join": len(feature_df),
            "unique_legacy_gene_id_count": len(set(genes)),
            "missing_label_count": int(label_manifest["label"].isna().sum()),
            "missing_expression_count": int((~node_manifest["has_expression"]).sum()),
            "missing_ortholog_count": int((~node_manifest["has_orthologs"]).sum()),
            "missing_subloc_count": int((~node_manifest["has_sublocalization"]).sum()),
            "duplicated_legacy_gene_id_count": int(node_manifest["legacy_gene_id"].duplicated().sum()),
        }
    ]
    pd.DataFrame(audit_rows).to_csv(os.path.join(out_dir, "dataset_alignment_audit.tsv"), sep="\t", index=False)

    enabled_blocks = []
    if orth_enabled:
        enabled_blocks.append("orthologs")
    if expr_enabled:
        enabled_blocks.append("expression")
    if sub_enabled:
        enabled_blocks.append("sublocalization")
    if degree_enabled:
        enabled_blocks.append("degree")

    summary_lines = [
        "# Legacy Dataset Summary",
        "",
        "- organism: {}".format(organism),
        "- total nodes: {}".format(len(genes)),
        "- total labeled genes: {}".format(len(labeled)),
        "- total edges: {}".format(len(edge_df)),
        "- feature dimension: {}".format(x.shape[1]),
        "- split counts: train={}, val={}, test={}".format(len(train_df), len(val_df), len(test_df)),
        "- enabled feature blocks: {}".format(" -> ".join(enabled_blocks) if enabled_blocks else "none"),
    ]
    with open(os.path.join(out_dir, "dataset_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines))
    return out_dir


def main():
    args = parse_args()
    config = load_yaml(args.config)
    out_dir = build_dataset(config)
    print("Legacy dataset built:", out_dir)


if __name__ == "__main__":
    main()
