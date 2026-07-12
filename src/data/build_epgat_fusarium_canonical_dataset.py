"""
Build canonicalized Fusarium replay dataset for Phase 1.5.
"""

import argparse
import os

import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import train_test_split

from src.data.epgat_fusarium_canonical_adapter import build_fusarium_canonical_lookup, map_legacy_ids_to_canonical
from src.features.epgat_legacy_features import append_degree_block, build_feature_schema, load_feature_table, zscore_matrix


def parse_args():
    parser = argparse.ArgumentParser(description="Build canonicalized Fusarium legacy replay dataset")
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
    if not use_weights:
        edge_weights = None
    return edge_df, edge_weights


def build_dataset(config):
    data_root = os.path.join(config["paths"]["legacy_epgat_root"], "data", "essential_genes", "fgraminearum")
    out_dir = os.path.join(config["paths"]["output_root"], config["run"]["name"])
    os.makedirs(out_dir, exist_ok=True)

    alias_lookup = build_fusarium_canonical_lookup(config["paths"]["registry_table"])

    labels = pd.read_csv(os.path.join(data_root, "EssentialGenes", "ogee.csv"), dtype=str).fillna("")
    labels["legacy_gene_id"] = labels["Gene"].astype(str)
    label_map = map_legacy_ids_to_canonical(labels["legacy_gene_id"].tolist(), alias_lookup)
    label_map = label_map.rename(columns={"legacy_id": "legacy_gene_id", "mapped_canonical_gene_id": "canonical_gene_id"})
    label_frame = labels.merge(label_map, on="legacy_gene_id", how="left")
    label_frame = label_frame[label_frame["mapping_status"] == "exact"].copy()
    label_frame["label"] = label_frame["Label"].astype(int)

    edge_df, edge_weights = load_legacy_ppi(
        os.path.join(data_root, "PPI", "STRING", "string.csv"),
        config["legacy"]["string_threshold"],
        bool(config["legacy"].get("use_weights", False)),
    )
    ppi_ids = sorted(set(edge_df["A"].astype(str)) | set(edge_df["B"].astype(str)))
    ppi_map = map_legacy_ids_to_canonical(ppi_ids, alias_lookup).rename(columns={"legacy_id": "legacy_gene_id", "mapped_canonical_gene_id": "canonical_gene_id"})
    ppi_lookup = dict(zip(ppi_map[ppi_map["mapping_status"] == "exact"]["legacy_gene_id"], ppi_map[ppi_map["mapping_status"] == "exact"]["canonical_gene_id"]))
    edge_df["source_canonical_gene_id"] = edge_df["A"].astype(str).map(ppi_lookup).fillna("")
    edge_df["target_canonical_gene_id"] = edge_df["B"].astype(str).map(ppi_lookup).fillna("")
    edge_df = edge_df[
        edge_df["source_canonical_gene_id"].astype(str).ne("") &
        edge_df["target_canonical_gene_id"].astype(str).ne("")
    ].copy()
    edge_df = edge_df.drop_duplicates(subset=["source_canonical_gene_id", "target_canonical_gene_id"], keep="first")

    def map_feature(path, prefix):
        df = load_feature_table(path, prefix)
        mapped = map_legacy_ids_to_canonical(df["legacy_gene_id"].tolist(), alias_lookup)
        mapped = mapped.rename(columns={"legacy_id": "legacy_gene_id", "mapped_canonical_gene_id": "canonical_gene_id"})
        out = df.merge(mapped, on="legacy_gene_id", how="left")
        out = out[out["mapping_status"] == "exact"].copy()
        cols = [c for c in out.columns if c.startswith(prefix + "_")]
        return out[["canonical_gene_id"] + cols].drop_duplicates(subset=["canonical_gene_id"], keep="first")

    orth_df = map_feature(os.path.join(data_root, "Orthologs", "orthologs.csv"), "ortholog")
    expr_df = map_feature(os.path.join(data_root, "Expression", "profile.csv"), "expression")
    sub_df = map_feature(os.path.join(data_root, "SubLocalizations", "subloc.csv"), "subloc")

    genes = sorted(set(edge_df["source_canonical_gene_id"]) | set(edge_df["target_canonical_gene_id"]) | set(label_frame["canonical_gene_id"]))
    feature_df = pd.DataFrame({"legacy_gene_id": genes, "canonical_gene_id": genes})
    feature_df = feature_df.merge(orth_df, on="canonical_gene_id", how="left")
    feature_df = feature_df.merge(expr_df, on="canonical_gene_id", how="left")
    feature_df = feature_df.merge(sub_df, on="canonical_gene_id", how="left")
    feature_df = feature_df.fillna(0.0)
    degree_edges = edge_df[["source_canonical_gene_id", "target_canonical_gene_id"]].rename(
        columns={"source_canonical_gene_id": "A", "target_canonical_gene_id": "B"}
    )
    feature_df = append_degree_block(feature_df, degree_edges)

    x_cols = [c for c in feature_df.columns if c not in ["legacy_gene_id", "canonical_gene_id"]]
    x = feature_df[x_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
    x = zscore_matrix(x).astype(np.float32)
    mapping = dict(zip(genes, range(len(genes))))
    edge_index = np.vectorize(mapping.__getitem__)(edge_df[["source_canonical_gene_id", "target_canonical_gene_id"]].to_numpy(dtype=object)).astype(np.int64)

    label_manifest = pd.DataFrame({"legacy_gene_id": genes, "canonical_gene_id": genes})
    label_manifest["label"] = label_manifest["canonical_gene_id"].map(dict(zip(label_frame["canonical_gene_id"], label_frame["label"])))
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
        for gene in frame["canonical_gene_id"].astype(str):
            split_map[gene] = split_name
    label_manifest["split"] = label_manifest["canonical_gene_id"].map(split_map).fillna("")

    node_manifest = pd.DataFrame({"legacy_gene_id": genes, "canonical_gene_id": genes})
    node_manifest["alignment_status"] = "aligned"
    node_manifest["in_ppi"] = node_manifest["canonical_gene_id"].isin(set(edge_df["source_canonical_gene_id"]) | set(edge_df["target_canonical_gene_id"]))
    node_manifest["legacy_species"] = "fgraminearum"

    schema = build_feature_schema(expr_df.rename(columns={"canonical_gene_id": "legacy_gene_id"}), orth_df.rename(columns={"canonical_gene_id": "legacy_gene_id"}), sub_df.rename(columns={"canonical_gene_id": "legacy_gene_id"}), include_degree=True)

    np.save(os.path.join(out_dir, "feature_matrix.npy"), x)
    np.save(os.path.join(out_dir, "edge_index.npy"), edge_index)
    node_manifest.to_csv(os.path.join(out_dir, "node_manifest.tsv"), sep="\t", index=False)
    label_manifest.to_csv(os.path.join(out_dir, "label_manifest.tsv"), sep="\t", index=False)
    edge_df.to_csv(os.path.join(out_dir, "edge_table.tsv"), sep="\t", index=False)
    schema.to_csv(os.path.join(out_dir, "feature_schema.tsv"), sep="\t", index=False)
    canonical_audit = pd.concat(
        [
            label_map.assign(source_block="labels"),
            ppi_map.assign(source_block="ppi"),
        ],
        ignore_index=True,
        sort=False,
    )
    canonical_audit.to_csv(os.path.join(out_dir, "canonical_mapping_audit.tsv"), sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "total_nodes_raw": len(genes),
                "canonical_nodes_final": len(node_manifest),
                "labeled_nodes": len(labeled),
                "positive_nodes": int((labeled["label"].astype(float) == 1).sum()),
                "negative_nodes": int((labeled["label"].astype(float) == 0).sum()),
                "feature_dim": x.shape[1],
                "ppi_edges_final": len(edge_df),
            }
        ]
    ).to_csv(os.path.join(out_dir, "dataset_alignment_audit.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "dataset_summary.md"), "w", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "# Fusarium Canonical Replay Summary",
                    "- final canonical nodes = {}".format(len(node_manifest)),
                    "- labeled nodes = {}".format(len(labeled)),
                    "- positive nodes = {}".format(int((labeled["label"].astype(float) == 1).sum())),
                    "- negative nodes = {}".format(int((labeled["label"].astype(float) == 0).sum())),
                    "- feature dimension = {}".format(x.shape[1]),
                    "- ppi edges final = {}".format(len(edge_df)),
                ]
            )
        )
    return out_dir


def main():
    args = parse_args()
    config = load_yaml(args.config)
    out_dir = build_dataset(config)
    print("Fusarium canonical legacy dataset built:", out_dir)


if __name__ == "__main__":
    main()
