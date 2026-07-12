"""
Build a slim Phase 2 extended dataset by appending one PLM block onto an
existing Phase 1/1.5 replay dataset.
"""

import argparse
import os

import numpy as np
import pandas as pd
import yaml

from src.features.epgat_extended_features import concatenate_feature_blocks, extend_feature_schema
from src.features.plm_loaders import get_plm_dim, load_h5_embeddings, zero_vector
from src.features.plm_manifest import build_plm_manifest


def parse_args():
    parser = argparse.ArgumentParser(description="Build slim extended EPGAT dataset")
    parser.add_argument("--config", required=True, type=str)
    parser.add_argument("--feature-mode", required=True, choices=["baseline", "baseline_plus_esm2", "baseline_plus_prott5"])
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def map_node_to_plm_key(row, species, mode):
    if species == "fgraminearum_canonical":
        canonical = str(row.get("canonical_gene_id", "")).strip()
        if "::" in canonical:
            return canonical.split("::", 1)[1]
        return canonical
    return str(row.get("legacy_gene_id", "")).strip()


def build_dataset(config, feature_mode):
    source_dir = config["source"]["dataset_dir"]
    out_dir = os.path.join(config["paths"]["output_root"], "{}__{}".format(config["run"]["name"], feature_mode))
    os.makedirs(out_dir, exist_ok=True)

    base_x = np.load(os.path.join(source_dir, "feature_matrix.npy"))
    edge_index = np.load(os.path.join(source_dir, "edge_index.npy"), allow_pickle=True)
    node_manifest = pd.read_csv(os.path.join(source_dir, "node_manifest.tsv"), sep="\t", dtype=str).fillna("")
    label_manifest = pd.read_csv(os.path.join(source_dir, "label_manifest.tsv"), sep="\t", dtype=str).fillna("")
    feature_schema = pd.read_csv(os.path.join(source_dir, "feature_schema.tsv"), sep="\t")
    dataset_audit = pd.read_csv(os.path.join(source_dir, "dataset_alignment_audit.tsv"), sep="\t")

    np.save(os.path.join(out_dir, "edge_index.npy"), edge_index)
    node_manifest.to_csv(os.path.join(out_dir, "node_manifest.tsv"), sep="\t", index=False)
    label_manifest.to_csv(os.path.join(out_dir, "label_manifest.tsv"), sep="\t", index=False)

    plm_rows = []
    coverage_rows = []
    species = config["species"]["name"]

    if feature_mode == "baseline":
        final_x = base_x.astype(np.float32)
        final_schema = feature_schema.copy()
        plm_name = ""
        plm_dim = 0
        coverage = 0.0
    else:
        plm_name = "esm2" if feature_mode.endswith("esm2") else "prott5"
        plm_path = config["plm"][plm_name]["path"]
        embedding_lookup = load_h5_embeddings(plm_path)
        sample_keys = sorted(list(embedding_lookup.keys()))[:5]
        manifest = build_plm_manifest(plm_name, plm_path, config["plm"][plm_name]["key_type"], len(embedding_lookup), sample_keys)
        manifest.to_csv(os.path.join(out_dir, "plm_manifest.tsv"), sep="\t", index=False)

        vectors = []
        present = 0
        for _, row in node_manifest.iterrows():
            raw_key = map_node_to_plm_key(row, species, feature_mode)
            canonical = row.get("canonical_gene_id", "")
            if raw_key in embedding_lookup:
                vec = embedding_lookup[raw_key]
                status = "matched"
                present += 1
            else:
                vec = zero_vector(plm_name)
                status = "absent"
            vectors.append(np.asarray(vec, dtype=np.float32))
            plm_rows.append(
                {
                    "raw_embedding_key": raw_key,
                    "mapped_internal_id": row.get("legacy_gene_id", ""),
                    "mapped_canonical_gene_id": canonical,
                    "mapping_status": status,
                    "embedding_present": status == "matched",
                }
            )
        extra_x = np.vstack(vectors).astype(np.float32)
        final_x = concatenate_feature_blocks(base_x, extra_x)
        final_schema = extend_feature_schema(
            feature_schema.copy(),
            plm_name,
            int(feature_schema["end_col"].max()) + 1,
            get_plm_dim(plm_name),
            plm_path,
        )
        plm_dim = get_plm_dim(plm_name)
        coverage = float(present) / float(len(node_manifest)) if len(node_manifest) else 0.0

    np.save(os.path.join(out_dir, "feature_matrix.npy"), final_x)
    final_schema.to_csv(os.path.join(out_dir, "feature_schema.tsv"), sep="\t", index=False)

    if plm_rows:
        pd.DataFrame(plm_rows).to_csv(os.path.join(out_dir, "plm_mapping_audit.tsv"), sep="\t", index=False)
    else:
        pd.DataFrame(
            [{"raw_embedding_key": "", "mapped_internal_id": "", "mapped_canonical_gene_id": "", "mapping_status": "baseline_no_plm", "embedding_present": False}]
        ).to_csv(os.path.join(out_dir, "plm_mapping_audit.tsv"), sep="\t", index=False)

    if feature_mode == "baseline":
        coverage_rows.append({"species": species, "feature_mode": feature_mode, "plm_name": "", "plm_dim": 0, "num_nodes": len(node_manifest), "num_present": 0, "coverage": 0.0, "key_type": ""})
    else:
        num_present = sum(1 for row in plm_rows if row["embedding_present"])
        coverage_rows.append({"species": species, "feature_mode": feature_mode, "plm_name": plm_name, "plm_dim": plm_dim, "num_nodes": len(node_manifest), "num_present": num_present, "coverage": coverage, "key_type": config["plm"][plm_name]["key_type"]})
    pd.DataFrame(coverage_rows).to_csv(os.path.join(out_dir, "plm_coverage.tsv"), sep="\t", index=False)

    audit = dataset_audit.copy()
    audit["feature_mode"] = feature_mode
    audit["plm_name"] = plm_name
    audit["plm_dim"] = plm_dim
    audit["plm_coverage"] = coverage
    audit.to_csv(os.path.join(out_dir, "dataset_alignment_audit.tsv"), sep="\t", index=False)

    summary_lines = [
        "# Slim Extended Dataset Summary",
        "",
        "- species: {}".format(species),
        "- feature_mode: {}".format(feature_mode),
        "- base_feature_dim: {}".format(base_x.shape[1]),
        "- final_feature_dim: {}".format(final_x.shape[1]),
        "- plm_name: {}".format(plm_name if plm_name else "baseline_only"),
        "- plm_coverage: {:.4f}".format(coverage),
    ]
    with open(os.path.join(out_dir, "dataset_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines))

    return out_dir


def main():
    args = parse_args()
    config = load_yaml(args.config)
    out_dir = build_dataset(config, args.feature_mode)
    print("Slim extended dataset built:", out_dir)


if __name__ == "__main__":
    main()
