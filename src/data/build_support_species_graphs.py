from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.data.gene_graph_adapter import _load_species_baseline_nodes, _valid_feature_ids, inspect_raw_gene_graph_asset
from src.features.load_embeddings import load_embedding_index
from src.train.train_baseline import load_config


def _read_manifest(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def main() -> None:
    config = load_config("configs/graph_ready.yaml")
    output_dir = Path("outputs/support_graphs")
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_manifest_df = _read_manifest(Path(config["paths"]["feature_manifest_path"]))
    feature_ids = _valid_feature_ids(feature_manifest_df)
    # Reuse the graph-ready manifest as the source of raw asset declarations.
    graph_manifest = _read_manifest(Path(config["paths"]["graph_ready_output_dir"]) / "graph_manifest.tsv")

    unmatched_rows = []
    summary_lines = [
        "# Support Graph Summary",
        "",
        "## Direct Answers",
    ]
    edge_variant_rows = []

    for species in config["graph"]["support_species"]:
        species_nodes = _load_species_baseline_nodes(config, species).copy()
        species_nodes = species_nodes.drop_duplicates(subset=["canonical_gene_id"], keep="first")
        nodes_df = species_nodes[["species", "canonical_gene_id", "gold_label", "label_status"]].copy()
        nodes_df["has_feature"] = nodes_df["canonical_gene_id"].isin(feature_ids)
        nodes_df["has_label"] = (
            nodes_df["gold_label"].astype(str).isin({"0", "1"})
            & nodes_df["label_status"].astype(str).eq("gold")
        )
        nodes_df["benchmark_membership"] = ""
        nodes_df["node_source"] = "baseline_support_gene_graph_placeholder"
        nodes_output = output_dir / f"{species}_gene_graph_nodes.tsv"
        edges_output = output_dir / f"{species}_gene_graph_edges.tsv"
        nodes_df[
            ["species", "canonical_gene_id", "has_feature", "has_label", "benchmark_membership", "node_source"]
        ].to_csv(nodes_output, sep="\t", index=False)
        pd.DataFrame(
            columns=[
                "source_canonical_gene_id",
                "target_canonical_gene_id",
                "edge_weight",
                "edge_type",
                "graph_id",
            ]
        ).to_csv(edges_output, sep="\t", index=False)

        suffixes = set(nodes_df["canonical_gene_id"].astype(str).str.split("::", n=1).str[1])
        ppi_rows = graph_manifest[
            (graph_manifest["species"].astype(str) == species)
            & (graph_manifest["graph_type"].astype(str) == "ppi_graph")
            & (graph_manifest["graph_file_path"].astype(str).str.strip() != "")
        ].copy()

        usable_ppi = False
        best_overlap = 0
        best_source_path = ""
        for _, row in ppi_rows.iterrows():
            inspection = inspect_raw_gene_graph_asset(row["graph_file_path"], suffixes)
            overlap = int(inspection["canonical_overlap_count"])
            if overlap > best_overlap:
                best_overlap = overlap
                best_source_path = row["graph_file_path"]
            reason = (
                "asset_node_without_safe_canonical_bridge"
                if overlap == 0
                else "partial_canonical_overlap_not_promoted_in_this_round"
            )
            for asset_node_id, normalized in inspection["normalized_lookup"].items():
                if overlap > 0 and any(candidate in suffixes for candidate in normalized):
                    continue
                unmatched_rows.append(
                    {
                        "species": species,
                        "graph_id": row["graph_id"],
                        "canonical_gene_id": "",
                        "asset_node_id": asset_node_id,
                        "normalized_node_id": "|".join(normalized),
                        "reason": reason,
                        "occurrence_count": 1,
                    }
                )
            usable_ppi = usable_ppi or overlap > 0

        summary_lines.append(
            f"- {species}: canonical placeholder nodes={len(nodes_df)}, pooled_feature_nodes={int(nodes_df['has_feature'].sum())}, "
            f"ppi_asset_overlap={best_overlap}, ppi_adapter_ready=false"
        )

        for variant_name in ["full_graph", "degree_capped_top10"]:
            edge_variant_rows.append(
                {
                    "species": species,
                    "graph_variant": variant_name,
                    "weighting_scheme": "unweighted" if variant_name == "full_graph" else "normalized_edge_weighted",
                    "edge_count": 0,
                    "node_count": len(nodes_df),
                    "status": "not_applicable",
                    "note": f"No safely canonicalized support-species PPI edge table is available yet for {species}.",
                }
            )

        edge_variant_rows.append(
            {
                "species": species,
                "graph_variant": "degree_capped_top10",
                "weighting_scheme": "topk_neighbor_weighted",
                "edge_count": 0,
                "node_count": len(nodes_df),
                "status": "not_applicable",
                "note": f"No safely canonicalized support-species PPI edge table is available yet for {species}.",
            }
        )

        (output_dir / f"{species}_graph_inventory.json").write_text(
            json.dumps(
                {
                    "species": species,
                    "placeholder_node_count": len(nodes_df),
                    "pooled_feature_nodes": int(nodes_df["has_feature"].sum()),
                    "best_ppi_overlap": best_overlap,
                    "best_ppi_source_path": best_source_path,
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    unmatched_df = pd.DataFrame(unmatched_rows)
    if unmatched_df.empty:
        unmatched_df = pd.DataFrame(
            columns=[
                "species",
                "graph_id",
                "canonical_gene_id",
                "asset_node_id",
                "normalized_node_id",
                "reason",
                "occurrence_count",
            ]
        )
    unmatched_df.to_csv(output_dir / "support_graph_unmatched_nodes.tsv", sep="\t", index=False)

    edge_variant_df = pd.DataFrame(edge_variant_rows)
    edge_variant_df.to_csv(output_dir / "support_edge_variant_summary.tsv", sep="\t", index=False)

    summary_lines.extend(
        [
            "",
            "## Interpretation",
            "- Support-species gene_graph placeholder tables are adapter-ready at the node layer because they are already canonical-gene-id aligned.",
            "- Support-species PPI edge assets remain inventory-ready but not edge-adapter-ready because their node ids still fail a safe bridge into the current v2 canonical contract.",
        ]
    )
    (output_dir / "support_graph_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")

    edge_md_lines = [
        "# Support Edge Variant Summary",
        "",
        "| species | graph_variant | weighting_scheme | edge_count | node_count | status | note |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in edge_variant_df.iterrows():
        edge_md_lines.append(
            f"| {row['species']} | {row['graph_variant']} | {row['weighting_scheme']} | {row['edge_count']} | "
            f"{row['node_count']} | {row['status']} | {row['note']} |"
        )
    (output_dir / "support_edge_variant_summary.md").write_text("\n".join(edge_md_lines), encoding="utf-8")


if __name__ == "__main__":
    main()
