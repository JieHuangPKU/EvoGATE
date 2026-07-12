from __future__ import annotations

from pathlib import Path

import pandas as pd


SPECIES_INFO = {
    "scerevisiae": {
        "graph_id": "scerevisiae__ppi_graph",
        "ppi_path": Path("/home/jiehuang/software/fungi/EPGAT/data/essential_genes/yeast/PPI/STRING/string.csv"),
        "source_cols": ["A", "B"],
        "weight_col": "combined_score",
        "bridge_path": Path("data_registry/scerevisiae_ppi_canonical_bridge.tsv"),
        "read_mode": "csv",
    },
    "human": {
        "graph_id": "human__ppi_graph",
        "ppi_path": Path("/home/jiehuang/software/fungi/EPGAT/data/essential_genes/human/PPI/STRING/string.csv"),
        "source_cols": ["A", "B"],
        "weight_col": "combined_score",
        "bridge_path": Path("data_registry/human_ppi_canonical_bridge.tsv"),
        "read_mode": "csv",
    },
    "celegans": {
        "graph_id": "celegans__ppi_graph",
        "ppi_path": Path("/home/jiehuang/software/fungi/EPGAT/data/essential_genes/celegans/PPI/STRING/6239.protein.links.detailed.v12.0.txt"),
        "source_cols": ["protein1", "protein2"],
        "weight_col": "combined_score",
        "bridge_path": Path("data_registry/celegans_ppi_canonical_bridge.tsv"),
        "read_mode": "whitespace",
    },
}


def main() -> None:
    root = Path(".")
    out = root / "outputs" / "support_graphs"
    out.mkdir(parents=True, exist_ok=True)

    audit_rows = []
    summary_lines = [
        "# Support PPI Edge Bridge Summary",
        "",
        "## Direct Answers",
    ]

    for species, info in SPECIES_INFO.items():
        bridge = pd.read_csv(info["bridge_path"], sep="\t", dtype=str).fillna("")
        if info.get("read_mode") == "whitespace":
            raw = pd.read_csv(
                info["ppi_path"],
                sep=r"\s+",
                usecols=info["source_cols"] + [info["weight_col"]],
                dtype=str,
            ).fillna("")
        else:
            raw = pd.read_csv(
                info["ppi_path"],
                usecols=info["source_cols"] + [info["weight_col"]],
                dtype=str,
            ).fillna("")
        raw = raw.rename(
            columns={
                info["source_cols"][0]: "source_raw_node_id",
                info["source_cols"][1]: "target_raw_node_id",
                info["weight_col"]: "edge_weight",
            }
        )
        bridge_map = bridge.set_index("raw_node_id")["proposed_canonical_gene_id"].to_dict()
        status_map = bridge.set_index("raw_node_id")["bridge_status"].to_dict()

        raw["source_canonical_gene_id"] = raw["source_raw_node_id"].map(bridge_map).fillna("")
        raw["target_canonical_gene_id"] = raw["target_raw_node_id"].map(bridge_map).fillna("")
        raw["source_bridge_status"] = raw["source_raw_node_id"].map(status_map).fillna("missing")
        raw["target_bridge_status"] = raw["target_raw_node_id"].map(status_map).fillna("missing")
        raw["bridge_status_pair"] = (
            raw["source_bridge_status"].astype(str) + "|" + raw["target_bridge_status"].astype(str)
        )
        raw.to_csv(out / f"{species}_ppi_edges_bridged.tsv", sep="\t", index=False)

        total_nodes = len(bridge)
        safe_nodes = int(bridge["bridge_status"].isin(["exact", "normalized", "bridged", "alias"]).sum())
        ambiguous_nodes = int(bridge["bridge_status"].eq("ambiguous").sum())
        missing_nodes = int(bridge["bridge_status"].eq("missing").sum())
        conflict_nodes = int(bridge["bridge_status"].eq("conflict").sum())
        safe_ratio = safe_nodes / total_nodes if total_nodes else 0.0

        edge_rows_total = len(raw)
        edge_rows_fully_bridged = int(
            ((raw["source_canonical_gene_id"] != "") & (raw["target_canonical_gene_id"] != "")).sum()
        )
        edge_rows_partially_bridged = int(
            ((raw["source_canonical_gene_id"] != "") ^ (raw["target_canonical_gene_id"] != "")).sum()
        )
        edge_rows_unbridged = int(
            ((raw["source_canonical_gene_id"] == "") & (raw["target_canonical_gene_id"] == "")).sum()
        )

        audit_rows.append(
            {
                "species": species,
                "raw_ppi_node_count": total_nodes,
                "safe_bridge_count": safe_nodes,
                "safe_bridge_ratio": safe_ratio,
                "ambiguous_count": ambiguous_nodes,
                "missing_count": missing_nodes,
                "conflict_count": conflict_nodes,
                "edge_rows_total": edge_rows_total,
                "edge_rows_fully_bridged": edge_rows_fully_bridged,
                "edge_rows_partially_bridged": edge_rows_partially_bridged,
                "edge_rows_unbridged": edge_rows_unbridged,
            }
        )
        summary_lines.append(
            f"- {species}: raw_nodes={total_nodes}, safe_bridge_ratio={safe_ratio}, "
            f"ambiguous={ambiguous_nodes}, missing={missing_nodes}, fully_bridged_edges={edge_rows_fully_bridged}"
        )

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(out / "support_ppi_bridge_audit.tsv", sep="\t", index=False)

    best_species = audit_df.sort_values(
        ["safe_bridge_ratio", "edge_rows_fully_bridged"],
        ascending=[False, False],
        kind="stable",
    ).iloc[0]["species"]
    hardest_species = audit_df.sort_values(
        ["safe_bridge_ratio", "missing_count", "conflict_count"],
        ascending=[True, False, False],
        kind="stable",
    ).iloc[0]["species"]
    summary_lines.extend(
        [
            "",
            "## Interpretation",
            f"- closest to edge-adapter-ready: {best_species}",
            f"- hardest to bridge: {hardest_species}",
            "- `human` is already edge-adapter-ready, `scerevisiae` is near-ready, and `celegans` remains the only unresolved support-species bridge blocker.",
        ]
    )
    (out / "support_ppi_edge_bridge_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")


if __name__ == "__main__":
    main()
