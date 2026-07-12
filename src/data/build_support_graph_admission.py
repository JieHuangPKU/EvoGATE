from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def main() -> None:
    config = load_config("configs/graph_ready.yaml")
    policy = config["support_graph_policy"]

    audit = pd.read_csv("outputs/support_graphs/support_ppi_bridge_audit.tsv", sep="\t")
    audit["safe_bridge_ratio"] = audit["safe_bridge_ratio"].astype(float)
    audit["edge_rows_total"] = audit["edge_rows_total"].astype(int)
    audit["edge_rows_fully_bridged"] = audit["edge_rows_fully_bridged"].astype(int)
    audit["edge_retention_rate"] = audit["edge_rows_fully_bridged"] / audit["edge_rows_total"].where(
        audit["edge_rows_total"] > 0, 1
    )

    # For celegans support-mode admission, use the processed PPI namespace (string.csv A/B)
    # rather than the larger raw STRING internal namespace. This preserves the intended
    # training-facing contract: gene-level graph support over the processed support PPI graph.
    ml = pd.read_csv("data_registry/master_label_table.preliminary.tsv", sep="\t", dtype=str).fillna("")
    ce = ml[ml["species"] == "celegans"].copy()
    wb_cols = ["taxid", "species_name", "WBGene_id", "gene_symbol", "sequence_id"]
    wb = pd.read_csv(
        "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/celegans/PPI/STRING/wormbase.WS240.gene_ids.txt",
        sep="\t",
        comment="#",
        names=wb_cols,
        dtype=str,
    ).fillna("")
    sym_counts = wb[wb["gene_symbol"].ne("")].groupby("gene_symbol")["WBGene_id"].nunique()
    safe_symbols = set(sym_counts[sym_counts == 1].index)
    symbol_to_wb = dict(
        zip(
            wb[wb["gene_symbol"].isin(safe_symbols)]["gene_symbol"].astype(str),
            wb[wb["gene_symbol"].isin(safe_symbols)]["WBGene_id"].astype(str),
        )
    )
    wbgene_to_canonical = dict(zip(ce["raw_gene_id"].astype(str), ce["canonical_gene_id"].astype(str)))
    ppi_ce = pd.read_csv(
        "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/celegans/PPI/STRING/string.csv",
        usecols=["A", "B"],
        dtype=str,
    ).fillna("")
    raw_nodes = sorted(set(ppi_ce["A"].astype(str)) | set(ppi_ce["B"].astype(str)))
    node_map = {}
    resolved = 0
    ambiguous = 0
    missing = 0
    ambiguous_symbols = set(sym_counts[sym_counts > 1].index)
    for raw in raw_nodes:
        if raw in ambiguous_symbols:
            node_map[raw] = ""
            ambiguous += 1
            continue
        wbg = symbol_to_wb.get(raw, "")
        canon = wbgene_to_canonical.get(wbg, "") if wbg else ""
        if canon:
            node_map[raw] = canon
            resolved += 1
        elif wbg:
            node_map[raw] = ""
            ambiguous += 1
        else:
            node_map[raw] = ""
            missing += 1
    src = ppi_ce["A"].map(node_map).fillna("").ne("")
    dst = ppi_ce["B"].map(node_map).fillna("").ne("")
    ce_mask = audit["species"].astype(str).eq("celegans")
    audit.loc[ce_mask, "raw_ppi_node_count"] = len(raw_nodes)
    audit.loc[ce_mask, "safe_bridge_count"] = resolved
    audit.loc[ce_mask, "safe_bridge_ratio"] = resolved / len(raw_nodes) if raw_nodes else 0.0
    audit.loc[ce_mask, "ambiguous_count"] = ambiguous
    audit.loc[ce_mask, "missing_count"] = missing
    audit.loc[ce_mask, "conflict_count"] = 0
    audit.loc[ce_mask, "edge_rows_total"] = len(ppi_ce)
    audit.loc[ce_mask, "edge_rows_fully_bridged"] = int((src & dst).sum())
    audit.loc[ce_mask, "edge_rows_partially_bridged"] = int((src ^ dst).sum())
    audit.loc[ce_mask, "edge_rows_unbridged"] = int((~src & ~dst).sum())
    audit.loc[ce_mask, "edge_retention_rate"] = float((src & dst).sum() / len(ppi_ce)) if len(ppi_ce) else 0.0

    def classify(row):
        if (
            row["safe_bridge_ratio"] >= float(policy["min_graph_complete_node_mapping_rate"])
            and row["edge_retention_rate"] >= float(policy["min_graph_complete_edge_retention"])
        ):
            return "graph-complete-support"
        if (
            row["safe_bridge_ratio"] >= float(policy["min_partial_support_node_mapping_rate"])
            and row["edge_retention_rate"] >= float(policy["min_partial_support_edge_retention"])
            and row["edge_rows_fully_bridged"] > 0
        ):
            return "partial-support-graph"
        return "unusable"

    audit["graph_completeness_class"] = audit.apply(classify, axis=1)
    audit["is_graph_complete_support"] = audit["graph_completeness_class"].eq("graph-complete-support")
    audit["is_partial_support_graph"] = audit["graph_completeness_class"].eq("partial-support-graph")
    audit["is_support_graph_usable"] = (
        audit["is_graph_complete_support"] | (bool(policy["allow_partial_support_graphs"]) & audit["is_partial_support_graph"])
    )

    complete_count = int(audit["is_graph_complete_support"].sum())
    usable_count = int(audit["is_support_graph_usable"].sum())
    graph_aware_training_allowed = complete_count >= 1 and usable_count >= 1

    out = Path("outputs/support_graphs")
    out.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out / "support_graph_admission.tsv", sep="\t", index=False)

    admitted_default = audit[
        audit["graph_completeness_class"].isin(["graph-complete-support", "partial-support-graph"])
    ]["species"].tolist()
    full_species = audit[audit["graph_completeness_class"] == "graph-complete-support"]["species"].tolist()

    lines = [
        "# Support Graph Admission",
        "",
        f"- graph_complete_thresholds: node_mapping>={policy['min_graph_complete_node_mapping_rate']}, edge_retention>={policy['min_graph_complete_edge_retention']}",
        f"- partial_support_thresholds: node_mapping>={policy['min_partial_support_node_mapping_rate']}, edge_retention>={policy['min_partial_support_edge_retention']}",
        f"- allow_partial_support_graphs: {str(policy['allow_partial_support_graphs']).lower()}",
        f"- graph_aware_support_training_allowed: {str(graph_aware_training_allowed).lower()}",
        f"- default_full_support_species: {', '.join(full_species) if full_species else 'none'}",
        f"- default_admitted_support_species: {', '.join(admitted_default) if admitted_default else 'none'}",
        "",
        "## Per-Species Classes",
    ]
    for _, row in audit.iterrows():
        lines.append(
            f"- {row['species']}: class={row['graph_completeness_class']}, "
            f"node_mapping_rate={row['safe_bridge_ratio']}, edge_retention_rate={row['edge_retention_rate']}"
        )
    (out / "support_graph_admission.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
