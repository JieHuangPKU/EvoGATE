import argparse
from pathlib import Path

import pandas as pd


SPECIES_PREFIX = "fgraminearum::"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Adapt Fusarium eFG PPIs into the frozen-protocol GraphSAGE graph format")
    parser.add_argument("--raw-efg", required=True, type=str)
    parser.add_argument("--split-manifest", required=True, type=str)
    parser.add_argument("--bridge", required=True, type=str)
    parser.add_argument(
        "--confidence-filter",
        default="ALL",
        type=str,
        help="One of ALL, HIGH_MEDIUM, HIGH. Applied before canonical ID mapping.",
    )
    parser.add_argument("--output-graph", required=True, type=str)
    parser.add_argument("--output-summary", required=True, type=str)
    parser.add_argument("--output-mapping", required=True, type=str)
    return parser.parse_args()


def classify_raw_id(token: str) -> str:
    text = str(token).strip()
    if not text:
        return "empty"
    if text.startswith(SPECIES_PREFIX):
        return "canonical_gene_id"
    if text.startswith("FGRAMPH1_"):
        return "ph1_gene_id"
    if text.startswith("FGSG_"):
        return "fgsg_id"
    if text.startswith("XP_"):
        return "protein_accession"
    return "other"


def build_mapping_tables(split_df: pd.DataFrame, bridge_df: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    exact_lookup: dict[str, str] = {}
    bridge_lookup: dict[str, str] = {}
    for row in split_df[["graph_gene_id", "canonical_gene_id"]].drop_duplicates().itertuples(index=False):
        graph_gene_id = str(row.graph_gene_id).strip()
        canonical_gene_id = str(row.canonical_gene_id).strip()
        canonical_stripped = canonical_gene_id.split("::", 1)[1] if "::" in canonical_gene_id else canonical_gene_id
        for key in [graph_gene_id, canonical_gene_id, canonical_stripped]:
            key = str(key).strip()
            if key:
                exact_lookup[key] = graph_gene_id

    resolved_bridge = bridge_df[bridge_df["bridge_status"].astype(str).eq("resolved")].copy()
    for row in resolved_bridge[["source_protein_id", "header_fgsg_id", "resolved_canonical_gene_id"]].itertuples(index=False):
        canonical_gene_id = str(row.resolved_canonical_gene_id).strip()
        if not canonical_gene_id:
            continue
        graph_gene_id = canonical_gene_id.split("::", 1)[1] if "::" in canonical_gene_id else canonical_gene_id
        for key in [row.source_protein_id, row.header_fgsg_id]:
            key = str(key).strip()
            if key:
                bridge_lookup[key] = graph_gene_id
    return exact_lookup, bridge_lookup


def normalize_confidence(value: str) -> str:
    return str(value).strip().upper()


def resolve_confidence_filter(filter_name: str) -> tuple[str, set[str]]:
    normalized = normalize_confidence(filter_name)
    allowed = {
        "ALL": {"HIGH", "MEDIUM", "LOW"},
        "HIGH_MEDIUM": {"HIGH", "MEDIUM"},
        "HIGH": {"HIGH"},
    }
    if normalized not in allowed:
        raise ValueError(f"Unsupported confidence filter '{filter_name}'. Expected one of {sorted(allowed)}")
    return normalized, allowed[normalized]


def map_raw_id(token: str, exact_lookup: dict[str, str], bridge_lookup: dict[str, str]) -> dict[str, str]:
    raw_id = str(token).strip()
    raw_id_type = classify_raw_id(raw_id)
    if not raw_id:
        return {
            "raw_id": raw_id,
            "raw_id_type": raw_id_type,
            "mapped_graph_gene_id": "",
            "mapped_canonical_gene_id": "",
            "mapping_status": "unresolved",
            "mapping_rule": "empty",
        }
    if raw_id in exact_lookup:
        mapped_graph = exact_lookup[raw_id]
        return {
            "raw_id": raw_id,
            "raw_id_type": raw_id_type,
            "mapped_graph_gene_id": mapped_graph,
            "mapped_canonical_gene_id": SPECIES_PREFIX + mapped_graph,
            "mapping_status": "resolved",
            "mapping_rule": "exact_mainline_graph_id",
        }
    if raw_id in bridge_lookup:
        mapped_graph = bridge_lookup[raw_id]
        return {
            "raw_id": raw_id,
            "raw_id_type": raw_id_type,
            "mapped_graph_gene_id": mapped_graph,
            "mapped_canonical_gene_id": SPECIES_PREFIX + mapped_graph,
            "mapping_status": "resolved",
            "mapping_rule": "protein_to_canonical_bridge",
        }
    return {
        "raw_id": raw_id,
        "raw_id_type": raw_id_type,
        "mapped_graph_gene_id": "",
        "mapped_canonical_gene_id": "",
        "mapping_status": "unresolved",
        "mapping_rule": "not_found",
    }


def main() -> None:
    args = parse_args()
    raw_path = Path(args.raw_efg)
    split_path = Path(args.split_manifest)
    bridge_path = Path(args.bridge)
    output_graph = Path(args.output_graph)
    output_summary = Path(args.output_summary)
    output_mapping = Path(args.output_mapping)

    output_graph.parent.mkdir(parents=True, exist_ok=True)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    output_mapping.parent.mkdir(parents=True, exist_ok=True)

    split_df = pd.read_csv(split_path, sep="\t", dtype=str).fillna("")
    bridge_df = pd.read_csv(bridge_path, sep="\t", dtype=str).fillna("")
    raw_df = pd.read_csv(raw_path, sep="\t", header=None, names=["A", "B", "source", "confidence"], dtype=str).fillna("")
    raw_df["confidence"] = raw_df["confidence"].astype(str).map(normalize_confidence)
    confidence_filter, allowed_confidence_values = resolve_confidence_filter(args.confidence_filter)
    observed_confidence_values = sorted({normalize_confidence(value) for value in raw_df["confidence"].astype(str) if str(value).strip()})
    filtered_raw_df = raw_df.loc[raw_df["confidence"].isin(allowed_confidence_values)].copy()

    exact_lookup, bridge_lookup = build_mapping_tables(split_df, bridge_df)

    raw_unique_ids = sorted(set(raw_df["A"].astype(str)) | set(raw_df["B"].astype(str)))
    unique_raw_ids = sorted(set(filtered_raw_df["A"].astype(str)) | set(filtered_raw_df["B"].astype(str)))
    mapping_records = [map_raw_id(raw_id, exact_lookup, bridge_lookup) for raw_id in unique_raw_ids]
    mapping_df = pd.DataFrame(mapping_records)
    mapping_lookup = mapping_df.set_index("raw_id").to_dict(orient="index")

    rows = []
    dropped_unresolved = 0
    dropped_self_loop = 0
    for row in filtered_raw_df.itertuples(index=False):
        map_a = mapping_lookup.get(str(row.A).strip(), {})
        map_b = mapping_lookup.get(str(row.B).strip(), {})
        gene_a = str(map_a.get("mapped_graph_gene_id", "")).strip()
        gene_b = str(map_b.get("mapped_graph_gene_id", "")).strip()
        if not gene_a or not gene_b:
            dropped_unresolved += 1
            continue
        if gene_a == gene_b:
            dropped_self_loop += 1
            continue
        src, dst = sorted([gene_a, gene_b])
        rows.append(
            {
                "A": src,
                "B": dst,
                "source": str(row.source).strip(),
                "confidence": str(row.confidence).strip(),
                "original_A": str(row.A).strip(),
                "original_B": str(row.B).strip(),
                "mapping_rule_A": map_a.get("mapping_rule", ""),
                "mapping_rule_B": map_b.get("mapping_rule", ""),
            }
        )

    edge_df = pd.DataFrame(rows)
    if edge_df.empty:
        raise ValueError("eFG graph adaptation produced zero mapped edges")

    edge_df = (
        edge_df.groupby(["A", "B"], as_index=False)
        .agg(
            source=("source", lambda values: "|".join(sorted({str(v).strip() for v in values if str(v).strip()}))),
            confidence=("confidence", lambda values: "|".join(sorted({str(v).strip() for v in values if str(v).strip()}))),
            original_edge_examples=("original_A", "first"),
            mapping_rule_A=("mapping_rule_A", "first"),
            mapping_rule_B=("mapping_rule_B", "first"),
        )
        .sort_values(["A", "B"], kind="stable")
        .reset_index(drop=True)
    )
    edge_df[["A", "B", "source", "confidence"]].to_csv(output_graph, sep="\t", index=False)
    mapping_df.to_csv(output_mapping, sep="\t", index=False)

    label_nodes = set(split_df["graph_gene_id"].astype(str))
    mapped_nodes = set(edge_df["A"].astype(str)) | set(edge_df["B"].astype(str))
    resolved_nodes = mapping_df[mapping_df["mapping_status"].astype(str).eq("resolved")]
    resolved_by_rule = resolved_nodes["mapping_rule"].value_counts().to_dict()
    summary_row = {
        "raw_efg_path": str(raw_path),
        "raw_edge_count": int(len(raw_df)),
        "filtered_raw_edge_count": int(len(filtered_raw_df)),
        "raw_unique_node_count": int(len(raw_unique_ids)),
        "filtered_raw_unique_node_count": int(len(unique_raw_ids)),
        "confidence_filter": confidence_filter,
        "allowed_confidence_values": "|".join(sorted(allowed_confidence_values)),
        "observed_confidence_values": "|".join(observed_confidence_values),
        "raw_id_type_top": "|".join(
            f"{key}:{value}"
            for key, value in mapping_df["raw_id_type"].value_counts().head(5).to_dict().items()
        ),
        "mapped_edge_count": int(len(edge_df)),
        "mapped_unique_node_count": int(len(mapped_nodes)),
        "dropped_unresolved_edge_count": int(dropped_unresolved),
        "dropped_self_loop_count": int(dropped_self_loop),
        "resolved_node_count": int(len(resolved_nodes)),
        "unresolved_node_count": int((mapping_df["mapping_status"].astype(str) != "resolved").sum()),
        "mapped_nodes_in_mainline_label_space": int(len(mapped_nodes & label_nodes)),
        "label_node_count": int(len(label_nodes)),
        "label_node_coverage_ratio": float(len(mapped_nodes & label_nodes) / len(label_nodes)),
        "mapping_rule_counts": "|".join(f"{key}:{value}" for key, value in resolved_by_rule.items()),
        "mapping_note": "Direct exact graph_gene_id mapping was attempted first; protein_to_canonical_bridge served as fallback.",
    }
    pd.DataFrame([summary_row]).to_csv(output_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
