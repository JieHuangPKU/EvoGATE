from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.schemas.graph_schema import validate_graph_manifest_file


DEFAULT_EDGE_WEIGHT_COLUMNS = [
    "combined_score",
    "edge_weight",
    "weight",
    "coexpression",
    "experimental",
    "database",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Adapt raw gene-level graph assets into canonical ProGATE_v2 tables")
    parser.add_argument("--config", type=str, required=True, help="Path to graph-ready YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict[str, Any]:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required input table not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")
    return pd.read_csv(path, dtype=str).fillna("")


def read_graph_asset_csv(path: str | Path) -> pd.DataFrame:
    return _read_csv(Path(path))


def _valid_feature_ids(feature_manifest_df: pd.DataFrame) -> set[str]:
    valid_features = feature_manifest_df[
        feature_manifest_df["feature_path"].astype(str).str.strip().ne("")
        & feature_manifest_df["exists"].astype(str).str.lower().isin({"true", "1", "yes"})
        & ~feature_manifest_df["needs_manual_review"].astype(str).str.lower().isin({"true", "1", "yes"})
    ].copy()
    return set(valid_features["canonical_gene_id"].astype(str))


def _load_species_baseline_nodes(config: dict[str, Any], species: str) -> pd.DataFrame:
    baseline_dir = Path(config["paths"]["baseline_dataset_dir"])
    support_df = _read_tsv(baseline_dir / "support_supervised_samples.tsv")
    inference_df = _read_tsv(baseline_dir / "fgraminearum_inference_pool.tsv")
    if species == config["graph"]["target_species"]:
        return inference_df[inference_df["species"] == species].copy()
    return support_df[support_df["species"] == species].copy()


def _load_benchmark_sets(config: dict[str, Any]) -> dict[str, set[str]]:
    baseline_dir = Path(config["paths"]["baseline_dataset_dir"])
    return {
        "broad79": set(_read_tsv(baseline_dir / "fgraminearum_broad79.tsv")["canonical_gene_id"].astype(str)),
        "strict29": set(_read_tsv(baseline_dir / "fgraminearum_strict29.tsv")["canonical_gene_id"].astype(str)),
        "conflict8": set(_read_tsv(baseline_dir / "fgraminearum_conflict8.tsv")["canonical_gene_id"].astype(str)),
    }


def _node_membership(canonical_gene_id: str, benchmark_sets: dict[str, set[str]]) -> str:
    memberships = [name for name, members in benchmark_sets.items() if canonical_gene_id in members]
    return "|".join(memberships)


def _canonical_suffix_map(canonical_ids: pd.Series) -> dict[str, str]:
    suffix_to_canonical: dict[str, str] = {}
    for canonical_id in canonical_ids.astype(str):
        suffix = canonical_id.split("::", 1)[1] if "::" in canonical_id else canonical_id
        suffix_to_canonical[suffix] = canonical_id
    return suffix_to_canonical


def _normalize_asset_node_id(value: str) -> list[str]:
    value = str(value).strip()
    if not value:
        return []
    candidates = [value]
    if "." in value:
        candidates.append(value.split(".")[-1])
    return list(dict.fromkeys(candidate for candidate in candidates if candidate))


def _pick_node_columns(raw_df: pd.DataFrame, suffix_to_canonical: dict[str, str]) -> tuple[str, str]:
    candidate_pairs = [
        ("protein1", "protein2"),
        ("A", "B"),
        ("source", "target"),
        ("gene_a", "gene_b"),
    ]
    best_pair: tuple[str, str] | None = None
    best_overlap = -1
    for source_col, target_col in candidate_pairs:
        if source_col not in raw_df.columns or target_col not in raw_df.columns:
            continue
        values = set(raw_df[source_col].astype(str)) | set(raw_df[target_col].astype(str))
        overlap = 0
        for value in values:
            normalized = _normalize_asset_node_id(value)
            if any(candidate in suffix_to_canonical for candidate in normalized):
                overlap += 1
        if overlap > best_overlap:
            best_overlap = overlap
            best_pair = (source_col, target_col)
    if best_pair is None:
        raise ValueError("No supported source/target node columns were found in the raw graph asset")
    return best_pair


def inspect_raw_gene_graph_asset(
    graph_file_path: str | Path,
    canonical_suffixes: set[str] | None = None,
) -> dict[str, Any]:
    raw_df = read_graph_asset_csv(graph_file_path)
    suffix_to_canonical = {suffix: suffix for suffix in sorted(canonical_suffixes or set())}
    source_col, target_col = _pick_node_columns(raw_df, suffix_to_canonical)
    values = pd.unique(pd.concat([raw_df[source_col], raw_df[target_col]], ignore_index=True))
    normalized = {str(value): _normalize_asset_node_id(str(value)) for value in values}
    overlap = 0
    if canonical_suffixes:
        for candidates in normalized.values():
            if any(candidate in canonical_suffixes for candidate in candidates):
                overlap += 1
    return {
        "raw_df": raw_df,
        "source_col": source_col,
        "target_col": target_col,
        "asset_node_count": int(len(values)),
        "asset_nodes": [str(value) for value in values],
        "normalized_lookup": normalized,
        "canonical_overlap_count": int(overlap),
    }


def _resolve_weight_column(raw_df: pd.DataFrame, config: dict[str, Any]) -> str | None:
    preferred = list(config.get("adapter", {}).get("edge_weight_columns", DEFAULT_EDGE_WEIGHT_COLUMNS))
    for column in preferred:
        if column in raw_df.columns:
            return column
    return None


def _bool(value: Any) -> bool:
    return bool(value)


def _build_placeholder_gene_graph_nodes(
    species_nodes_df: pd.DataFrame,
    feature_ids: set[str],
    benchmark_sets: dict[str, set[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    nodes_df = species_nodes_df[["species", "canonical_gene_id", "gold_label", "label_status"]].drop_duplicates().copy()
    nodes_df["has_feature"] = nodes_df["canonical_gene_id"].isin(feature_ids)
    nodes_df["has_label"] = (
        nodes_df["gold_label"].astype(str).isin({"0", "1"}) & nodes_df["label_status"].astype(str).eq("gold")
    )
    nodes_df["benchmark_membership"] = nodes_df["canonical_gene_id"].map(lambda value: _node_membership(value, benchmark_sets))
    nodes_df["node_source"] = "baseline_gene_graph_placeholder"

    unmatched_rows = nodes_df[~nodes_df["has_feature"]].copy()
    unmatched_rows = unmatched_rows.assign(
        asset_node_id="",
        normalized_node_id="",
        reason="baseline_node_without_feature",
        occurrence_count=1,
    )[
        ["species", "canonical_gene_id", "asset_node_id", "normalized_node_id", "reason", "occurrence_count"]
    ]

    edges_df = pd.DataFrame(
        columns=[
            "source_canonical_gene_id",
            "target_canonical_gene_id",
            "edge_weight",
            "edge_type",
            "graph_id",
        ]
    )
    metadata = {
        "raw_asset_node_count": 0,
        "mapped_asset_node_count": int(len(nodes_df)),
        "final_edge_count": 0,
        "drop_reason_counts": {"placeholder_graph_has_no_edges": 1},
        "selected_node_columns": [],
        "selected_weight_column": "",
    }
    return nodes_df, edges_df, unmatched_rows, metadata


def _build_ppi_graph_dataset(
    graph_row: pd.Series,
    species_nodes_df: pd.DataFrame,
    feature_ids: set[str],
    benchmark_sets: dict[str, set[str]],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    raw_graph_path = Path(graph_row["graph_file_path"])
    print(f"[gene_graph_adapter] load raw ppi asset: {raw_graph_path}")
    raw_df = _read_csv(raw_graph_path)
    print(f"[gene_graph_adapter] raw ppi rows={len(raw_df)}")
    source_col, target_col = _pick_node_columns(raw_df, _canonical_suffix_map(species_nodes_df["canonical_gene_id"]))
    weight_col = _resolve_weight_column(raw_df, config)
    print(f"[gene_graph_adapter] selected node columns=({source_col}, {target_col}) weight_column={weight_col or 'NONE'}")

    suffix_to_canonical = _canonical_suffix_map(species_nodes_df["canonical_gene_id"])
    baseline_canonical_ids = set(species_nodes_df["canonical_gene_id"].astype(str))
    raw_edges = raw_df[[source_col, target_col] + ([weight_col] if weight_col else [])].copy()
    raw_edges = raw_edges.rename(columns={source_col: "source_asset_node_id", target_col: "target_asset_node_id"})
    raw_edges["edge_weight"] = (
        pd.to_numeric(raw_edges[weight_col], errors="coerce").fillna(1.0) if weight_col else 1.0
    )

    def map_asset_id(asset_id: str) -> str:
        for candidate in _normalize_asset_node_id(asset_id):
            if candidate in suffix_to_canonical:
                return suffix_to_canonical[candidate]
        return ""

    raw_edges["source_canonical_gene_id"] = raw_edges["source_asset_node_id"].map(map_asset_id)
    raw_edges["target_canonical_gene_id"] = raw_edges["target_asset_node_id"].map(map_asset_id)
    raw_edges["drop_reason"] = ""

    missing_endpoint = raw_edges["source_canonical_gene_id"].eq("") | raw_edges["target_canonical_gene_id"].eq("")
    raw_edges.loc[missing_endpoint, "drop_reason"] = "unmatched_source_or_target"

    drop_self_loops = _bool(config.get("adapter", {}).get("drop_self_loops", True))
    if drop_self_loops:
        self_loops = (
            raw_edges["drop_reason"].eq("")
            & raw_edges["source_canonical_gene_id"].eq(raw_edges["target_canonical_gene_id"])
        )
        raw_edges.loc[self_loops, "drop_reason"] = "self_loop"

    mapped_nodes = pd.unique(
        pd.concat(
            [
                raw_edges.loc[raw_edges["source_canonical_gene_id"].ne(""), "source_canonical_gene_id"],
                raw_edges.loc[raw_edges["target_canonical_gene_id"].ne(""), "target_canonical_gene_id"],
            ],
            ignore_index=True,
        )
    )
    mapped_nodes_df = pd.DataFrame({"canonical_gene_id": mapped_nodes})
    mapped_nodes_df["species"] = graph_row["species"]
    mapped_nodes_df["has_feature"] = mapped_nodes_df["canonical_gene_id"].isin(feature_ids)

    label_lookup = species_nodes_df.set_index("canonical_gene_id")[["gold_label", "label_status"]].copy()
    mapped_nodes_df["gold_label"] = mapped_nodes_df["canonical_gene_id"].map(label_lookup["gold_label"]).fillna("")
    mapped_nodes_df["label_status"] = mapped_nodes_df["canonical_gene_id"].map(label_lookup["label_status"]).fillna("")
    mapped_nodes_df["has_label"] = (
        mapped_nodes_df["gold_label"].astype(str).isin({"0", "1"})
        & mapped_nodes_df["label_status"].astype(str).eq("gold")
    )
    mapped_nodes_df["benchmark_membership"] = mapped_nodes_df["canonical_gene_id"].map(
        lambda value: _node_membership(value, benchmark_sets)
    )
    mapped_nodes_df["node_source"] = "raw_ppi_graph_asset"

    require_feature_for_edges = _bool(config.get("adapter", {}).get("require_feature_for_edges", True))
    feature_lookup = mapped_nodes_df.set_index("canonical_gene_id")["has_feature"].to_dict()
    if require_feature_for_edges:
        missing_feature_edge = raw_edges["drop_reason"].eq("") & (
            ~raw_edges["source_canonical_gene_id"].map(feature_lookup).fillna(False)
            | ~raw_edges["target_canonical_gene_id"].map(feature_lookup).fillna(False)
        )
        raw_edges.loc[missing_feature_edge, "drop_reason"] = "missing_feature_source_or_target"

    kept_edges = raw_edges[raw_edges["drop_reason"].eq("")].copy()
    deduplicate_undirected = _bool(config.get("adapter", {}).get("deduplicate_undirected_edges", True))
    if deduplicate_undirected and not kept_edges.empty:
        ordered = kept_edges.apply(
            lambda row: sorted([row["source_canonical_gene_id"], row["target_canonical_gene_id"]]),
            axis=1,
            result_type="expand",
        )
        kept_edges["source_canonical_gene_id"] = ordered[0]
        kept_edges["target_canonical_gene_id"] = ordered[1]
        kept_edges = (
            kept_edges.groupby(["source_canonical_gene_id", "target_canonical_gene_id"], as_index=False)["edge_weight"]
            .max()
        )
    else:
        kept_edges = kept_edges[["source_canonical_gene_id", "target_canonical_gene_id", "edge_weight"]].copy()

    kept_edges["edge_type"] = graph_row["graph_type"]
    kept_edges["graph_id"] = graph_row["graph_id"]

    raw_asset_nodes = pd.unique(
        pd.concat([raw_edges["source_asset_node_id"], raw_edges["target_asset_node_id"]], ignore_index=True)
    )
    unmatched_rows: list[dict[str, Any]] = []
    for asset_node_id in raw_asset_nodes:
        normalized = _normalize_asset_node_id(asset_node_id)
        canonical = ""
        for candidate in normalized:
            if candidate in suffix_to_canonical:
                canonical = suffix_to_canonical[candidate]
                break
        if not canonical:
            unmatched_rows.append(
                {
                    "species": graph_row["species"],
                    "canonical_gene_id": "",
                    "asset_node_id": asset_node_id,
                    "normalized_node_id": "|".join(normalized),
                    "reason": "asset_node_without_canonical_match",
                    "occurrence_count": 1,
                }
            )

    asset_canonical_set = set(mapped_nodes_df["canonical_gene_id"].astype(str))
    for canonical_gene_id in sorted(baseline_canonical_ids - asset_canonical_set):
        unmatched_rows.append(
            {
                "species": graph_row["species"],
                "canonical_gene_id": canonical_gene_id,
                "asset_node_id": "",
                "normalized_node_id": canonical_gene_id.split("::", 1)[1] if "::" in canonical_gene_id else canonical_gene_id,
                "reason": "baseline_canonical_not_present_in_graph_asset",
                "occurrence_count": 1,
            }
        )

    feature_missing_nodes = mapped_nodes_df[~mapped_nodes_df["has_feature"]]["canonical_gene_id"].astype(str).tolist()
    for canonical_gene_id in feature_missing_nodes:
        unmatched_rows.append(
            {
                "species": graph_row["species"],
                "canonical_gene_id": canonical_gene_id,
                "asset_node_id": "",
                "normalized_node_id": canonical_gene_id.split("::", 1)[1] if "::" in canonical_gene_id else canonical_gene_id,
                "reason": "mapped_canonical_without_feature",
                "occurrence_count": 1,
            }
        )

    unmatched_df = pd.DataFrame(unmatched_rows)
    if unmatched_df.empty:
        unmatched_df = pd.DataFrame(
            columns=["species", "canonical_gene_id", "asset_node_id", "normalized_node_id", "reason", "occurrence_count"]
        )

    metadata = {
        "raw_asset_node_count": int(len(raw_asset_nodes)),
        "mapped_asset_node_count": int(len(mapped_nodes_df)),
        "feature_ready_mapped_node_count": int(mapped_nodes_df["has_feature"].sum()),
        "final_edge_count": int(len(kept_edges)),
        "drop_reason_counts": {
            key: int(value)
            for key, value in raw_edges.loc[raw_edges["drop_reason"].ne(""), "drop_reason"].value_counts().sort_index().items()
        },
        "selected_node_columns": [source_col, target_col],
        "selected_weight_column": weight_col or "",
    }
    return mapped_nodes_df, kept_edges, unmatched_df, metadata


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    output_dir = Path(config["paths"]["graph_ready_output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_df = validate_graph_manifest_file(output_dir / "graph_manifest.tsv")

    adapter_cfg = dict(config.get("adapter", {}))
    primary_graph_id = str(adapter_cfg.get("primary_graph_id", "fgraminearum__ppi_graph"))
    graph_rows = manifest_df[manifest_df["graph_id"].astype(str) == primary_graph_id].copy()
    if graph_rows.empty:
        raise ValueError(f"primary_graph_id not found in graph manifest: {primary_graph_id}")
    graph_row = graph_rows.iloc[0]
    print(f"[gene_graph_adapter] graph_id={graph_row['graph_id']} graph_type={graph_row['graph_type']} species={graph_row['species']}")
    if graph_row["graph_type"] not in {"gene_graph", "ppi_graph"}:
        raise ValueError(
            f"gene_graph_adapter only supports gene_graph or ppi_graph, got: {graph_row['graph_type']}"
        )

    species = str(graph_row["species"])
    species_nodes_df = _load_species_baseline_nodes(config, species)
    feature_manifest_df = _read_tsv(Path(config["paths"]["feature_manifest_path"]))
    feature_ids = _valid_feature_ids(feature_manifest_df)
    benchmark_sets = _load_benchmark_sets(config)

    if graph_row["graph_type"] == "gene_graph":
        nodes_df, edges_df, unmatched_df, metadata = _build_placeholder_gene_graph_nodes(
            species_nodes_df,
            feature_ids,
            benchmark_sets,
        )
    else:
        nodes_df, edges_df, unmatched_df, metadata = _build_ppi_graph_dataset(
            graph_row,
            species_nodes_df,
            feature_ids,
            benchmark_sets,
            config,
        )
    print(
        f"[gene_graph_adapter] adapted nodes={len(nodes_df)} edges={len(edges_df)} "
        f"unmatched={len(unmatched_df)}"
    )

    nodes_output = Path(adapter_cfg["output_nodes_path"])
    edges_output = Path(adapter_cfg["output_edges_path"])
    unmatched_output = Path(adapter_cfg["output_unmatched_nodes_path"])
    summary_output = Path(adapter_cfg["output_summary_path"])
    metadata_output = Path(adapter_cfg["output_metadata_path"])

    nodes_df = nodes_df[
        ["species", "canonical_gene_id", "has_feature", "has_label", "benchmark_membership", "node_source"]
    ].copy()
    nodes_df.to_csv(nodes_output, sep="\t", index=False)
    edges_df.to_csv(edges_output, sep="\t", index=False)
    unmatched_df.to_csv(unmatched_output, sep="\t", index=False)

    broad_count = int(nodes_df["benchmark_membership"].astype(str).str.contains("broad79", regex=False).sum())
    strict_count = int(nodes_df["benchmark_membership"].astype(str).str.contains("strict29", regex=False).sum())
    conflict_count = int(nodes_df["benchmark_membership"].astype(str).str.contains("conflict8", regex=False).sum())
    unmatched_reason_counts = unmatched_df["reason"].astype(str).value_counts().sort_index().to_dict()

    summary_lines = [
        "# Gene Graph Summary",
        "",
        "## Direct Answers",
        f"- graph_id: {graph_row['graph_id']}",
        f"- graph_type: {graph_row['graph_type']}",
        f"- species: {species}",
        "",
        "## Node Coverage",
        f"- total nodes: {len(nodes_df)}",
        f"- canonical-aligned nodes: {metadata.get('mapped_asset_node_count', len(nodes_df))}",
        f"- nodes with pooled features: {int(nodes_df['has_feature'].astype(bool).sum())}",
        f"- nodes with gold labels: {int(nodes_df['has_label'].astype(bool).sum())}",
        f"- broad79 nodes: {broad_count}",
        f"- strict29 nodes: {strict_count}",
        f"- conflict8 nodes: {conflict_count}",
        "",
        "## Edge Coverage",
        f"- final edges: {len(edges_df)}",
        f"- selected node columns: {metadata.get('selected_node_columns', [])}",
        f"- selected weight column: {metadata.get('selected_weight_column', '')}",
        "",
        "## Unmatched / Dropped Reasons",
    ]
    if not unmatched_reason_counts and not metadata.get("drop_reason_counts"):
        summary_lines.append("- none")
    else:
        for reason, count in unmatched_reason_counts.items():
            summary_lines.append(f"- unmatched::{reason}: {count}")
        for reason, count in metadata.get("drop_reason_counts", {}).items():
            summary_lines.append(f"- dropped_edge::{reason}: {count}")

    summary_output.write_text("\n".join(summary_lines), encoding="utf-8")
    metadata_output.write_text(
        json.dumps(
            {
                "graph_id": graph_row["graph_id"],
                "graph_type": graph_row["graph_type"],
                "species": species,
                "node_count": int(len(nodes_df)),
                "edge_count": int(len(edges_df)),
                "metadata": metadata,
                "unmatched_reason_counts": unmatched_reason_counts,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    print(f"Wrote adapted gene graph nodes to: {nodes_output}")
    print(f"Wrote adapted gene graph edges to: {edges_output}")
    print(f"Wrote unmatched node report to: {unmatched_output}")
    print(f"Wrote gene graph summary to: {summary_output}")


if __name__ == "__main__":
    main()
