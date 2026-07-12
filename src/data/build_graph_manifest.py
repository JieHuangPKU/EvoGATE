from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.registry.load_registry import load_registry_bundle
from src.schemas.graph_schema import SUPPORTED_GRAPH_TYPES, validate_graph_manifest_frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build graph-ready manifest tables for ProGATE_v2")
    parser.add_argument("--config", type=str, required=True, help="Path to graph-ready YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict[str, Any]:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required input table not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _path_from_config(value: Any) -> Path | None:
    if value in [None, ""]:
        return None
    text = str(value).strip()
    return Path(text) if text else None


def _bool_string(value: bool) -> str:
    return "true" if bool(value) else "false"


def _canonical_suffix_set(series: pd.Series) -> set[str]:
    suffixes = set()
    for value in series.astype(str):
        if "::" in value:
            suffixes.add(value.split("::", 1)[1])
        else:
            suffixes.add(value)
    return suffixes


def _normalize_asset_id(value: str) -> set[str]:
    value = str(value).strip()
    if not value:
        return set()
    candidates = {value}
    if "::" in value:
        candidates.add(value.split("::", 1)[1])
    if "." in value:
        candidates.add(value.split(".")[-1])
    return {candidate for candidate in candidates if candidate}


def _load_asset_node_ids(graph_type: str, graph_file_path: Path) -> tuple[set[str], str]:
    if not graph_file_path.exists():
        return set(), "graph asset path does not exist"

    if graph_type == "residue_graph":
        if not graph_file_path.is_dir():
            return set(), "residue graph asset must be a directory of raw protein files"
        return {path.stem for path in graph_file_path.glob("*.pt")}, ""

    if graph_type == "ppi_graph":
        df = pd.read_csv(graph_file_path, usecols=lambda column: column in {"A", "B", "protein1", "protein2"}, dtype=str).fillna("")
        if df.empty:
            return set(), "ppi asset is empty"
        values = set()
        for column in df.columns:
            values.update(df[column].astype(str).tolist())
        return values, ""

    if graph_type == "orthology_graph":
        df = pd.read_csv(graph_file_path, usecols=lambda column: column == "Gene", dtype=str).fillna("")
        if "Gene" not in df.columns:
            return set(), "orthology asset lacks Gene column"
        return set(df["Gene"].astype(str).tolist()), ""

    if graph_type == "coexpression_graph":
        df = pd.read_csv(graph_file_path, usecols=lambda column: column in {"A", "B", "gene_a", "gene_b", "source", "target"}, dtype=str).fillna("")
        if df.empty:
            return set(), "coexpression asset has no recognized node columns"
        values = set()
        for column in df.columns:
            values.update(df[column].astype(str).tolist())
        return values, ""

    return set(), ""


def _asset_overlap(asset_ids: set[str], baseline_suffixes: set[str]) -> tuple[int, float]:
    matched = set()
    for asset_id in asset_ids:
        normalized = _normalize_asset_id(asset_id)
        if normalized & baseline_suffixes:
            matched.add(asset_id)
    overlap = len(matched)
    fraction = overlap / len(baseline_suffixes) if baseline_suffixes else 0.0
    return overlap, fraction


def _write_metadata(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _adapted_graph_status(config: dict[str, Any]) -> dict[str, Any]:
    nodes_path = _path_from_config(config["paths"].get("adapted_gene_graph_nodes_path"))
    edges_path = _path_from_config(config["paths"].get("adapted_gene_graph_edges_path"))
    unmatched_path = _path_from_config(config["paths"].get("adapted_gene_graph_unmatched_nodes_path"))
    summary_path = _path_from_config(config["paths"].get("adapted_gene_graph_summary_path"))
    metadata_path = _path_from_config(config["paths"].get("adapted_gene_graph_metadata_path"))

    status = {
        "adapter_graph_id": str(config.get("adapter", {}).get("primary_graph_id", "")),
        "adapted_nodes_path": str(nodes_path) if nodes_path and nodes_path.exists() else "",
        "adapted_edges_path": str(edges_path) if edges_path and edges_path.exists() else "",
        "adapted_unmatched_nodes_path": str(unmatched_path) if unmatched_path and unmatched_path.exists() else "",
        "adapted_summary_path": str(summary_path) if summary_path and summary_path.exists() else "",
        "adapted_metadata_path": str(metadata_path) if metadata_path and metadata_path.exists() else "",
        "adapter_ready": False,
        "adapter_node_count": 0,
        "adapter_edge_count": 0,
        "adapter_unmatched_count": 0,
    }
    if status["adapted_nodes_path"] and status["adapted_edges_path"]:
        status["adapter_ready"] = True
        status["adapter_node_count"] = len(pd.read_csv(status["adapted_nodes_path"], sep="\t", dtype=str).fillna(""))
        status["adapter_edge_count"] = len(pd.read_csv(status["adapted_edges_path"], sep="\t", dtype=str).fillna(""))
        if status["adapted_unmatched_nodes_path"]:
            status["adapter_unmatched_count"] = len(
                pd.read_csv(status["adapted_unmatched_nodes_path"], sep="\t", dtype=str).fillna("")
            )
    return status


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    output_dir = Path(config["paths"]["graph_ready_output_dir"])
    metadata_dir = output_dir / "metadata"
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    registry_dir = Path(config["paths"]["registry_dir"])
    baseline_dir = Path(config["paths"]["baseline_dataset_dir"])
    feature_manifest_path = Path(config["paths"]["feature_manifest_path"])

    bundle = load_registry_bundle(registry_dir)
    support_df = _read_tsv(baseline_dir / "support_supervised_samples.tsv")
    inference_df = _read_tsv(baseline_dir / "fgraminearum_inference_pool.tsv")
    feature_manifest_df = _read_tsv(feature_manifest_path)

    valid_features = feature_manifest_df[
        feature_manifest_df["feature_path"].astype(str).str.strip().ne("")
        & feature_manifest_df["exists"].astype(str).str.lower().isin({"true", "1", "yes"})
        & ~feature_manifest_df["needs_manual_review"].astype(str).str.lower().isin({"true", "1", "yes"})
    ].copy()

    broad_ids = set(bundle.fg_broad79["canonical_gene_id"].astype(str))
    strict_ids = set(bundle.fg_strict29["canonical_gene_id"].astype(str))
    conflict_ids = set(bundle.fg_conflict8["canonical_gene_id"].astype(str))

    target_graph_types = list(config["graph"]["target_graph_types"])
    unknown_graph_types = sorted(set(target_graph_types) - SUPPORTED_GRAPH_TYPES)
    if unknown_graph_types:
        raise ValueError(f"Unsupported target graph types in config: {unknown_graph_types}")

    rows: list[dict[str, Any]] = []
    species_order = list(config["graph"]["support_species"]) + [config["graph"]["target_species"]]
    adapted_status = _adapted_graph_status(config)

    for species in species_order:
        if species == config["graph"]["target_species"]:
            species_nodes_df = inference_df[inference_df["species"] == species].copy()
        else:
            species_nodes_df = support_df[support_df["species"] == species].copy()
        species_feature_df = valid_features[valid_features["species"] == species].copy()

        canonical_ids = set(species_nodes_df["canonical_gene_id"].astype(str))
        baseline_suffixes = _canonical_suffix_set(species_nodes_df["canonical_gene_id"])
        feature_ids = set(species_feature_df["canonical_gene_id"].astype(str))
        feature_ready_count = len(canonical_ids & feature_ids)
        benchmark_overlap = {
            "broad79": len(canonical_ids & broad_ids) if species == "fgraminearum" else 0,
            "strict29": len(canonical_ids & strict_ids) if species == "fgraminearum" else 0,
            "conflict8": len(canonical_ids & conflict_ids) if species == "fgraminearum" else 0,
        }

        ppi_asset_path = _path_from_config(config["graph_assets"].get("ppi", {}).get(species, ""))
        orthology_asset_path = _path_from_config(config["graph_assets"].get("orthology", {}).get(species, ""))
        coexpression_asset_path = _path_from_config(config["graph_assets"].get("coexpression", {}).get(species, ""))
        residue_asset_path = _path_from_config(config["graph_assets"].get("residue", {}).get(species, ""))

        for graph_type in target_graph_types:
            graph_id = f"{species}__{graph_type}"
            metadata_path = metadata_dir / f"{graph_id}.json"
            graph_file_path: Path | None = None
            edge_source = "placeholder"
            node_feature_source = config["graph"]["node_feature_source"]
            node_key_type = "canonical_gene_id"
            is_directed = False
            is_weighted = False
            notes = ""
            asset_ids: set[str] = set()
            missing_reason = ""

            if graph_type == "gene_graph":
                notes = (
                    "Placeholder gene-level graph node universe derived from the baseline dataset; "
                    "no edge asset is required at graph-ready stage."
                )
            elif graph_type == "ppi_graph":
                graph_file_path = ppi_asset_path
                edge_source = "epgat_ppi_candidate"
                is_weighted = True
                notes = "PPI candidate asset discovered from EPGAT; canonical-id compatibility must be checked explicitly."
            elif graph_type == "coexpression_graph":
                graph_file_path = coexpression_asset_path
                edge_source = "coexpression_candidate"
                is_weighted = True
                notes = "Coexpression graph type is schema-ready only in this round."
            elif graph_type == "orthology_graph":
                graph_file_path = orthology_asset_path
                edge_source = "epgat_orthology_candidate"
                notes = "Orthology graph type is cross-species by definition and needs dedicated adapter logic before training."
            elif graph_type == "residue_graph":
                graph_file_path = residue_asset_path
                edge_source = "bingo_contact_map_candidate"
                node_feature_source = "bingo_raw_pt_feature_representation"
                node_key_type = "residue_index"
                notes = (
                    "Residue-level graph assets come from Bingo raw protein `.pt` files and must remain isolated "
                    "from gene-level graph runtime."
                )

            asset_exists = False
            asset_overlap_count = 0
            asset_overlap_fraction = 0.0
            asset_node_count = 0

            if graph_type == "gene_graph":
                asset_exists = True
                asset_overlap_count = len(canonical_ids)
                asset_overlap_fraction = 1.0 if canonical_ids else 0.0
                asset_node_count = len(canonical_ids)
            elif graph_file_path is not None:
                asset_exists = graph_file_path.exists()
                if asset_exists:
                    asset_ids, missing_reason = _load_asset_node_ids(graph_type, graph_file_path)
                    asset_node_count = len(asset_ids)
                    asset_overlap_count, asset_overlap_fraction = _asset_overlap(asset_ids, baseline_suffixes)
                else:
                    missing_reason = "configured graph asset path does not exist"
            else:
                missing_reason = "no graph asset configured"

            schema_ready = True
            compatible_with_baseline = graph_type == "gene_graph" or asset_overlap_count > 0
            runtime_ready = graph_type == "gene_graph" or (
                asset_exists and compatible_with_baseline and graph_type != "residue_graph"
            )

            if graph_type == "residue_graph":
                runtime_ready = False
                if not missing_reason:
                    missing_reason = "residue graph is intentionally schema-ready only in this round"

            if graph_type == "orthology_graph":
                runtime_ready = False
                if not missing_reason:
                    missing_reason = "orthology graph still needs a dedicated cross-species adapter before runtime use"

            if graph_type != "gene_graph" and compatible_with_baseline and graph_type in {"ppi_graph", "orthology_graph"} and species != "fgraminearum":
                notes += " Asset exists, but node ids are only partially canonical-compatible and still require mapping."
            elif graph_type != "gene_graph" and not compatible_with_baseline and not missing_reason:
                missing_reason = "graph asset node ids are not yet canonical-id compatible with the baseline dataset"

            metadata_payload = {
                "species": species,
                "graph_id": graph_id,
                "graph_type": graph_type,
                "baseline_node_count": len(canonical_ids),
                "feature_ready_node_count": feature_ready_count,
                "asset_node_count": asset_node_count,
                "asset_overlap_count": asset_overlap_count,
                "asset_overlap_fraction": asset_overlap_fraction,
                "asset_exists": asset_exists,
                "compatible_with_baseline": compatible_with_baseline,
                "runtime_ready": runtime_ready,
                "benchmark_overlap": benchmark_overlap,
                "graph_file_path": str(graph_file_path) if graph_file_path is not None else "",
                "notes": notes,
                "missing_reason": missing_reason,
            }
            _write_metadata(metadata_path, metadata_payload)

            rows.append(
                {
                    "species": species,
                    "graph_id": graph_id,
                    "graph_type": graph_type,
                    "node_key_type": node_key_type,
                    "edge_source": edge_source,
                    "node_feature_source": node_feature_source,
                    "graph_file_path": str(graph_file_path) if graph_file_path is not None else "",
                    "metadata_path": str(metadata_path),
                    "is_directed": _bool_string(is_directed),
                    "is_weighted": _bool_string(is_weighted),
                    "notes": notes,
                    "asset_exists": _bool_string(asset_exists),
                    "schema_ready": _bool_string(schema_ready),
                    "runtime_ready": _bool_string(runtime_ready),
                    "compatible_with_baseline": _bool_string(compatible_with_baseline),
                    "baseline_node_count": len(canonical_ids),
                    "feature_ready_node_count": feature_ready_count,
                    "asset_node_count": asset_node_count,
                    "asset_overlap_count": asset_overlap_count,
                    "asset_overlap_fraction": asset_overlap_fraction,
                    "benchmark_overlap_broad79": benchmark_overlap["broad79"],
                    "benchmark_overlap_strict29": benchmark_overlap["strict29"],
                    "benchmark_overlap_conflict8": benchmark_overlap["conflict8"],
                    "missing_reason": missing_reason,
                    "adapter_ready": "false",
                    "adapted_nodes_path": "",
                    "adapted_edges_path": "",
                    "adapted_unmatched_nodes_path": "",
                    "adapted_summary_path": "",
                    "adapted_metadata_path": "",
                    "adapter_node_count": 0,
                    "adapter_edge_count": 0,
                    "adapter_unmatched_count": 0,
                }
            )

    manifest_df = pd.DataFrame(rows)
    if adapted_status["adapter_ready"] and adapted_status["adapter_graph_id"] in set(manifest_df["graph_id"].astype(str)):
        mask = manifest_df["graph_id"].astype(str).eq(adapted_status["adapter_graph_id"])
        manifest_df.loc[mask, "adapter_ready"] = "true"
        manifest_df.loc[mask, "adapted_nodes_path"] = adapted_status["adapted_nodes_path"]
        manifest_df.loc[mask, "adapted_edges_path"] = adapted_status["adapted_edges_path"]
        manifest_df.loc[mask, "adapted_unmatched_nodes_path"] = adapted_status["adapted_unmatched_nodes_path"]
        manifest_df.loc[mask, "adapted_summary_path"] = adapted_status["adapted_summary_path"]
        manifest_df.loc[mask, "adapted_metadata_path"] = adapted_status["adapted_metadata_path"]
        manifest_df.loc[mask, "adapter_node_count"] = adapted_status["adapter_node_count"]
        manifest_df.loc[mask, "adapter_edge_count"] = adapted_status["adapter_edge_count"]
        manifest_df.loc[mask, "adapter_unmatched_count"] = adapted_status["adapter_unmatched_count"]
    validate_graph_manifest_frame(manifest_df)
    missing_df = manifest_df[
        ~manifest_df["runtime_ready"].astype(str).str.lower().isin({"true", "1", "yes"})
        | manifest_df["missing_reason"].astype(str).str.strip().ne("")
    ].copy()

    manifest_df.to_csv(output_dir / "graph_manifest.tsv", sep="\t", index=False)
    missing_df.to_csv(output_dir / "graph_manifest_missing.tsv", sep="\t", index=False)

    print(f"Wrote graph manifest to: {output_dir / 'graph_manifest.tsv'}")
    print(f"Wrote graph missing report to: {output_dir / 'graph_manifest_missing.tsv'}")
    print(f"Graph manifest rows: {len(manifest_df)}")
    print(f"Runtime-ready rows: {int(manifest_df['runtime_ready'].astype(str).str.lower().isin(['true', '1', 'yes']).sum())}")


if __name__ == "__main__":
    main()
