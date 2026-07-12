from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from src.schemas.sample_schema import ALLOWED_SPECIES


SUPPORTED_GRAPH_TYPES = {
    "gene_graph",
    "ppi_graph",
    "coexpression_graph",
    "orthology_graph",
    "residue_graph",
}

SUPPORTED_NODE_KEY_TYPES = {
    "canonical_gene_id",
    "residue_index",
}

GENE_LEVEL_GRAPH_TYPES = {
    "gene_graph",
    "ppi_graph",
    "coexpression_graph",
    "orthology_graph",
}

REQUIRED_GRAPH_MANIFEST_COLUMNS = [
    "species",
    "graph_id",
    "graph_type",
    "node_key_type",
    "edge_source",
    "node_feature_source",
    "graph_file_path",
    "metadata_path",
    "is_directed",
    "is_weighted",
    "notes",
]


@dataclass
class GraphManifestRecord:
    species: str
    graph_id: str
    graph_type: str
    node_key_type: str
    edge_source: str
    node_feature_source: str
    graph_file_path: str
    metadata_path: str
    is_directed: bool
    is_weighted: bool
    notes: str

    def to_dict(self) -> dict:
        return asdict(self)


def is_gene_level_graph(graph_type: str) -> bool:
    return str(graph_type) in GENE_LEVEL_GRAPH_TYPES


def _missing_columns(df: pd.DataFrame, required: Iterable[str]) -> list[str]:
    return [column for column in required if column not in df.columns]


def _normalize_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().map(
        {
            "true": True,
            "1": True,
            "yes": True,
            "false": False,
            "0": False,
            "no": False,
        }
    )


def validate_graph_manifest_frame(
    df: pd.DataFrame,
    required_columns: Iterable[str] | None = None,
) -> None:
    required_columns = list(required_columns or REQUIRED_GRAPH_MANIFEST_COLUMNS)
    missing = _missing_columns(df, required_columns)
    if missing:
        raise ValueError(f"Graph manifest is missing required columns: {missing}")

    if df["graph_id"].astype(str).str.strip().eq("").any():
        raise ValueError("Graph manifest contains empty graph_id values")

    bad_species = sorted(set(df.loc[~df["species"].isin(ALLOWED_SPECIES), "species"].astype(str)))
    if bad_species:
        raise ValueError(f"Unsupported graph manifest species values detected: {bad_species}")

    bad_graph_types = sorted(set(df.loc[~df["graph_type"].isin(SUPPORTED_GRAPH_TYPES), "graph_type"].astype(str)))
    if bad_graph_types:
        raise ValueError(f"Unsupported graph_type values detected: {bad_graph_types}")

    bad_node_key_types = sorted(
        set(df.loc[~df["node_key_type"].isin(SUPPORTED_NODE_KEY_TYPES), "node_key_type"].astype(str))
    )
    if bad_node_key_types:
        raise ValueError(f"Unsupported node_key_type values detected: {bad_node_key_types}")

    duplicated_ids = df["graph_id"].astype(str).duplicated()
    if duplicated_ids.any():
        duplicated = sorted(df.loc[duplicated_ids, "graph_id"].astype(str).unique().tolist())
        raise ValueError(f"Duplicate graph_id values detected: {duplicated[:20]}")

    for bool_column in ["is_directed", "is_weighted"]:
        normalized = _normalize_bool_series(df[bool_column])
        if normalized.isna().any():
            bad_values = sorted(df.loc[normalized.isna(), bool_column].astype(str).unique().tolist())
            raise ValueError(f"Invalid boolean-like values detected in {bool_column}: {bad_values}")

    gene_level = df["graph_type"].astype(str).isin(GENE_LEVEL_GRAPH_TYPES)
    bad_gene_level = df.loc[gene_level & df["node_key_type"].astype(str).ne("canonical_gene_id"), "graph_id"]
    if not bad_gene_level.empty:
        raise ValueError(
            "All gene-level graphs must use node_key_type=canonical_gene_id. "
            f"Violations: {bad_gene_level.astype(str).tolist()[:20]}"
        )

    bad_residue = df.loc[df["graph_type"].astype(str).eq("residue_graph") & df["node_key_type"].astype(str).ne("residue_index"), "graph_id"]
    if not bad_residue.empty:
        raise ValueError(
            "All residue_graph rows must use node_key_type=residue_index. "
            f"Violations: {bad_residue.astype(str).tolist()[:20]}"
        )

    missing_metadata = df["metadata_path"].astype(str).str.strip().eq("")
    if missing_metadata.any():
        raise ValueError("Graph manifest contains empty metadata_path values")


def validate_graph_manifest_file(path: str | Path) -> pd.DataFrame:
    manifest_path = Path(path)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Graph manifest not found: {manifest_path}")
    df = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    validate_graph_manifest_frame(df)
    return df
