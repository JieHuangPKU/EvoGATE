from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch


@dataclass
class EmbeddingIndex:
    source_name: str
    manifest_path: Path
    table: pd.DataFrame


def discover_embedding_manifest(path_value: str | Path | None) -> Path | None:
    if not path_value:
        return None
    path = Path(path_value)
    return path if path.exists() else None


def load_embedding_index(path_value: str | Path | None, source_name: str = "embedding_source") -> EmbeddingIndex:
    manifest_path = discover_embedding_manifest(path_value)
    if manifest_path is None:
        raise FileNotFoundError(
            "Embedding manifest is not available. Set configs/baseline.yaml -> embeddings.manifest_path "
            "to a TSV with columns canonical_gene_id and feature_path."
        )

    table = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    required = {"canonical_gene_id", "feature_path"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Embedding manifest is missing required columns: {missing}")

    if "feature_format" not in table.columns:
        table["feature_format"] = table["feature_path"].map(_infer_feature_format)
    if "exists" not in table.columns:
        table["exists"] = table["feature_path"].map(lambda value: str(Path(value).exists()).lower() if value else "false")
    if "needs_manual_review" not in table.columns:
        table["needs_manual_review"] = "false"
    if "species" not in table.columns:
        table["species"] = table["canonical_gene_id"].astype(str).map(_species_from_canonical_id)
    if "pooled" not in table.columns:
        table["pooled"] = table["feature_format"].astype(str).isin(["npy", "npz"]).map(lambda value: str(value).lower())
    if "effective_canonical_gene_id" not in table.columns:
        table["effective_canonical_gene_id"] = table["canonical_gene_id"]
    if "alignment_patch_applied" not in table.columns:
        table["alignment_patch_applied"] = "false"
    if "patch_rule" not in table.columns:
        table["patch_rule"] = ""

    duplicates = table["canonical_gene_id"].astype(str).duplicated()
    if duplicates.any():
        duplicated_ids = sorted(table.loc[duplicates, "canonical_gene_id"].astype(str).unique().tolist())
        raise ValueError(f"Embedding manifest contains duplicate canonical_gene_id values: {duplicated_ids[:20]}")

    return EmbeddingIndex(source_name=source_name, manifest_path=manifest_path, table=table)


def attach_embedding_status(samples: pd.DataFrame, embedding_index: EmbeddingIndex | None) -> pd.DataFrame:
    output = samples.copy()
    if embedding_index is None:
        output["feature_path"] = ""
        output["embedding_available"] = False
        output["embedding_source"] = ""
        output["effective_canonical_gene_id"] = output["canonical_gene_id"]
        output["alignment_patch_applied"] = "false"
        return output

    merged = output.merge(
        embedding_index.table[
            ["canonical_gene_id", "feature_path", "effective_canonical_gene_id", "alignment_patch_applied"]
        ],
        on="canonical_gene_id",
        how="left",
        suffixes=("", "_manifest"),
    )
    if "feature_path_manifest" in merged.columns:
        merged["feature_path"] = merged["feature_path_manifest"]
        merged = merged.drop(columns=["feature_path_manifest"])
    merged["embedding_available"] = merged["feature_path"].astype(str).str.strip().ne("")
    merged["embedding_source"] = embedding_index.source_name
    merged["effective_canonical_gene_id"] = merged["effective_canonical_gene_id"].where(
        merged["effective_canonical_gene_id"].astype(str).str.strip().ne(""),
        merged["canonical_gene_id"],
    )
    merged["alignment_patch_applied"] = merged["alignment_patch_applied"].where(
        merged["alignment_patch_applied"].astype(str).str.strip().ne(""),
        "false",
    )
    return merged


def load_vector(feature_path: str | Path) -> np.ndarray:
    feature_path = Path(feature_path)
    if not feature_path.exists():
        raise FileNotFoundError(f"Feature file not found: {feature_path}")

    if feature_path.suffix == ".npy":
        return np.load(feature_path)

    if feature_path.suffix == ".npz":
        data = np.load(feature_path)
        if "embedding" in data:
            return data["embedding"]
        first_key = list(data.keys())[0]
        return data[first_key]

    if feature_path.suffix == ".pt":
        data = torch.load(feature_path, map_location="cpu")
        if isinstance(data, dict):
            if "feature_representation" in data:
                feature = data["feature_representation"]
                if isinstance(feature, torch.Tensor):
                    feature = feature.detach().cpu().numpy()
                feature = np.asarray(feature)
                if feature.ndim == 1:
                    return feature.astype(np.float32, copy=False)
                return feature.mean(axis=0).astype(np.float32, copy=False)
            if "embedding" in data:
                feature = data["embedding"]
                if isinstance(feature, torch.Tensor):
                    feature = feature.detach().cpu().numpy()
                return np.asarray(feature, dtype=np.float32)
        if isinstance(data, torch.Tensor):
            return data.detach().cpu().numpy().astype(np.float32, copy=False)
        raise ValueError(f"Unsupported .pt embedding payload structure in {feature_path}")

    raise ValueError(f"Unsupported feature file suffix: {feature_path.suffix}")


def _infer_feature_format(feature_path: str) -> str:
    if not feature_path:
        return ""
    suffix = Path(feature_path).suffix.lower()
    if suffix in {".npy", ".npz", ".pt"}:
        return suffix.lstrip(".")
    return "unsupported"


def _species_from_canonical_id(canonical_gene_id: str) -> str:
    if "::" not in canonical_gene_id:
        return ""
    return canonical_gene_id.split("::", 1)[0]


def validate_manifest_rows(table: pd.DataFrame) -> None:
    duplicated_existing = table[table["feature_path"].astype(str).str.strip().ne("")]["canonical_gene_id"].duplicated()
    if duplicated_existing.any():
        duplicated_ids = (
            table.loc[duplicated_existing, "canonical_gene_id"].astype(str).drop_duplicates().tolist()
        )
        raise ValueError(f"Embedding manifest contains duplicated mapped canonical_gene_id values: {duplicated_ids[:20]}")


def load_feature_matrix(
    samples: pd.DataFrame,
    embedding_index: EmbeddingIndex,
    require_all: bool = True,
    require_pooled_features: bool = False,
) -> tuple[np.ndarray, pd.DataFrame]:
    merged = samples.merge(
        embedding_index.table[
            [
                "canonical_gene_id",
                "effective_canonical_gene_id",
                "feature_path",
                "feature_format",
                "pooled",
                "exists",
                "needs_manual_review",
                "alignment_patch_applied",
                "patch_rule",
                "species",
            ]
        ],
        on="canonical_gene_id",
        how="left",
        suffixes=("", "_manifest"),
    )
    for column in [
        "effective_canonical_gene_id",
        "feature_path",
        "feature_format",
        "pooled",
        "exists",
        "needs_manual_review",
        "alignment_patch_applied",
        "patch_rule",
        "species",
    ]:
        manifest_column = f"{column}_manifest"
        if manifest_column in merged.columns:
            merged[column] = merged[manifest_column]
            merged = merged.drop(columns=[manifest_column])
    merged["feature_path"] = merged["feature_path"].astype(str)
    merged["exists"] = merged["exists"].astype(str).str.lower()
    merged["needs_manual_review"] = merged["needs_manual_review"].astype(str).str.lower()
    merged["pooled"] = merged["pooled"].astype(str).str.lower()
    merged["alignment_patch_applied"] = merged["alignment_patch_applied"].astype(str).str.lower()
    merged["effective_canonical_gene_id"] = merged["effective_canonical_gene_id"].where(
        merged["effective_canonical_gene_id"].astype(str).str.strip().ne(""),
        merged["canonical_gene_id"],
    )

    valid_mask = (
        merged["feature_path"].str.strip().ne("")
        & merged["exists"].isin({"true", "1", "yes"})
        & ~merged["needs_manual_review"].isin({"true", "1", "yes"})
    )

    if require_pooled_features:
        non_pooled_mask = valid_mask & ~merged["pooled"].isin({"true", "1", "yes"})
        if non_pooled_mask.any():
            example_rows = merged.loc[non_pooled_mask, ["canonical_gene_id", "feature_path", "feature_format"]].head(20)
            raise ValueError(
                "Pooled features are required, but non-pooled manifest rows were selected. "
                f"Examples: {example_rows.to_dict(orient='records')}"
            )

    if require_all and not valid_mask.all():
        missing_ids = sorted(merged.loc[~valid_mask, "canonical_gene_id"].astype(str).tolist())
        raise FileNotFoundError(
            "Embeddings are missing or not ready for some required samples. "
            f"First missing canonical_gene_id values: {missing_ids[:20]}"
        )

    aligned = merged.loc[valid_mask].copy()
    if aligned.empty:
        raise FileNotFoundError("No usable embeddings were found after applying the manifest and validation filters.")

    vectors = [np.ravel(load_vector(path)) for path in aligned["feature_path"]]
    dimensions = sorted({int(vector.shape[0]) for vector in vectors})
    if len(dimensions) != 1:
        raise ValueError(f"Embedding dimensionality mismatch detected: {dimensions}")

    x = np.vstack(vectors)
    return x, aligned.reset_index(drop=True)
