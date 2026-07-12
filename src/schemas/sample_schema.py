from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Iterable

import pandas as pd


ALLOWED_SPECIES = {"scerevisiae", "human", "celegans", "fgraminearum"}
ALLOWED_LABEL_STATUS = {"gold", "weak_support", "excluded", "unresolved", "inference"}
ALLOWED_GOLD_LABELS = {"0", "1", ""}

REQUIRED_SAMPLE_COLUMNS = [
    "species",
    "canonical_gene_id",
    "gold_label",
    "label_status",
    "label_confidence",
    "split",
    "embedding_key",
    "notes",
]


@dataclass
class SampleRecord:
    species: str
    canonical_gene_id: str
    gold_label: str
    label_status: str
    label_confidence: str
    split: str
    embedding_key: str
    notes: str = ""
    feature_path: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def _missing_columns(df: pd.DataFrame, required: Iterable[str]) -> list[str]:
    return [column for column in required if column not in df.columns]


def validate_sample_frame(df: pd.DataFrame, required_columns: Iterable[str] | None = None) -> None:
    required_columns = list(required_columns or REQUIRED_SAMPLE_COLUMNS)
    missing = _missing_columns(df, required_columns)
    if missing:
        raise ValueError(f"Sample frame is missing required columns: {missing}")

    if df["canonical_gene_id"].isna().any() or (df["canonical_gene_id"].astype(str).str.strip() == "").any():
        raise ValueError("Sample frame contains empty canonical_gene_id values")

    bad_species = sorted(set(df.loc[~df["species"].isin(ALLOWED_SPECIES), "species"].astype(str)))
    if bad_species:
        raise ValueError(f"Unsupported species values detected: {bad_species}")

    if not df["canonical_gene_id"].astype(str).str.contains("::", regex=False).all():
        raise ValueError("All canonical_gene_id values must use 'species::gene_id' format")

    bad_label_status = sorted(
        set(df.loc[~df["label_status"].isin(ALLOWED_LABEL_STATUS), "label_status"].astype(str))
    )
    if bad_label_status:
        raise ValueError(f"Invalid label_status values detected: {bad_label_status}")

    bad_labels = sorted(set(df.loc[~df["gold_label"].astype(str).isin(ALLOWED_GOLD_LABELS), "gold_label"].astype(str)))
    if bad_labels:
        raise ValueError(f"Invalid gold_label values detected: {bad_labels}")

    duplicates = df["canonical_gene_id"].astype(str).duplicated()
    if duplicates.any():
        duplicated_ids = sorted(df.loc[duplicates, "canonical_gene_id"].astype(str).unique().tolist())
        raise ValueError(f"Duplicate canonical_gene_id values detected: {duplicated_ids[:20]}")
