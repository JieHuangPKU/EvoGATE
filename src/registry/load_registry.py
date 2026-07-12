from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from src.schemas.sample_schema import ALLOWED_GOLD_LABELS, ALLOWED_LABEL_STATUS


MASTER_LABEL_REQUIRED = {
    "species",
    "canonical_gene_id",
    "gold_label",
    "label_status",
    "label_confidence",
    "needs_manual_review",
}

FG_BENCHMARK_REQUIRED = {
    "species",
    "canonical_gene_id",
    "gold_label",
    "label_status",
    "mapping_confidence",
}


@dataclass
class RegistryBundle:
    registry_dir: Path
    master_label_table: pd.DataFrame
    master_evidence_table: pd.DataFrame
    fg_broad79: pd.DataFrame
    fg_strict29: pd.DataFrame
    fg_conflict8: pd.DataFrame


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_registry_dir() -> Path:
    return repo_root() / "data_registry"


def _require_columns(df: pd.DataFrame, required: Iterable[str], table_name: str) -> None:
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def _read_tsv(path: Path, table_name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required registry file not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _validate_label_table(df: pd.DataFrame, table_name: str) -> None:
    _require_columns(df, MASTER_LABEL_REQUIRED, table_name)

    if df["canonical_gene_id"].astype(str).str.strip().eq("").any():
        raise ValueError(f"{table_name} contains empty canonical_gene_id values")

    duplicates = df["canonical_gene_id"].astype(str).duplicated()
    if duplicates.any():
        duplicated_ids = sorted(df.loc[duplicates, "canonical_gene_id"].astype(str).unique().tolist())
        raise ValueError(f"{table_name} contains duplicate canonical_gene_id values: {duplicated_ids[:20]}")

    invalid_labels = sorted(set(df.loc[~df["gold_label"].astype(str).isin(ALLOWED_GOLD_LABELS), "gold_label"]))
    if invalid_labels:
        raise ValueError(f"{table_name} contains invalid gold_label values: {invalid_labels}")

    invalid_status = sorted(set(df.loc[~df["label_status"].astype(str).isin(ALLOWED_LABEL_STATUS), "label_status"]))
    if invalid_status:
        raise ValueError(f"{table_name} contains invalid label_status values: {invalid_status}")


def _validate_fg_table(df: pd.DataFrame, table_name: str) -> None:
    _require_columns(df, FG_BENCHMARK_REQUIRED, table_name)

    if not (df["species"] == "fgraminearum").all():
        raise ValueError(f"{table_name} contains non-fgraminearum rows")

    duplicates = df["canonical_gene_id"].astype(str).duplicated()
    if duplicates.any():
        duplicated_ids = sorted(df.loc[duplicates, "canonical_gene_id"].astype(str).unique().tolist())
        raise ValueError(f"{table_name} contains duplicate canonical_gene_id values: {duplicated_ids[:20]}")


def load_registry_bundle(registry_dir: str | Path | None = None) -> RegistryBundle:
    registry_dir = Path(registry_dir) if registry_dir else default_registry_dir()

    master_label_table = _read_tsv(registry_dir / "master_label_table.preliminary.tsv", "master_label_table")
    master_evidence_table = _read_tsv(registry_dir / "master_evidence_table.preliminary.tsv", "master_evidence_table")
    fg_broad79 = _read_tsv(registry_dir / "fgraminearum_gold_positive.broad79.tsv", "fgraminearum_gold_positive.broad79")
    fg_strict29 = _read_tsv(registry_dir / "fgraminearum_gold_positive.strict29.tsv", "fgraminearum_gold_positive.strict29")
    fg_conflict8 = _read_tsv(registry_dir / "fgraminearum_gold_positive.conflict.tsv", "fgraminearum_gold_positive.conflict")

    _validate_label_table(master_label_table, "master_label_table")
    _validate_fg_table(fg_broad79, "fgraminearum_gold_positive.broad79")
    _validate_fg_table(fg_strict29, "fgraminearum_gold_positive.strict29")
    _validate_fg_table(fg_conflict8, "fgraminearum_gold_positive.conflict")

    return RegistryBundle(
        registry_dir=registry_dir,
        master_label_table=master_label_table,
        master_evidence_table=master_evidence_table,
        fg_broad79=fg_broad79,
        fg_strict29=fg_strict29,
        fg_conflict8=fg_conflict8,
    )
