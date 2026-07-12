from __future__ import annotations

import argparse
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import torch
import yaml


POOLING_STRATEGY = "mean_pool_token_axis0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare pooled PLM embedding cache assets")
    parser.add_argument("--config", type=str, required=True, help="Path to baseline YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _species_from_canonical_id(canonical_gene_id: str) -> str:
    return canonical_gene_id.split("::", 1)[0] if "::" in canonical_gene_id else ""


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required TSV file not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _load_requested_samples(dataset_dir: Path) -> pd.DataFrame:
    frames = []
    for file_name in ["support_supervised_samples.tsv", "fgraminearum_inference_pool.tsv"]:
        path = dataset_dir / file_name
        if not path.exists():
            raise FileNotFoundError(f"Required baseline dataset file is missing: {path}")
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        frames.append(frame[["species", "canonical_gene_id"]])
    requested = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["canonical_gene_id"], keep="first")
    return requested.sort_values(["species", "canonical_gene_id"], kind="stable").reset_index(drop=True)


def _load_patch_lookup(config: dict) -> dict[str, dict[str, str]]:
    if not config["embeddings"].get("use_alignment_patch", False):
        return {}
    path = Path(config["embeddings"]["alignment_patch_path"])
    if not path.exists():
        return {}
    patch_df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    lookup = {}
    for record in patch_df.to_dict(orient="records"):
        if record["needs_manual_review"].lower() == "false" and record["patched_canonical_gene_id"].strip():
            lookup[record["original_canonical_gene_id"]] = record
    return lookup


def _pool_raw_pt(raw_path: Path) -> np.ndarray:
    data = torch.load(raw_path, map_location="cpu")
    if not isinstance(data, dict) or "feature_representation" not in data:
        raise ValueError(f"{raw_path} does not contain a supported feature_representation payload")
    feature = data["feature_representation"]
    if isinstance(feature, torch.Tensor):
        feature = feature.detach().cpu().numpy()
    feature = np.asarray(feature)
    if feature.ndim == 1:
        pooled = feature.astype(np.float32, copy=False)
    else:
        pooled = feature.mean(axis=0).astype(np.float32, copy=False)
    if pooled.ndim != 1:
        raise ValueError(f"{raw_path} pooled to non-vector shape: {pooled.shape}")
    return pooled


def _convert_one_raw_asset(raw_path: Path, output_dir: Path, species: str) -> dict[str, str]:
    gene_id = raw_path.stem
    canonical_gene_id = f"{species}::{gene_id}"
    feature_path = output_dir / f"{gene_id}.npy"
    status = "cached"
    notes = ""
    vector_dim = ""
    try:
        if feature_path.exists():
            vector = np.load(feature_path)
        else:
            vector = _pool_raw_pt(raw_path)
            np.save(feature_path, vector.astype(np.float32, copy=False))
            status = "converted"
        vector_dim = str(int(np.asarray(vector).shape[0]))
    except Exception as exc:
        status = "failed"
        notes = str(exc)

    return {
        "species": species,
        "canonical_gene_id": canonical_gene_id,
        "source_gene_id": gene_id,
        "embedding_source": "bingo_pooled_embedding",
        "raw_feature_path": str(raw_path),
        "feature_path": str(feature_path) if status != "failed" else "",
        "feature_format": "npy" if status != "failed" else "",
        "pooled": "true" if status != "failed" else "false",
        "pooling_strategy": POOLING_STRATEGY,
        "vector_dim": vector_dim,
        "exists": str(status != "failed").lower(),
        "conversion_status": status,
        "notes": notes,
    }


def build_pooled_assets(config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    dataset_dir = Path(config["paths"]["baseline_dataset_dir"])
    pooled_root = Path(config["paths"]["pooled_embeddings_dir"])
    pooled_root.mkdir(parents=True, exist_ok=True)

    requested = _load_requested_samples(dataset_dir)
    patch_lookup = _load_patch_lookup(config)

    raw_dirs = config["embeddings"]["raw_candidate_dirs"]
    manifest_rows = []
    asset_lookup: dict[str, str] = {}
    max_workers = int(config["embeddings"].get("pooling_workers", max(4, min(16, (os.cpu_count() or 8)))))

    for raw_dir_value in raw_dirs:
        raw_dir = Path(raw_dir_value)
        if not raw_dir.exists():
            continue
        species = raw_dir.name if raw_dir.name in {"scerevisiae", "human", "celegans", "fgraminearum"} else raw_dir.parent.name
        if species not in {"scerevisiae", "human", "celegans", "fgraminearum"}:
            continue
        output_dir = pooled_root / species
        output_dir.mkdir(parents=True, exist_ok=True)
        raw_files = sorted(raw_dir.glob("*.pt"))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(_convert_one_raw_asset, raw_path, output_dir, species): raw_path
                for raw_path in raw_files
            }
            for index, future in enumerate(as_completed(future_map), start=1):
                row = future.result()
                if row["conversion_status"] != "failed":
                    asset_lookup[row["canonical_gene_id"]] = row["feature_path"]
                manifest_rows.append(row)
                if index % 2000 == 0:
                    print(f"[plm_prepare] {species}: processed {index}/{len(raw_files)} raw assets", flush=True)

    manifest_df = pd.DataFrame(manifest_rows).sort_values(["species", "canonical_gene_id"], kind="stable").reset_index(drop=True)

    missing_rows = []
    for record in requested.to_dict(orient="records"):
        original_canonical_gene_id = record["canonical_gene_id"]
        species = record["species"] or _species_from_canonical_id(original_canonical_gene_id)
        patch_record = patch_lookup.get(original_canonical_gene_id)
        effective_canonical_gene_id = (
            patch_record["patched_canonical_gene_id"] if patch_record else original_canonical_gene_id
        )
        feature_path = asset_lookup.get(effective_canonical_gene_id, "")
        if feature_path:
            continue
        missing_rows.append(
            {
                "species": species,
                "original_canonical_gene_id": original_canonical_gene_id,
                "effective_canonical_gene_id": effective_canonical_gene_id,
                "alignment_patch_applied": str(bool(patch_record)).lower(),
                "patch_rule": patch_record["patch_rule"] if patch_record else "",
                "feature_path": "",
                "missing_reason": "no pooled asset available for requested canonical id after optional alignment patch",
            }
        )

    failed_rows = manifest_df[manifest_df["conversion_status"] == "failed"].copy()
    if not failed_rows.empty:
        failed_rows = failed_rows.assign(
            original_canonical_gene_id=failed_rows["canonical_gene_id"],
            effective_canonical_gene_id=failed_rows["canonical_gene_id"],
            alignment_patch_applied="false",
            patch_rule="",
            missing_reason=failed_rows["notes"],
        )[
            [
                "species",
                "original_canonical_gene_id",
                "effective_canonical_gene_id",
                "alignment_patch_applied",
                "patch_rule",
                "feature_path",
                "missing_reason",
            ]
        ]
        missing_rows.extend(failed_rows.to_dict(orient="records"))

    missing_df = pd.DataFrame(missing_rows)
    if missing_df.empty:
        missing_df = pd.DataFrame(
            columns=[
                "species",
                "original_canonical_gene_id",
                "effective_canonical_gene_id",
                "alignment_patch_applied",
                "patch_rule",
                "feature_path",
                "missing_reason",
            ]
        )

    return manifest_df, missing_df.sort_values(["species", "original_canonical_gene_id"], kind="stable").reset_index(drop=True)


def write_summary(summary_path: Path, manifest_df: pd.DataFrame, missing_df: pd.DataFrame) -> None:
    converted = manifest_df[manifest_df["conversion_status"].isin(["converted", "cached"])].copy()
    failed = manifest_df[manifest_df["conversion_status"] == "failed"].copy()

    lines = [
        "# 33 PLM Pooling Spec",
        "",
        "## Pooling Strategy",
        f"- strategy: {POOLING_STRATEGY}",
        "- implementation: load Bingo raw `.pt` payload, read `feature_representation`, apply mean pooling over token axis 0",
        "- output dtype: float32",
        "",
        "## Why This Strategy",
        "- Bingo upstream `runners/esm/data_loader.py` already uses `.mean(0)` on token-level representations",
        "- mean pooling is deterministic, reproducible, and produces one fixed-length vector per gene without sample-order assumptions",
        "",
        "## Input / Output Format",
        "- input: gene-keyed Bingo raw `.pt` assets under config-driven raw candidate directories",
        "- output asset: `outputs/pooled_embeddings/<species>/<gene_id>.npy`",
        "- output manifest: `outputs/pooled_embeddings/pooled_embedding_manifest.tsv`",
        "- output missing report: `outputs/pooled_embeddings/pooled_embedding_missing.tsv`",
        "",
        "## Conversion Summary",
        f"- successful pooled assets: {len(converted)}",
        f"- failed raw asset conversions: {len(failed)}",
        f"- requested samples still missing after optional patch: {len(missing_df)}",
        "",
        "## Unconverted / Missing Reasons",
    ]

    if missing_df.empty:
        lines.append("- none")
    else:
        for _, row in missing_df.head(30).iterrows():
            lines.append(f"- {row['original_canonical_gene_id']}: {row['missing_reason']}")
        if len(missing_df) > 30:
            lines.append(f"- ... {len(missing_df) - 30} more rows omitted from summary; see pooled_embedding_missing.tsv")

    summary_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    pooled_dir = Path(config["paths"]["pooled_embeddings_dir"])
    pooled_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = pooled_dir / "pooled_embedding_manifest.tsv"
    missing_path = pooled_dir / "pooled_embedding_missing.tsv"

    manifest_df, missing_df = build_pooled_assets(config)
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    missing_df.to_csv(missing_path, sep="\t", index=False)

    summary_path = Path("33_plm_pooling_spec.md")
    write_summary(summary_path, manifest_df, missing_df)

    converted = int(manifest_df["conversion_status"].isin(["converted", "cached"]).sum())
    failed = int((manifest_df["conversion_status"] == "failed").sum())
    print(f"Wrote pooled embedding manifest to: {manifest_path}")
    print(f"Wrote pooled embedding missing report to: {missing_path}")
    print(f"Successful pooled assets: {converted}")
    print(f"Failed raw asset conversions: {failed}")
    print(f"Wrote pooling summary to: {summary_path}")


if __name__ == "__main__":
    main()
