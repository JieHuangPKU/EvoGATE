from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from src.features.load_embeddings import validate_manifest_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build embedding manifest aligned by canonical_gene_id")
    parser.add_argument("--config", type=str, required=True, help="Path to baseline YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _species_from_canonical_id(canonical_gene_id: str) -> str:
    if "::" not in canonical_gene_id:
        return ""
    return canonical_gene_id.split("::", 1)[0]


def _collect_requested_samples(dataset_dir: Path) -> pd.DataFrame:
    frames = []
    for file_name in ["support_supervised_samples.tsv", "fgraminearum_inference_pool.tsv"]:
        path = dataset_dir / file_name
        if not path.exists():
            raise FileNotFoundError(f"Baseline dataset file is missing: {path}. Run build_baseline_dataset.py first.")
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        required = ["species", "canonical_gene_id"]
        optional = [column for column in ["raw_gene_id", "sample_id"] if column in frame.columns]
        frames.append(frame[required + optional])

    requested = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["canonical_gene_id"], keep="first")
    requested = requested.sort_values(["species", "canonical_gene_id"], kind="stable").reset_index(drop=True)
    if "raw_gene_id" not in requested.columns:
        requested["raw_gene_id"] = ""
    requested["canonical_gene_suffix"] = requested["canonical_gene_id"].astype(str).map(_canonical_gene_suffix)
    return requested


def _canonical_gene_suffix(canonical_gene_id: str) -> str:
    if "::" not in canonical_gene_id:
        return canonical_gene_id
    return canonical_gene_id.split("::", 1)[1]


def _load_candidate_manifest_tables(candidate_paths: list[str]) -> pd.DataFrame:
    frames = []
    for raw_path in candidate_paths:
        if not raw_path:
            continue
        path = Path(raw_path)
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        required = {"canonical_gene_id", "feature_path"}
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"Candidate manifest {path} is missing required columns: {missing}")
        if "embedding_source" not in frame.columns:
            frame["embedding_source"] = "candidate_manifest"
        if "feature_format" not in frame.columns:
            frame["feature_format"] = frame["feature_path"].map(lambda value: Path(value).suffix.lstrip(".").lower())
        if "pooled" not in frame.columns:
            frame["pooled"] = "false"
        if "effective_canonical_gene_id" not in frame.columns:
            frame["effective_canonical_gene_id"] = frame["canonical_gene_id"]
        if "alignment_patch_applied" not in frame.columns:
            frame["alignment_patch_applied"] = "false"
        if "patch_rule" not in frame.columns:
            frame["patch_rule"] = ""
        if "notes" not in frame.columns:
            frame["notes"] = ""
        frame["notes"] = frame["notes"].astype(str)
        frames.append(
            frame[
                [
                    "canonical_gene_id",
                    "effective_canonical_gene_id",
                    "feature_path",
                    "embedding_source",
                    "feature_format",
                    "pooled",
                    "alignment_patch_applied",
                    "patch_rule",
                    "notes",
                ]
            ]
        )

    if not frames:
        return pd.DataFrame(
            columns=[
                "canonical_gene_id",
                "effective_canonical_gene_id",
                "feature_path",
                "embedding_source",
                "feature_format",
                "pooled",
                "alignment_patch_applied",
                "patch_rule",
                "notes",
            ]
        )

    merged = pd.concat(frames, ignore_index=True)
    return merged.drop_duplicates(subset=["canonical_gene_id", "feature_path"], keep="first")


def _infer_species_from_path(path: Path) -> str:
    known_species = {"scerevisiae", "human", "celegans", "fgraminearum"}
    for part in path.parts:
        if part in known_species:
            return part
    return ""


def _build_scan_index(candidate_dirs: list[str], allowed_suffixes: set[str]) -> tuple[dict[str, list[Path]], dict[tuple[str, str], list[Path]]]:
    exact_index: dict[str, list[Path]] = {}
    species_gene_index: dict[tuple[str, str], list[Path]] = {}
    for raw_dir in candidate_dirs:
        if not raw_dir:
            continue
        directory = Path(raw_dir)
        if not directory.exists():
            continue
        for path in directory.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() not in allowed_suffixes:
                continue
            stem = path.stem
            for key in {stem, stem.replace("__", "::")}:
                if "::" in key:
                    exact_index.setdefault(key, []).append(path)
            species = _infer_species_from_path(path)
            if species:
                species_gene_index.setdefault((species, stem), []).append(path)
    return exact_index, species_gene_index


def _load_alignment_patch(config: dict) -> dict[str, dict[str, str]]:
    embeddings_cfg = config.get("embeddings", {})
    if not embeddings_cfg.get("use_alignment_patch", False):
        return {}
    path = Path(embeddings_cfg.get("alignment_patch_path", ""))
    if not path.exists():
        return {}
    patch_df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    lookup = {}
    for record in patch_df.to_dict(orient="records"):
        if record["needs_manual_review"].lower() == "false" and record["patched_canonical_gene_id"].strip():
            lookup[record["original_canonical_gene_id"]] = record
    return lookup


def build_manifest_rows(requested: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    embeddings_cfg = config.get("embeddings", {})
    candidate_manifest_paths = embeddings_cfg.get("candidate_manifest_paths", []) or []
    candidate_dirs = embeddings_cfg.get("candidate_dirs", []) or []
    allowed_suffixes = {suffix.lower() for suffix in embeddings_cfg.get("allowed_suffixes", [".npy", ".npz"])}
    manifest_candidates = _load_candidate_manifest_tables(candidate_manifest_paths)
    exact_scan_index, species_gene_scan_index = _build_scan_index(candidate_dirs, allowed_suffixes)
    alignment_patch_lookup = _load_alignment_patch(config)

    candidate_lookup: dict[str, list[dict[str, str]]] = {}
    for record in manifest_candidates.to_dict(orient="records"):
        candidate_lookup.setdefault(record["canonical_gene_id"], []).append(record)
    for canonical_gene_id, paths in exact_scan_index.items():
        for path in paths:
            candidate_lookup.setdefault(canonical_gene_id, []).append(
                {
                    "canonical_gene_id": canonical_gene_id,
                    "effective_canonical_gene_id": canonical_gene_id,
                    "feature_path": str(path),
                    "embedding_source": embeddings_cfg.get("source_name", "embedding_source"),
                    "feature_format": path.suffix.lstrip(".").lower(),
                    "pooled": "true" if path.suffix.lower() in {".npy", ".npz"} else "false",
                    "alignment_patch_applied": "false",
                    "patch_rule": "",
                    "notes": "matched by exact canonical_gene_id filename stem",
                }
            )

    manifest_rows = []
    missing_rows = []
    for record in requested.to_dict(orient="records"):
        canonical_gene_id = record["canonical_gene_id"]
        species = record["species"] or _species_from_canonical_id(canonical_gene_id)
        raw_gene_id = record.get("raw_gene_id", "")
        canonical_gene_suffix = record.get("canonical_gene_suffix", _canonical_gene_suffix(canonical_gene_id))
        patch_record = alignment_patch_lookup.get(canonical_gene_id)
        effective_canonical_gene_id = patch_record["patched_canonical_gene_id"] if patch_record else canonical_gene_id
        effective_suffix = _canonical_gene_suffix(effective_canonical_gene_id)
        alignment_patch_applied = str(bool(patch_record)).lower()
        patch_rule = patch_record["patch_rule"] if patch_record else ""

        candidates = list(candidate_lookup.get(effective_canonical_gene_id, []))
        for gene_key, note in [
            (raw_gene_id, "matched by species + raw_gene_id filename stem"),
            (effective_suffix, "matched by species + effective canonical gene suffix filename stem"),
            (canonical_gene_suffix, "matched by species + original canonical gene suffix filename stem"),
        ]:
            if not gene_key:
                continue
            for path in species_gene_scan_index.get((species, gene_key), []):
                candidates.append(
                    {
                        "canonical_gene_id": canonical_gene_id,
                        "effective_canonical_gene_id": effective_canonical_gene_id,
                        "feature_path": str(path),
                        "embedding_source": embeddings_cfg.get("source_name", "embedding_source"),
                        "feature_format": path.suffix.lstrip(".").lower(),
                        "pooled": "true" if path.suffix.lower() in {".npy", ".npz"} else "false",
                        "alignment_patch_applied": alignment_patch_applied,
                        "patch_rule": patch_rule,
                        "notes": note,
                    }
                )

        deduped_candidates = []
        seen_paths = set()
        for candidate in candidates:
            feature_path = candidate["feature_path"]
            if feature_path in seen_paths:
                continue
            seen_paths.add(feature_path)
            deduped_candidates.append(candidate)
        candidates = deduped_candidates

        if len(candidates) == 1:
            candidate = candidates[0]
            feature_path = candidate["feature_path"]
            exists = Path(feature_path).exists()
            row = {
                "species": species,
                "canonical_gene_id": canonical_gene_id,
                "effective_canonical_gene_id": candidate.get("effective_canonical_gene_id", effective_canonical_gene_id),
                "embedding_source": candidate["embedding_source"],
                "feature_path": feature_path,
                "feature_format": candidate["feature_format"],
                "pooled": str(candidate.get("pooled", "false")).lower(),
                "alignment_patch_applied": alignment_patch_applied,
                "patch_rule": patch_rule,
                "exists": str(exists).lower(),
                "needs_manual_review": "false",
                "notes": candidate.get("notes", ""),
            }
            manifest_rows.append(row)
            if not exists:
                missing_rows.append({**row, "notes": "feature_path listed but file does not exist"})
            continue

        if len(candidates) > 1:
            row = {
                "species": species,
                "canonical_gene_id": canonical_gene_id,
                "effective_canonical_gene_id": effective_canonical_gene_id,
                "embedding_source": embeddings_cfg.get("source_name", "embedding_source"),
                "feature_path": "",
                "feature_format": "",
                "pooled": "false",
                "alignment_patch_applied": alignment_patch_applied,
                "patch_rule": patch_rule,
                "exists": "false",
                "needs_manual_review": "true",
                "notes": f"multiple candidate embeddings found: {len(candidates)}",
            }
            manifest_rows.append(row)
            missing_rows.append(row.copy())
            continue

        row = {
            "species": species,
            "canonical_gene_id": canonical_gene_id,
            "effective_canonical_gene_id": effective_canonical_gene_id,
            "embedding_source": embeddings_cfg.get("source_name", "embedding_source"),
            "feature_path": "",
            "feature_format": "",
            "pooled": "false",
            "alignment_patch_applied": alignment_patch_applied,
            "patch_rule": patch_rule,
            "exists": "false",
            "needs_manual_review": "false",
            "notes": "no aligned embedding asset found by canonical_gene_id, optional patched canonical_gene_id, raw_gene_id, or canonical suffix",
        }
        manifest_rows.append(row)
        missing_rows.append(row.copy())

    manifest_df = pd.DataFrame(manifest_rows)
    missing_df = pd.DataFrame(missing_rows)
    validate_manifest_rows(manifest_df)
    return manifest_df, missing_df


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    dataset_dir = Path(config["paths"]["baseline_dataset_dir"])
    output_manifest = Path(config["embeddings"]["manifest_path"])
    missing_report = Path(config["embeddings"]["missing_report_path"])
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    missing_report.parent.mkdir(parents=True, exist_ok=True)

    requested = _collect_requested_samples(dataset_dir)
    manifest_df, missing_df = build_manifest_rows(requested, config)

    manifest_df.to_csv(output_manifest, sep="\t", index=False)
    missing_df.to_csv(missing_report, sep="\t", index=False)

    resolved = int((manifest_df["feature_path"].astype(str).str.strip().ne("")).sum())
    missing = int((manifest_df["feature_path"].astype(str).str.strip().eq("")).sum())
    review = int((manifest_df["needs_manual_review"].astype(str).str.lower() == "true").sum())

    print(f"Wrote embedding manifest to: {output_manifest}")
    print(f"Wrote embedding missing report to: {missing_report}")
    print(f"Requested canonical_gene_id rows: {len(requested)}")
    print(f"Resolved embeddings: {resolved}")
    print(f"Missing embeddings: {missing}")
    print(f"Needs manual review: {review}")


if __name__ == "__main__":
    main()
