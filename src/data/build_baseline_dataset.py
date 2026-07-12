from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml
from sklearn.model_selection import train_test_split

from src.features.load_embeddings import attach_embedding_status, discover_embedding_manifest, load_embedding_index
from src.registry.load_registry import load_registry_bundle
from src.schemas.sample_schema import validate_sample_frame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build baseline dataset tables for ProGATE_v2")
    parser.add_argument("--config", type=str, required=True, help="Path to baseline YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def assign_splits(df: pd.DataFrame, seed: int, val_fraction: float, test_fraction: float) -> pd.DataFrame:
    output_parts = []
    for species, species_df in df.groupby("species", sort=True):
        working = species_df.copy()
        labels = working["gold_label"].astype(int)

        train_val_idx, test_idx = train_test_split(
            working.index,
            test_size=test_fraction,
            random_state=seed,
            stratify=labels,
        )

        train_idx, val_idx = train_test_split(
            train_val_idx,
            test_size=val_fraction / (1.0 - test_fraction),
            random_state=seed,
            stratify=labels.loc[train_val_idx],
        )

        working["split"] = ""
        working.loc[train_idx, "split"] = "train"
        working.loc[val_idx, "split"] = "val"
        working.loc[test_idx, "split"] = "test"
        output_parts.append(working)

    output = pd.concat(output_parts, ignore_index=True)
    return output


def build_support_supervised_samples(bundle, config: dict) -> pd.DataFrame:
    train_species = set(config["train"]["species"])
    labels = bundle.master_label_table.copy()

    samples = labels[
        (labels["species"].isin(train_species))
        & (labels["label_status"] == "gold")
        & (labels["gold_label"].isin(["0", "1"]))
        & (labels["needs_manual_review"].str.lower() == "false")
    ].copy()

    samples = samples[
        [
            "species",
            "canonical_gene_id",
            "gold_label",
            "label_status",
            "label_confidence",
            "notes",
        ]
    ].copy()
    samples["embedding_key"] = samples["canonical_gene_id"]
    samples["feature_path"] = ""
    samples["split"] = ""
    return samples


def build_fusarium_inference_pool(bundle) -> pd.DataFrame:
    labels = bundle.master_label_table.copy()
    fg = labels[labels["species"] == "fgraminearum"].copy()
    fg = fg.drop_duplicates(subset=["canonical_gene_id"], keep="first")
    fg = fg[
        [
            "species",
            "canonical_gene_id",
            "gold_label",
            "label_status",
            "label_confidence",
            "notes",
        ]
    ].copy()
    fg["split"] = "inference"
    fg["embedding_key"] = fg["canonical_gene_id"]
    fg["feature_path"] = ""
    return fg


def build_embedding_gap_report(samples: pd.DataFrame) -> pd.DataFrame:
    gap = samples.copy()
    gap["missing_embedding"] = ~gap["embedding_available"].astype(bool)
    return gap[
        [
            "species",
            "canonical_gene_id",
            "split",
            "gold_label",
            "label_status",
            "feature_path",
            "embedding_available",
            "missing_embedding",
        ]
    ]


def write_summary(output_dir: Path, supervised: pd.DataFrame, fusarium_pool: pd.DataFrame, fg_broad: pd.DataFrame, fg_strict: pd.DataFrame, fg_conflict: pd.DataFrame, embedding_manifest: Path | None) -> None:
    summary_lines = [
        "# Baseline Dataset Summary",
        "",
        "## Support Species Supervised Samples",
        f"- total rows: {len(supervised)}",
    ]

    for species, species_df in supervised.groupby("species", sort=True):
        summary_lines.append(f"- {species}: {len(species_df)}")

    summary_lines.extend(
        [
            "",
            "## Fusarium Pools",
            f"- inference pool rows: {len(fusarium_pool)}",
            f"- broad79 rows: {len(fg_broad)}",
            f"- strict29 rows: {len(fg_strict)}",
            f"- conflict8 rows: {len(fg_conflict)}",
            "",
            "## Embeddings",
            f"- embedding manifest: {embedding_manifest if embedding_manifest else 'NOT PROVIDED'}",
        ]
    )

    (output_dir / "dataset_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    registry_dir = Path(config["paths"]["registry_dir"])
    output_dir = Path(config["paths"]["baseline_dataset_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    bundle = load_registry_bundle(registry_dir)

    supervised = build_support_supervised_samples(bundle, config)
    supervised = assign_splits(
        supervised,
        seed=int(config["train"]["random_seed"]),
        val_fraction=float(config["train"]["val_fraction"]),
        test_fraction=float(config["train"]["test_fraction"]),
    )

    fusarium_pool = build_fusarium_inference_pool(bundle)
    fg_broad = bundle.fg_broad79.copy()
    fg_strict = bundle.fg_strict29.copy()
    fg_conflict = bundle.fg_conflict8.copy()

    embedding_manifest = discover_embedding_manifest(config["embeddings"].get("manifest_path"))
    embedding_index = load_embedding_index(embedding_manifest, config["embeddings"]["source_name"]) if embedding_manifest else None

    supervised = attach_embedding_status(supervised, embedding_index)
    fusarium_pool = attach_embedding_status(fusarium_pool, embedding_index)
    fg_broad = attach_embedding_status(fg_broad, embedding_index)
    fg_strict = attach_embedding_status(fg_strict, embedding_index)
    fg_conflict = attach_embedding_status(fg_conflict, embedding_index)

    validate_sample_frame(
        supervised[
            [
                "species",
                "canonical_gene_id",
                "gold_label",
                "label_status",
                "label_confidence",
                "split",
                "embedding_key",
                "notes",
            ]
        ]
    )

    supervised.to_csv(output_dir / "support_supervised_samples.tsv", sep="\t", index=False)
    fusarium_pool.to_csv(output_dir / "fgraminearum_inference_pool.tsv", sep="\t", index=False)
    fg_broad.to_csv(output_dir / "fgraminearum_broad79.tsv", sep="\t", index=False)
    fg_strict.to_csv(output_dir / "fgraminearum_strict29.tsv", sep="\t", index=False)
    fg_conflict.to_csv(output_dir / "fgraminearum_conflict8.tsv", sep="\t", index=False)

    combined_gap = pd.concat(
        [
            build_embedding_gap_report(supervised),
            build_embedding_gap_report(fusarium_pool),
        ],
        ignore_index=True,
    )
    combined_gap.to_csv(output_dir / "embedding_gap_report.tsv", sep="\t", index=False)

    write_summary(output_dir, supervised, fusarium_pool, fg_broad, fg_strict, fg_conflict, embedding_manifest)

    print(f"Wrote baseline dataset tables to: {output_dir}")
    print(f"Support supervised rows: {len(supervised)}")
    print(f"Fusarium inference pool rows: {len(fusarium_pool)}")
    print(f"Embedding manifest: {embedding_manifest if embedding_manifest else 'NOT PROVIDED'}")


if __name__ == "__main__":
    main()
