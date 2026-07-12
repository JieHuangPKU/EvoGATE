import argparse
from pathlib import Path

import pandas as pd
from sklearn.model_selection import train_test_split


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build label-scarcity split manifests with fixed test set.")
    parser.add_argument("--base-split-manifest", required=True, type=str)
    parser.add_argument("--output-dir", required=True, type=str)
    parser.add_argument("--train-fractions", required=True, type=str, help="Comma-separated list, e.g. 0.10,0.20,...,0.90")
    parser.add_argument("--seeds", required=True, type=str, help="Comma-separated list of integer seeds")
    return parser.parse_args()


def parse_fraction_csv(text: str) -> list[float]:
    values = [float(token.strip()) for token in str(text).split(",") if token.strip()]
    if not values:
        raise ValueError("No train fractions were provided")
    return values


def parse_seed_csv(text: str) -> list[int]:
    values = [int(token.strip()) for token in str(text).split(",") if token.strip()]
    if not values:
        raise ValueError("No seeds were provided")
    return values


def format_fraction_tag(fraction: float) -> str:
    return f"{int(round(float(fraction) * 100)):02d}"


def compute_val_relative_fraction(base_split: pd.DataFrame) -> float:
    test_fraction = float(pd.to_numeric(base_split["test_fraction"], errors="raise").iloc[0])
    val_fraction = float(pd.to_numeric(base_split["val_fraction"], errors="raise").iloc[0])
    if not 0.0 < test_fraction < 1.0:
        raise ValueError(f"Unexpected test_fraction={test_fraction}")
    val_relative = val_fraction / (1.0 - test_fraction)
    if not 0.0 < val_relative < 1.0:
        raise ValueError(f"Derived invalid within-train validation fraction={val_relative}")
    return float(val_relative)


def build_fraction_split(base_split: pd.DataFrame, fraction: float, seed: int, val_relative_fraction: float) -> pd.DataFrame:
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"train_fraction must be in (0, 1], got {fraction}")

    working = base_split.copy().reset_index(drop=True)
    test_mask = working["split"].astype(str).eq("test")
    pool_mask = ~test_mask
    pool = working.loc[pool_mask].copy()
    pool_labels = pool["label"].astype(int)

    selected_idx, _ = train_test_split(
        pool.index.to_numpy(),
        train_size=float(fraction),
        random_state=int(seed),
        stratify=pool_labels,
    )
    selected_idx = pd.Index(sorted(int(idx) for idx in selected_idx))
    selected = working.loc[selected_idx].copy()

    train_idx, val_idx = train_test_split(
        selected.index.to_numpy(),
        test_size=float(val_relative_fraction),
        random_state=int(seed),
        stratify=selected["label"].astype(int),
    )
    train_idx = pd.Index(sorted(int(idx) for idx in train_idx))
    val_idx = pd.Index(sorted(int(idx) for idx in val_idx))

    working["split"] = "unused"
    working.loc[test_mask, "split"] = "test"
    working.loc[train_idx, "split"] = "train"
    working.loc[val_idx, "split"] = "val"
    working["split_seed"] = int(seed)
    working["split_strategy"] = "label_scarcity_fixed_test_stratified_subsample"
    working["train_fraction"] = float(fraction)
    working["train_fraction_tag"] = format_fraction_tag(fraction)
    working["base_split_version"] = str(base_split["split_version"].iloc[0])
    working["split_version"] = (
        f"{working['base_split_version'].iloc[0]}_label_scarcity_train{format_fraction_tag(fraction)}_seed{int(seed)}"
    )
    working = working.loc[working["split"].isin(["train", "val", "test"])].copy()
    return working.sort_values(["split", "canonical_gene_id"], kind="stable").reset_index(drop=True)


def write_manifest(output_dir: Path, manifest: pd.DataFrame, fraction: float, seed: int) -> Path:
    fraction_tag = format_fraction_tag(fraction)
    output_path = output_dir / f"train_fraction_{fraction_tag}" / f"split_seed_{int(seed)}.tsv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_path, sep="\t", index=False)
    return output_path


def main() -> None:
    args = parse_args()
    base_split = pd.read_csv(args.base_split_manifest, sep="\t", dtype=str).fillna("")
    val_relative_fraction = compute_val_relative_fraction(base_split)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    for fraction in parse_fraction_csv(args.train_fractions):
        for seed in parse_seed_csv(args.seeds):
            manifest = build_fraction_split(base_split, fraction=fraction, seed=seed, val_relative_fraction=val_relative_fraction)
            output_path = write_manifest(output_dir, manifest, fraction=fraction, seed=seed)
            split_counts = manifest["split"].value_counts().to_dict()
            manifest_rows.append(
                {
                    "train_fraction": float(fraction),
                    "train_fraction_tag": format_fraction_tag(fraction),
                    "seed": int(seed),
                    "split_manifest_path": str(output_path),
                    "split_version": str(manifest["split_version"].iloc[0]),
                    "train_count": int(split_counts.get("train", 0)),
                    "val_count": int(split_counts.get("val", 0)),
                    "test_count": int(split_counts.get("test", 0)),
                    "train_positive_count": int(((manifest["split"] == "train") & (manifest["label"].astype(int) == 1)).sum()),
                    "val_positive_count": int(((manifest["split"] == "val") & (manifest["label"].astype(int) == 1)).sum()),
                    "test_positive_count": int(((manifest["split"] == "test") & (manifest["label"].astype(int) == 1)).sum()),
                    "base_split_version": str(manifest["base_split_version"].iloc[0]),
                }
            )

    pd.DataFrame(manifest_rows).sort_values(["train_fraction", "seed"]).to_csv(
        output_dir / "label_scarcity_split_manifest_index.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
