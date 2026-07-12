from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.eval.publication_summary import build_publication_summary, publication_markdown


METRICS = ["auroc", "auprc", "mcc", "f1", "precision", "recall", "accuracy", "specificity"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate one Phase 1 benchmark target across fixed-split seeds")
    parser.add_argument("--run-root", required=True, help="Directory containing run_<seed>/metrics.tsv files")
    parser.add_argument("--target", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--regime", required=True)
    parser.add_argument("--feature-setting", required=True)
    parser.add_argument("--model-name", required=True)
    parser.add_argument("--evaluation-contract", required=True)
    parser.add_argument("--embedding-contract-group", required=True)
    parser.add_argument("--split-manifest", required=True)
    parser.add_argument("--config-used", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def _format_seed(seed_value: object) -> str:
    if pd.isna(seed_value):
        return ""
    return str(seed_value).strip()


def _load_per_run_rows(run_root: Path, args: argparse.Namespace) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    split_manifest_df = pd.read_csv(args.split_manifest, sep="\t")
    split_version = str(split_manifest_df["split_version"].iloc[0]) if "split_version" in split_manifest_df.columns and not split_manifest_df.empty else "frozen_fixed_split"
    for metrics_path in sorted(run_root.glob("run_*/metrics.tsv")):
        metrics_df = pd.read_csv(metrics_path, sep="\t")
        split_lookup = {str(row["eval_split"]): row for _, row in metrics_df.iterrows()}
        val_row = split_lookup.get("val", {})
        test_row = split_lookup.get("test", {})
        seed_token = metrics_path.parent.name.replace("run_", "", 1)
        row = {
            "target": args.target,
            "protocol": args.protocol,
            "species": args.species,
            "regime": args.regime,
            "model": args.model_name,
            "feature_setting": args.feature_setting,
            "label_regime": args.regime,
            "run_id": f"seed_{seed_token}",
            "seed": seed_token,
            "is_deterministic": "false",
            "split_version": split_version,
            "evaluation_contract": args.evaluation_contract,
            "feature_contract_group": args.embedding_contract_group,
            "split_manifest": args.split_manifest,
            "config_used": args.config_used,
            "embedding_cache_scope": "species",
        }
        for split_name, source_row in [("val", val_row), ("test", test_row)]:
            for metric in METRICS:
                row[f"{split_name}_{metric}"] = source_row.get(metric, pd.NA)
        rows.append(row)
    return pd.DataFrame(rows)


def _aggregate(per_run: pd.DataFrame) -> pd.DataFrame:
    group_columns = [
        "target",
        "protocol",
        "species",
        "regime",
        "model",
        "feature_setting",
        "label_regime",
        "split_version",
        "evaluation_contract",
        "feature_contract_group",
        "split_manifest",
        "config_used",
        "embedding_cache_scope",
        "is_deterministic",
    ]
    metric_columns = [column for column in per_run.columns if column.startswith("val_") or column.startswith("test_")]
    if per_run.empty:
        return pd.DataFrame(columns=group_columns + ["n_runs", "seed_list", "run_ids"] + [f"{column}_{suffix}" for column in metric_columns for suffix in ["mean", "std"]])

    rows: list[dict[str, object]] = []
    for group_key, group_df in per_run.groupby(group_columns, dropna=False, sort=True):
        row = {column: value for column, value in zip(group_columns, group_key)}
        row["n_runs"] = int(len(group_df))
        row["seed_list"] = ",".join(seed for seed in (_format_seed(value) for value in group_df["seed"]) if seed)
        row["run_ids"] = ",".join(group_df["run_id"].astype(str).tolist())
        for metric_column in metric_columns:
            values = pd.to_numeric(group_df[metric_column], errors="coerce")
            row[f"{metric_column}_mean"] = float(values.mean())
            row[f"{metric_column}_std"] = float(values.std(ddof=0))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["target", "model", "feature_setting"], kind="stable").reset_index(drop=True)


def main() -> None:
    args = parse_args()
    run_root = Path(args.run_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    per_run = _load_per_run_rows(run_root, args)
    aggregated = _aggregate(per_run)
    publication = build_publication_summary(aggregated)

    per_run.to_csv(output_dir / "per_run_metrics.tsv", sep="\t", index=False)
    aggregated.to_csv(output_dir / "aggregated_metrics.tsv", sep="\t", index=False)
    publication.to_csv(output_dir / "final_summary.tsv", sep="\t", index=False)
    (output_dir / "per_run_metrics.md").write_text(
        "# Per-Run Metrics\n\n" + (per_run.to_markdown(index=False) if not per_run.empty else "No per-run rows available.") + "\n",
        encoding="utf-8",
    )
    (output_dir / "aggregated_metrics.md").write_text(
        "# Aggregated Metrics\n\n" + (aggregated.to_markdown(index=False) if not aggregated.empty else "No aggregated rows available.") + "\n",
        encoding="utf-8",
    )
    (output_dir / "final_summary.md").write_text(
        publication_markdown(
            publication,
            "Final Summary",
            "Publication-facing Phase 1 benchmark summary with frozen split manifests, five seeds, and mean ± std test metrics.",
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
