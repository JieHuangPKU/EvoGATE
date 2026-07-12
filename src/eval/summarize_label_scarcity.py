import argparse
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


METRICS = ["test_auprc", "test_auroc", "test_mcc", "test_f1", "test_precision", "test_recall"]
MAIN_MODELS = ["GraphSAGE", "MLP", "SVM", "RF", "N2V_MLP", "DC", "CC"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize label-scarcity benchmark outputs.")
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--split-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--feature-setting", required=True, type=str)
    parser.add_argument("--seeds", required=True, type=str)
    parser.add_argument("--train-fractions", required=True, type=str)
    return parser.parse_args()


def parse_fraction_csv(text: str) -> list[float]:
    return [float(token.strip()) for token in str(text).split(",") if token.strip()]


def parse_seed_csv(text: str) -> list[int]:
    return [int(token.strip()) for token in str(text).split(",") if token.strip()]


def fraction_tag(value: float) -> str:
    return f"{int(round(float(value) * 100)):02d}"


def collect_metrics(output_root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    pattern = re.compile(r"train_fraction_(\d{2})")
    for root, _, files in os.walk(output_root):
        if "metrics.tsv" not in files:
            continue
        metrics_path = Path(root) / "metrics.tsv"
        row = pd.read_csv(metrics_path, sep="\t").iloc[0].to_dict()
        fraction_match = pattern.search(str(metrics_path))
        if not fraction_match:
            raise ValueError(f"Could not parse train_fraction from path: {metrics_path}")
        row["train_fraction_tag"] = fraction_match.group(1)
        row["train_fraction"] = int(fraction_match.group(1)) / 100.0
        row["metrics_path"] = str(metrics_path)
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["seed"] = pd.to_numeric(df["seed"], errors="coerce").astype("Int64")
    df["train_fraction"] = pd.to_numeric(df["train_fraction"], errors="raise")
    return df.sort_values(["model", "train_fraction", "seed"], kind="stable").reset_index(drop=True)


def aggregate_metrics(per_run: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = ["protocol", "species", "regime", "model", "feature_setting", "train_fraction", "train_fraction_tag"]
    for keys, group_df in per_run.groupby(group_cols, dropna=False, sort=True):
        row = {column: value for column, value in zip(group_cols, keys)}
        row["n_runs"] = int(len(group_df))
        row["seed_list"] = ",".join(str(int(seed)) for seed in pd.to_numeric(group_df["seed"], errors="coerce").dropna().astype(int).tolist())
        row["split_version_count"] = int(group_df["split_version"].astype(str).nunique())
        row["split_manifest_count"] = int(group_df["split_manifest"].astype(str).nunique())
        for metric in METRICS:
            values = pd.to_numeric(group_df[metric], errors="coerce")
            row[f"mean_{metric.replace('test_', '').upper()}"] = float(values.mean())
            row[f"sd_{metric.replace('test_', '').upper()}"] = float(values.std(ddof=0))
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values(["model", "train_fraction"], kind="stable").reset_index(drop=True)
    baseline = (
        summary.loc[summary["train_fraction"].eq(0.90), ["model", "mean_AUPRC"]]
        .rename(columns={"mean_AUPRC": "baseline_AUPRC_at_0.90"})
        .reset_index(drop=True)
    )
    summary = summary.merge(baseline, on="model", how="left")
    summary["performance_retention_AUPRC"] = summary["mean_AUPRC"] / summary["baseline_AUPRC_at_0.90"]
    summary["performance_drop_AUPRC"] = summary["baseline_AUPRC_at_0.90"] - summary["mean_AUPRC"]
    return summary


def build_ranking_table(summary: pd.DataFrame) -> pd.DataFrame:
    low = summary.loc[summary["train_fraction"].eq(0.10), ["model", "mean_AUPRC"]].rename(columns={"mean_AUPRC": "AUPRC@10%"})
    high = summary.loc[summary["train_fraction"].eq(0.90), ["model", "mean_AUPRC"]].rename(columns={"mean_AUPRC": "AUPRC@90%"})
    retention = summary.loc[summary["train_fraction"].eq(0.10), ["model", "performance_retention_AUPRC"]].rename(
        columns={"performance_retention_AUPRC": "retention_AUPRC"}
    )
    ranking = low.merge(high, on="model", how="outer").merge(retention, on="model", how="outer")
    ranking = ranking[ranking["model"].isin(MAIN_MODELS)].copy()
    ranking["rank_low_label"] = ranking["AUPRC@10%"].rank(method="min", ascending=False).astype("Int64")
    ranking["rank_retention"] = ranking["retention_AUPRC"].rank(method="min", ascending=False).astype("Int64")
    order = {model: idx for idx, model in enumerate(MAIN_MODELS)}
    ranking["plot_order"] = ranking["model"].map(order)
    return ranking.sort_values(["rank_low_label", "plot_order"], kind="stable").drop(columns=["plot_order"]).reset_index(drop=True)


def build_coverage_audit(per_run: pd.DataFrame, split_root: Path, expected_models: list[str], expected_fractions: list[float], expected_seeds: list[int]) -> pd.DataFrame:
    split_index = pd.read_csv(split_root / "label_scarcity_split_manifest_index.tsv", sep="\t", dtype=str).fillna("")
    split_manifest_paths = [Path(value) for value in split_index["split_manifest_path"].astype(str).tolist()]

    test_sets: dict[str, tuple[str, ...]] = {}
    for path in split_manifest_paths:
        split_df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        test_genes = tuple(sorted(split_df.loc[split_df["split"].eq("test"), "graph_gene_id"].astype(str).tolist()))
        test_sets[str(path)] = test_genes

    unique_test_sets = {genes for genes in test_sets.values()}
    expected_fraction_tags = [fraction_tag(value) for value in expected_fractions]
    observed_grid = (
        per_run[["model", "train_fraction_tag", "seed"]]
        .drop_duplicates()
        .assign(seed=lambda df: pd.to_numeric(df["seed"], errors="coerce").astype("Int64"))
    )

    audit_rows: list[dict[str, object]] = []
    audit_rows.append(
        {
            "audit_check": "RF_present",
            "status": "PASS" if "RF" in set(per_run["model"].astype(str)) else "FAIL",
            "detail": "RF detected in benchmark outputs" if "RF" in set(per_run["model"].astype(str)) else "RF missing from benchmark outputs",
        }
    )
    audit_rows.append(
        {
            "audit_check": "train_fraction_complete",
            "status": "PASS" if set(split_index["train_fraction_tag"].astype(str)) == set(expected_fraction_tags) else "FAIL",
            "detail": ",".join(sorted(split_index["train_fraction_tag"].astype(str).unique().tolist())),
        }
    )
    audit_rows.append(
        {
            "audit_check": "seed_complete",
            "status": "PASS" if set(pd.to_numeric(split_index["seed"], errors="coerce").dropna().astype(int).tolist()) == set(expected_seeds) else "FAIL",
            "detail": ",".join(str(seed) for seed in sorted(pd.to_numeric(split_index["seed"], errors="coerce").dropna().astype(int).unique().tolist())),
        }
    )
    audit_rows.append(
        {
            "audit_check": "test_set_fixed_across_all_manifests",
            "status": "PASS" if len(unique_test_sets) == 1 else "FAIL",
            "detail": f"unique_test_sets={len(unique_test_sets)}",
        }
    )

    for model in expected_models:
        model_df = observed_grid.loc[observed_grid["model"].astype(str).eq(model)].copy()
        observed_pairs = {
            (str(row["train_fraction_tag"]), int(row["seed"]))
            for _, row in model_df.dropna(subset=["seed"]).iterrows()
        }
        missing = [
            f"{tag}:seed_{seed}"
            for tag in expected_fraction_tags
            for seed in expected_seeds
            if (tag, int(seed)) not in observed_pairs
        ]
        audit_rows.append(
            {
                "audit_check": f"grid_complete::{model}",
                "status": "PASS" if not missing else "FAIL",
                "detail": "" if not missing else ",".join(missing),
            }
        )

    return pd.DataFrame(audit_rows)


def build_report(summary: pd.DataFrame, ranking: pd.DataFrame, output_path: Path, protocol: str, feature_setting: str) -> None:
    low_label = ranking.sort_values(["rank_low_label", "rank_retention"], kind="stable").reset_index(drop=True)
    best_low = str(low_label.iloc[0]["model"]) if not low_label.empty else "NA"
    strongest_retention = ranking.sort_values("rank_retention", kind="stable").reset_index(drop=True)
    best_retention = str(strongest_retention.iloc[0]["model"]) if not strongest_retention.empty else "NA"
    collapse = ranking.loc[ranking["retention_AUPRC"].lt(0.50), "model"].astype(str).tolist()
    lines = [
        "# Label Scarcity Benchmark Report",
        "",
        f"- protocol: `{protocol}`",
        f"- feature_setting: `{feature_setting}`",
        f"- main robustness answer at 10% labels: `{best_low}`",
        f"- strongest retention answer: `{best_retention}`",
        f"- models with retention_AUPRC < 0.50 at 10%: `{', '.join(collapse) if collapse else 'none'}`",
        "",
        "## How To Run",
        "- Run `scripts/run_label_scarcity_benchmark.sh` to build the benchmark.",
        "- Run `Rscript src/plot/plot_label_scarcity.R --summary-dir results/Figure2_label_scarcity/summary --output-dir results/Figure2_label_scarcity/plots` to regenerate the figures.",
        "",
        "## Interpretation",
        "- Judge robustness primarily from `performance_retention_AUPRC` and secondarily from `AUPRC@10%`.",
        "- If `GraphSAGE` is top-ranked in both low-label AUPRC and retention, the benchmark supports the claim that graph-based learning is the most robust under label scarcity.",
        "- Classical ML degradation is reflected by lower `AUPRC@10%` and larger `performance_drop_AUPRC` relative to `GraphSAGE`.",
        "- `N2V_MLP` should be interpreted as a shallow topology-aware buffer model between GraphSAGE and classical/tabular baselines.",
        "- `DC` and `CC` provide topology-only reference floors rather than competitive predictors.",
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    split_root = Path(args.split_root)
    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)

    per_run = collect_metrics(output_root)
    if per_run.empty:
        raise FileNotFoundError(f"No metrics.tsv files found under {output_root}")
    summary = aggregate_metrics(per_run)
    ranking = build_ranking_table(summary)
    audit = build_coverage_audit(
        per_run,
        split_root=split_root,
        expected_models=MAIN_MODELS,
        expected_fractions=parse_fraction_csv(args.train_fractions),
        expected_seeds=parse_seed_csv(args.seeds),
    )

    per_run.to_csv(summary_dir / "label_scarcity_per_run_metrics.tsv", sep="\t", index=False)
    summary.to_csv(summary_dir / "label_scarcity_summary.tsv", sep="\t", index=False)
    ranking.to_csv(summary_dir / "label_scarcity_ranking_table.tsv", sep="\t", index=False)
    audit.to_csv(summary_dir / "label_scarcity_coverage_audit.tsv", sep="\t", index=False)
    build_report(
        summary=summary,
        ranking=ranking,
        output_path=summary_dir / "label_scarcity_report.md",
        protocol=args.protocol,
        feature_setting=args.feature_setting,
    )


if __name__ == "__main__":
    main()
