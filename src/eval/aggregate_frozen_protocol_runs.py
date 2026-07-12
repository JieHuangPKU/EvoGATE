import argparse
import os
from pathlib import Path

import pandas as pd

from src.eval.publication_summary import build_publication_summary, publication_markdown


GROUP_COLUMNS = [
    "protocol",
    "species",
    "regime",
    "model",
    "feature_setting",
    "feature_contract_group",
    "label_regime",
    "split_version",
    "graph_source",
    "graph_contract",
    "threshold_strategy",
    "evaluation_contract",
    "label_manifest",
    "split_manifest",
    "config_used",
    "is_deterministic",
    "require_true_node2vec",
    "embedding_method",
    "embedding_backend",
    "fallback_used",
    "esm2_dim",
]
METRIC_COLUMNS = [
    "test_auroc",
    "test_auprc",
    "test_mcc",
    "test_f1",
    "test_precision",
    "test_recall",
    "test_accuracy",
    "test_specificity",
    "val_auroc",
    "val_auprc",
    "val_mcc",
    "val_f1",
    "val_precision",
    "val_recall",
    "val_accuracy",
    "val_specificity",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate frozen protocol benchmark runs")
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    parser.add_argument("--output-prefix", default="", type=str)
    parser.add_argument(
        "--include-model-feature-combos",
        default="",
        type=str,
        help="Optional comma-separated model:feature_setting whitelist used to ignore stale outputs.",
    )
    return parser.parse_args()


def collect_metrics(output_root):
    rows = []
    metric_paths = []
    for root, _, files in os.walk(output_root, followlinks=True):
        if "metrics.tsv" in files:
            metric_paths.append(Path(root) / "metrics.tsv")
    for metrics_path in sorted(metric_paths):
        if not metrics_path.exists():
            continue
        row = pd.read_csv(metrics_path, sep="\t").iloc[0].to_dict()
        row["metrics_path"] = str(metrics_path)
        rows.append(row)
    expected_columns = GROUP_COLUMNS + ["run_id", "seed", "metrics_path"] + METRIC_COLUMNS
    if not rows:
        return pd.DataFrame(columns=expected_columns)
    per_run = pd.DataFrame(rows)
    for column in expected_columns:
        if column not in per_run.columns:
            per_run[column] = pd.NA
    return per_run


def parse_model_feature_combos(raw_value: str) -> set[tuple[str, str]]:
    combos: set[tuple[str, str]] = set()
    for item in raw_value.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid model/feature combo '{item}'. Expected MODEL:FEATURE_SETTING.")
        model, feature_setting = item.split(":", 1)
        combos.add((model.strip(), feature_setting.strip()))
    return combos


def filter_model_feature_combos(per_run: pd.DataFrame, combos: set[tuple[str, str]]) -> pd.DataFrame:
    if not combos or per_run.empty:
        return per_run
    df = per_run.copy()
    combo_index = pd.MultiIndex.from_frame(df[["model", "feature_setting"]].astype(str))
    allowed_index = pd.MultiIndex.from_tuples(sorted(combos), names=["model", "feature_setting"])
    return df[combo_index.isin(allowed_index)].reset_index(drop=True)


def aggregate_runs(per_run):
    per_run = per_run.copy()
    for metric in METRIC_COLUMNS:
        if metric not in per_run.columns:
            per_run[metric] = pd.NA
    if per_run.empty:
        columns = GROUP_COLUMNS + ["n_runs", "seed_list", "run_ids"] + [f"{metric}_{suffix}" for metric in METRIC_COLUMNS for suffix in ["mean", "std"]]
        return pd.DataFrame(columns=columns)

    rows = []
    for group_key, group_df in per_run.groupby(GROUP_COLUMNS, dropna=False, sort=True):
        row = {column: value for column, value in zip(GROUP_COLUMNS, group_key)}
        row["n_runs"] = int(len(group_df))
        seed_values = [str(value) for value in group_df["seed"].tolist() if str(value).strip()]
        row["seed_list"] = ",".join(seed_values)
        row["run_ids"] = ",".join(group_df["run_id"].astype(str).tolist())
        for metric in METRIC_COLUMNS:
            values = pd.to_numeric(group_df[metric], errors="coerce")
            row[f"{metric}_mean"] = float(values.mean())
            row[f"{metric}_std"] = float(values.std(ddof=0))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["protocol", "model", "feature_setting"], kind="stable").reset_index(drop=True)


def dataframe_to_markdown(df: pd.DataFrame, title: str, intro: str) -> str:
    lines = [f"# {title}", "", intro, ""]
    if df.empty:
        lines.extend(["No rows available.", ""])
    else:
        lines.append(df.to_markdown(index=False))
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)

    include_combos = parse_model_feature_combos(args.include_model_feature_combos)
    per_run = filter_model_feature_combos(collect_metrics(output_root), include_combos)
    aggregated = aggregate_runs(per_run)
    publication = build_publication_summary(aggregated)

    prefix = args.output_prefix.strip()
    per_run_name = f"{prefix}_per_run_metrics" if prefix else "per_run_metrics"
    aggregated_name = f"{prefix}_aggregated_metrics" if prefix else "aggregated_metrics"
    publication_name = f"{prefix}_publication_summary" if prefix else "final_summary"

    per_run.to_csv(summary_dir / f"{per_run_name}.tsv", sep="\t", index=False)
    aggregated.to_csv(summary_dir / f"{aggregated_name}.tsv", sep="\t", index=False)
    publication.to_csv(summary_dir / f"{publication_name}.tsv", sep="\t", index=False)
    (summary_dir / f"{per_run_name}.md").write_text(
        dataframe_to_markdown(per_run, "Per-Run Metrics", "Per-run benchmark rows collected from the output root."),
        encoding="utf-8",
    )
    (summary_dir / f"{aggregated_name}.md").write_text(
        dataframe_to_markdown(aggregated, "Aggregated Metrics", "Aggregated benchmark rows grouped by protocol/model/feature contract."),
        encoding="utf-8",
    )
    (summary_dir / f"{publication_name}.md").write_text(
        publication_markdown(
            publication,
            "Publication Summary",
            "Publication-facing benchmark summary with fixed split, five seeds, and mean ± std test metrics.",
        ),
        encoding="utf-8",
    )

    if prefix:
        per_run.to_csv(summary_dir / "per_run_metrics.tsv", sep="\t", index=False)
        aggregated.to_csv(summary_dir / "aggregated_metrics.tsv", sep="\t", index=False)
        publication.to_csv(summary_dir / "final_summary.tsv", sep="\t", index=False)
        (summary_dir / "per_run_metrics.md").write_text(
            dataframe_to_markdown(per_run, "Per-Run Metrics", "Per-run benchmark rows collected from the output root."),
            encoding="utf-8",
        )
        (summary_dir / "aggregated_metrics.md").write_text(
            dataframe_to_markdown(
                aggregated,
                "Aggregated Metrics",
                "Aggregated benchmark rows grouped by protocol/model/feature contract.",
            ),
            encoding="utf-8",
        )
        (summary_dir / "final_summary.md").write_text(
            publication_markdown(
                publication,
                "Final Summary",
                "Legacy compatibility summary; prefer the prefixed publication summary files.",
            ),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
