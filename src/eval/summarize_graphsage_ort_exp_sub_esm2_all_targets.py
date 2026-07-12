import argparse
from pathlib import Path

import pandas as pd

from src.eval.aggregate_frozen_protocol_runs import aggregate_runs, collect_metrics
from src.eval.publication_summary import build_publication_summary, publication_markdown


TARGET_PROTOCOLS = [
    "human",
    "celegans",
    "scerevisiae",
    "dmelanogaster",
    "fgraminearum_oldlabel",
    "fgraminearum_newlabel",
]
MODEL = "GraphSAGE"
MODEL_VARIANT = "GraphSAGE_ORT_EXP_SUB_ESM2"
FEATURE_SETTING = "ORT_EXP_SUB_ESM2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize GraphSAGE + ORT_EXP_SUB_ESM2 across all frozen benchmark targets"
    )
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    return parser.parse_args()


def filter_scope(per_run: pd.DataFrame) -> pd.DataFrame:
    df = per_run.copy()
    if df.empty:
        return df
    if "model_variant" not in df.columns:
        raise ValueError("Expected column 'model_variant' in metrics rows, but it was not found")
    filtered = df[
        (df["protocol"].astype(str).isin(TARGET_PROTOCOLS))
        & (df["model"].astype(str) == MODEL)
        & (df["model_variant"].astype(str) == MODEL_VARIANT)
        & (df["feature_setting"].astype(str) == FEATURE_SETTING)
    ].copy()
    filtered["protocol"] = pd.Categorical(filtered["protocol"], categories=TARGET_PROTOCOLS, ordered=True)
    sort_columns = [column for column in ["protocol", "seed", "run_id"] if column in filtered.columns]
    return filtered.sort_values(sort_columns, kind="stable").reset_index(drop=True)


def validate_scope(aggregated: pd.DataFrame) -> None:
    found = sorted(aggregated["protocol"].astype(str).tolist()) if not aggregated.empty else []
    expected = sorted(TARGET_PROTOCOLS)
    if found != expected:
        raise ValueError(
            "Expected exactly these protocols in the GraphSAGE ORT_EXP_SUB_ESM2 summary scope: "
            f"{expected}; found: {found}"
        )


def mode_intro(publication: pd.DataFrame) -> str:
    if publication.empty:
        return "No benchmark rows were found for the requested GraphSAGE ORT_EXP_SUB_ESM2 scope."
    run_counts = sorted({int(value) for value in publication["Runs"].tolist()})
    if run_counts == [1]:
        return (
            "Single-seed summary for GraphSAGE + ORT_EXP_SUB_ESM2 across all six frozen benchmark targets. "
            "All `_std` columns are reported as `0.000` because each target has one run."
        )
    if run_counts == [5]:
        return (
            "Five-seed summary for GraphSAGE + ORT_EXP_SUB_ESM2 across all six frozen benchmark targets. "
            "Each row reports mean ± std across the five frozen seeds."
        )
    return (
        "Summary for GraphSAGE + ORT_EXP_SUB_ESM2 across all six frozen benchmark targets. "
        "Rows are restricted to the requested model-feature subset and report the observed run counts."
    )


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)

    per_run = filter_scope(collect_metrics(output_root))
    aggregated = aggregate_runs(per_run)
    if not aggregated.empty:
        aggregated["protocol"] = pd.Categorical(aggregated["protocol"], categories=TARGET_PROTOCOLS, ordered=True)
        aggregated = aggregated.sort_values(["protocol"], kind="stable").reset_index(drop=True)
    validate_scope(aggregated)
    publication = build_publication_summary(aggregated)

    per_run.to_csv(summary_dir / "per_run_metrics.tsv", sep="\t", index=False)
    aggregated.to_csv(summary_dir / "aggregated_metrics.tsv", sep="\t", index=False)
    publication.to_csv(summary_dir / "final_summary.tsv", sep="\t", index=False)
    (summary_dir / "final_summary.md").write_text(
        publication_markdown(
            publication,
            "Final Summary",
            mode_intro(publication),
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
