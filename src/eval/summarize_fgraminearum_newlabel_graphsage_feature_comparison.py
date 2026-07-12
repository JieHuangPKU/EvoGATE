import argparse
from pathlib import Path

import pandas as pd

from src.eval.aggregate_frozen_protocol_runs import aggregate_runs, collect_metrics
from src.eval.publication_summary import format_float


PROTOCOL = "fgraminearum_newlabel"
MODEL = "GraphSAGE"
MODEL_VARIANTS = ["GraphSAGE_ORT_EXP_SUB", "GraphSAGE_ESM2", "GraphSAGE_ORT_EXP_SUB_ESM2"]
FEATURE_SETTINGS = ["ORT_EXP_SUB", "ESM2", "ORT_EXP_SUB_ESM2"]
METRIC_LABELS = [
    ("test_auroc", "AUROC"),
    ("test_auprc", "AUPRC"),
    ("test_mcc", "MCC"),
    ("test_f1", "F1"),
    ("test_precision", "Precision"),
    ("test_recall", "Recall"),
    ("test_specificity", "Specificity"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize the fgraminearum_newlabel GraphSAGE feature comparison benchmark"
    )
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    return parser.parse_args()


def dataframe_to_markdown(df: pd.DataFrame, title: str, intro: str) -> str:
    lines = [f"# {title}", "", intro, ""]
    if df.empty:
        lines.extend(["No rows available.", ""])
    else:
        lines.append(df.to_markdown(index=False))
        lines.append("")
    return "\n".join(lines)


def filter_scope(per_run: pd.DataFrame) -> pd.DataFrame:
    df = per_run.copy()
    if df.empty:
        return df
    if "model_variant" not in df.columns:
        raise ValueError("Expected column 'model_variant' in metrics rows, but it was not found")
    filtered = df[
        (df["protocol"].astype(str) == PROTOCOL)
        & (df["model"].astype(str) == MODEL)
        & (df["model_variant"].astype(str).isin(MODEL_VARIANTS))
        & (df["feature_setting"].astype(str).isin(FEATURE_SETTINGS))
    ].copy()
    filtered["feature_setting"] = pd.Categorical(
        filtered["feature_setting"], categories=FEATURE_SETTINGS, ordered=True
    )
    sort_columns = [column for column in ["feature_setting", "seed", "run_id"] if column in filtered.columns]
    return filtered.sort_values(sort_columns, kind="stable").reset_index(drop=True)


def validate_scope(aggregated: pd.DataFrame) -> None:
    found = sorted(aggregated["feature_setting"].astype(str).tolist()) if not aggregated.empty else []
    expected = sorted(FEATURE_SETTINGS)
    if found != expected:
        raise ValueError(
            "Expected exactly these GraphSAGE feature settings in the summary scope: "
            f"{expected}; found: {found}"
        )


def normalize_seed_list(seed_list: object) -> str:
    tokens = []
    for token in str(seed_list).split(","):
        value = str(token).strip()
        if not value:
            continue
        try:
            numeric = float(value)
        except ValueError:
            tokens.append(value)
            continue
        if numeric.is_integer():
            tokens.append(str(int(numeric)))
        else:
            tokens.append(value)
    return ",".join(tokens)


def build_final_summary(aggregated: pd.DataFrame) -> pd.DataFrame:
    summary = aggregated.copy()
    summary["feature_setting"] = pd.Categorical(
        summary["feature_setting"], categories=FEATURE_SETTINGS, ordered=True
    )
    summary = summary.sort_values(["feature_setting"], kind="stable").reset_index(drop=True)

    rows = []
    for _, row in summary.iterrows():
        out = {
            "Comparison": f"{MODEL} + {row['feature_setting']}",
            "Target": str(row.get("protocol", PROTOCOL)),
            "Model": str(row.get("model", MODEL)),
            "Feature_Setting": str(row["feature_setting"]),
            "Runs": int(row.get("n_runs", 0)),
            "Seed_List": normalize_seed_list(row.get("seed_list", "")),
        }
        for technical_name, public_name in METRIC_LABELS:
            out[f"{public_name}_mean"] = float(row[f"{technical_name}_mean"])
            out[f"{public_name}_std"] = float(row[f"{technical_name}_std"])
        rows.append(out)
    return pd.DataFrame(rows)


def mode_intro(final_summary: pd.DataFrame) -> str:
    if final_summary.empty:
        return "No benchmark rows were found for the requested GraphSAGE comparison scope."
    run_counts = sorted({int(value) for value in final_summary["Runs"].tolist()})
    if run_counts == [1]:
        return (
            "Single-seed mode summary for `fgraminearum_newlabel` under the GraphSAGE feature comparison subset. "
            "All `_std` columns are reported as `0.000` because each feature setting has one run."
        )
    if run_counts == [5]:
        return (
            "Five-seed summary for `fgraminearum_newlabel` under the GraphSAGE feature comparison subset. "
            "Each row reports mean ± std across the five frozen seeds."
        )
    return (
        "Summary for `fgraminearum_newlabel` under the GraphSAGE feature comparison subset. "
        "Rows are restricted to the three requested feature settings and report the observed run counts."
    )


def build_final_markdown(final_summary: pd.DataFrame) -> str:
    lines = ["# Final Summary", "", mode_intro(final_summary), ""]
    if final_summary.empty:
        lines.extend(["No benchmark rows available.", ""])
        return "\n".join(lines)

    display_df = final_summary.copy()
    for _, public_name in METRIC_LABELS:
        display_df[public_name] = display_df.apply(
            lambda row: f"{format_float(row[f'{public_name}_mean'])} ± {format_float(row[f'{public_name}_std'])}",
            axis=1,
        )
    metric_columns = [public_name for _, public_name in METRIC_LABELS]
    display_df = display_df[
        ["Comparison", "Runs"] + metric_columns + ["Seed_List", "Feature_Setting", "Model", "Target"]
    ].copy()
    lines.append(display_df.to_markdown(index=False))
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)

    per_run = filter_scope(collect_metrics(output_root))
    aggregated = aggregate_runs(per_run)
    if not aggregated.empty:
        aggregated["feature_setting"] = pd.Categorical(
            aggregated["feature_setting"], categories=FEATURE_SETTINGS, ordered=True
        )
        aggregated = aggregated.sort_values(["feature_setting"], kind="stable").reset_index(drop=True)
    validate_scope(aggregated)
    final_summary = build_final_summary(aggregated)

    per_run.to_csv(summary_dir / "per_run_metrics.tsv", sep="\t", index=False)
    aggregated.to_csv(summary_dir / "aggregated_metrics.tsv", sep="\t", index=False)
    final_summary.to_csv(summary_dir / "final_summary.tsv", sep="\t", index=False)
    (summary_dir / "final_summary.md").write_text(build_final_markdown(final_summary), encoding="utf-8")


if __name__ == "__main__":
    main()
