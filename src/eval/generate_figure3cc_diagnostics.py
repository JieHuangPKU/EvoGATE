from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.eval.aggregate_frozen_protocol_runs import dataframe_to_markdown
from src.eval.summarize_figure3c_threshold_tuned_metrics import (
    compute_score_metrics,
    compute_threshold_metrics,
    find_best_threshold,
    parse_bool_series,
    pr_curve_points,
)


FEATURE_SPECS = [
    ("baseline concat", "ORT_EXP_SUB_ESM2", "comparison"),
    ("old gated", "ORT_EXP_SUB_ESM2_GATED", "comparison"),
    ("residual gated", "ORT_EXP_SUB_ESM2_GATED_RESIDUAL", "comparison"),
    ("residual gated + WBCE", "ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE", "comparison"),
    ("old gated + WBCE", "ORT_EXP_SUB_ESM2_OLD_GATED_WBCE", "new"),
]
BOOL_TRUE = {"true", "t", "1", "yes", "y"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Figure3cC comparison and diagnostic summaries.")
    parser.add_argument("--comparison-root", required=True, type=str)
    parser.add_argument("--new-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    return parser.parse_args()


def load_predictions(predictions_path: Path) -> pd.DataFrame:
    df = pd.read_csv(predictions_path, sep="\t")
    df["is_labeled_flag"] = parse_bool_series(df["is_labeled"])
    df["label"] = pd.to_numeric(df["label"], errors="coerce")
    df["pred_score"] = pd.to_numeric(df["pred_score"], errors="coerce")
    filtered = df.loc[df["is_labeled_flag"]].copy()
    filtered = filtered.loc[filtered["split"].isin(["val", "test"])].copy()
    return filtered.dropna(subset=["label", "pred_score"]).copy()


def expected_setting_root(comparison_root: Path, new_root: Path, feature_setting: str, location: str) -> Path:
    base_root = comparison_root if location == "comparison" else new_root
    return base_root / "fgraminearum_newlabel" / "GraphSAGE" / feature_setting


def probability_summary(scores: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    positives = scores[labels == 1]
    negatives = scores[labels == 0]
    return {
        "mean": float(np.mean(scores)),
        "std": float(np.std(scores, ddof=0)),
        "min": float(np.min(scores)),
        "max": float(np.max(scores)),
        "positive_class_mean_score": float(np.mean(positives)) if positives.size else float("nan"),
        "negative_class_mean_score": float(np.mean(negatives)) if negatives.size else float("nan"),
    }


def calibration_summary(scores: np.ndarray, labels: np.ndarray, bins: int = 10) -> dict[str, float]:
    labels = labels.astype(float)
    scores = scores.astype(float)
    brier = float(np.mean((scores - labels) ** 2))
    bin_edges = np.linspace(0.0, 1.0, bins + 1)
    ece = 0.0
    rows = 0
    for idx in range(bins):
        left = bin_edges[idx]
        right = bin_edges[idx + 1]
        if idx == bins - 1:
            mask = (scores >= left) & (scores <= right)
        else:
            mask = (scores >= left) & (scores < right)
        if not np.any(mask):
            continue
        rows += 1
        bin_scores = scores[mask]
        bin_labels = labels[mask]
        ece += abs(float(np.mean(bin_scores)) - float(np.mean(bin_labels))) * (float(bin_scores.size) / float(scores.size))
    return {
        "brier_score": brier,
        "ece_10bin": float(ece),
        "nonempty_bins": int(rows),
    }


def build_rows(comparison_root: Path, new_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    threshold_rows = []
    probability_rows = []
    calibration_rows = []
    pr_rows = []

    for comparison_name, feature_setting, location in FEATURE_SPECS:
        setting_root = expected_setting_root(comparison_root, new_root, feature_setting, location)
        metrics_paths = sorted(setting_root.glob("run_*/metrics.tsv"))
        if not metrics_paths:
            raise FileNotFoundError(f"No metrics.tsv files found under {setting_root}")
        for metrics_path in metrics_paths:
            run_dir = metrics_path.parent
            predictions_path = run_dir / "predictions.tsv"
            if not predictions_path.exists():
                raise FileNotFoundError(f"Missing predictions.tsv for {run_dir}")
            metrics_row = pd.read_csv(metrics_path, sep="\t").iloc[0].to_dict()
            predictions = load_predictions(predictions_path)
            val_df = predictions.loc[predictions["split"] == "val"].copy()
            test_df = predictions.loc[predictions["split"] == "test"].copy()
            if val_df.empty or test_df.empty:
                raise ValueError(f"Expected non-empty val/test rows in {predictions_path}")

            val_labels = val_df["label"].to_numpy(dtype=int)
            val_scores = val_df["pred_score"].to_numpy(dtype=float)
            test_labels = test_df["label"].to_numpy(dtype=int)
            test_scores = test_df["pred_score"].to_numpy(dtype=float)

            val_score_metrics = compute_score_metrics(val_labels, val_scores)
            test_score_metrics = compute_score_metrics(test_labels, test_scores)
            best_threshold_by_f1, best_f1_metrics = find_best_threshold(val_labels, val_scores, "f1", lambda m: m["f1"])
            best_threshold_by_mcc, best_mcc_metrics = find_best_threshold(val_labels, val_scores, "mcc", lambda m: m["mcc"])
            test_metrics_at_best_mcc = compute_threshold_metrics(test_labels, test_scores, float(best_threshold_by_mcc))

            common = {
                "comparison": comparison_name,
                "feature_setting": feature_setting,
                "protocol": metrics_row.get("protocol", ""),
                "run_id": metrics_row.get("run_id", ""),
                "seed": metrics_row.get("seed", ""),
                "metrics_path": str(metrics_path.resolve()),
                "predictions_path": str(predictions_path.resolve()),
            }

            threshold_rows.append(
                {
                    **common,
                    "best_threshold": float(best_threshold_by_mcc),
                    "best_threshold_by_mcc": float(best_threshold_by_mcc),
                    "best_threshold_by_f1": float(best_threshold_by_f1),
                    "val_auroc": val_score_metrics["auroc"],
                    "val_auprc": val_score_metrics["auprc"],
                    "test_auroc": test_score_metrics["auroc"],
                    "test_auprc": test_score_metrics["auprc"],
                    "val_mcc_at_best_threshold": float(best_mcc_metrics["mcc"]),
                    "val_f1_at_best_threshold": float(best_mcc_metrics["f1"]),
                    "val_precision_at_best_threshold": float(best_mcc_metrics["precision"]),
                    "val_recall_at_best_threshold": float(best_mcc_metrics["recall"]),
                    "test_mcc_at_best_threshold": float(test_metrics_at_best_mcc["mcc"]),
                    "test_f1_at_best_threshold": float(test_metrics_at_best_mcc["f1"]),
                    "test_precision_at_best_threshold": float(test_metrics_at_best_mcc["precision"]),
                    "test_recall_at_best_threshold": float(test_metrics_at_best_mcc["recall"]),
                    "test_specificity_at_best_threshold": float(test_metrics_at_best_mcc["specificity"]),
                }
            )

            probability_rows.append(
                {
                    **common,
                    "split": "test",
                    **probability_summary(test_scores, test_labels),
                }
            )

            calibration_rows.append(
                {
                    **common,
                    "split": "test",
                    **calibration_summary(test_scores, test_labels),
                }
            )

            recalls, precisions, auprc = pr_curve_points(test_labels, test_scores)
            for idx, (recall, precision) in enumerate(zip(recalls, precisions)):
                pr_rows.append(
                    {
                        **common,
                        "curve_index": int(idx),
                        "recall": float(recall),
                        "precision": float(precision),
                        "auprc": float(auprc),
                    }
                )

    return (
        pd.DataFrame(threshold_rows),
        pd.DataFrame(probability_rows),
        pd.DataFrame(calibration_rows),
        pd.DataFrame(pr_rows),
    )


def aggregate_numeric(df: pd.DataFrame, group_cols: list[str], metric_cols: list[str]) -> pd.DataFrame:
    rows = []
    for group_key, group_df in df.groupby(group_cols, sort=True, dropna=False):
        row = {column: value for column, value in zip(group_cols, group_key)}
        row["runs"] = int(len(group_df))
        for metric in metric_cols:
            values = pd.to_numeric(group_df[metric], errors="coerce")
            row[f"{metric}_mean"] = float(values.mean())
            row[f"{metric}_std"] = float(values.std(ddof=0))
        rows.append(row)
    return pd.DataFrame(rows)


def write_table(path_tsv: Path, path_md: Path, dataframe: pd.DataFrame, title: str, intro: str) -> None:
    dataframe.to_csv(path_tsv, sep="\t", index=False)
    path_md.write_text(dataframe_to_markdown(dataframe, title, intro), encoding="utf-8")


def main() -> None:
    args = parse_args()
    comparison_root = Path(args.comparison_root).resolve()
    new_root = Path(args.new_root).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    summary_dir.mkdir(parents=True, exist_ok=True)

    threshold_rows, probability_rows, calibration_rows, pr_rows = build_rows(comparison_root, new_root)

    group_cols = ["comparison", "feature_setting"]
    threshold_summary = aggregate_numeric(
        threshold_rows,
        group_cols,
        [
            "best_threshold",
            "best_threshold_by_f1",
            "best_threshold_by_mcc",
            "val_auroc",
            "val_auprc",
            "test_auroc",
            "test_auprc",
            "val_mcc_at_best_threshold",
            "val_f1_at_best_threshold",
            "val_precision_at_best_threshold",
            "val_recall_at_best_threshold",
            "test_mcc_at_best_threshold",
            "test_f1_at_best_threshold",
            "test_precision_at_best_threshold",
            "test_recall_at_best_threshold",
            "test_specificity_at_best_threshold",
        ],
    )
    probability_summary_df = aggregate_numeric(
        probability_rows,
        group_cols,
        [
            "mean",
            "std",
            "min",
            "max",
            "positive_class_mean_score",
            "negative_class_mean_score",
        ],
    )
    calibration_summary_df = aggregate_numeric(
        calibration_rows,
        group_cols,
        ["brier_score", "ece_10bin", "nonempty_bins"],
    )

    comparison_summary = threshold_summary.merge(
        probability_summary_df,
        on=["comparison", "feature_setting", "runs"],
        how="left",
    )
    comparison_summary = comparison_summary.rename(
        columns={
            "test_auroc_mean": "AUROC",
            "test_auprc_mean": "AUPRC",
            "test_mcc_at_best_threshold_mean": "MCC",
            "test_f1_at_best_threshold_mean": "F1",
            "test_precision_at_best_threshold_mean": "Precision",
            "test_recall_at_best_threshold_mean": "Recall",
            "best_threshold_mean": "best_threshold",
        }
    )

    write_table(
        summary_dir / "Figure3cC_comparison_summary.tsv",
        summary_dir / "Figure3cC_comparison_summary.md",
        comparison_summary,
        "Figure3cC Comparison Summary",
        "Five-way comparison with validation-selected thresholds and test-set threshold metrics.",
    )
    write_table(
        summary_dir / "Figure3cC_threshold_diagnostics.tsv",
        summary_dir / "Figure3cC_threshold_diagnostics.md",
        threshold_rows.sort_values(["comparison", "seed"], kind="stable").reset_index(drop=True),
        "Figure3cC Threshold Diagnostics",
        "Per-seed validation-selected threshold diagnostics. `best_threshold` uses the validation MCC optimum.",
    )
    write_table(
        summary_dir / "Figure3cC_probability_summary.tsv",
        summary_dir / "Figure3cC_probability_summary.md",
        probability_rows.sort_values(["comparison", "seed"], kind="stable").reset_index(drop=True),
        "Figure3cC Probability Summary",
        "Per-seed test prediction probability summary.",
    )
    write_table(
        summary_dir / "Figure3cC_calibration_summary.tsv",
        summary_dir / "Figure3cC_calibration_summary.md",
        calibration_rows.sort_values(["comparison", "seed"], kind="stable").reset_index(drop=True),
        "Figure3cC Calibration Summary",
        "Per-seed test calibration summary with Brier score and 10-bin ECE.",
    )
    write_table(
        summary_dir / "Figure3cC_pr_curve_data.tsv",
        summary_dir / "Figure3cC_pr_curve_data.md",
        pr_rows.sort_values(["comparison", "seed", "curve_index"], kind="stable").reset_index(drop=True),
        "Figure3cC PR Curve Data",
        "Per-seed test precision-recall curve points.",
    )


if __name__ == "__main__":
    main()
