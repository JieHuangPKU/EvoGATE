from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from src.eval.aggregate_frozen_protocol_runs import dataframe_to_markdown


BOOL_TRUE = {"true", "t", "1", "yes", "y"}
THRESHOLD_MODES = [
    ("fixed_0.5", "fixed", 0.5),
    ("val_f1_opt", "validation", None),
    ("val_mcc_opt", "validation", None),
]
FINAL_METRICS = [
    ("test_auroc_mean", "AUROC_mean"),
    ("test_auroc_std", "AUROC_std"),
    ("test_auprc_mean", "AUPRC_mean"),
    ("test_auprc_std", "AUPRC_std"),
    ("test_mcc_mean", "MCC_mean"),
    ("test_mcc_std", "MCC_std"),
    ("test_f1_mean", "F1_mean"),
    ("test_f1_std", "F1_std"),
    ("test_precision_mean", "Precision_mean"),
    ("test_precision_std", "Precision_std"),
    ("test_recall_mean", "Recall_mean"),
    ("test_recall_std", "Recall_std"),
    ("test_specificity_mean", "Specificity_mean"),
    ("test_specificity_std", "Specificity_std"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Post-process Figure3c predictions with validation-based threshold tuning.")
    parser.add_argument(
        "--figure3c-root",
        default="outputs/Figure3c/fgraminearum_newlabel/GraphSAGE",
        type=str,
        help="Root directory containing Figure3c feature-setting subdirectories.",
    )
    parser.add_argument(
        "--summary-dir",
        default="results/Figure3c_threshold_tuned",
        type=str,
        help="Directory for threshold-tuned summaries.",
    )
    parser.add_argument(
        "--feature-settings",
        default="ORT_EXP_SUB_ESM2,ORT_EXP_SUB_ESM2_GATED",
        type=str,
        help="Comma-separated Figure3c feature settings to compare.",
    )
    return parser.parse_args()


def parse_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(BOOL_TRUE)


def safe_divide(num: float, den: float) -> float:
    return float(num) / float(den) if den else 0.0


def roc_curve_points(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    order = np.argsort(-scores, kind="mergesort")
    labels_sorted = labels[order]
    positives = int(np.sum(labels_sorted == 1))
    negatives = int(np.sum(labels_sorted == 0))
    if positives == 0 or negatives == 0:
        raise ValueError("ROC requires both positive and negative labels.")
    tp = np.cumsum(labels_sorted == 1)
    fp = np.cumsum(labels_sorted == 0)
    tpr = np.concatenate(([0.0], tp / positives, [1.0]))
    fpr = np.concatenate(([0.0], fp / negatives, [1.0]))
    auc = float(np.trapz(tpr, fpr))
    return fpr, tpr, auc


def pr_curve_points(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    order = np.argsort(-scores, kind="mergesort")
    labels_sorted = labels[order]
    positives = int(np.sum(labels_sorted == 1))
    if positives == 0:
        raise ValueError("PR requires at least one positive label.")
    tp = np.cumsum(labels_sorted == 1)
    fp = np.cumsum(labels_sorted == 0)
    recall = tp / positives
    precision = tp / np.maximum(tp + fp, 1)
    recall_curve = np.concatenate(([0.0], recall))
    precision_curve = np.concatenate(([1.0], precision))
    auprc = float(np.sum(np.diff(recall_curve) * precision_curve[1:]))
    return recall_curve, precision_curve, auprc


def compute_score_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    _, _, auroc = roc_curve_points(labels, scores)
    _, _, auprc = pr_curve_points(labels, scores)
    return {"auroc": auroc, "auprc": auprc}


def compute_threshold_metrics(labels: np.ndarray, scores: np.ndarray, threshold: float) -> dict[str, float]:
    preds = (scores >= threshold).astype(int)
    tp = int(np.sum((preds == 1) & (labels == 1)))
    tn = int(np.sum((preds == 0) & (labels == 0)))
    fp = int(np.sum((preds == 1) & (labels == 0)))
    fn = int(np.sum((preds == 0) & (labels == 1)))
    precision = safe_divide(tp, tp + fp)
    recall = safe_divide(tp, tp + fn)
    specificity = safe_divide(tn, tn + fp)
    accuracy = safe_divide(tp + tn, tp + tn + fp + fn)
    f1 = safe_divide(2 * precision * recall, precision + recall)
    denom = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = safe_divide((tp * tn) - (fp * fn), np.sqrt(denom)) if denom > 0 else 0.0
    return {
        "mcc": mcc,
        "f1": f1,
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "accuracy": accuracy,
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }


def threshold_grid(scores: np.ndarray) -> np.ndarray:
    grid = np.concatenate((np.array([0.0, 0.5, 1.0]), np.unique(scores)))
    grid = grid[(grid >= 0.0) & (grid <= 1.0)]
    return np.unique(np.round(grid, 12))


def find_best_threshold(
    labels: np.ndarray,
    scores: np.ndarray,
    objective_name: str,
    objective_fn: Callable[[dict[str, float]], float],
) -> tuple[float, dict[str, float]]:
    best_threshold = 0.5
    best_metrics = compute_threshold_metrics(labels, scores, best_threshold)
    best_value = objective_fn(best_metrics)
    best_distance = abs(best_threshold - 0.5)

    for threshold in threshold_grid(scores):
        metrics = compute_threshold_metrics(labels, scores, float(threshold))
        value = objective_fn(metrics)
        distance = abs(float(threshold) - 0.5)
        if (
            value > best_value + 1e-12
            or (
                abs(value - best_value) <= 1e-12
                and (distance < best_distance - 1e-12 or (abs(distance - best_distance) <= 1e-12 and float(threshold) < best_threshold))
            )
        ):
            best_threshold = float(threshold)
            best_metrics = metrics
            best_value = value
            best_distance = distance
    best_metrics = dict(best_metrics)
    best_metrics[f"best_{objective_name}"] = best_value
    return best_threshold, best_metrics


def load_predictions(predictions_path: Path) -> pd.DataFrame:
    df = pd.read_csv(predictions_path, sep="\t")
    df["is_labeled_flag"] = parse_bool_series(df["is_labeled"])
    df["label"] = pd.to_numeric(df["label"], errors="coerce")
    df["pred_score"] = pd.to_numeric(df["pred_score"], errors="coerce")
    filtered = df.loc[df["is_labeled_flag"]].copy()
    filtered = filtered.loc[filtered["split"].isin(["val", "test"])].copy()
    filtered = filtered.dropna(subset=["label", "pred_score"])
    return filtered


def build_per_run_rows(figure_root: Path, feature_settings: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    split_versions_seen: set[str] = set()

    for feature_setting in feature_settings:
        setting_root = figure_root / feature_setting
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
                raise ValueError(f"Expected non-empty val/test labeled rows in {predictions_path}")

            val_labels = val_df["label"].to_numpy(dtype=int)
            val_scores = val_df["pred_score"].to_numpy(dtype=float)
            test_labels = test_df["label"].to_numpy(dtype=int)
            test_scores = test_df["pred_score"].to_numpy(dtype=float)

            base_score_metrics_val = compute_score_metrics(val_labels, val_scores)
            base_score_metrics_test = compute_score_metrics(test_labels, test_scores)
            best_threshold_by_f1, _ = find_best_threshold(val_labels, val_scores, "f1", lambda m: m["f1"])
            best_threshold_by_mcc, _ = find_best_threshold(val_labels, val_scores, "mcc", lambda m: m["mcc"])
            split_versions_seen.update(test_df["split_version"].dropna().astype(str).unique().tolist())

            for threshold_mode, threshold_source, fixed_threshold in THRESHOLD_MODES:
                threshold_value = fixed_threshold
                if threshold_mode == "val_f1_opt":
                    threshold_value = best_threshold_by_f1
                elif threshold_mode == "val_mcc_opt":
                    threshold_value = best_threshold_by_mcc

                val_threshold_metrics = compute_threshold_metrics(val_labels, val_scores, float(threshold_value))
                test_threshold_metrics = compute_threshold_metrics(test_labels, test_scores, float(threshold_value))

                row = dict(metrics_row)
                row.update(
                    {
                        "target": metrics_row.get("protocol", metrics_row.get("target")),
                        "Threshold_Mode": threshold_mode,
                        "Threshold_Source": threshold_source,
                        "Threshold_Value": float(threshold_value),
                        "best_threshold_by_f1": float(best_threshold_by_f1),
                        "best_threshold_by_mcc": float(best_threshold_by_mcc),
                        "threshold_strategy": threshold_mode,
                        "evaluation_contract": f"val_threshold_tuned_{threshold_mode}",
                        "metrics_path": str(metrics_path.resolve()),
                        "predictions_path": str(predictions_path.resolve()),
                    }
                )
                row.update(
                    {
                        "val_auroc": base_score_metrics_val["auroc"],
                        "val_auprc": base_score_metrics_val["auprc"],
                        "val_mcc": val_threshold_metrics["mcc"],
                        "val_f1": val_threshold_metrics["f1"],
                        "val_precision": val_threshold_metrics["precision"],
                        "val_recall": val_threshold_metrics["recall"],
                        "val_specificity": val_threshold_metrics["specificity"],
                        "val_accuracy": val_threshold_metrics["accuracy"],
                        "test_auroc": base_score_metrics_test["auroc"],
                        "test_auprc": base_score_metrics_test["auprc"],
                        "test_mcc": test_threshold_metrics["mcc"],
                        "test_f1": test_threshold_metrics["f1"],
                        "test_precision": test_threshold_metrics["precision"],
                        "test_recall": test_threshold_metrics["recall"],
                        "test_specificity": test_threshold_metrics["specificity"],
                        "test_accuracy": test_threshold_metrics["accuracy"],
                        "val_count": int(len(val_df)),
                        "test_count": int(len(test_df)),
                    }
                )
                rows.append(row)

    if len(split_versions_seen) != 1:
        raise ValueError(f"Expected exactly one split_version across Figure3c runs, found: {sorted(split_versions_seen)}")
    return pd.DataFrame(rows)


def aggregate_per_run(per_run: pd.DataFrame) -> pd.DataFrame:
    group_columns = [
        "target",
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
        "Threshold_Mode",
        "Threshold_Source",
        "label_manifest",
        "split_manifest",
        "config_used",
        "is_deterministic",
        "require_true_node2vec",
        "embedding_method",
        "embedding_backend",
        "fallback_used",
        "esm2_dim",
        "feature_combo",
        "fusion_mode",
        "fusion_hidden_dim",
        "fusion_dropout",
    ]
    metric_bases = [
        "Threshold_Value",
        "best_threshold_by_f1",
        "best_threshold_by_mcc",
        "test_auroc",
        "test_auprc",
        "test_mcc",
        "test_f1",
        "test_precision",
        "test_recall",
        "test_specificity",
        "test_accuracy",
        "val_auroc",
        "val_auprc",
        "val_mcc",
        "val_f1",
        "val_precision",
        "val_recall",
        "val_specificity",
        "val_accuracy",
    ]

    rows: list[dict[str, object]] = []
    for column in group_columns:
        if column not in per_run.columns:
            per_run[column] = pd.NA
    for group_key, group_df in per_run.groupby(group_columns, dropna=False, sort=True):
        row = {column: value for column, value in zip(group_columns, group_key)}
        seeds = pd.to_numeric(group_df["seed"], errors="coerce").dropna().astype(int).tolist()
        run_ids = group_df["run_id"].astype(str).tolist()
        row["Runs"] = int(len(group_df))
        row["Seed_List"] = ",".join(str(seed) for seed in sorted(seeds))
        row["Run_IDs"] = ",".join(run_ids)
        for metric in metric_bases:
            values = pd.to_numeric(group_df[metric], errors="coerce")
            row[f"{metric}_mean"] = float(values.mean())
            row[f"{metric}_std"] = float(values.std(ddof=0))
        rows.append(row)

    aggregated = pd.DataFrame(rows)
    return aggregated.sort_values(["feature_setting", "Threshold_Mode"], kind="stable").reset_index(drop=True)


def build_final_summary(aggregated: pd.DataFrame) -> pd.DataFrame:
    keep_columns = {
        "target": "Target",
        "species": "Species",
        "regime": "Regime",
        "model": "Model",
        "feature_setting": "Feature_Setting",
        "label_regime": "Label_Regime",
        "split_version": "Split_Version",
        "Runs": "Runs",
        "Seed_List": "Seed_List",
        "esm2_dim": "ESM2_Dim",
        "Threshold_Mode": "Threshold_Mode",
        "Threshold_Source": "Threshold_Source",
        "Threshold_Value_mean": "Threshold_Value_Mean",
        "Threshold_Value_std": "Threshold_Value_Std",
        "best_threshold_by_f1_mean": "Best_Threshold_By_F1_Mean",
        "best_threshold_by_f1_std": "Best_Threshold_By_F1_Std",
        "best_threshold_by_mcc_mean": "Best_Threshold_By_MCC_Mean",
        "best_threshold_by_mcc_std": "Best_Threshold_By_MCC_Std",
    }
    for source_name, display_name in FINAL_METRICS:
        keep_columns[source_name] = display_name

    selected = aggregated[[column for column in keep_columns if column in aggregated.columns]].rename(columns=keep_columns)
    return selected.sort_values(["Feature_Setting", "Threshold_Mode"], kind="stable").reset_index(drop=True)


def final_summary_markdown(final_summary: pd.DataFrame) -> str:
    lines = [
        "# Figure3c Threshold Tuned Final Summary",
        "",
        "Validation-tuned Figure3c comparison for GraphSAGE baseline vs gated fusion on fgraminearum newlabel.",
        "",
    ]
    if final_summary.empty:
        lines.extend(["No rows available.", ""])
        return "\n".join(lines)

    display = final_summary.copy()
    metric_pairs = [
        ("AUROC_mean", "AUROC_std", "AUROC"),
        ("AUPRC_mean", "AUPRC_std", "AUPRC"),
        ("MCC_mean", "MCC_std", "MCC"),
        ("F1_mean", "F1_std", "F1"),
        ("Precision_mean", "Precision_std", "Precision"),
        ("Recall_mean", "Recall_std", "Recall"),
        ("Specificity_mean", "Specificity_std", "Specificity"),
    ]
    for mean_col, std_col, label in metric_pairs:
        if mean_col in display.columns and std_col in display.columns:
            display[label] = display.apply(lambda row: f"{row[mean_col]:.3f} ± {row[std_col]:.3f}", axis=1)
    drop_cols = [column for pair in metric_pairs for column in pair[:2] if column in display.columns]
    display = display.drop(columns=drop_cols)
    lines.append(display.to_markdown(index=False))
    lines.append("")
    return "\n".join(lines)


def write_outputs(summary_dir: Path, prefix: str, per_run: pd.DataFrame, aggregated: pd.DataFrame, final_summary: pd.DataFrame) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)
    per_run_tsv = summary_dir / f"{prefix}_per_run_metrics.tsv"
    per_run_md = summary_dir / f"{prefix}_per_run_metrics.md"
    aggregated_tsv = summary_dir / f"{prefix}_aggregated_metrics.tsv"
    aggregated_md = summary_dir / f"{prefix}_aggregated_metrics.md"
    final_tsv = summary_dir / f"{prefix}_final_summary.tsv"
    final_md = summary_dir / f"{prefix}_final_summary.md"

    per_run.to_csv(per_run_tsv, sep="\t", index=False)
    aggregated.to_csv(aggregated_tsv, sep="\t", index=False)
    final_summary.to_csv(final_tsv, sep="\t", index=False)

    per_run_md.write_text(
        dataframe_to_markdown(per_run, f"{prefix} Per-Run Metrics", "Per-seed threshold-tuned Figure3c metrics computed from predictions.tsv."),
        encoding="utf-8",
    )
    aggregated_md.write_text(
        dataframe_to_markdown(aggregated, f"{prefix} Aggregated Metrics", "Mean ± std aggregates grouped by feature setting and threshold mode."),
        encoding="utf-8",
    )
    final_md.write_text(final_summary_markdown(final_summary), encoding="utf-8")


def main() -> None:
    args = parse_args()
    figure_root = Path(args.figure3c_root).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    feature_settings = [item.strip() for item in args.feature_settings.split(",") if item.strip()]
    per_run = build_per_run_rows(figure_root, feature_settings)
    aggregated = aggregate_per_run(per_run)
    final_summary = build_final_summary(aggregated)
    write_outputs(summary_dir, "Figure3c_threshold_tuned", per_run, aggregated, final_summary)


if __name__ == "__main__":
    main()
