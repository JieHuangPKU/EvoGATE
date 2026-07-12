from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


OUTPUT_ROOT = Path("results/FigureS4C_threshold_sweep")
TABLE_DIR = OUTPUT_ROOT / "tables"
PLOT_DIR = OUTPUT_ROOT / "plots"
SUMMARY_DIR = OUTPUT_ROOT / "summary"

THRESHOLDS = np.round(np.linspace(0.0, 1.0, 101), 2)
BOOL_TRUE = {"true", "t", "1", "yes", "y"}

VARIANT_SPECS = [
    {
        "fusion_method": "Concat",
        "feature_setting": "ORT_EXP_SUB_ESM2",
        "search_roots": [
            Path("outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE"),
            Path("outputs/Figure3c/fgraminearum_newlabel/GraphSAGE"),
            Path("outputs/Figure3cA/fgraminearum_newlabel/GraphSAGE"),
        ],
    },
    {
        "fusion_method": "Gated",
        "feature_setting": "ORT_EXP_SUB_ESM2_GATED",
        "search_roots": [
            Path("outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE"),
            Path("outputs/Figure3c/fgraminearum_newlabel/GraphSAGE"),
            Path("outputs/Figure3cA/fgraminearum_newlabel/GraphSAGE"),
        ],
    },
    {
        "fusion_method": "Gated+WBCE",
        "feature_setting": "ORT_EXP_SUB_ESM2_OLD_GATED_WBCE",
        "search_roots": [
            Path("outputs/Figure3cC/fgraminearum_newlabel/GraphSAGE"),
        ],
    },
    {
        "fusion_method": "Residual gated",
        "feature_setting": "ORT_EXP_SUB_ESM2_GATED_RESIDUAL",
        "search_roots": [
            Path("outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE"),
            Path("outputs/Figure3cA/fgraminearum_newlabel/GraphSAGE"),
        ],
    },
    {
        "fusion_method": "Residual gated+WBCE",
        "feature_setting": "ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE",
        "search_roots": [
            Path("outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE"),
        ],
    },
]

COLUMN_ALIASES = {
    "gene_id": [
        "gene_id",
        "canonical_gene_id",
        "graph_gene_id",
        "node_id",
        "protein_id",
        "gene",
    ],
    "label": [
        "y_true",
        "label",
        "true_label",
        "target",
        "gold_label",
    ],
    "score": [
        "y_score",
        "pred_score",
        "prob",
        "probability",
        "sigmoid_score",
        "essentiality_score",
        "score",
    ],
    "split": [
        "split",
        "train_val_test",
        "dataset_split",
        "partition",
    ],
    "is_labeled": [
        "is_labeled",
        "labeled",
        "is_labelled",
    ],
    "seed": [
        "seed",
    ],
    "run_id": [
        "run_id",
        "run",
    ],
    "feature_setting": [
        "feature_setting",
        "feature_combo",
        "model_variant",
    ],
}


def ensure_dirs() -> None:
    for path in (TABLE_DIR, PLOT_DIR, SUMMARY_DIR):
        path.mkdir(parents=True, exist_ok=True)


def detect_column(columns: list[str], aliases: list[str], required: bool = True) -> str | None:
    lowered = {col.lower(): col for col in columns}
    for alias in aliases:
        if alias.lower() in lowered:
            return lowered[alias.lower()]
    if required:
        raise KeyError(f"Could not detect required column from aliases: {aliases}")
    return None


def parse_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(BOOL_TRUE)


def safe_div(num: float, den: float) -> float:
    return float(num) / float(den) if den else 0.0


def compute_confusion_metrics(labels: np.ndarray, scores: np.ndarray, threshold: float) -> dict[str, float]:
    preds = (scores >= threshold).astype(int)
    tp = int(np.sum((preds == 1) & (labels == 1)))
    fp = int(np.sum((preds == 1) & (labels == 0)))
    tn = int(np.sum((preds == 0) & (labels == 0)))
    fn = int(np.sum((preds == 0) & (labels == 1)))
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    f1 = safe_div(2 * precision * recall, precision + recall)
    denom = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = safe_div((tp * tn) - (fp * fn), math.sqrt(denom)) if denom > 0 else 0.0
    return {
        "MCC": mcc,
        "F1": f1,
        "precision": precision,
        "recall": recall,
        "TP": tp,
        "FP": fp,
        "TN": tn,
        "FN": fn,
    }


def summarize_metric(values: pd.Series) -> tuple[float, float, int]:
    arr = values.astype(float).to_numpy()
    n = int(arr.size)
    mean = float(np.mean(arr)) if n else float("nan")
    sd = float(np.std(arr, ddof=1)) if n > 1 else 0.0
    return mean, sd, n


def load_prediction_frame(predictions_path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    df = pd.read_csv(predictions_path, sep="\t")
    cols = list(df.columns)
    detected = {
        "gene_id": detect_column(cols, COLUMN_ALIASES["gene_id"]),
        "label": detect_column(cols, COLUMN_ALIASES["label"]),
        "score": detect_column(cols, COLUMN_ALIASES["score"]),
        "split": detect_column(cols, COLUMN_ALIASES["split"], required=False),
        "is_labeled": detect_column(cols, COLUMN_ALIASES["is_labeled"], required=False),
        "seed": detect_column(cols, COLUMN_ALIASES["seed"], required=False),
        "run_id": detect_column(cols, COLUMN_ALIASES["run_id"], required=False),
        "feature_setting": detect_column(cols, COLUMN_ALIASES["feature_setting"], required=False),
    }

    working = df.copy()
    working["_gene_id"] = working[detected["gene_id"]].astype(str)
    working["_label"] = pd.to_numeric(working[detected["label"]], errors="coerce")
    working["_score"] = pd.to_numeric(working[detected["score"]], errors="coerce")
    if detected["split"] is not None:
        working["_split"] = working[detected["split"]].astype(str).str.strip().str.lower()
    else:
        working["_split"] = "all"
    if detected["is_labeled"] is not None:
        working["_is_labeled"] = parse_bool(working[detected["is_labeled"]])
    else:
        working["_is_labeled"] = working["_label"].notna()

    working = working.loc[working["_is_labeled"]].copy()
    working = working.dropna(subset=["_label", "_score"])
    working["_label"] = working["_label"].astype(int)
    return working, detected


def locate_run_files(feature_setting: str, roots: list[Path]) -> list[dict[str, object]]:
    for root in roots:
        runs: list[dict[str, object]] = []
        run_dirs = sorted((root / feature_setting).glob("run_*"))
        for run_dir in run_dirs:
            predictions_path = run_dir / "predictions.tsv"
            metrics_path = run_dir / "metrics.tsv"
            if not predictions_path.exists() or not metrics_path.exists():
                continue
            runs.append(
                {
                    "run_dir": run_dir,
                    "predictions_path": predictions_path.resolve(),
                    "metrics_path": metrics_path.resolve(),
                }
            )
        if runs:
            return runs
    if not roots:
        raise FileNotFoundError(f"No prediction runs found for feature setting {feature_setting}")
    raise FileNotFoundError(f"No prediction runs found for feature setting {feature_setting}")


def write_input_audit(audit_rows: list[dict[str, object]], val_available: bool) -> None:
    path = SUMMARY_DIR / "input_audit.md"
    lines = [
        "# Figure S4C threshold sweep input audit",
        "",
        f"Validation split available across all included runs: {'yes' if val_available else 'no'}",
        "",
        "| fusion_method | internal_feature_setting | run_id | seed | predictions_path | metrics_path | gene_id_col | label_col | score_col | split_col | is_labeled_col | val_rows | test_rows | labeled_rows |",
        "|:--|:--|:--|--:|:--|:--|:--|:--|:--|:--|:--|--:|--:|--:|",
    ]
    for row in audit_rows:
        lines.append(
            "| {fusion_method} | {feature_setting} | {run_id} | {seed} | {predictions_path} | {metrics_path} | {gene_id_col} | {label_col} | {score_col} | {split_col} | {is_labeled_col} | {val_rows} | {test_rows} | {labeled_rows} |".format(
                **row
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary(
    file_rows: list[dict[str, object]],
    split_used_for_sweep: str,
    optimization_split: str,
    warning: str | None,
    optima_df: pd.DataFrame,
    genes_per_split: dict[str, int],
    runs_per_variant: dict[str, int],
) -> None:
    lines = [
        "# Figure S4C threshold sweep summary",
        "",
        "Decision metrics such as MCC and F1 vary with the classification threshold applied to the model output score, whereas AUROC and AUPRC remain threshold-invariant ranking metrics.",
        "",
        "## Prediction files used",
        "",
    ]
    for row in file_rows:
        lines.append(f"- `{row['fusion_method']}`: `{row['predictions_path']}`")
    lines.extend(
        [
            "",
            f"## Split used for threshold optimization",
            "",
            f"- `{optimization_split}`",
            "",
            "## Split used for threshold sweep curves",
            "",
            f"- `{split_used_for_sweep}`",
            "",
            "## Fusion variants included",
            "",
        ]
    )
    for fusion_method in runs_per_variant:
        lines.append(f"- `{fusion_method}`")
    lines.extend(
        [
            "",
            "## Number of genes",
            "",
        ]
    )
    for split_name, count in genes_per_split.items():
        lines.append(f"- `{split_name}`: {count}")
    lines.extend(
        [
            "",
            "## Number of seeds/runs",
            "",
        ]
    )
    for fusion_method, count in runs_per_variant.items():
        lines.append(f"- `{fusion_method}`: {count}")

    f1_rows = optima_df.loc[optima_df["threshold_type"] == "val_f1_opt"].copy()
    mcc_rows = optima_df.loc[optima_df["threshold_type"] == "val_mcc_opt"].copy()
    lines.extend(
        [
            "",
            "## Optimal F1 thresholds per fusion method",
            "",
            "| fusion_method | threshold | metric_value |",
            "|:--|--:|--:|",
        ]
    )
    for _, row in f1_rows.iterrows():
        lines.append(f"| {row['fusion_method']} | {row['threshold']:.2f} | {row['metric_value']:.6f} |")
    lines.extend(
        [
            "",
            "## Optimal MCC thresholds per fusion method",
            "",
            "| fusion_method | threshold | metric_value |",
            "|:--|--:|--:|",
        ]
    )
    for _, row in mcc_rows.iterrows():
        lines.append(f"| {row['fusion_method']} | {row['threshold']:.2f} | {row['metric_value']:.6f} |")

    if warning:
        lines.extend(
            [
                "",
                "## Warning",
                "",
                f"- {warning}",
            ]
        )

    (SUMMARY_DIR / "FigureS4C_threshold_sweep_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    all_run_rows: list[dict[str, object]] = []
    sweep_rows: list[dict[str, object]] = []
    optima_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []

    validation_available_for_all = True
    split_used_for_sweep = "test"
    optimization_split = "validation"
    warning: str | None = None
    genes_per_split: dict[str, int] = {}
    runs_per_variant: dict[str, int] = {}

    for spec in VARIANT_SPECS:
        runs = locate_run_files(spec["feature_setting"], spec["search_roots"])
        runs_per_variant[spec["fusion_method"]] = len(runs)

        validation_metric_rows: list[dict[str, object]] = []

        for run in runs:
            metrics_df = pd.read_csv(run["metrics_path"], sep="\t")
            metrics_row = metrics_df.iloc[0].to_dict()
            frame, detected = load_prediction_frame(Path(run["predictions_path"]))

            run_id = str(metrics_row.get("run_id", Path(run["run_dir"]).name))
            seed_value = metrics_row.get("seed", pd.NA)

            split_counts = frame["_split"].value_counts().to_dict()
            val_df = frame.loc[frame["_split"] == "val"].copy()
            test_df = frame.loc[frame["_split"] == "test"].copy()
            use_validation = not val_df.empty
            eval_df = test_df if not test_df.empty else val_df
            if eval_df.empty:
                raise ValueError(f"No labeled val/test rows available in {run['predictions_path']}")
            if not use_validation:
                validation_available_for_all = False

            labels_eval = eval_df["_label"].to_numpy(dtype=int)
            scores_eval = eval_df["_score"].to_numpy(dtype=float)
            labels_opt = (val_df if use_validation else eval_df)["_label"].to_numpy(dtype=int)
            scores_opt = (val_df if use_validation else eval_df)["_score"].to_numpy(dtype=float)

            if use_validation:
                genes_per_split.setdefault("validation", int(val_df["_gene_id"].nunique()))
            if not test_df.empty:
                genes_per_split.setdefault("test", int(test_df["_gene_id"].nunique()))
            if not use_validation and "available_split" not in genes_per_split:
                genes_per_split["available_split"] = int(eval_df["_gene_id"].nunique())

            audit_rows.append(
                {
                    "fusion_method": spec["fusion_method"],
                    "feature_setting": spec["feature_setting"],
                    "run_id": run_id,
                    "seed": seed_value,
                    "predictions_path": run["predictions_path"],
                    "metrics_path": run["metrics_path"],
                    "gene_id_col": detected["gene_id"],
                    "label_col": detected["label"],
                    "score_col": detected["score"],
                    "split_col": detected["split"] or "not_present",
                    "is_labeled_col": detected["is_labeled"] or "not_present",
                    "val_rows": int(len(val_df)),
                    "test_rows": int(len(test_df)),
                    "labeled_rows": int(len(frame)),
                }
            )

            for threshold in THRESHOLDS:
                metrics = compute_confusion_metrics(labels_eval, scores_eval, float(threshold))
                for metric_name in ("MCC", "F1", "precision", "recall"):
                    sweep_rows.append(
                        {
                            "fusion_method": spec["fusion_method"],
                            "run_id": run_id,
                            "seed": seed_value,
                            "threshold": float(threshold),
                            "metric": metric_name,
                            "value": float(metrics[metric_name]),
                            "split_used": "test" if not test_df.empty else "available_split",
                        }
                    )

                validation_metric_rows.append(
                    {
                        "fusion_method": spec["fusion_method"],
                        "run_id": run_id,
                        "seed": seed_value,
                        "threshold": float(threshold),
                        "MCC": float(compute_confusion_metrics(labels_opt, scores_opt, float(threshold))["MCC"]),
                        "F1": float(compute_confusion_metrics(labels_opt, scores_opt, float(threshold))["F1"]),
                    }
                )

            all_run_rows.append(
                {
                    "fusion_method": spec["fusion_method"],
                    "feature_setting": spec["feature_setting"],
                    "run_id": run_id,
                    "seed": seed_value,
                    "predictions_path": str(run["predictions_path"]),
                    "metrics_path": str(run["metrics_path"]),
                }
            )

        val_metric_df = pd.DataFrame(validation_metric_rows)
        summary_by_threshold = (
            val_metric_df.groupby("threshold", as_index=False)[["MCC", "F1"]]
            .mean()
            .sort_values("threshold")
            .reset_index(drop=True)
        )

        for threshold_type, metric_name in (("val_f1_opt", "F1"), ("val_mcc_opt", "MCC")):
            metric_series = summary_by_threshold[metric_name]
            best_value = float(metric_series.max())
            tied = summary_by_threshold.loc[np.isclose(metric_series.to_numpy(), best_value)]
            best_threshold = float(tied["threshold"].min())
            optima_rows.append(
                {
                    "fusion_method": spec["fusion_method"],
                    "threshold_type": threshold_type,
                    "threshold": best_threshold,
                    "metric_optimized": metric_name,
                    "metric_value": best_value,
                    "split_used": "validation" if validation_available_for_all else "available_split",
                }
            )

        optima_rows.append(
            {
                "fusion_method": spec["fusion_method"],
                "threshold_type": "fixed_0.5",
                "threshold": 0.5,
                "metric_optimized": "fixed",
                "metric_value": float("nan"),
                "split_used": "test",
            }
        )

    if not validation_available_for_all:
        optimization_split = "available_split"
        split_used_for_sweep = "available_split"
        warning = (
            "Validation split was unavailable for at least one included run. "
            "Threshold optima are descriptive on the available evaluation split, not validation-selected."
        )
    else:
        warning = None

    write_input_audit(audit_rows, validation_available_for_all)

    sweep_df = pd.DataFrame(sweep_rows)
    aggregated_rows: list[dict[str, object]] = []
    for (fusion_method, threshold, metric), group in sweep_df.groupby(["fusion_method", "threshold", "metric"], sort=True):
        mean, sd, n = summarize_metric(group["value"])
        aggregated_rows.append(
            {
                "fusion_method": fusion_method,
                "threshold": float(threshold),
                "metric": metric,
                "mean": mean,
                "sd": sd,
                "n": n,
                "split_used": group["split_used"].iloc[0],
            }
        )
    aggregated_df = pd.DataFrame(aggregated_rows).sort_values(["metric", "fusion_method", "threshold"]).reset_index(drop=True)
    aggregated_df.to_csv(TABLE_DIR / "threshold_sweep_metrics.tsv", sep="\t", index=False)

    optima_df = pd.DataFrame(optima_rows).sort_values(["fusion_method", "threshold_type"]).reset_index(drop=True)
    optima_df.to_csv(TABLE_DIR / "threshold_optima.tsv", sep="\t", index=False)

    write_summary(
        file_rows=all_run_rows,
        split_used_for_sweep=split_used_for_sweep,
        optimization_split=optimization_split,
        warning=warning,
        optima_df=optima_df,
        genes_per_split=genes_per_split,
        runs_per_variant=runs_per_variant,
    )

    print(str((TABLE_DIR / "threshold_sweep_metrics.tsv").resolve()))
    print(str((TABLE_DIR / "threshold_optima.tsv").resolve()))
    print(str((SUMMARY_DIR / "input_audit.md").resolve()))
    print(str((SUMMARY_DIR / "FigureS4C_threshold_sweep_summary.md").resolve()))


if __name__ == "__main__":
    main()
