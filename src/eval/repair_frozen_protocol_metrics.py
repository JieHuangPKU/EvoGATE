import argparse
from pathlib import Path

import pandas as pd

from src.train.run_frozen_protocol_model import compute_binary_metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Repair or complete frozen-protocol metrics.tsv columns for one run directory")
    parser.add_argument("--run-dir", required=True, type=str)
    return parser.parse_args()


def _compute_split_metrics(predictions: pd.DataFrame, split_name: str) -> dict:
    split_df = predictions[predictions["split"].astype(str) == split_name].copy()
    split_df = split_df[pd.notna(split_df["label"])].copy()
    if split_df.empty:
        return compute_binary_metrics([], [], [])
    y_true = split_df["label"].astype(int).to_numpy()
    y_score = pd.to_numeric(split_df["pred_score"], errors="coerce").to_numpy()
    y_pred = pd.to_numeric(split_df["pred_label"], errors="coerce").fillna(0).astype(int).to_numpy()
    return compute_binary_metrics(y_true, y_score, y_pred)


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    metrics_path = run_dir / "metrics.tsv"
    predictions_path = run_dir / "predictions.tsv"
    esm2_metadata_path = run_dir / "esm2_feature_metadata.tsv"
    if not metrics_path.exists():
        raise FileNotFoundError(f"Missing metrics file: {metrics_path}")
    if not predictions_path.exists():
        raise FileNotFoundError(f"Missing predictions file: {predictions_path}")

    metrics_df = pd.read_csv(metrics_path, sep="\t")
    if metrics_df.empty:
        raise ValueError(f"Empty metrics file: {metrics_path}")
    metrics_row = metrics_df.iloc[0].to_dict()

    predictions = pd.read_csv(predictions_path, sep="\t")
    val_metrics = _compute_split_metrics(predictions, "val")
    test_metrics = _compute_split_metrics(predictions, "test")

    for split_name, split_metrics in [("val", val_metrics), ("test", test_metrics)]:
        for metric_name, metric_value in split_metrics.items():
            metrics_row[f"{split_name}_{metric_name}"] = metric_value

    if "target" not in metrics_row or pd.isna(metrics_row["target"]):
        metrics_row["target"] = metrics_row.get("species")
    if ("esm2_dim" not in metrics_row or pd.isna(metrics_row["esm2_dim"]) or str(metrics_row["esm2_dim"]).strip() == "") and esm2_metadata_path.exists():
        esm2_metadata = pd.read_csv(esm2_metadata_path, sep="\t")
        if not esm2_metadata.empty and "embedding_dim" in esm2_metadata.columns:
            metrics_row["esm2_dim"] = esm2_metadata.iloc[0]["embedding_dim"]

    pd.DataFrame([metrics_row]).to_csv(metrics_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
