import argparse
from pathlib import Path

import pandas as pd
import yaml

from src.train.run_frozen_protocol_model import compute_binary_metrics, is_gated_feature_setting, resolve_feature_contract_group


def parse_args():
    parser = argparse.ArgumentParser(description="Repair/upgrade a run metrics.tsv using predictions.tsv")
    parser.add_argument("--run-dir", required=True, type=str)
    return parser.parse_args()


def _subset_metrics(predictions: pd.DataFrame, split_name: str):
    subset = predictions[predictions["split"].astype(str) == split_name].copy()
    y_true = pd.to_numeric(subset["label"], errors="coerce").astype(int).to_numpy()
    y_score = pd.to_numeric(subset["pred_score"], errors="coerce").astype(float).to_numpy()
    y_pred = pd.to_numeric(subset["pred_label"], errors="coerce").astype(int).to_numpy()
    return compute_binary_metrics(y_true, y_score, y_pred)


def main():
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    metrics_path = run_dir / "metrics.tsv"
    predictions_path = run_dir / "predictions.tsv"
    resolved_path = run_dir / "resolved_config.yaml"
    esm2_metadata_path = run_dir / "esm2_feature_metadata.tsv"

    metrics_row = pd.read_csv(metrics_path, sep="\t").iloc[0].to_dict()
    predictions = pd.read_csv(predictions_path, sep="\t")
    val_metrics = _subset_metrics(predictions, "val")
    test_metrics = _subset_metrics(predictions, "test")

    feature_setting = str(metrics_row.get("feature_setting", "")).strip().upper()
    model = str(metrics_row.get("model", "")).strip()

    metrics_row.update(
        {
            "model_variant": metrics_row.get("model_variant", model),
            "graph_contract": metrics_row.get("graph_contract", "undirected_symmetrized"),
            "threshold_strategy": metrics_row.get("threshold_strategy", "fixed_0.5"),
            "evaluation_contract": metrics_row.get("evaluation_contract", "auroc_auprc_mcc_specificity_fixed_0.5"),
            "feature_contract_group": metrics_row.get("feature_contract_group", resolve_feature_contract_group(model, feature_setting)),
            "val_precision": val_metrics["precision"],
            "val_recall": val_metrics["recall"],
            "val_specificity": val_metrics["specificity"],
            "test_precision": test_metrics["precision"],
            "test_recall": test_metrics["recall"],
            "test_specificity": test_metrics["specificity"],
            "fusion_mode": metrics_row.get("fusion_mode", "") if is_gated_feature_setting(feature_setting) else "",
            "fusion_hidden_dim": metrics_row.get("fusion_hidden_dim", "") if is_gated_feature_setting(feature_setting) else "",
            "fusion_dropout": metrics_row.get("fusion_dropout", "") if is_gated_feature_setting(feature_setting) else "",
            "best_checkpoint": str(run_dir / "best_model.pt"),
        }
    )
    if esm2_metadata_path.exists():
        esm2_meta = pd.read_csv(esm2_metadata_path, sep="\t").iloc[0].to_dict()
        metrics_row["esm2_cache_path"] = str(esm2_meta.get("cache_path", metrics_row.get("esm2_cache_path", "")))

    pd.DataFrame([metrics_row]).to_csv(metrics_path, sep="\t", index=False)

    if resolved_path.exists():
        resolved = yaml.safe_load(resolved_path.read_text(encoding="utf-8"))
        resolved["output_dir"] = str(run_dir)
        resolved["metrics"] = metrics_row
        resolved_path.write_text(yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8")


if __name__ == "__main__":
    main()
