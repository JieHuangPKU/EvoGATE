import argparse
import pickle
from pathlib import Path

import pandas as pd
import torch
import yaml

from src.data.frozen_protocol_loader import load_protocol_dataset
from src.train.run_frozen_protocol_feature_combo_model import resolve_graph_model_config, save_pickle
from src.train.run_frozen_protocol_model import (
    GRAPH_MODELS,
    build_prediction_table,
    compute_binary_metrics,
    is_gated_feature_setting,
    normalize_model_name,
    resolve_feature_contract_group,
    save_yaml,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Materialize frozen-protocol GraphSAGE outputs from an existing best_model.pt checkpoint")
    parser.add_argument("--config", required=True, type=str)
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--model", required=True, type=str)
    parser.add_argument("--feature-setting", required=True, type=str)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=str)
    parser.add_argument("--graph-source", default=None, type=str)
    parser.add_argument("--graph-source-name", default=None, type=str)
    parser.add_argument("--string-threshold", default=None, type=float)
    return parser.parse_args()


def build_graph_model_for_inference(bundle, model_name, model_cfg):
    from src.models.epgat_gcn import EPGATOriginalGCN
    from src.models.epgat_gin import EPGATOriginalGIN
    from src.models.epgat_original import EPGATOriginalGAT
    from src.models.epgat_sage import EPGATOriginalSAGE, EPGATOriginalSAGEWithFusion

    x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
    normalized_model = normalize_model_name(model_name)
    if normalized_model == "GAT":
        return EPGATOriginalGAT(
            in_feats=x.shape[1],
            h_feats=[12, 1] if bundle["species"] != "human" else [16, 1],
            heads=[8, 1],
            dropout=float(model_cfg["dropout"]),
            negative_slope=0.2,
        )
    if normalized_model == "GCN":
        return EPGATOriginalGCN(
            in_feats=x.shape[1],
            h_layers=[int(model_cfg["hidden_dim"]), 1],
            dropout=float(model_cfg["dropout"]),
        )
    if normalized_model == "GIN":
        return EPGATOriginalGIN(
            in_feats=x.shape[1],
            dim_h=int(model_cfg["dim_h"]),
            dropout=float(model_cfg["dropout"]),
        )
    if normalized_model == "GraphSAGE":
        if is_gated_feature_setting(bundle["feature_setting"]):
            fusion_partition = dict(bundle.get("feature_metadata", {}).get("fusion_partition", {}))
            return EPGATOriginalSAGEWithFusion(
                in_feats=x.shape[1],
                fusion_partition=fusion_partition,
                fusion_hidden_dim=int(model_cfg.get("fusion_hidden_dim", 256)),
                fusion_dropout=float(model_cfg.get("fusion_dropout", 0.2)),
                fusion_mode=str(model_cfg.get("fusion_mode", "gated")),
                n_hidden=int(model_cfg["n_hidden"]),
                n_layers=int(model_cfg["n_layers"]),
                dropout=float(model_cfg["dropout"]),
                aggregator_type=str(model_cfg["aggregator_type"]),
            )
        return EPGATOriginalSAGE(
            in_feats=x.shape[1],
            n_hidden=int(model_cfg["n_hidden"]),
            n_layers=int(model_cfg["n_layers"]),
            dropout=float(model_cfg["dropout"]),
            aggregator_type=str(model_cfg["aggregator_type"]),
        )
    raise ValueError("Unsupported graph model: {}".format(model_name))


def main():
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output_dir / "best_model.pt"
    if not checkpoint_path.exists():
        raise FileNotFoundError("Checkpoint not found: {}".format(checkpoint_path))

    base_config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    model_key = normalize_model_name(args.model)
    if model_key not in GRAPH_MODELS:
        raise ValueError("Only graph models are supported, got '{}'".format(model_key))
    model_cfg = resolve_graph_model_config(base_config, model_key)
    feature_setting = str(args.feature_setting).strip().upper()

    bundle = load_protocol_dataset(
        args.config,
        args.protocol,
        feature_setting,
        graph_source_path=args.graph_source,
        string_threshold=args.string_threshold,
        graph_source_name=args.graph_source_name,
    )
    model = build_graph_model_for_inference(bundle, model_key, model_cfg)
    checkpoint = torch.load(str(checkpoint_path), map_location="cpu")
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()

    x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
    edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
    with torch.no_grad():
        pred_score = torch.sigmoid(model(x, edge_index).view(-1)).cpu().numpy()
    pred_label = (pred_score >= 0.5).astype(int)

    test_idx = bundle["test_idx"]
    val_idx = bundle["val_idx"]
    test_metrics = compute_binary_metrics(bundle["y_all"][test_idx].astype(int), pred_score[test_idx], pred_label[test_idx])
    val_metrics = compute_binary_metrics(bundle["y_all"][val_idx].astype(int), pred_score[val_idx], pred_label[val_idx])

    prediction_table = build_prediction_table(bundle, pred_score, pred_label, model_key)
    prediction_table["feature_combo"] = feature_setting
    prediction_table.to_csv(output_dir / "predictions.tsv", sep="\t", index=False)
    bundle["feature_schema"].to_csv(output_dir / "feature_schema.tsv", sep="\t", index=False)
    pd.DataFrame(bundle["edge_table"]).to_csv(output_dir / "edge_table.tsv", sep="\t", index=False)
    bundle["split_manifest"].to_csv(output_dir / "split_manifest.tsv", sep="\t", index=False)
    feature_metadata = dict(bundle.get("feature_metadata", {}))
    graph_metadata = dict(bundle.get("graph_metadata", {}))
    if "esm2_alignment_audit" in feature_metadata:
        feature_metadata["esm2_alignment_audit"].to_csv(output_dir / "esm2_alignment_audit.tsv", sep="\t", index=False)
    if "esm2_metadata" in feature_metadata:
        pd.DataFrame([feature_metadata["esm2_metadata"]]).to_csv(output_dir / "esm2_feature_metadata.tsv", sep="\t", index=False)

    metrics_row = {
        "protocol": bundle["protocol_name"],
        "species": bundle["species"],
        "regime": bundle["regime"],
        "model": model_key,
        "model_variant": args.model,
        "feature_setting": feature_setting,
        "feature_combo": feature_setting,
        "label_regime": bundle["label_regime"],
        "run_id": "seed_{}".format(int(args.seed)),
        "seed": int(args.seed),
        "is_deterministic": "false",
        "split_version": bundle["split_version"],
        "graph_source": bundle["graph_source"],
        "graph_source_name": bundle.get("graph_source_name", ""),
        "graph_threshold": bundle.get("string_threshold", ""),
        "graph_contract": bundle["graph_contract"],
        "graph_score_column": str(graph_metadata.get("score_column", "")),
        "graph_source_columns": "|".join(str(column) for column in graph_metadata.get("source_columns", [])),
        "graph_has_edge_weights": str(bool(graph_metadata.get("has_edge_weights", False))).lower(),
        "threshold_strategy": "fixed_0.5",
        "evaluation_contract": "auroc_auprc_mcc_specificity_fixed_0.5",
        "feature_contract_group": resolve_feature_contract_group(model_key, feature_setting),
        "label_manifest": bundle["label_manifest_path"],
        "split_manifest": bundle["split_manifest_path"],
        "config_used": str(Path(args.config)),
        "val_auroc": val_metrics["auroc"],
        "val_auprc": val_metrics["auprc"],
        "val_mcc": val_metrics["mcc"],
        "val_f1": val_metrics["f1"],
        "val_precision": val_metrics["precision"],
        "val_recall": val_metrics["recall"],
        "val_accuracy": val_metrics["accuracy"],
        "val_specificity": val_metrics["specificity"],
        "test_auroc": test_metrics["auroc"],
        "test_auprc": test_metrics["auprc"],
        "test_mcc": test_metrics["mcc"],
        "test_f1": test_metrics["f1"],
        "test_precision": test_metrics["precision"],
        "test_recall": test_metrics["recall"],
        "test_accuracy": test_metrics["accuracy"],
        "test_specificity": test_metrics["specificity"],
        "train_count": int(len(bundle["train_idx"])),
        "val_count": int(len(bundle["val_idx"])),
        "test_count": int(len(bundle["test_idx"])),
        "feature_dim": int(bundle["feature_matrix"].shape[1]),
        "esm2_dim": feature_metadata.get("esm2_metadata", {}).get("embedding_dim", ""),
        "node_count": int(len(bundle["node_manifest"])),
        "edge_count": int(bundle["edge_index"].shape[0]),
        "esm2_cache_path": str(feature_metadata.get("esm2_metadata", {}).get("cache_path", "")),
        "best_checkpoint": str(checkpoint_path),
        "best_epoch": checkpoint.get("epoch", ""),
        "best_val_score": checkpoint.get("best_metric", ""),
        "materialized_from_existing_checkpoint": "true",
    }
    pd.DataFrame([metrics_row]).to_csv(output_dir / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "epoch": checkpoint.get("epoch", ""),
                "best_val_score": checkpoint.get("best_metric", ""),
                "materialized_from_existing_checkpoint": True,
            }
        ]
    ).to_csv(output_dir / "training_log.tsv", sep="\t", index=False)
    save_pickle({"model": model_key, "state_dict": model.state_dict(), "config": model_cfg}, output_dir / "model.pkl")

    resolved = {
        "protocol": bundle["protocol_name"],
        "species": bundle["species"],
        "regime": bundle["regime"],
        "model": model_key,
        "model_variant": args.model,
        "feature_setting": feature_setting,
        "seed": int(args.seed),
        "graph_source": bundle["graph_source"],
        "graph_source_name": bundle.get("graph_source_name", ""),
        "graph_threshold": bundle.get("string_threshold", ""),
        "graph_score_column": str(graph_metadata.get("score_column", "")),
        "graph_source_columns": list(graph_metadata.get("source_columns", [])),
        "graph_has_edge_weights": bool(graph_metadata.get("has_edge_weights", False)),
        "output_dir": str(output_dir),
        "checkpoint_materialization": True,
    }
    save_yaml(resolved, output_dir / "resolved_config.yaml")


if __name__ == "__main__":
    main()
