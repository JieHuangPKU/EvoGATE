import argparse
import pickle
from pathlib import Path

import pandas as pd
import torch
import yaml

from src.data.frozen_protocol_loader import load_protocol_dataset
from src.train.run_frozen_protocol_model import (
    GRAPH_MODELS,
    build_prediction_table,
    compute_binary_metrics,
    fit_graph_model,
    is_gated_feature_setting,
    normalize_model_name,
    resolve_feature_contract_group,
    save_yaml,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one frozen-protocol graph model with an explicit feature combo override")
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


def save_pickle(payload, output_path: Path) -> None:
    with output_path.open("wb") as handle:
        pickle.dump(payload, handle)


def resolve_graph_model_config(base_config: dict, model_key: str) -> dict:
    if model_key in base_config["models"]:
        return dict(base_config["models"][model_key])
    fallback_keys = {
        "GraphSAGE": ["GraphSAGE_ORT_EXP_SUB", "GraphSAGE_ESM2", "GraphSAGE_ORT_EXP_SUB_ESM2"],
        "GAT": ["GAT"],
        "GCN": ["GCN"],
        "GIN": ["GIN"],
    }.get(model_key, [])
    for fallback_key in fallback_keys:
        if fallback_key in base_config["models"]:
            return dict(base_config["models"][fallback_key])
    raise ValueError(f"Unknown model '{model_key}' in config")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    base_config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    model_key = normalize_model_name(args.model)
    if model_key not in GRAPH_MODELS:
        raise ValueError(f"Feature-combo benchmark only supports graph models; got '{model_key}'")
    model_cfg = resolve_graph_model_config(base_config, model_key)
    feature_setting = str(args.feature_setting).strip().upper()
    seed = int(args.seed)

    bundle = load_protocol_dataset(
        args.config,
        args.protocol,
        feature_setting,
        graph_source_path=args.graph_source,
        string_threshold=args.string_threshold,
        graph_source_name=args.graph_source_name,
    )
    print(f"[run_frozen_protocol_feature_combo_model] loaded dataset nodes={len(bundle['node_manifest'])} edges={bundle['edge_index'].shape[0]}", flush=True)
    checkpoint_path = output_dir / "best_model.pt"
    model, pred_score, training_log, best_info = fit_graph_model(bundle, model_key, model_cfg, seed, checkpoint_path)
    print("[run_frozen_protocol_feature_combo_model] finished fit_graph_model", flush=True)

    pred_label = (pred_score >= 0.5).astype(int)
    test_idx = bundle["test_idx"]
    val_idx = bundle["val_idx"]
    test_metrics = compute_binary_metrics(
        bundle["y_all"][test_idx].astype(int),
        pred_score[test_idx],
        pred_label[test_idx],
    )
    val_metrics = compute_binary_metrics(
        bundle["y_all"][val_idx].astype(int),
        pred_score[val_idx],
        pred_label[val_idx],
    )
    print("[run_frozen_protocol_feature_combo_model] computed validation/test metrics", flush=True)

    prediction_table = build_prediction_table(bundle, pred_score, pred_label, model_key)
    prediction_table["feature_combo"] = feature_setting
    prediction_table.to_csv(output_dir / "predictions.tsv", sep="\t", index=False)
    bundle["feature_schema"].to_csv(output_dir / "feature_schema.tsv", sep="\t", index=False)
    pd.DataFrame(bundle["edge_table"]).to_csv(output_dir / "edge_table.tsv", sep="\t", index=False)
    bundle["split_manifest"].to_csv(output_dir / "split_manifest.tsv", sep="\t", index=False)
    print("[run_frozen_protocol_feature_combo_model] wrote predictions/schema/edge/split tables", flush=True)
    feature_metadata = dict(bundle.get("feature_metadata", {}))
    graph_metadata = dict(bundle.get("graph_metadata", {}))
    if "esm2_alignment_audit" in feature_metadata:
        feature_metadata["esm2_alignment_audit"].to_csv(output_dir / "esm2_alignment_audit.tsv", sep="\t", index=False)
    if "esm2_metadata" in feature_metadata:
        pd.DataFrame([feature_metadata["esm2_metadata"]]).to_csv(output_dir / "esm2_feature_metadata.tsv", sep="\t", index=False)
    if hasattr(model, "gate_statistics"):
        model.gate_statistics(
            torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32),
            bundle["node_manifest"],
            bundle["split_manifest"],
        ).to_csv(output_dir / "gate_statistics.tsv", sep="\t", index=False)

    is_gated = is_gated_feature_setting(feature_setting)

    metrics_row = {
        "protocol": bundle["protocol_name"],
        "species": bundle["species"],
        "regime": bundle["regime"],
        "model": model_key,
        "model_variant": args.model,
        "feature_setting": feature_setting,
        "feature_combo": feature_setting,
        "label_regime": bundle["label_regime"],
        "run_id": f"seed_{seed}",
        "seed": seed,
        "is_deterministic": "false",
        "split_version": bundle["split_version"],
        "graph_source": bundle["graph_source"],
        "graph_source_name": bundle.get("graph_source_name", ""),
        "graph_contract": bundle["graph_contract"],
        "graph_threshold": bundle.get("string_threshold", ""),
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
        "fusion_mode": str(model_cfg.get("fusion_mode", "")) if is_gated else "",
        "fusion_hidden_dim": model_cfg.get("fusion_hidden_dim", "") if is_gated else "",
        "fusion_dropout": model_cfg.get("fusion_dropout", "") if is_gated else "",
        "loss_type": str(model_cfg.get("loss_type", "balanced_bce")),
        "pos_weight_mode": str(model_cfg.get("pos_weight_mode", "")),
        "pos_weight_scale": model_cfg.get("pos_weight_scale", ""),
        "best_checkpoint": str(checkpoint_path) if checkpoint_path.exists() else "",
        **best_info,
    }
    pd.DataFrame([metrics_row]).to_csv(output_dir / "metrics.tsv", sep="\t", index=False)
    training_log.to_csv(output_dir / "training_log.tsv", sep="\t", index=False)
    save_pickle({"model": model_key, "state_dict": model.state_dict(), "config": model_cfg}, output_dir / "model.pkl")
    print("[run_frozen_protocol_feature_combo_model] wrote metrics/training/model artifacts", flush=True)

    resolved = {
        "protocol": bundle["protocol_name"],
        "species": bundle["species"],
        "regime": bundle["regime"],
        "model": model_key,
        "model_variant": args.model,
        "feature_setting": feature_setting,
        "feature_combo": feature_setting,
        "seed": seed,
        "label_regime": bundle["label_regime"],
        "split_version": bundle["split_version"],
        "graph_source": bundle["graph_source"],
        "graph_source_name": bundle.get("graph_source_name", ""),
        "graph_threshold": bundle.get("string_threshold", ""),
        "graph_contract": bundle["graph_contract"],
        "graph_score_column": str(graph_metadata.get("score_column", "")),
        "graph_source_columns": list(graph_metadata.get("source_columns", [])),
        "graph_has_edge_weights": bool(graph_metadata.get("has_edge_weights", False)),
        "label_manifest": bundle["label_manifest_path"],
        "split_manifest": bundle["split_manifest_path"],
        "output_dir": str(output_dir),
        "config_path": str(Path(args.config).resolve()),
        "metrics": metrics_row,
    }
    save_yaml(resolved, output_dir / "resolved_config.yaml")
    print("[run_frozen_protocol_feature_combo_model] completed successfully", flush=True)


if __name__ == "__main__":
    main()
