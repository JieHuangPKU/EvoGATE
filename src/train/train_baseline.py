from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from src.features.load_embeddings import load_embedding_index, load_feature_matrix
from src.models.baseline_models import create_baseline_model, normalize_model_name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train embedding-first baseline model")
    parser.add_argument("--config", type=str, required=True, help="Path to baseline YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def resolve_model_name(config: dict[str, Any]) -> str:
    model_cfg = dict(config.get("model", {}))
    raw_name = model_cfg.get("name", model_cfg.get("type", "logistic_regression"))
    return normalize_model_name(raw_name)


def resolve_training_output_dir(config: dict[str, Any]) -> Path:
    root = Path(config["paths"]["training_output_dir"])
    return root / resolve_model_name(config)


def _safe_binary_metrics(y_true: np.ndarray, y_score: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    metrics = {"auroc": float("nan"), "auprc": float("nan"), "mcc": float("nan")}
    unique = np.unique(y_true)
    if unique.size >= 2:
        metrics["auroc"] = float(roc_auc_score(y_true, y_score))
        metrics["auprc"] = float(average_precision_score(y_true, y_score))
    metrics["mcc"] = float(matthews_corrcoef(y_true, y_pred))
    return metrics


def _load_split_with_logging(
    split_name: str,
    split_df: pd.DataFrame,
    embedding_index,
    require_all: bool,
    require_pooled_features: bool,
) -> tuple[np.ndarray, pd.DataFrame]:
    print(f"[train_baseline] load split={split_name} rows={len(split_df)}")
    x, ready = load_feature_matrix(
        split_df,
        embedding_index,
        require_all=require_all,
        require_pooled_features=require_pooled_features,
    )
    patched = int(ready["alignment_patch_applied"].astype(str).str.lower().isin(["true", "1", "yes"]).sum())
    print(
        f"[train_baseline] loaded split={split_name} features={len(ready)} "
        f"missing={len(split_df) - len(ready)} patched={patched}"
    )
    return x, ready


def _write_train_summary(
    output_path: Path,
    metrics_df: pd.DataFrame,
    prediction_table: pd.DataFrame,
    config: dict,
    train_ready: pd.DataFrame,
    val_ready: pd.DataFrame,
    test_ready: pd.DataFrame,
    inference_ready: pd.DataFrame,
    inference_requested: int,
) -> None:
    model_name = resolve_model_name(config)
    patch_enabled = bool(config["embeddings"].get("use_alignment_patch", False))
    patch_used = int(
        prediction_table.get("alignment_patch_applied", pd.Series(dtype=str)).astype(str).str.lower().isin(["true", "1", "yes"]).sum()
    )
    lines = [
        "# Baseline Train Summary",
        "",
        "## Run",
        f"- model: {model_name}",
        f"- output_dir: {output_path.parent}",
        "",
        "## Loading",
        f"- loaded samples: train={len(train_ready)}, val={len(val_ready)}, test={len(test_ready)}, fusarium={len(inference_ready)}",
        f"- loaded features: train={len(train_ready)}, val={len(val_ready)}, test={len(test_ready)}, fusarium={len(inference_ready)}",
        f"- missing features: train=0, val=0, test=0, fusarium={inference_requested - len(inference_ready)}",
        f"- alignment patch enabled: {str(patch_enabled).lower()}",
        f"- alignment patch used rows: {patch_used}",
        "",
        "## Metrics",
    ]
    for _, row in metrics_df.iterrows():
        lines.append(
            f"- {row['subset']}: n={int(row['n_rows'])}, AUROC={row['auroc']}, "
            f"AUPRC={row['auprc']}, MCC={row['mcc']}"
        )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    model_name = resolve_model_name(config)

    dataset_dir = Path(config["paths"]["baseline_dataset_dir"])
    output_dir = resolve_training_output_dir(config)
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_path = dataset_dir / "support_supervised_samples.tsv"
    if not dataset_path.exists():
        raise FileNotFoundError(
            f"Baseline dataset not found: {dataset_path}. Run build_baseline_dataset.py first."
        )

    manifest_path = config["embeddings"].get("manifest_path")
    embedding_index = load_embedding_index(manifest_path, config["embeddings"]["source_name"])
    require_pooled_features = bool(config["embeddings"].get("require_pooled_features", False))

    samples = pd.read_csv(dataset_path, sep="\t", dtype=str).fillna("")
    train_df = samples[samples["split"] == "train"].copy()
    val_df = samples[samples["split"] == "val"].copy()
    test_df = samples[samples["split"] == "test"].copy()
    inference_df = pd.read_csv(dataset_dir / "fgraminearum_inference_pool.tsv", sep="\t", dtype=str).fillna("")

    model = create_baseline_model(
        model_name=model_name,
        config=config["model"].get(model_name, {}),
        random_state=int(config["train"]["random_seed"]),
    )

    x_train, train_ready = _load_split_with_logging(
        "train",
        train_df,
        embedding_index,
        require_all=True,
        require_pooled_features=require_pooled_features,
    )
    y_train = train_ready["gold_label"].astype(int).to_numpy()
    x_val, val_ready = _load_split_with_logging(
        "val",
        val_df,
        embedding_index,
        require_all=True,
        require_pooled_features=require_pooled_features,
    )
    y_val = val_ready["gold_label"].astype(int).to_numpy()
    x_test, test_ready = _load_split_with_logging(
        "test",
        test_df,
        embedding_index,
        require_all=True,
        require_pooled_features=require_pooled_features,
    )
    y_test = test_ready["gold_label"].astype(int).to_numpy()

    print(f"[train_baseline] model={model_name} train start")
    model.fit(x_train, y_train)
    print(f"[train_baseline] model={model_name} train finish")
    model.save(output_dir / "baseline_model.pkl")

    predictions = []
    for split_name, split_df, x_split in [
        ("train", train_ready, x_train),
        ("val", val_ready, x_val),
        ("test", test_ready, x_test),
    ]:
        probs = model.predict_proba(x_split)
        preds = (probs >= 0.5).astype(int)
        temp = split_df.copy()
        temp["pred_score"] = probs
        temp["pred_label"] = preds
        temp["eval_group"] = "support_species"
        temp["model_name"] = model_name
        predictions.append(temp)

    try:
        x_fg, inference_with_embeddings = _load_split_with_logging(
            "fgraminearum_inference",
            inference_df,
            embedding_index,
            require_all=False,
            require_pooled_features=require_pooled_features,
        )
        fg_probs = model.predict_proba(x_fg)
        inference_with_embeddings["pred_score"] = fg_probs
        inference_with_embeddings["pred_label"] = (fg_probs >= 0.5).astype(int)
        inference_with_embeddings["eval_group"] = "fgraminearum_inference"
        inference_with_embeddings["model_name"] = model_name
        predictions.append(inference_with_embeddings)
    except FileNotFoundError:
        inference_with_embeddings = pd.DataFrame(
            columns=inference_df.columns.tolist() + ["pred_score", "pred_label", "eval_group", "model_name"]
        )

    prediction_table = pd.concat(predictions, ignore_index=True)
    prediction_table.to_csv(output_dir / "predictions.tsv", sep="\t", index=False)
    prediction_table.to_csv(output_dir / "prediction_table.tsv", sep="\t", index=False)

    metric_rows = []
    for split_name, x_split, y_split in [
        ("train", x_train, y_train),
        ("val", x_val, y_val),
        ("test", x_test, y_test),
    ]:
        y_score = model.predict_proba(x_split)
        y_pred = model.predict(x_split)
        base_metrics = _safe_binary_metrics(y_split, y_score, y_pred)
        metric_rows.append(
            {
                "model_name": model_name,
                "metric_scope": "training_split",
                "subset": split_name,
                "n_rows": int(len(y_split)),
                "auroc": base_metrics["auroc"],
                "auprc": base_metrics["auprc"],
                "mcc": base_metrics["mcc"],
            }
        )

    metrics_df = pd.DataFrame(metric_rows)
    metrics_df.to_csv(output_dir / "metrics.tsv", sep="\t", index=False)
    metrics_df.to_csv(output_dir / "training_fit_metrics.tsv", sep="\t", index=False)
    (output_dir / "training_metrics.json").write_text(
        json.dumps(metric_rows, indent=2, allow_nan=True),
        encoding="utf-8",
    )
    _write_train_summary(
        output_dir / "train_summary.md",
        metrics_df,
        prediction_table,
        config,
        train_ready,
        val_ready,
        test_ready,
        inference_with_embeddings,
        len(inference_df),
    )

    print(f"[train_baseline] model={model_name} artifact={output_dir / 'baseline_model.pkl'}")
    print(f"[train_baseline] model={model_name} predictions={output_dir / 'predictions.tsv'}")
    print(f"[train_baseline] model={model_name} metrics={output_dir / 'metrics.tsv'}")
    print(f"[train_baseline] model={model_name} fusarium_predictions={len(inference_with_embeddings)}")


if __name__ == "__main__":
    main()
