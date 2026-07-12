from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, train_test_split
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train a minimal MLP baseline on pooled embedding features")
    parser.add_argument("--dataset-pt", required=True, help="Input dataset artifact created by build_feature_matrix.py")
    parser.add_argument("--output-dir", required=True, help="Directory for model artifacts")
    parser.add_argument("--metrics-json", required=True, help="Output metrics JSON path")
    parser.add_argument("--metrics-tsv", required=True, help="Output metrics TSV path")
    parser.add_argument("--split-mode", choices=["predefined", "stratified", "kfold"], default="predefined")
    parser.add_argument("--train-split", default="train")
    parser.add_argument("--val-split", default="val")
    parser.add_argument("--test-split", default="test")
    parser.add_argument("--val-fraction", type=float, default=0.2)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--k-folds", type=int, default=0)
    parser.add_argument("--random-seed", type=int, default=20260411)
    parser.add_argument("--hidden-dims", nargs="+", type=int, default=[64, 16])
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--max-epochs", type=int, default=20)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--device", default="cpu")
    return parser.parse_args()


class MLPBinaryClassifier(nn.Module):
    def __init__(self, input_dim: int, hidden_dims: list[int], dropout: float) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        last_dim = input_dim
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(last_dim, hidden_dim))
            layers.append(nn.ReLU())
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            last_dim = hidden_dim
        layers.append(nn.Linear(last_dim, 1))
        self.network = nn.Sequential(*layers)

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        return self.network(features).squeeze(-1)


def safe_metric(default: float, fn, *args) -> float:
    try:
        value = fn(*args)
        if isinstance(value, (np.floating, float)) and math.isnan(float(value)):
            return default
        return float(value)
    except Exception:
        return default


def compute_metrics(y_true: np.ndarray, y_score: np.ndarray) -> dict[str, float]:
    y_pred = (y_score >= 0.5).astype(int)
    specificity = float("nan")
    negatives = y_true == 0
    if negatives.any():
        specificity = float(((y_pred == 0) & negatives).sum() / negatives.sum())

    return {
        "n_samples": int(len(y_true)),
        "positive_rate": float(np.mean(y_true)),
        "loss": float("nan"),
        "accuracy": safe_metric(float("nan"), accuracy_score, y_true, y_pred),
        "precision": safe_metric(float("nan"), lambda a, b: precision_score(a, b, zero_division=0), y_true, y_pred),
        "recall": safe_metric(float("nan"), lambda a, b: recall_score(a, b, zero_division=0), y_true, y_pred),
        "f1": safe_metric(float("nan"), lambda a, b: f1_score(a, b, zero_division=0), y_true, y_pred),
        "mcc": safe_metric(float("nan"), matthews_corrcoef, y_true, y_pred),
        "auroc": safe_metric(float("nan"), roc_auc_score, y_true, y_score),
        "auprc": safe_metric(float("nan"), average_precision_score, y_true, y_score),
        "specificity": specificity,
    }


def build_dataloader(features: torch.Tensor, labels: torch.Tensor, batch_size: int, shuffle: bool) -> DataLoader:
    dataset = TensorDataset(features, labels)
    return DataLoader(dataset, batch_size=min(batch_size, len(dataset)), shuffle=shuffle)


def standardize(
    train_x: torch.Tensor,
    other_tensors: list[torch.Tensor],
) -> tuple[torch.Tensor, list[torch.Tensor], torch.Tensor, torch.Tensor]:
    mean = train_x.mean(dim=0, keepdim=True)
    std = train_x.std(dim=0, keepdim=True, unbiased=False).clamp(min=1e-6)
    standardized_train = (train_x - mean) / std
    standardized_others = [(tensor - mean) / std for tensor in other_tensors]
    return standardized_train, standardized_others, mean.squeeze(0), std.squeeze(0)


def train_one_model(
    train_x: torch.Tensor,
    train_y: torch.Tensor,
    val_x: torch.Tensor,
    val_y: torch.Tensor,
    input_dim: int,
    hidden_dims: list[int],
    dropout: float,
    learning_rate: float,
    weight_decay: float,
    batch_size: int,
    max_epochs: int,
    patience: int,
    device: str,
) -> tuple[MLPBinaryClassifier, dict[str, list[float]], float]:
    resolved_device = torch.device(device)
    model = MLPBinaryClassifier(input_dim=input_dim, hidden_dims=hidden_dims, dropout=dropout).to(resolved_device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    criterion = nn.BCEWithLogitsLoss()

    train_loader = build_dataloader(train_x, train_y, batch_size=batch_size, shuffle=True)
    history = {"train_loss": [], "val_loss": []}

    best_state = None
    best_val_loss = float("inf")
    stale_epochs = 0

    for _epoch in range(max_epochs):
        model.train()
        epoch_losses: list[float] = []
        for batch_x, batch_y in train_loader:
            batch_x = batch_x.to(resolved_device)
            batch_y = batch_y.to(resolved_device)
            optimizer.zero_grad()
            logits = model(batch_x)
            loss = criterion(logits, batch_y)
            loss.backward()
            optimizer.step()
            epoch_losses.append(float(loss.item()))

        model.eval()
        with torch.inference_mode():
            val_logits = model(val_x.to(resolved_device))
            val_loss = float(criterion(val_logits, val_y.to(resolved_device)).item())

        history["train_loss"].append(float(np.mean(epoch_losses)))
        history["val_loss"].append(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            stale_epochs = 0
        else:
            stale_epochs += 1
            if stale_epochs >= patience:
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model, history, best_val_loss


def predict_probabilities(model: nn.Module, features: torch.Tensor, device: str) -> np.ndarray:
    model.eval()
    with torch.inference_mode():
        logits = model(features.to(torch.device(device)))
        probs = torch.sigmoid(logits).detach().cpu().numpy()
    return probs.astype(np.float64)


def save_model_artifact(
    output_path: Path,
    model: nn.Module,
    mean: torch.Tensor,
    std: torch.Tensor,
    metadata: dict[str, object],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "scaler_mean": mean,
            "scaler_std": std,
            "metadata": metadata,
        },
        output_path,
    )


def _select_predefined_indices(
    table: pd.DataFrame,
    train_split: str,
    val_split: str,
    test_split: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if "split" not in table.columns:
        raise ValueError("Dataset artifact does not include a split column required by split_mode=predefined")
    train_idx = np.flatnonzero(table["split"].astype(str).to_numpy() == train_split)
    val_idx = np.flatnonzero(table["split"].astype(str).to_numpy() == val_split)
    test_idx = np.flatnonzero(table["split"].astype(str).to_numpy() == test_split)
    if len(train_idx) == 0 or len(val_idx) == 0 or len(test_idx) == 0:
        raise ValueError("Predefined split requires non-empty train/val/test partitions")
    return train_idx, val_idx, test_idx


def _select_stratified_indices(
    y: np.ndarray,
    test_fraction: float,
    val_fraction: float,
    random_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    all_idx = np.arange(len(y))
    train_val_idx, test_idx = train_test_split(
        all_idx,
        test_size=test_fraction,
        random_state=random_seed,
        stratify=y,
    )
    relative_val_fraction = val_fraction / max(1e-8, 1.0 - test_fraction)
    train_idx, val_idx = train_test_split(
        train_val_idx,
        test_size=relative_val_fraction,
        random_state=random_seed,
        stratify=y[train_val_idx],
    )
    return train_idx, val_idx, test_idx


def run_predefined_or_stratified(args: argparse.Namespace, payload: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, object]]:
    output_dir = Path(args.output_dir)
    features = payload["X"].detach().cpu().to(torch.float32)
    labels = payload["y"].detach().cpu().to(torch.float32)
    table = pd.DataFrame(payload.get("table", []))

    if args.split_mode == "predefined":
        if table.empty and payload.get("splits") is not None:
            table = pd.DataFrame({"split": payload["splits"]})
        train_idx, val_idx, test_idx = _select_predefined_indices(table, args.train_split, args.val_split, args.test_split)
    else:
        train_idx, val_idx, test_idx = _select_stratified_indices(
            y=payload["y"].detach().cpu().numpy(),
            test_fraction=float(args.test_fraction),
            val_fraction=float(args.val_fraction),
            random_seed=int(args.random_seed),
        )

    train_x = features[train_idx]
    val_x = features[val_idx]
    test_x = features[test_idx]
    train_y = labels[train_idx]
    val_y = labels[val_idx]
    test_y = labels[test_idx]

    train_x, [val_x, test_x], mean, std = standardize(train_x, [val_x, test_x])
    model, history, best_val_loss = train_one_model(
        train_x=train_x,
        train_y=train_y,
        val_x=val_x,
        val_y=val_y,
        input_dim=int(features.shape[1]),
        hidden_dims=list(args.hidden_dims),
        dropout=float(args.dropout),
        learning_rate=float(args.learning_rate),
        weight_decay=float(args.weight_decay),
        batch_size=int(args.batch_size),
        max_epochs=int(args.max_epochs),
        patience=int(args.patience),
        device=args.device,
    )

    save_model_artifact(
        output_dir / "best_model.pt",
        model=model,
        mean=mean,
        std=std,
        metadata={
            "mode": args.split_mode,
            "input_dim": int(features.shape[1]),
            "hidden_dims": list(args.hidden_dims),
            "dropout": float(args.dropout),
        },
    )
    (output_dir / "history.json").write_text(json.dumps(history, indent=2), encoding="utf-8")

    metric_rows: list[dict[str, object]] = []
    for split_name, split_x, split_y in [
        ("train", train_x, train_y),
        ("val", val_x, val_y),
        ("test", test_x, test_y),
    ]:
        y_true = split_y.detach().cpu().numpy().astype(int)
        y_score = predict_probabilities(model, split_x, args.device)
        metrics = compute_metrics(y_true, y_score)
        metrics["loss"] = float(best_val_loss) if split_name == "val" else float("nan")
        metric_rows.append({"fold": "holdout", "eval_split": split_name, **metrics})

    run_metadata = {
        "mode": args.split_mode,
        "n_samples": int(features.shape[0]),
        "feature_dim": int(features.shape[1]),
        "history": history,
    }
    return metric_rows, run_metadata


def run_kfold(args: argparse.Namespace, payload: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, object]]:
    output_dir = Path(args.output_dir)
    features = payload["X"].detach().cpu().to(torch.float32)
    labels_np = payload["y"].detach().cpu().numpy().astype(int)
    labels = payload["y"].detach().cpu().to(torch.float32)

    splitter = StratifiedKFold(
        n_splits=int(args.k_folds),
        shuffle=True,
        random_state=int(args.random_seed),
    )

    metric_rows: list[dict[str, object]] = []
    histories: dict[str, object] = {}

    for fold_id, (train_val_idx, test_idx) in enumerate(splitter.split(np.zeros(len(labels_np)), labels_np)):
        train_idx, val_idx = train_test_split(
            train_val_idx,
            test_size=float(args.val_fraction),
            random_state=int(args.random_seed) + fold_id,
            stratify=labels_np[train_val_idx],
        )

        train_x = features[train_idx]
        val_x = features[val_idx]
        test_x = features[test_idx]
        train_y = labels[train_idx]
        val_y = labels[val_idx]
        test_y = labels[test_idx]

        train_x, [val_x, test_x], mean, std = standardize(train_x, [val_x, test_x])
        model, history, best_val_loss = train_one_model(
            train_x=train_x,
            train_y=train_y,
            val_x=val_x,
            val_y=val_y,
            input_dim=int(features.shape[1]),
            hidden_dims=list(args.hidden_dims),
            dropout=float(args.dropout),
            learning_rate=float(args.learning_rate),
            weight_decay=float(args.weight_decay),
            batch_size=int(args.batch_size),
            max_epochs=int(args.max_epochs),
            patience=int(args.patience),
            device=args.device,
        )

        fold_dir = output_dir / f"fold_{fold_id}"
        save_model_artifact(
            fold_dir / "model.pt",
            model=model,
            mean=mean,
            std=std,
            metadata={
                "mode": "kfold",
                "fold": fold_id,
                "input_dim": int(features.shape[1]),
                "hidden_dims": list(args.hidden_dims),
                "dropout": float(args.dropout),
            },
        )
        histories[f"fold_{fold_id}"] = history

        for split_name, split_x, split_y in [
            ("val", val_x, val_y),
            ("test", test_x, test_y),
        ]:
            y_true = split_y.detach().cpu().numpy().astype(int)
            y_score = predict_probabilities(model, split_x, args.device)
            metrics = compute_metrics(y_true, y_score)
            metrics["loss"] = float(best_val_loss) if split_name == "val" else float("nan")
            metric_rows.append({"fold": f"fold_{fold_id}", "eval_split": split_name, **metrics})

    run_metadata = {
        "mode": "kfold",
        "k_folds": int(args.k_folds),
        "n_samples": int(features.shape[0]),
        "feature_dim": int(features.shape[1]),
        "histories": histories,
    }
    return metric_rows, run_metadata


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_path = Path(args.dataset_pt)
    try:
        payload = torch.load(dataset_path, map_location="cpu", weights_only=False)
    except TypeError:
        payload = torch.load(dataset_path, map_location="cpu")
    if not isinstance(payload, dict) or "X" not in payload or "y" not in payload:
        raise ValueError(f"Unsupported dataset artifact: {args.dataset_pt}")

    if args.k_folds > 1 or args.split_mode == "kfold":
        if args.k_folds < 2:
            raise ValueError("k-fold mode requires --k-folds >= 2")
        metric_rows, run_metadata = run_kfold(args, payload)
    else:
        metric_rows, run_metadata = run_predefined_or_stratified(args, payload)

    metrics_df = pd.DataFrame(metric_rows)
    metrics_df.to_csv(args.metrics_tsv, sep="\t", index=False)

    metrics_json_payload = {
        "dataset_pt": args.dataset_pt,
        "output_dir": str(output_dir),
        "hyperparameters": {
            "hidden_dims": list(args.hidden_dims),
            "dropout": float(args.dropout),
            "learning_rate": float(args.learning_rate),
            "weight_decay": float(args.weight_decay),
            "batch_size": int(args.batch_size),
            "max_epochs": int(args.max_epochs),
            "patience": int(args.patience),
            "device": args.device,
        },
        "run_metadata": run_metadata,
        "metrics": metric_rows,
    }
    Path(args.metrics_json).write_text(json.dumps(metrics_json_payload, indent=2), encoding="utf-8")
    print(json.dumps({"metrics_tsv": args.metrics_tsv, "metrics_json": args.metrics_json}, indent=2))


if __name__ == "__main__":
    main()
