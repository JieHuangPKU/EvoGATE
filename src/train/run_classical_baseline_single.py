import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

from src.classical_baselines.common import (
    TRAINABLE_METHODS,
    build_dataset_for_benchmark,
    build_prediction_table,
    compute_binary_metrics,
    dump_yaml,
    load_dataset_bundle,
    load_yaml,
    parse_feature_setting,
    set_seed,
)
from src.graph.run_node2vec_embedding import train_node2vec_embeddings, write_node2vec_outputs


class BalancedBCE:
    def __init__(self, y):
        self.y = y
        self.pos_mask = y == 1
        self.neg_mask = y == 0

    def __call__(self, logits):
        pos_loss = F.binary_cross_entropy_with_logits(logits[self.pos_mask].view(-1), self.y[self.pos_mask].view(-1))
        neg_loss = F.binary_cross_entropy_with_logits(logits[self.neg_mask].view(-1), self.y[self.neg_mask].view(-1))
        return pos_loss + neg_loss


class TorchMLP(nn.Module):
    def __init__(self, input_dim, hidden_dim, dropout):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x):
        return self.network(x)


def parse_args():
    parser = argparse.ArgumentParser(description="Run one classical baseline benchmark task")
    parser.add_argument("--config", required=True, type=str)
    parser.add_argument("--species", required=True, type=str)
    parser.add_argument("--method", required=True, choices=sorted(TRAINABLE_METHODS))
    parser.add_argument("--feature-setting", required=True, type=str)
    parser.add_argument("--run-id", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=str)
    return parser.parse_args()


def _to_tensor(matrix):
    return torch.as_tensor(np.asarray(matrix, dtype=np.float32), dtype=torch.float32)


def fit_torch_mlp(x_train, y_train, x_val, y_val, params, seed):
    set_seed(seed)
    model = TorchMLP(input_dim=int(x_train.shape[1]), hidden_dim=int(params["hidden_dim"]), dropout=float(params["dropout"]))
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(params["learning_rate"]),
        weight_decay=float(params["weight_decay"]),
    )
    train_loss_fn = BalancedBCE(y_train)
    val_loss_fn = BalancedBCE(y_val)

    best_state = None
    best_val_loss = float("inf")
    patience_counter = 0
    history_rows = []
    for epoch in range(int(params["epochs"])):
        model.train()
        train_logits = model(x_train)
        train_loss = train_loss_fn(train_logits)
        optimizer.zero_grad()
        train_loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            val_logits = model(x_val)
            val_loss = val_loss_fn(val_logits)
        train_value = float(train_loss.item())
        val_value = float(val_loss.item())
        history_rows.append({"epoch": epoch + 1, "train_loss": train_value, "val_loss": val_value})
        if val_value < best_val_loss:
            best_val_loss = val_value
            best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= int(params["patience"]):
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model, pd.DataFrame(history_rows)


def predict_torch_model(model, matrix):
    model.eval()
    with torch.no_grad():
        logits = model(_to_tensor(matrix)).view(-1)
    return torch.sigmoid(logits).cpu().numpy()


def fit_sklearn_model(method, params, seed):
    if method == "RF":
        return RandomForestClassifier(
            n_estimators=int(params["n_estimators"]),
            class_weight=params.get("class_weight", "balanced"),
            n_jobs=int(params.get("n_jobs", -1)),
            random_state=seed,
        )
    if method == "SVM":
        return SVC(
            probability=True,
            C=float(params["c"]),
            kernel=str(params["kernel"]),
            gamma=params.get("gamma", "scale"),
            class_weight=params.get("class_weight", "balanced"),
            random_state=seed,
        )
    if method == "NB":
        return GaussianNB()
    raise ValueError(f"Unsupported sklearn method: {method}")


def save_pickle(payload, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as handle:
        pickle.dump(payload, handle)


def write_feature_summary(output_dir, species, method, feature_setting, label_regime, run_id, seed, bundle):
    feature_schema = bundle["feature_schema"]
    flags = parse_feature_setting(feature_setting)
    row = {
        "species": species,
        "method": method,
        "feature_setting": feature_setting,
        "run_id": int(run_id),
        "seed": int(seed),
        "label_regime": label_regime,
        "feature_dim": int(np.load(output_dir / "feature_matrix.npy").shape[1]),
        "feature_blocks": "|".join(feature_schema["feature_block"].astype(str).tolist()) if not feature_schema.empty else "",
        "orthologs_enabled": bool(flags["orthologs"]),
        "expression_enabled": bool(flags["expression"]),
        "sublocalization_enabled": bool(flags["sublocalization"]),
        "node_count": int(len(bundle["node_manifest"])),
        "labeled_count": int(len(bundle["labeled"])),
        "essential_count": int((bundle["labeled"]["label_numeric"].astype(float) == 1).sum()),
        "test_count": int(len(bundle["test_idx"])),
    }
    pd.DataFrame([row]).to_csv(output_dir / "feature_summary.tsv", sep="\t", index=False)


def main():
    args = parse_args()
    config = load_yaml(args.config)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    method_cfg = dict(config["methods"][args.method])
    seed = int(config["runtime"]["base_seed"]) + int(args.run_id)

    dataset_meta = build_dataset_for_benchmark(
        config=config,
        species=args.species,
        feature_setting=args.feature_setting,
        output_dir=output_dir,
        seed=seed,
    )
    bundle = load_dataset_bundle(output_dir)

    if args.method == "N2V_MLP":
        node2vec_params = {
            "embedding_dim": int(method_cfg["embedding_dim"]),
            "walk_length": int(method_cfg["walk_length"]),
            "context_size": int(method_cfg["context_size"]),
            "walks_per_node": int(method_cfg["walks_per_node"]),
            "num_negative_samples": int(method_cfg["num_negative_samples"]),
            "epochs": int(method_cfg["epochs"]),
            "batch_size": int(method_cfg["batch_size"]),
            "learning_rate": float(method_cfg["learning_rate"]),
        }
        embeddings, history_df = train_node2vec_embeddings(
            edge_array=bundle["edge_index"],
            num_nodes=len(bundle["node_manifest"]),
            params=node2vec_params,
            seed=seed,
            device=str(config["runtime"].get("device", "cpu")),
        )
        np.save(output_dir / "feature_matrix.npy", embeddings.astype(np.float32, copy=False))
        write_node2vec_outputs(embeddings, node2vec_params, output_dir, history_df)
        bundle = load_dataset_bundle(output_dir)
        mlp_params = {
            "hidden_dim": int(method_cfg["mlp_hidden_dim"]),
            "dropout": float(method_cfg["mlp_dropout"]),
            "learning_rate": float(method_cfg["mlp_learning_rate"]),
            "weight_decay": float(method_cfg["mlp_weight_decay"]),
            "epochs": int(method_cfg["mlp_epochs"]),
            "patience": int(method_cfg["mlp_patience"]),
        }
        model, training_log = fit_torch_mlp(
            _to_tensor(bundle["feature_matrix"][bundle["train_idx"]]),
            torch.as_tensor(bundle["y_all"][bundle["train_idx"]], dtype=torch.float32),
            _to_tensor(bundle["feature_matrix"][bundle["val_idx"]]),
            torch.as_tensor(bundle["y_all"][bundle["val_idx"]], dtype=torch.float32),
            mlp_params,
            seed,
        )
        pred_score = predict_torch_model(model, bundle["feature_matrix"])
        save_pickle({"method": args.method, "state_dict": model.state_dict(), "params": mlp_params}, output_dir / "model.pkl")
        training_log.to_csv(output_dir / "training_log.tsv", sep="\t", index=False)
        method_params = {"node2vec": node2vec_params, "mlp": mlp_params}
    elif args.method == "MLP":
        mlp_params = {
            "hidden_dim": int(method_cfg["hidden_dim"]),
            "dropout": float(method_cfg["dropout"]),
            "learning_rate": float(method_cfg["learning_rate"]),
            "weight_decay": float(method_cfg["weight_decay"]),
            "epochs": int(method_cfg["epochs"]),
            "patience": int(method_cfg["patience"]),
        }
        model, training_log = fit_torch_mlp(
            _to_tensor(bundle["feature_matrix"][bundle["train_idx"]]),
            torch.as_tensor(bundle["y_all"][bundle["train_idx"]], dtype=torch.float32),
            _to_tensor(bundle["feature_matrix"][bundle["val_idx"]]),
            torch.as_tensor(bundle["y_all"][bundle["val_idx"]], dtype=torch.float32),
            mlp_params,
            seed,
        )
        pred_score = predict_torch_model(model, bundle["feature_matrix"])
        save_pickle({"method": args.method, "state_dict": model.state_dict(), "params": mlp_params}, output_dir / "model.pkl")
        training_log.to_csv(output_dir / "training_log.tsv", sep="\t", index=False)
        method_params = mlp_params
    else:
        sklearn_model = fit_sklearn_model(args.method, method_cfg, seed)
        sklearn_model.fit(bundle["feature_matrix"][bundle["train_idx"]], bundle["y_all"][bundle["train_idx"]].astype(int))
        pred_score = sklearn_model.predict_proba(bundle["feature_matrix"])[:, 1]
        save_pickle({"method": args.method, "model": sklearn_model, "params": method_cfg}, output_dir / "model.pkl")
        method_params = method_cfg

    pred_label = (pred_score >= 0.5).astype(int)
    predictions = build_prediction_table(
        bundle["node_manifest"],
        bundle["label_manifest"],
        pred_score=pred_score,
        pred_label=pred_label,
        method=args.method,
        feature_setting=args.feature_setting,
    )
    predictions.to_csv(output_dir / "predictions.tsv", sep="\t", index=False)

    test_idx = bundle["test_idx"]
    y_test = bundle["y_all"][test_idx].astype(int)
    test_score = pred_score[test_idx]
    test_pred = pred_label[test_idx]
    metrics = compute_binary_metrics(y_test, test_score, test_pred)

    val_idx = bundle["val_idx"]
    val_metrics = compute_binary_metrics(bundle["y_all"][val_idx].astype(int), pred_score[val_idx], pred_label[val_idx])
    metrics_row = {
        "species": args.species,
        "method": args.method,
        "feature_setting": args.feature_setting,
        "label_regime": dataset_meta["label_regime"],
        "run_id": f"run_{args.run_id}",
        "seed": int(seed),
        "test_count": int(len(test_idx)),
        "val_auroc": val_metrics["auroc"],
        "val_auprc": val_metrics["auprc"],
        "val_mcc": val_metrics["mcc"],
        **metrics,
    }
    pd.DataFrame([metrics_row]).to_csv(output_dir / "metrics.tsv", sep="\t", index=False)

    write_feature_summary(
        output_dir=output_dir,
        species=args.species,
        method=args.method,
        feature_setting=args.feature_setting,
        label_regime=dataset_meta["label_regime"],
        run_id=args.run_id,
        seed=seed,
        bundle=bundle,
    )
    resolved_config = {
        "species": args.species,
        "method": args.method,
        "feature_setting": args.feature_setting,
        "run_id": int(args.run_id),
        "seed": int(seed),
        "output_dir": str(output_dir),
        "label_regime": dataset_meta["label_regime"],
        "positive_set_path": dataset_meta["positive_set_path"],
        "negative_set_path": dataset_meta["negative_set_path"],
        "dataset_config": dataset_meta["builder_config"],
        "method_params": method_params,
        "metrics": metrics_row,
    }
    dump_yaml(resolved_config, output_dir / "resolved_config.yaml")


if __name__ == "__main__":
    main()
