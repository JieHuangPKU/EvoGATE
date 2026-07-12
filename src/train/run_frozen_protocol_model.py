import argparse
import math
import pickle
import random
from pathlib import Path
from typing import Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, matthews_corrcoef, precision_score, recall_score, roc_auc_score
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC

from src.data.frozen_protocol_loader import load_protocol_dataset


GRAPH_MODELS = {"GAT", "GCN", "GIN", "GraphSAGE"}
TRAINABLE_MODELS = {
    "MLP",
    "RF",
    "SVM",
    "NB",
    "N2V_MLP",
    "GAT",
    "GCN",
    "GIN",
    "GRAPHSAGE",
    "GRAPHSAGE_ORT_EXP_SUB",
    "GRAPHSAGE_ESM2",
    "GRAPHSAGE_ORT_EXP_SUB_ESM2",
    "GRAPHSAGE_ORT_EXP_SUB_ESM2_GATED",
}
DETERMINISTIC_MODELS = {"DC", "CC"}


class BalancedBCE:
    def __init__(self, y: torch.Tensor):
        self.y = y.view(-1)
        self.pos_mask = self.y == 1
        self.neg_mask = self.y == 0

    def __call__(self, logits: torch.Tensor) -> torch.Tensor:
        logits = logits.view(-1)
        pos_loss = F.binary_cross_entropy_with_logits(logits[self.pos_mask], self.y[self.pos_mask])
        neg_loss = F.binary_cross_entropy_with_logits(logits[self.neg_mask], self.y[self.neg_mask])
        return pos_loss + neg_loss


class WeightedBCE:
    def __init__(self, y: torch.Tensor, pos_weight_value: float):
        self.y = y.view(-1)
        pos_weight = torch.tensor([float(pos_weight_value)], dtype=torch.float32, device=self.y.device)
        self.loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_weight, reduction="mean")

    def __call__(self, logits: torch.Tensor) -> torch.Tensor:
        logits = logits.view(-1)
        return self.loss_fn(logits, self.y)


def is_gated_feature_setting(feature_setting: str) -> bool:
    return "_GATED" in str(feature_setting).strip().upper()


def build_loss_metadata(labels: np.ndarray, model_cfg: dict) -> Tuple[str, Optional[float]]:
    loss_type = str(model_cfg.get("loss_type", "balanced_bce")).strip().lower()
    if loss_type == "balanced_bce":
        return loss_type, None
    if loss_type != "weighted_bce":
        raise ValueError(f"Unsupported loss_type '{loss_type}'")

    labels = np.asarray(labels).astype(int)
    num_pos = int((labels == 1).sum())
    num_neg = int((labels == 0).sum())
    if num_pos == 0:
        return loss_type, 1.0
    raw_ratio = float(num_neg) / float(num_pos)
    mode = str(model_cfg.get("pos_weight_mode", "sqrt_ratio")).strip().lower()
    scale = float(model_cfg.get("pos_weight_scale", 1.0))
    if mode == "sqrt_ratio":
        base_weight = math.sqrt(raw_ratio)
    elif mode == "ratio":
        base_weight = raw_ratio
    elif mode == "none":
        base_weight = 1.0
    else:
        raise ValueError(f"Unsupported pos_weight_mode '{mode}'")
    return loss_type, float(base_weight * scale)


def create_loss_fn(labels: torch.Tensor, model_cfg: dict, pos_weight_value: Optional[float] = None):
    loss_type = str(model_cfg.get("loss_type", "balanced_bce")).strip().lower()
    if loss_type == "balanced_bce":
        return BalancedBCE(labels)
    if loss_type != "weighted_bce":
        raise ValueError(f"Unsupported loss_type '{loss_type}'")
    if pos_weight_value is None:
        _, pos_weight_value = build_loss_metadata(labels.detach().cpu().numpy(), model_cfg)
    return WeightedBCE(labels, float(pos_weight_value))


class TorchMLP(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, dropout: float):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one frozen-protocol benchmark task")
    parser.add_argument("--config", required=True, type=str)
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--model", required=True, type=str)
    parser.add_argument("--feature-setting", default=None, type=str)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--output-dir", required=True, type=str)
    parser.add_argument("--graph-contract", default=None, type=str)
    parser.add_argument("--split-manifest", default=None, type=str)
    return parser.parse_args()


def normalize_model_name(model_name: str) -> str:
    raw = str(model_name).strip().lower()
    aliases = {
        "mlp": "MLP",
        "rf": "RF",
        "svm": "SVM",
        "nb": "NB",
        "n2v_mlp": "N2V_MLP",
        "dc": "DC",
        "cc": "CC",
        "gat": "GAT",
        "gcn": "GCN",
        "gin": "GIN",
        "graphsage": "GraphSAGE",
        "sage": "GraphSAGE",
        "graphsage_ort_exp_sub": "GraphSAGE",
        "graphsage_esm2": "GraphSAGE",
        "graphsage_ort_exp_sub_esm2": "GraphSAGE",
        "graphsage_ort_exp_sub_esm2_gated": "GraphSAGE",
    }
    if raw not in aliases:
        raise ValueError(f"Unsupported model '{model_name}'")
    return aliases[raw]


def resolve_feature_contract_group(model_key: str, feature_setting: str) -> str:
    normalized_feature = str(feature_setting).strip().upper()
    if model_key in {"N2V_MLP", "DC", "CC"}:
        return "topology_embedding_contract"
    if normalized_feature == "ORT_EXP_SUB":
        return "ort_exp_sub_contract"
    if normalized_feature == "ORT_ESM2":
        return "ort_plus_esm2_species_level_pooled_contract"
    if normalized_feature == "ESM2":
        return "esm2_species_level_pooled_contract"
    if normalized_feature == "ORT_EXP_SUB_ESM2":
        return "ort_exp_sub_plus_esm2_species_level_pooled_contract"
    if normalized_feature == "ORT_EXP_SUB_ESM2_GATED":
        return "ort_exp_sub_plus_esm2_species_level_pooled_contract"
    if normalized_feature in {
        "ORT_EXP_SUB_ESM2_OLD_GATED_WBCE",
        "ORT_EXP_SUB_ESM2_GATED_RESIDUAL",
        "ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE",
    }:
        return "ort_exp_sub_plus_esm2_species_level_pooled_contract"
    return "feature_contract"


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def save_yaml(payload: dict, output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def compute_binary_metrics(y_true, y_score, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    y_pred = np.asarray(y_pred).astype(int)
    out = {
        "auroc": float("nan"),
        "auprc": float("nan"),
        "mcc": float("nan"),
        "f1": float("nan"),
        "precision": float("nan"),
        "recall": float("nan"),
        "accuracy": float("nan"),
        "specificity": float("nan"),
    }
    if y_true.size == 0:
        return out
    if np.unique(y_true).size > 1:
        out["auroc"] = float(roc_auc_score(y_true, y_score))
        out["auprc"] = float(average_precision_score(y_true, y_score))
        out["mcc"] = float(matthews_corrcoef(y_true, y_pred))
    out["f1"] = float(f1_score(y_true, y_pred, zero_division=0))
    out["precision"] = float(precision_score(y_true, y_pred, zero_division=0))
    out["recall"] = float(recall_score(y_true, y_pred, zero_division=0))
    out["accuracy"] = float(accuracy_score(y_true, y_pred))
    negatives = y_true == 0
    if int(negatives.sum()) > 0:
        out["specificity"] = float(((y_pred == 0) & negatives).sum() / negatives.sum())
    return out


def build_prediction_table(bundle: dict, pred_score: np.ndarray, pred_label: np.ndarray, model_name: str) -> pd.DataFrame:
    table = bundle["node_manifest"][["canonical_gene_id", "graph_gene_id", "split", "label", "is_labeled", "in_graph"]].copy()
    table["pred_score"] = np.asarray(pred_score, dtype=float)
    table["pred_label"] = np.asarray(pred_label).astype(int)
    table["protocol"] = bundle["protocol_name"]
    table["species"] = bundle["species"]
    table["regime"] = bundle["regime"]
    table["model"] = model_name
    table["feature_setting"] = bundle["feature_setting"]
    table["split_version"] = bundle["split_version"]
    table["graph_contract"] = bundle["graph_contract"]
    return table


def _selection_metric(metric_row):
    auprc = metric_row.get("auprc", float("nan"))
    if not math.isnan(auprc):
        return float(auprc)
    auroc = metric_row.get("auroc", float("nan"))
    if not math.isnan(auroc):
        return float(auroc)
    return float(metric_row.get("accuracy", float("-inf")))


def _save_best_checkpoint(model, checkpoint_path, epoch, best_metric):
    torch.save(
        {
            "epoch": int(epoch),
            "best_metric": float(best_metric),
            "state_dict": model.state_dict(),
        },
        checkpoint_path,
    )


def fit_torch_mlp(
    bundle: dict,
    model_cfg: dict,
    seed: int,
    checkpoint_path,
):
    set_seed(seed)
    x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
    y_all = torch.as_tensor(bundle["y_all"], dtype=torch.float32)
    train_idx = bundle["train_idx"]
    val_idx = bundle["val_idx"]
    model = TorchMLP(x.shape[1], int(model_cfg["hidden_dim"]), float(model_cfg["dropout"]))
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(model_cfg["learning_rate"]),
        weight_decay=float(model_cfg["weight_decay"]),
    )
    _, pos_weight_value = build_loss_metadata(bundle["y_all"][train_idx], model_cfg)
    train_loss_fn = create_loss_fn(y_all[train_idx], model_cfg, pos_weight_value=pos_weight_value)
    val_loss_fn = create_loss_fn(y_all[val_idx], model_cfg, pos_weight_value=pos_weight_value)
    best_val = float("-inf")
    best_epoch = 0
    patience_counter = 0
    history = []
    early_stopping = bool(model_cfg.get("early_stopping", True))

    for epoch in range(1, int(model_cfg["epochs"]) + 1):
        model.train()
        logits = model(x).view(-1)
        train_loss = train_loss_fn(logits[train_idx])
        optimizer.zero_grad()
        train_loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            logits_eval = model(x).view(-1)
            val_scores = torch.sigmoid(logits_eval[val_idx]).cpu().numpy()
            val_metrics = compute_binary_metrics(
                bundle["y_all"][val_idx].astype(int),
                val_scores,
                (val_scores >= 0.5).astype(int),
            )
            val_loss = float(val_loss_fn(logits_eval[val_idx]).item())
        selection_value = _selection_metric(val_metrics)
        history.append(
            {
                "epoch": epoch,
                "train_loss": float(train_loss.item()),
                "val_loss": val_loss,
                "val_auroc": val_metrics["auroc"],
                "val_auprc": val_metrics["auprc"],
                "val_mcc": val_metrics["mcc"],
                "val_f1": val_metrics["f1"],
                "val_accuracy": val_metrics["accuracy"],
                "best_checkpoint_updated": selection_value > best_val,
            }
        )
        if selection_value > best_val:
            best_val = selection_value
            best_epoch = epoch
            _save_best_checkpoint(model, checkpoint_path, epoch, selection_value)
            patience_counter = 0
        else:
            patience_counter += 1
            if early_stopping and patience_counter >= int(model_cfg["patience"]):
                break

    checkpoint = torch.load(checkpoint_path, map_location="cpu")
    model.load_state_dict(checkpoint["state_dict"])
    with torch.no_grad():
        pred_score = torch.sigmoid(model(x).view(-1)).cpu().numpy()
    return model, pred_score, pd.DataFrame(history), {"best_epoch": int(best_epoch), "best_val_score": float(best_val)}


def fit_graph_model(
    bundle: dict,
    model_name: str,
    model_cfg: dict,
    seed: int,
    checkpoint_path,
):
    from src.models.epgat_gcn import EPGATOriginalGCN
    from src.models.epgat_gin import EPGATOriginalGIN
    from src.models.epgat_original import EPGATOriginalGAT
    from src.models.epgat_sage import EPGATOriginalSAGE, EPGATOriginalSAGEWithFusion

    set_seed(seed)
    x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
    edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
    y_all = torch.as_tensor(bundle["y_all"], dtype=torch.float32)
    train_idx = bundle["train_idx"]
    val_idx = bundle["val_idx"]

    normalized_model = normalize_model_name(model_name)
    if normalized_model == "GAT":
        model = EPGATOriginalGAT(
            in_feats=x.shape[1],
            h_feats=[12, 1] if bundle["species"] != "human" else [16, 1],
            heads=[8, 1],
            dropout=float(model_cfg["dropout"]),
            negative_slope=0.2,
        )
    elif normalized_model == "GCN":
        model = EPGATOriginalGCN(
            in_feats=x.shape[1],
            h_layers=[int(model_cfg["hidden_dim"]), 1],
            dropout=float(model_cfg["dropout"]),
        )
    elif normalized_model == "GIN":
        model = EPGATOriginalGIN(
            in_feats=x.shape[1],
            dim_h=int(model_cfg["dim_h"]),
            dropout=float(model_cfg["dropout"]),
        )
    elif normalized_model == "GraphSAGE":
        if is_gated_feature_setting(bundle["feature_setting"]):
            fusion_partition = dict(bundle.get("feature_metadata", {}).get("fusion_partition", {}))
            if not fusion_partition:
                raise ValueError(f"Missing fusion_partition metadata for gated feature_setting '{bundle['feature_setting']}'")
            model = EPGATOriginalSAGEWithFusion(
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
        else:
            model = EPGATOriginalSAGE(
                in_feats=x.shape[1],
                n_hidden=int(model_cfg["n_hidden"]),
                n_layers=int(model_cfg["n_layers"]),
                dropout=float(model_cfg["dropout"]),
                aggregator_type=str(model_cfg["aggregator_type"]),
            )
    else:
        raise ValueError(f"Unsupported graph model: {model_name}")

    optimizer = torch.optim.Adam(model.parameters(), lr=float(model_cfg["lr"]), weight_decay=float(model_cfg["weight_decay"]))
    loss_type, pos_weight_value = build_loss_metadata(bundle["y_all"][train_idx], model_cfg)
    train_loss_fn = create_loss_fn(y_all[train_idx], model_cfg, pos_weight_value=pos_weight_value)
    val_loss_fn = create_loss_fn(y_all[val_idx], model_cfg, pos_weight_value=pos_weight_value)
    best_val = float("-inf")
    best_epoch = 0
    patience_counter = 0
    history = []
    early_stopping = bool(model_cfg.get("early_stopping", True))

    for epoch in range(1, int(model_cfg["epochs"]) + 1):
        model.train()
        logits = model(x, edge_index).view(-1)
        train_loss = train_loss_fn(logits[train_idx])
        optimizer.zero_grad()
        train_loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            logits_eval = model(x, edge_index).view(-1)
            val_scores = torch.sigmoid(logits_eval[val_idx]).cpu().numpy()
            val_metrics = compute_binary_metrics(
                bundle["y_all"][val_idx].astype(int),
                val_scores,
                (val_scores >= 0.5).astype(int),
            )
            val_loss = float(val_loss_fn(logits_eval[val_idx]).item())
        selection_value = _selection_metric(val_metrics)
        history.append(
            {
                "epoch": epoch,
                "train_loss": float(train_loss.item()),
                "val_loss": val_loss,
                "val_auroc": val_metrics["auroc"],
                "val_auprc": val_metrics["auprc"],
                "val_mcc": val_metrics["mcc"],
                "val_f1": val_metrics["f1"],
                "val_accuracy": val_metrics["accuracy"],
                "best_checkpoint_updated": selection_value > best_val,
            }
        )
        if selection_value > best_val:
            best_val = selection_value
            best_epoch = epoch
            _save_best_checkpoint(model, checkpoint_path, epoch, selection_value)
            patience_counter = 0
        else:
            patience_counter += 1
            if early_stopping and patience_counter >= int(model_cfg["patience"]):
                break

    checkpoint = torch.load(checkpoint_path, map_location="cpu")
    model.load_state_dict(checkpoint["state_dict"])
    with torch.no_grad():
        pred_score = torch.sigmoid(model(x, edge_index).view(-1)).cpu().numpy()
    extra_info = {"best_epoch": int(best_epoch), "best_val_score": float(best_val)}
    extra_info["loss_type"] = loss_type
    extra_info["pos_weight"] = pos_weight_value if pos_weight_value is not None else ""
    if is_gated_feature_setting(bundle["feature_setting"]):
        extra_info["fusion_hidden_dim"] = int(model_cfg.get("fusion_hidden_dim", 256))
        extra_info["fusion_dropout"] = float(model_cfg.get("fusion_dropout", 0.2))
        extra_info["fusion_mode"] = str(model_cfg.get("fusion_mode", "gated"))
    return model, pred_score, pd.DataFrame(history), extra_info


def fit_sklearn_model(bundle, model_name, model_cfg, seed):
    x = bundle["feature_matrix"]
    y = bundle["y_all"].astype(int)
    train_idx = bundle["train_idx"]
    if model_name == "RF":
        model = RandomForestClassifier(
            n_estimators=int(model_cfg["n_estimators"]),
            class_weight=model_cfg.get("class_weight", "balanced"),
            n_jobs=int(model_cfg.get("n_jobs", -1)),
            random_state=seed,
        )
    elif model_name == "SVM":
        model = SVC(
            probability=True,
            C=float(model_cfg["c"]),
            kernel=str(model_cfg["kernel"]),
            gamma=model_cfg.get("gamma", "scale"),
            class_weight=model_cfg.get("class_weight", "balanced"),
            random_state=seed,
        )
    elif model_name == "NB":
        base_model = GaussianNB()
        train_y = y[train_idx]
        class_counts = np.bincount(train_y, minlength=2)
        min_class_count = int(class_counts.min()) if class_counts.size else 0
        if min_class_count >= 2:
            n_splits = min(3, min_class_count)
            cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
            model = CalibratedClassifierCV(base_model, method="sigmoid", cv=cv)
        else:
            model = base_model
    else:
        raise ValueError(f"Unsupported sklearn model: {model_name}")
    model.fit(x[train_idx], y[train_idx])
    pred_score = model.predict_proba(x)[:, 1]
    return model, pred_score


def fit_n2v_mlp(
    bundle: dict,
    model_cfg: dict,
    seed: int,
    output_dir,
    checkpoint_path,
):
    from src.graph.run_node2vec_embedding import train_node2vec_embeddings, write_node2vec_outputs

    node2vec_params = {
        "embedding_dim": int(model_cfg["embedding_dim"]),
        "walk_length": int(model_cfg["walk_length"]),
        "context_size": int(model_cfg["context_size"]),
        "walks_per_node": int(model_cfg["walks_per_node"]),
        "num_negative_samples": int(model_cfg["num_negative_samples"]),
        "epochs": int(model_cfg["epochs"]),
        "batch_size": int(model_cfg["batch_size"]),
        "learning_rate": float(model_cfg["learning_rate"]),
        "p": float(model_cfg.get("p", 1.0)),
        "q": float(model_cfg.get("q", 1.0)),
    }
    embeddings, n2v_history = train_node2vec_embeddings(
        edge_array=bundle["edge_index"],
        num_nodes=len(bundle["node_manifest"]),
        params=node2vec_params,
        seed=seed,
        device=str(bundle["config"]["runtime"].get("device", "cpu")),
        backend=str(bundle["config"]["runtime"].get("node2vec_backend", "auto")),
        require_true_node2vec=bool(bundle["config"]["runtime"].get("require_true_node2vec", False)),
    )
    backend_name = str(n2v_history["backend"].iloc[0]) if not n2v_history.empty and "backend" in n2v_history.columns else "unknown"
    write_node2vec_outputs(
        embeddings,
        node2vec_params,
        output_dir,
        n2v_history,
        metadata={
            "embedding_method": "node2vec_walk_based",
            "embedding_backend": backend_name,
            "fallback_used": False,
            "require_true_node2vec": bool(bundle["config"]["runtime"].get("require_true_node2vec", False)),
            "graph_contract": bundle["graph_contract"],
        },
    )
    mlp_bundle = dict(bundle)
    mlp_bundle["feature_matrix"] = embeddings.astype(np.float32, copy=False)
    mlp_cfg = {
        "hidden_dim": int(model_cfg["mlp_hidden_dim"]),
        "dropout": float(model_cfg["mlp_dropout"]),
        "learning_rate": float(model_cfg["mlp_learning_rate"]),
        "weight_decay": float(model_cfg["mlp_weight_decay"]),
        "epochs": int(model_cfg["mlp_epochs"]),
        "patience": int(model_cfg["mlp_patience"]),
        "early_stopping": bool(model_cfg.get("early_stopping", True)),
    }
    model, pred_score, training_log, best_info = fit_torch_mlp(mlp_bundle, mlp_cfg, seed, checkpoint_path)
    payload = {
        "node2vec_params": node2vec_params,
        "node2vec_backend": backend_name,
        "mlp_state_dict": model.state_dict(),
        "mlp_params": mlp_cfg,
    }
    return payload, pred_score, training_log, best_info


def run_network_heuristic(bundle, model_name):
    graph = nx.Graph()
    graph.add_nodes_from(range(len(bundle["node_manifest"])))
    graph.add_edges_from((int(src), int(dst)) for src, dst in bundle["edge_index"].tolist())
    if model_name == "DC":
        score_map = nx.degree_centrality(graph)
        score_name = "degree_centrality"
    elif model_name == "CC":
        score_map = nx.clustering(graph)
        score_name = "clustering_coefficient"
    else:
        raise ValueError(model_name)
    scores = np.array([float(score_map.get(idx, 0.0)) for idx in range(len(bundle["node_manifest"]))], dtype=np.float64)
    labeled = bundle["split_manifest"].copy()
    positive_total = int((labeled["label"].astype(int) == 1).sum())
    labeled_indices = np.array([bundle["mapping"][gene] for gene in labeled["graph_gene_id"].astype(str)], dtype=np.int64)
    pred_all = np.zeros(scores.shape[0], dtype=np.int64)
    ranked = labeled_indices[np.argsort(-scores[labeled_indices], kind="mergesort")]
    pred_all[ranked[:positive_total]] = 1
    return score_name, scores, pred_all


def save_pickle(payload, output_path):
    with output_path.open("wb") as handle:
        pickle.dump(payload, handle)


def run_benchmark_task(
    config_path,
    protocol_name,
    model_name,
    output_dir,
    seed=None,
    graph_contract=None,
    split_manifest=None,
    feature_setting_override=None,
):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    base_config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    requested_model_key = str(model_name).strip()
    model_key = normalize_model_name(requested_model_key)
    if requested_model_key not in base_config["models"]:
        raise ValueError(f"Unknown model config key '{requested_model_key}'")
    model_cfg = dict(base_config["models"][requested_model_key])
    feature_setting = str(feature_setting_override if feature_setting_override not in {None, ""} else model_cfg["feature_setting"]).strip()
    seed = int(seed) if seed is not None else int(base_config["runtime"]["split_seed"])
    threshold_strategy = "top_k_labeled_positive_count" if model_key in DETERMINISTIC_MODELS else "fixed_0.5"
    feature_contract_group = resolve_feature_contract_group(model_key, feature_setting)
    evaluation_contract = "auroc_auprc_mcc_specificity_fixed_0.5" if model_key not in DETERMINISTIC_MODELS else "auroc_auprc_mcc_specificity_top_k"
    runtime_graph_contract = graph_contract or base_config["runtime"].get("graph_contract", "undirected_symmetrized")

    bundle = load_protocol_dataset(
        config_path,
        protocol_name,
        feature_setting,
        graph_contract=runtime_graph_contract,
        split_manifest_override=split_manifest,
    )

    training_log = pd.DataFrame()
    best_info = {}
    model_payload = None
    feature_schema_to_write = bundle["feature_schema"]
    feature_dim = int(bundle["feature_matrix"].shape[1])
    checkpoint_path = output_dir / "best_model.pt"
    node2vec_backend = ""
    embedding_method = ""
    fallback_used = ""
    if model_key in GRAPH_MODELS:
        model, pred_score, training_log, best_info = fit_graph_model(bundle, model_key, model_cfg, seed, checkpoint_path)
        model_payload = {"model": model_key, "state_dict": model.state_dict(), "config": model_cfg}
    elif model_key == "MLP":
        model, pred_score, training_log, best_info = fit_torch_mlp(bundle, model_cfg, seed, checkpoint_path)
        model_payload = {"model": model_key, "state_dict": model.state_dict(), "config": model_cfg}
    elif model_key in {"RF", "SVM", "NB"}:
        model, pred_score = fit_sklearn_model(bundle, model_key, model_cfg, seed)
        model_payload = {"model": model_key, "model_object": model, "config": model_cfg}
    elif model_key == "N2V_MLP":
        model_payload, pred_score, training_log, best_info = fit_n2v_mlp(bundle, model_cfg, seed, output_dir, checkpoint_path)
        feature_schema_to_write = pd.read_csv(output_dir / "feature_schema.tsv", sep="\t")
        feature_dim = int(model_cfg["embedding_dim"])
        node2vec_summary = pd.read_csv(output_dir / "node2vec_summary.tsv", sep="\t").iloc[0]
        node2vec_backend = str(node2vec_summary["embedding_backend"])
        embedding_method = str(node2vec_summary["embedding_method"])
        fallback_used = str(node2vec_summary["fallback_used"])
    elif model_key in DETERMINISTIC_MODELS:
        score_name, pred_score, pred_label = run_network_heuristic(bundle, model_key)
        model_payload = {"model": model_key, "score_name": score_name}
    else:
        raise ValueError(model_key)

    if model_key not in DETERMINISTIC_MODELS:
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

    prediction_table = build_prediction_table(bundle, pred_score, pred_label, model_key)
    prediction_table.to_csv(output_dir / "predictions.tsv", sep="\t", index=False)
    feature_schema_to_write.to_csv(output_dir / "feature_schema.tsv", sep="\t", index=False)
    pd.DataFrame(bundle["edge_table"]).to_csv(output_dir / "edge_table.tsv", sep="\t", index=False)
    bundle["split_manifest"].to_csv(output_dir / "split_manifest.tsv", sep="\t", index=False)
    feature_metadata = dict(bundle.get("feature_metadata", {}))
    if "esm2_alignment_audit" in feature_metadata:
        feature_metadata["esm2_alignment_audit"].to_csv(output_dir / "esm2_alignment_audit.tsv", sep="\t", index=False)
    if "esm2_metadata" in feature_metadata:
        pd.DataFrame([feature_metadata["esm2_metadata"]]).to_csv(output_dir / "esm2_feature_metadata.tsv", sep="\t", index=False)
    if "model" in locals() and hasattr(model, "gate_statistics"):
        model.gate_statistics(
            torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32),
            bundle["node_manifest"],
            bundle["split_manifest"],
        ).to_csv(output_dir / "gate_statistics.tsv", sep="\t", index=False)

    metrics_row = {
        "protocol": bundle["protocol_name"],
        "species": bundle["species"],
        "regime": bundle["regime"],
        "model": model_key,
        "model_variant": requested_model_key,
        "feature_setting": feature_setting,
        "label_regime": bundle["label_regime"],
        "run_id": f"seed_{seed}" if seed is not None else ("deterministic" if model_key in DETERMINISTIC_MODELS else ""),
        "seed": seed if seed is not None else "",
        "is_deterministic": str(model_key in DETERMINISTIC_MODELS).lower(),
        "split_version": bundle["split_version"],
        "graph_source": bundle["graph_source"],
        "graph_contract": bundle["graph_contract"],
        "threshold_strategy": threshold_strategy,
        "evaluation_contract": evaluation_contract,
        "feature_contract_group": feature_contract_group,
        "label_manifest": bundle["label_manifest_path"],
        "split_manifest": bundle["split_manifest_path"],
        "config_used": str(Path(config_path)),
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
        "feature_dim": int(feature_dim),
        "node_count": int(len(bundle["node_manifest"])),
        "edge_count": int(bundle["edge_index"].shape[0]),
        "require_true_node2vec": str(bool(base_config["runtime"].get("require_true_node2vec", False))).lower(),
        "embedding_method": embedding_method,
        "embedding_backend": node2vec_backend,
        "fallback_used": fallback_used,
        "esm2_cache_path": str(feature_metadata.get("esm2_metadata", {}).get("cache_path", "")),
        "fusion_mode": str(model_cfg.get("fusion_mode", "")) if is_gated_feature_setting(feature_setting) else "",
        "fusion_hidden_dim": model_cfg.get("fusion_hidden_dim", "") if is_gated_feature_setting(feature_setting) else "",
        "fusion_dropout": model_cfg.get("fusion_dropout", "") if is_gated_feature_setting(feature_setting) else "",
        "loss_type": str(model_cfg.get("loss_type", "balanced_bce")),
        "pos_weight_mode": str(model_cfg.get("pos_weight_mode", "")),
        "pos_weight_scale": model_cfg.get("pos_weight_scale", ""),
        "best_checkpoint": str(checkpoint_path) if checkpoint_path.exists() else "",
        **best_info,
    }
    pd.DataFrame([metrics_row]).to_csv(output_dir / "metrics.tsv", sep="\t", index=False)
    if not training_log.empty:
        training_log.to_csv(output_dir / "training_log.tsv", sep="\t", index=False)
    save_pickle(model_payload, output_dir / "model.pkl")

    resolved = {
        "protocol": bundle["protocol_name"],
        "species": bundle["species"],
        "regime": bundle["regime"],
        "model": model_key,
        "model_variant": requested_model_key,
        "feature_setting": feature_setting,
        "seed": seed if seed is not None else None,
        "label_regime": bundle["label_regime"],
        "split_version": bundle["split_version"],
        "graph_source": bundle["graph_source"],
        "graph_contract": bundle["graph_contract"],
        "threshold_strategy": threshold_strategy,
        "evaluation_contract": evaluation_contract,
        "feature_contract_group": feature_contract_group,
        "label_manifest": bundle["label_manifest_path"],
        "split_manifest": bundle["split_manifest_path"],
        "output_dir": str(output_dir),
        "config_path": str(Path(config_path).resolve()),
        "metrics": metrics_row,
    }
    save_yaml(resolved, output_dir / "resolved_config.yaml")
    return metrics_row


def main() -> None:
    args = parse_args()
    run_benchmark_task(
        config_path=args.config,
        protocol_name=args.protocol,
        model_name=args.model,
        output_dir=args.output_dir,
        feature_setting_override=args.feature_setting,
        seed=args.seed,
        graph_contract=args.graph_contract,
        split_manifest=args.split_manifest,
    )


if __name__ == "__main__":
    main()
