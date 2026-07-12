"""
Train the minimal EPGAT original-compatible GAT inside ProGATE_v2.

This is a Phase 1 compatibility entrypoint and intentionally serves only the
legacy dataset builder outputs.
"""

import argparse
import os

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import yaml
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, matthews_corrcoef, roc_auc_score

from src.models.epgat_original import EPGATOriginalGAT


LEGACY_PARAM_PRESETS = {
    "fgraminearum": {"lr": 0.005, "weight_decay": 2e-4, "h_feats": [12, 1], "heads": [8, 1], "dropout": 0.3, "negative_slope": 0.2},
    "scerevisiae": {"lr": 0.005, "weight_decay": 2e-4, "h_feats": [12, 1], "heads": [8, 1], "dropout": 0.3, "negative_slope": 0.2},
    "human": {"lr": 0.005, "weight_decay": 5e-4, "h_feats": [16, 1], "heads": [8, 1], "dropout": 0.4, "negative_slope": 0.2},
    "celegans": {"lr": 0.005, "weight_decay": 5e-4, "h_feats": [12, 1], "heads": [8, 1], "dropout": 0.3, "negative_slope": 0.2},
}


def parse_args():
    parser = argparse.ArgumentParser(description="Train EPGAT original-compatible legacy GAT")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)


class LegacyBalancedLoss(object):
    def __init__(self, y, idx):
        idx = np.asarray(idx)
        self.y = y
        self.pos = idx[y.cpu() == 1]
        self.neg = idx[y.cpu() == 0]
        self.y_pos = y[y == 1]
        self.y_neg = y[y == 0]

    def __call__(self, out):
        loss_p = F.binary_cross_entropy_with_logits(out[self.pos].squeeze(), self.y_pos)
        loss_n = F.binary_cross_entropy_with_logits(out[self.neg].squeeze(), self.y_neg)
        return loss_p + loss_n


def legacy_auc(y_true, scores):
    return roc_auc_score(y_true.cpu().numpy(), torch.sigmoid(scores).detach().cpu().numpy())


def train_from_config(config):
    dataset_dir = os.path.join(config["paths"]["output_root"], config["run"]["name"])
    run_out = dataset_dir
    feature_matrix = np.load(os.path.join(dataset_dir, "feature_matrix.npy"))
    edge_pairs = np.load(os.path.join(dataset_dir, "edge_index.npy"), allow_pickle=True)
    label_manifest = pd.read_csv(os.path.join(dataset_dir, "label_manifest.tsv"), sep="\t", dtype=str).fillna("")
    node_manifest = pd.read_csv(os.path.join(dataset_dir, "node_manifest.tsv"), sep="\t", dtype=str).fillna("")

    mapping = dict(zip(node_manifest["legacy_gene_id"].astype(str), range(len(node_manifest))))
    if np.issubdtype(edge_pairs.dtype, np.integer):
        edge_index = edge_pairs.astype(np.int64, copy=False)
    else:
        edge_index = np.vectorize(mapping.__getitem__)(edge_pairs).astype(np.int64)
    edge_index = torch.from_numpy(edge_index.T).to(torch.long).contiguous()

    x = torch.from_numpy(feature_matrix).to(torch.float32)
    labeled = label_manifest[label_manifest["is_labeled"].astype(str).str.lower().isin(["true", "1", "yes"])].copy()
    train_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "train"]["legacy_gene_id"].astype(str)])
    val_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "val"]["legacy_gene_id"].astype(str)])
    test_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "test"]["legacy_gene_id"].astype(str)])

    label_map = dict(
        zip(
            labeled["legacy_gene_id"].astype(str),
            labeled["label"].astype(float).astype(int),
        )
    )
    y_all = np.zeros(len(node_manifest), dtype=np.float32)
    for gene_id, idx in mapping.items():
        if gene_id in label_map:
            y_all[idx] = float(label_map[gene_id])
    train_y = torch.tensor(y_all[train_idx], dtype=torch.float32)
    val_y = torch.tensor(y_all[val_idx], dtype=torch.float32)
    test_y = torch.tensor(y_all[test_idx], dtype=torch.float32)

    params = dict(LEGACY_PARAM_PRESETS[config["legacy"]["organism"]])
    params["dropout"] = float(config["train"]["dropout"])
    set_seed(int(config["legacy"]["seed"]))

    model = EPGATOriginalGAT(in_feats=x.shape[1], h_feats=params["h_feats"], heads=params["heads"], dropout=params["dropout"], negative_slope=params["negative_slope"])
    optimizer = torch.optim.Adam(model.parameters(), lr=params["lr"], weight_decay=params["weight_decay"])
    train_loss_fn = LegacyBalancedLoss(train_y, train_idx)
    val_loss_fn = LegacyBalancedLoss(val_y, val_idx)

    best_state = None
    best_val_auc = -1.0
    epochs = int(config["train"]["epochs"])
    epoch_rows = []

    for epoch in range(epochs):
        model.train()
        logits = model(x, edge_index)
        optimizer.zero_grad()
        loss = train_loss_fn(logits)
        loss.backward()
        optimizer.step()

        logits_detached = logits.detach()
        train_auc = legacy_auc(train_y, logits_detached[train_idx])
        val_auc = legacy_auc(val_y, logits_detached[val_idx])
        val_loss = val_loss_fn(logits_detached)
        epoch_rows.append(
            {
                "epoch": epoch + 1,
                "train_loss": float(loss.item()),
                "val_loss": float(val_loss.item()),
                "train_auc": float(train_auc),
                "val_auc": float(val_auc),
            }
        )
        if val_auc > best_val_auc:
            best_val_auc = val_auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

    reload_best_state = bool(config.get("train", {}).get("reload_best_state", True))
    if best_state is not None and reload_best_state:
        model.load_state_dict(best_state)

    model.eval()
    with torch.no_grad():
        logits = model(x, edge_index)
        probs = torch.sigmoid(logits).cpu().numpy().reshape(-1)

    test_probs = probs[test_idx]
    test_pred = (test_probs >= 0.5).astype(int)
    test_true = test_y.cpu().numpy().astype(int)
    metrics = {
        "auroc": float(roc_auc_score(test_true, test_probs)),
        "auprc": float(average_precision_score(test_true, test_probs)),
        "accuracy": float(accuracy_score(test_true, test_pred)),
        "f1": float(f1_score(test_true, test_pred)),
        "mcc": float(matthews_corrcoef(test_true, test_pred)),
        "best_val_auc": float(best_val_auc),
        "test_count": int(len(test_idx)),
    }

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "legacy_params": params,
            "best_val_auc": best_val_auc,
        },
        os.path.join(run_out, "checkpoint.pt"),
    )
    pd.DataFrame(epoch_rows).to_csv(os.path.join(run_out, "training_log.tsv"), sep="\t", index=False)
    pd.DataFrame([metrics]).to_csv(os.path.join(run_out, "metrics.tsv"), sep="\t", index=False)

    predictions = node_manifest[["legacy_gene_id"]].copy()
    predictions["pred_score"] = probs
    predictions["split"] = predictions["legacy_gene_id"].map(dict(zip(label_manifest["legacy_gene_id"], label_manifest["split"]))).fillna("")
    predictions["label"] = predictions["legacy_gene_id"].map(dict(zip(label_manifest["legacy_gene_id"], label_manifest["label"]))).fillna("")
    predictions.to_csv(os.path.join(run_out, "predictions.tsv"), sep="\t", index=False)

    with open(os.path.join(run_out, "run_config_frozen.yaml"), "w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)

    legacy_export = pd.DataFrame(
        [
            {
                "Method": "GAT",
                "Organism": config["legacy"]["organism"],
                "PPI": config["legacy"]["ppi"],
                "Expression": bool(config["legacy"]["expression"]),
                "Orthologs": bool(config["legacy"]["orthologs"]),
                "Sublocalization": bool(config["legacy"]["sublocalization"]),
                "AUC": metrics["auroc"],
                "AUPR": metrics["auprc"],
                "F1": metrics["f1"],
                "Accuracy": metrics["accuracy"],
                "MCC": metrics["mcc"],
            }
        ]
    )
    legacy_export.to_csv(os.path.join(run_out, "legacy_GAT_results_compatible.csv"), index=False)

    return {
        "run_out": run_out,
        "metrics": metrics,
        "predictions": predictions,
        "epoch_rows": epoch_rows,
    }


def main():
    args = parse_args()
    config = load_yaml(args.config)
    result = train_from_config(config)

    print("Legacy GAT training complete:", result["run_out"])
    print(result["metrics"])


if __name__ == "__main__":
    main()
