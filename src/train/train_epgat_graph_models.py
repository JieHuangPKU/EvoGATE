"""
Archived replay graph-model benchmark trainer.
"""

import argparse
import os

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import yaml
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, matthews_corrcoef, roc_auc_score

from src.models.epgat_gcn import EPGATOriginalGCN
from src.models.epgat_gin import EPGATOriginalGIN
from src.models.epgat_original import EPGATOriginalGAT
from src.models.epgat_sage import EPGATOriginalSAGE


GAT_PARAMS = {
    "fgraminearum": {"lr": 0.005, "weight_decay": 2e-4, "dropout": 0.3, "h_feats": [12, 1], "heads": [8, 1], "negative_slope": 0.2},
    "scerevisiae": {"lr": 0.005, "weight_decay": 2e-4, "dropout": 0.3, "h_feats": [12, 1], "heads": [8, 1], "negative_slope": 0.2},
    "human": {"lr": 0.005, "weight_decay": 5e-4, "dropout": 0.4, "h_feats": [16, 1], "heads": [8, 1], "negative_slope": 0.2},
    "celegans": {"lr": 0.005, "weight_decay": 5e-4, "dropout": 0.3, "h_feats": [12, 1], "heads": [8, 1], "negative_slope": 0.2},
}
GCN_PARAMS = {
    "default": {"lr": 1e-3, "weight_decay": 1e-4, "hidden_dim": 64, "dropout": 0.5},
}
GIN_PARAMS = {
    "default": {"lr": 0.005, "weight_decay": 5e-4, "dim_h": 32, "dropout": 0.5},
    "human": {"lr": 0.001, "weight_decay": 5e-4, "dim_h": 64, "dropout": 0.5},
    "scerevisiae": {"lr": 0.005, "weight_decay": 1e-4, "dim_h": 64, "dropout": 0.4},
    "fgraminearum": {"lr": 0.005, "weight_decay": 1e-4, "dim_h": 64, "dropout": 0.4},
}
SAGE_PARAMS = {
    "default": {"lr": 0.005, "weight_decay": 5e-4, "n_hidden": 64, "n_layers": 2, "dropout": 0.5, "aggregator_type": "mean"},
    "scerevisiae": {"lr": 0.005, "weight_decay": 1e-4, "n_hidden": 64, "n_layers": 2, "dropout": 0.4, "aggregator_type": "pool"},
    "fgraminearum": {"lr": 0.005, "weight_decay": 1e-4, "n_hidden": 64, "n_layers": 2, "dropout": 0.4, "aggregator_type": "pool"},
}


class LegacyBalancedLoss(object):
    def __init__(self, y, idx):
        idx = np.asarray(idx)
        self.y = y
        self.pos = idx[y.cpu() == 1]
        self.neg = idx[y.cpu() == 0]
        self.y_pos = y[y == 1]
        self.y_neg = y[y == 0]

    def __call__(self, out):
        out = out.squeeze()
        return F.binary_cross_entropy_with_logits(out[self.pos], self.y_pos) + F.binary_cross_entropy_with_logits(out[self.neg], self.y_neg)


def parse_args():
    parser = argparse.ArgumentParser(description="Train graph benchmark model")
    parser.add_argument("--config", required=True, type=str)
    parser.add_argument("--model", required=True, choices=["gat", "gcn", "gin", "sage"])
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)


def build_model(model_name, species, in_feats):
    if model_name == "gat":
        params = dict(GAT_PARAMS[species])
        model = EPGATOriginalGAT(in_feats=in_feats, h_feats=params["h_feats"], heads=params["heads"], dropout=params["dropout"], negative_slope=params["negative_slope"])
        optim_params = {"lr": params["lr"], "weight_decay": params["weight_decay"]}
    elif model_name == "gcn":
        params = dict(GCN_PARAMS["default"])
        model = EPGATOriginalGCN(in_feats=in_feats, h_layers=[params["hidden_dim"], 1], dropout=params["dropout"])
        optim_params = {"lr": params["lr"], "weight_decay": params["weight_decay"]}
    elif model_name == "gin":
        params = dict(GIN_PARAMS.get(species, GIN_PARAMS["default"]))
        model = EPGATOriginalGIN(in_feats=in_feats, dim_h=params["dim_h"], dropout=params["dropout"])
        optim_params = {"lr": params["lr"], "weight_decay": params["weight_decay"]}
    elif model_name == "sage":
        params = dict(SAGE_PARAMS.get(species, SAGE_PARAMS["default"]))
        model = EPGATOriginalSAGE(in_feats=in_feats, n_hidden=params["n_hidden"], n_layers=params["n_layers"], dropout=params["dropout"], aggregator_type=params["aggregator_type"])
        optim_params = {"lr": params["lr"], "weight_decay": params["weight_decay"]}
    else:
        raise ValueError(model_name)
    return model, params, optim_params


def main():
    args = parse_args()
    config = load_yaml(args.config)
    dataset_dir = config["dataset"]["source_dir"]
    run_dir = os.path.join(config["paths"]["output_root"], config["benchmark"]["species"], args.model)
    os.makedirs(run_dir, exist_ok=True)

    feature_matrix = np.load(os.path.join(dataset_dir, "feature_matrix.npy"))
    edge_index = np.load(os.path.join(dataset_dir, "edge_index.npy"), allow_pickle=True)
    node_manifest = pd.read_csv(os.path.join(dataset_dir, "node_manifest.tsv"), sep="\t", dtype=str).fillna("")
    label_manifest = pd.read_csv(os.path.join(dataset_dir, "label_manifest.tsv"), sep="\t", dtype=str).fillna("")
    feature_schema = pd.read_csv(os.path.join(dataset_dir, "feature_schema.tsv"), sep="\t")
    dataset_audit = pd.read_csv(os.path.join(dataset_dir, "dataset_alignment_audit.tsv"), sep="\t")

    node_manifest.to_csv(os.path.join(run_dir, "node_manifest.tsv"), sep="\t", index=False)
    label_manifest.to_csv(os.path.join(run_dir, "label_manifest.tsv"), sep="\t", index=False)
    feature_schema.to_csv(os.path.join(run_dir, "feature_schema.tsv"), sep="\t", index=False)
    dataset_audit.to_csv(os.path.join(run_dir, "dataset_alignment_audit.tsv"), sep="\t", index=False)

    mapping = dict(zip(node_manifest["legacy_gene_id"].astype(str), range(len(node_manifest))))
    if np.issubdtype(edge_index.dtype, np.integer):
        edge_index = edge_index.astype(np.int64, copy=False)
    else:
        edge_index = np.vectorize(mapping.__getitem__)(edge_index).astype(np.int64)
    edge_index = torch.from_numpy(edge_index.T).to(torch.long).contiguous()
    x = torch.from_numpy(feature_matrix).to(torch.float32)

    labeled = label_manifest[label_manifest["is_labeled"].astype(str).str.lower().isin(["true", "1", "yes"])].copy()
    train_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "train"]["legacy_gene_id"].astype(str)])
    val_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "val"]["legacy_gene_id"].astype(str)])
    test_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "test"]["legacy_gene_id"].astype(str)])
    label_map = dict(zip(labeled["legacy_gene_id"].astype(str), labeled["label"].astype(float).astype(int)))
    y_all = np.zeros(len(node_manifest), dtype=np.float32)
    for gene_id, idx in mapping.items():
        if gene_id in label_map:
            y_all[idx] = float(label_map[gene_id])
    train_y = torch.tensor(y_all[train_idx], dtype=torch.float32)
    val_y = torch.tensor(y_all[val_idx], dtype=torch.float32)
    test_y = torch.tensor(y_all[test_idx], dtype=torch.float32)

    set_seed(int(config["benchmark"]["seed"]))
    model, model_params, optim_params = build_model(args.model, config["benchmark"]["species"], x.shape[1])
    optimizer = torch.optim.Adam(model.parameters(), lr=optim_params["lr"], weight_decay=optim_params["weight_decay"])
    train_loss_fn = LegacyBalancedLoss(train_y, train_idx)
    val_loss_fn = LegacyBalancedLoss(val_y, val_idx)

    best_val_auc = -1.0
    best_state = None
    epoch_rows = []
    for epoch in range(int(config["benchmark"]["epochs"])):
        model.train()
        logits = model(x, edge_index)
        optimizer.zero_grad()
        loss = train_loss_fn(logits)
        loss.backward()
        optimizer.step()

        logits_detached = logits.detach().squeeze()
        train_auc = roc_auc_score(train_y.cpu().numpy(), torch.sigmoid(logits_detached[train_idx]).cpu().numpy())
        val_auc = roc_auc_score(val_y.cpu().numpy(), torch.sigmoid(logits_detached[val_idx]).cpu().numpy())
        val_loss = val_loss_fn(logits_detached)
        epoch_rows.append({"epoch": epoch + 1, "train_loss": float(loss.item()), "val_loss": float(val_loss.item()), "train_auc": float(train_auc), "val_auc": float(val_auc)})
        if val_auc > best_val_auc:
            best_val_auc = val_auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)

    model.eval()
    with torch.no_grad():
        logits = model(x, edge_index).squeeze()
        probs = torch.sigmoid(logits).cpu().numpy()
    test_probs = probs[test_idx]
    test_true = test_y.cpu().numpy().astype(int)
    test_pred = (test_probs >= 0.5).astype(int)
    metrics = {
        "auroc": float(roc_auc_score(test_true, test_probs)),
        "auprc": float(average_precision_score(test_true, test_probs)),
        "accuracy": float(accuracy_score(test_true, test_pred)),
        "f1": float(f1_score(test_true, test_pred)),
        "mcc": float(matthews_corrcoef(test_true, test_pred)),
        "best_val_auc": float(best_val_auc),
        "test_count": int(len(test_idx)),
    }

    pd.DataFrame(epoch_rows).to_csv(os.path.join(run_dir, "training_log.tsv"), sep="\t", index=False)
    pd.DataFrame([metrics]).to_csv(os.path.join(run_dir, "metrics.tsv"), sep="\t", index=False)
    predictions = node_manifest[["legacy_gene_id"]].copy()
    predictions["pred_score"] = probs
    predictions["split"] = predictions["legacy_gene_id"].map(dict(zip(label_manifest["legacy_gene_id"], label_manifest["split"]))).fillna("")
    predictions["label"] = predictions["legacy_gene_id"].map(dict(zip(label_manifest["legacy_gene_id"], label_manifest["label"]))).fillna("")
    predictions.to_csv(os.path.join(run_dir, "predictions.tsv"), sep="\t", index=False)

    frozen = {"config": config, "model": args.model, "model_params": model_params, "optim_params": optim_params}
    with open(os.path.join(run_dir, "run_config_frozen.yaml"), "w", encoding="utf-8") as handle:
        yaml.safe_dump(frozen, handle, sort_keys=False)

    num_positive = int((labeled["label"].astype(float).astype(int) == 1).sum())
    num_negative = int((labeled["label"].astype(float).astype(int) == 0).sum())
    lines = [
        "# Graph Benchmark Run Summary",
        "",
        "- species: {}".format(config["benchmark"]["species"]),
        "- model: {}".format(args.model),
        "- num_nodes_final: {}".format(len(node_manifest)),
        "- num_labeled: {}".format(len(labeled)),
        "- num_positive: {}".format(num_positive),
        "- num_negative: {}".format(num_negative),
        "- auroc: {:.4f}".format(metrics["auroc"]),
        "- auprc: {:.4f}".format(metrics["auprc"]),
        "- accuracy: {:.4f}".format(metrics["accuracy"]),
        "- f1: {:.4f}".format(metrics["f1"]),
        "- mcc: {:.4f}".format(metrics["mcc"]),
    ]
    with open(os.path.join(run_dir, "run_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    print("Graph benchmark run complete:", run_dir)
    print(metrics)


if __name__ == "__main__":
    main()
