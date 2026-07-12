"""
Run one fixed Phase 2A species x model x run task.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import yaml
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, matthews_corrcoef, roc_auc_score

from src.data.build_epgat_legacy_dataset import build_dataset
from src.models.epgat_gcn import EPGATOriginalGCN
from src.models.epgat_gin import EPGATOriginalGIN
from src.models.epgat_original import EPGATOriginalGAT
from src.models.epgat_sage import EPGATOriginalSAGE


SPECIES_DIR = {
    "human": "human",
    "celegans": "celegans",
    "scerevisiae": "yeast",
    "fgraminearum": "fgraminearum",
}

MODEL_NAME_MAP = {
    "gat": "gat",
    "gcn": "gcn",
    "gin": "gin",
    "sage": "sage",
    "graphsage": "sage",
}

DISPLAY_MODEL_MAP = {
    "gat": "GAT",
    "gcn": "GCN",
    "gin": "GIN",
    "sage": "GraphSAGE",
}
FEATURE_COMBOS = ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB"]

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
        self.pos = idx[y.cpu() == 1]
        self.neg = idx[y.cpu() == 0]
        self.y_pos = y[y == 1]
        self.y_neg = y[y == 0]

    def __call__(self, out):
        out = out.squeeze()
        return F.binary_cross_entropy_with_logits(out[self.pos], self.y_pos) + F.binary_cross_entropy_with_logits(out[self.neg], self.y_neg)


def parse_args():
    parser = argparse.ArgumentParser(description="Run one fixed Phase 2A task")
    parser.add_argument("--species", required=True, choices=sorted(SPECIES_DIR.keys()))
    parser.add_argument("--model", required=True, type=str)
    parser.add_argument("--feature_combo", required=True, choices=FEATURE_COMBOS)
    parser.add_argument("--string_thr", required=True, type=int)
    parser.add_argument("--run_id", required=True, type=int)
    parser.add_argument("--include_degree", required=True, type=str)
    parser.add_argument("--output_dir", required=True, type=str)
    parser.add_argument("--positive_set_path", default="", type=str)
    parser.add_argument("--negative_set_path", default="", type=str)
    parser.add_argument("--label_regime", default="old440", type=str)
    parser.add_argument("--legacy_epgat_root", default="/home/jiehuang/software/fungi/EPGAT", type=str)
    parser.add_argument("--base_seed", default=1029, type=int)
    parser.add_argument("--epochs", default=100, type=int)
    return parser.parse_args()


def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)


def _normalize_model(model_name):
    normalized = str(model_name).strip().lower()
    if normalized not in MODEL_NAME_MAP:
        raise ValueError("Unsupported model: {}".format(model_name))
    return MODEL_NAME_MAP[normalized]


def _bool_text(value):
    return str(value).strip().lower() in ["true", "1", "yes"]


def _feature_combo_flags(feature_combo):
    tokens = set([token for token in str(feature_combo).strip().upper().split("_") if token])
    if not tokens:
        raise ValueError("Empty feature_combo is not allowed")
    unknown = sorted(tokens.difference({"ORT", "EXP", "SUB"}))
    if unknown:
        raise ValueError("Unsupported tokens in feature_combo {}: {}".format(feature_combo, ",".join(unknown)))
    normalized = "_".join([token for token in ["ORT", "EXP", "SUB"] if token in tokens])
    if normalized not in FEATURE_COMBOS:
        raise ValueError("Unsupported feature_combo: {}".format(feature_combo))
    return normalized, {
        "orthologs": "ORT" in tokens,
        "expression": "EXP" in tokens,
        "sublocalization": "SUB" in tokens,
    }


def _build_model(model_name, species, in_feats):
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


def _train_one(run_dir, species, model_name, seed, epochs):
    feature_matrix = np.load(os.path.join(run_dir, "feature_matrix.npy"))
    edge_pairs = np.load(os.path.join(run_dir, "edge_index.npy"), allow_pickle=True)
    label_manifest = pd.read_csv(os.path.join(run_dir, "label_manifest.tsv"), sep="\t", dtype=str).fillna("")
    node_manifest = pd.read_csv(os.path.join(run_dir, "node_manifest.tsv"), sep="\t", dtype=str).fillna("")

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

    label_map = dict(zip(labeled["legacy_gene_id"].astype(str), labeled["label"].astype(float).astype(int)))
    y_all = np.zeros(len(node_manifest), dtype=np.float32)
    for gene_id, idx in mapping.items():
        if gene_id in label_map:
            y_all[idx] = float(label_map[gene_id])
    train_y = torch.tensor(y_all[train_idx], dtype=torch.float32)
    val_y = torch.tensor(y_all[val_idx], dtype=torch.float32)
    test_y = torch.tensor(y_all[test_idx], dtype=torch.float32)

    set_seed(seed)
    model, model_params, optim_params = _build_model(model_name, species, x.shape[1])
    optimizer = torch.optim.Adam(model.parameters(), lr=optim_params["lr"], weight_decay=optim_params["weight_decay"])
    train_loss_fn = LegacyBalancedLoss(train_y, train_idx)
    val_loss_fn = LegacyBalancedLoss(val_y, val_idx)

    best_state = None
    best_val_auc = -1.0
    epoch_rows = []
    for epoch in range(int(epochs)):
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

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "model_name": model_name,
            "model_params": model_params,
            "optim_params": optim_params,
            "seed": int(seed),
            "best_val_auc": float(best_val_auc),
        },
        os.path.join(run_dir, "checkpoint.pt"),
    )
    return metrics, model_params, optim_params


def _require_files(base_dir, filenames, stage_name):
    missing = []
    for filename in filenames:
        path = os.path.join(base_dir, filename)
        if not os.path.exists(path):
            missing.append(path)
    if missing:
        raise RuntimeError(
            "{} did not produce required files under output_dir={}: {}".format(
                stage_name,
                base_dir,
                ", ".join(missing),
            )
        )


def main():
    args = parse_args()
    model_name = _normalize_model(args.model)
    include_degree = _bool_text(args.include_degree)
    feature_combo, feature_flags = _feature_combo_flags(args.feature_combo)
    output_dir = os.path.abspath(args.output_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    seed = int(args.base_seed) + int(args.run_id)
    positive_set_path = str(args.positive_set_path).strip()
    negative_set_path = str(args.negative_set_path).strip()

    builder_config = {
        "paths": {
            "legacy_epgat_root": args.legacy_epgat_root,
            "output_root": os.path.dirname(output_dir),
        },
        "run": {
            "name": os.path.basename(output_dir),
        },
        "legacy": {
            "organism": args.species,
            "legacy_species_dir": SPECIES_DIR[args.species],
            "ppi": "string",
            "expression": bool(feature_flags["expression"]),
            "orthologs": bool(feature_flags["orthologs"]),
            "sublocalization": bool(feature_flags["sublocalization"]),
            "include_degree": bool(include_degree),
            "string_threshold": int(args.string_thr),
            "use_weights": False,
            "positive_set_path": positive_set_path,
            "negative_set_path": negative_set_path,
            "label_regime": str(args.label_regime),
            "seed": seed,
            "test_fraction": 0.20,
            "val_fraction": 0.04,
        },
        "train": {
            "epochs": int(args.epochs),
            "dropout": GAT_PARAMS.get(args.species, GAT_PARAMS["human"])["dropout"],
            "reload_best_state": True,
        },
    }
    dataset_dir = os.path.abspath(build_dataset(builder_config))
    if dataset_dir != output_dir:
        raise RuntimeError(
            "Dataset builder wrote to an unexpected directory: expected {}, got {}".format(
                output_dir,
                dataset_dir,
            )
        )
    _require_files(
        output_dir,
        [
            "feature_matrix.npy",
            "edge_index.npy",
            "label_manifest.tsv",
            "node_manifest.tsv",
            "feature_schema.tsv",
        ],
        "build_dataset",
    )

    metrics, model_params, optim_params = _train_one(output_dir, args.species, model_name, seed, args.epochs)

    feature_schema = pd.read_csv(os.path.join(output_dir, "feature_schema.tsv"), sep="\t")
    feature_summary = pd.DataFrame(
        [
            {
                "species": args.species,
                "model": DISPLAY_MODEL_MAP[model_name],
                "feature_combo": feature_combo,
                "string_thr": int(args.string_thr),
                "include_degree": bool(include_degree),
                "run_id": int(args.run_id),
                "seed": int(seed),
                "feature_dim": int(np.load(os.path.join(output_dir, "feature_matrix.npy")).shape[1]),
                "feature_blocks": "|".join(feature_schema["feature_block"].astype(str).tolist()) if not feature_schema.empty else "",
                "orthologs_enabled": bool(feature_flags["orthologs"]),
                "expression_enabled": bool(feature_flags["expression"]),
                "sublocalization_enabled": bool(feature_flags["sublocalization"]),
                "label_regime": str(args.label_regime),
            }
        ]
    )
    feature_summary.to_csv(os.path.join(output_dir, "feature_summary.tsv"), sep="\t", index=False)

    resolved_config = {
        "species": args.species,
        "model": DISPLAY_MODEL_MAP[model_name],
        "model_internal": model_name,
        "feature_combo": feature_combo,
        "string_thr": int(args.string_thr),
        "include_degree": bool(include_degree),
        "feature_flags": feature_flags,
        "label_regime": str(args.label_regime),
        "positive_set_path": positive_set_path,
        "negative_set_path": negative_set_path,
        "run_id": int(args.run_id),
        "seed": int(seed),
        "output_dir": output_dir,
        "legacy_epgat_root": args.legacy_epgat_root,
        "dataset_config": builder_config,
        "model_params": model_params,
        "optim_params": optim_params,
        "metrics": metrics,
    }
    with open(os.path.join(output_dir, "resolved_config.yaml"), "w", encoding="utf-8") as handle:
        yaml.safe_dump(resolved_config, handle, sort_keys=False)
    _require_files(
        output_dir,
        [
            "metrics.tsv",
            "resolved_config.yaml",
            "feature_schema.tsv",
            "feature_summary.tsv",
        ],
        "run_epgat_phase2a_single",
    )


if __name__ == "__main__":
    main()
