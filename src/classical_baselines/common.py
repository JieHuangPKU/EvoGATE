import random
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, matthews_corrcoef, roc_auc_score

from src.data.build_epgat_legacy_dataset import build_dataset


FEATURE_SETTINGS = ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB", "N2V", "network"]
TABULAR_FEATURE_SETTINGS = FEATURE_SETTINGS[:7]
TRAINABLE_METHODS = {"MLP", "RF", "SVM", "NB", "N2V_MLP"}
HEURISTIC_METHODS = {"DC", "CC"}


def load_yaml(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def dump_yaml(payload, path):
    with Path(path).open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def parse_feature_setting(feature_setting):
    feature_setting = str(feature_setting).strip().upper()
    if feature_setting in {"N2V", "NETWORK"}:
        return {"orthologs": False, "expression": False, "sublocalization": False}
    tokens = {token for token in feature_setting.split("_") if token}
    unknown = sorted(tokens.difference({"ORT", "EXP", "SUB"}))
    if unknown:
        raise ValueError(f"Unsupported feature_setting '{feature_setting}'; unknown tokens: {unknown}")
    normalized = "_".join(token for token in ["ORT", "EXP", "SUB"] if token in tokens)
    if normalized not in TABULAR_FEATURE_SETTINGS:
        raise ValueError(f"Unsupported feature_setting '{feature_setting}'")
    return {
        "orthologs": "ORT" in tokens,
        "expression": "EXP" in tokens,
        "sublocalization": "SUB" in tokens,
    }


def resolve_species_config(config, species):
    species_cfg = dict(config["species"][species])
    species_cfg["species"] = species
    return species_cfg


def build_dataset_for_benchmark(
    config,
    species,
    feature_setting,
    output_dir,
    seed,
):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    species_cfg = resolve_species_config(config, species)
    feature_flags = parse_feature_setting(feature_setting)
    positive_set_path = ""
    negative_set_path = ""
    label_regime = species_cfg["label_regime"]
    if label_regime == "new_label":
        positive_set_path = str(config["paths"]["new_label_positive_path"])
        negative_set_path = str(config["paths"]["new_label_negative_path"])

    builder_config = {
        "paths": {
            "legacy_epgat_root": config["paths"]["legacy_epgat_root"],
            "output_root": str(output_dir.parent),
        },
        "run": {"name": output_dir.name},
        "legacy": {
            "organism": species,
            "legacy_species_dir": species_cfg["legacy_species_dir"],
            "ppi": "string",
            "expression": bool(feature_flags["expression"]),
            "orthologs": bool(feature_flags["orthologs"]),
            "sublocalization": bool(feature_flags["sublocalization"]),
            "include_degree": bool(config["runtime"].get("include_degree", False)),
            "string_threshold": int(config["runtime"]["string_threshold"]),
            "use_weights": False,
            "positive_set_path": positive_set_path,
            "negative_set_path": negative_set_path,
            "label_regime": label_regime,
            "seed": int(seed),
            "test_fraction": float(config["runtime"]["test_fraction"]),
            "val_fraction": float(config["runtime"]["val_fraction"]),
        },
        "train": {
            "reload_best_state": True,
        },
    }
    dataset_dir = Path(build_dataset(builder_config)).resolve()
    if dataset_dir != output_dir:
        raise RuntimeError(f"Dataset builder wrote to {dataset_dir}, expected {output_dir}")

    return {
        "dataset_dir": str(dataset_dir),
        "builder_config": builder_config,
        "feature_flags": feature_flags,
        "label_regime": label_regime,
        "positive_set_path": positive_set_path,
        "negative_set_path": negative_set_path,
    }


def resolve_edge_pairs(edge_pairs, node_manifest):
    if edge_pairs.ndim != 2:
        raise ValueError(f"Expected 2D edge array, got shape {edge_pairs.shape}")
    if edge_pairs.shape[1] != 2 and edge_pairs.shape[0] == 2:
        edge_pairs = edge_pairs.T
    if edge_pairs.shape[1] != 2:
        raise ValueError(f"Expected edge array with two columns, got shape {edge_pairs.shape}")
    if np.issubdtype(edge_pairs.dtype, np.integer):
        return edge_pairs.astype(np.int64, copy=False)
    mapping = dict(zip(node_manifest["legacy_gene_id"].astype(str), range(len(node_manifest))))
    source = np.vectorize(mapping.__getitem__)(edge_pairs[:, 0]).astype(np.int64)
    target = np.vectorize(mapping.__getitem__)(edge_pairs[:, 1]).astype(np.int64)
    return np.column_stack([source, target])


def load_dataset_bundle(dataset_dir):
    dataset_dir = Path(dataset_dir)
    node_manifest = pd.read_csv(dataset_dir / "node_manifest.tsv", sep="\t", dtype=str).fillna("")
    label_manifest = pd.read_csv(dataset_dir / "label_manifest.tsv", sep="\t", dtype=str).fillna("")
    feature_schema_path = dataset_dir / "feature_schema.tsv"
    if feature_schema_path.exists() and feature_schema_path.stat().st_size > 0:
        try:
            feature_schema = pd.read_csv(feature_schema_path, sep="\t")
        except pd.errors.EmptyDataError:
            feature_schema = pd.DataFrame(columns=["feature_block", "start_col", "end_col", "dimension", "data_source", "missing_strategy"])
    else:
        feature_schema = pd.DataFrame(columns=["feature_block", "start_col", "end_col", "dimension", "data_source", "missing_strategy"])
    feature_matrix = np.load(dataset_dir / "feature_matrix.npy")
    edge_pairs = np.load(dataset_dir / "edge_index.npy", allow_pickle=True)
    edge_index = resolve_edge_pairs(edge_pairs, node_manifest)

    label_manifest["is_labeled_flag"] = label_manifest["is_labeled"].astype(str).str.lower().isin(["true", "1", "yes"])
    label_manifest["label_numeric"] = pd.to_numeric(label_manifest["label"], errors="coerce")

    mapping = dict(zip(node_manifest["legacy_gene_id"].astype(str), range(len(node_manifest))))
    labeled = label_manifest[label_manifest["is_labeled_flag"]].copy()
    train_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "train"]["legacy_gene_id"].astype(str)], dtype=np.int64)
    val_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "val"]["legacy_gene_id"].astype(str)], dtype=np.int64)
    test_idx = np.array([mapping[g] for g in labeled[labeled["split"] == "test"]["legacy_gene_id"].astype(str)], dtype=np.int64)

    y_all = np.full(len(node_manifest), np.nan, dtype=np.float32)
    for _, row in labeled.iterrows():
        y_all[mapping[str(row["legacy_gene_id"])]] = float(row["label_numeric"])

    return {
        "dataset_dir": str(dataset_dir),
        "node_manifest": node_manifest,
        "label_manifest": label_manifest,
        "feature_schema": feature_schema,
        "feature_matrix": feature_matrix,
        "edge_pairs": edge_pairs,
        "edge_index": edge_index,
        "mapping": mapping,
        "labeled": labeled,
        "train_idx": train_idx,
        "val_idx": val_idx,
        "test_idx": test_idx,
        "y_all": y_all,
    }


def compute_specificity(y_true, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    negatives = y_true == 0
    if negatives.sum() == 0:
        return float("nan")
    true_negative = int(((y_pred == 0) & negatives).sum())
    false_positive = int(((y_pred == 1) & negatives).sum())
    denom = true_negative + false_positive
    if denom == 0:
        return float("nan")
    return float(true_negative / denom)


def compute_binary_metrics(y_true, y_score, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    y_pred = np.asarray(y_pred).astype(int)
    metrics = {
        "auroc": float("nan"),
        "auprc": float("nan"),
        "mcc": float("nan"),
        "f1": float("nan"),
        "accuracy": float("nan"),
        "specificity": float("nan"),
    }
    if y_true.size == 0:
        return metrics
    if np.unique(y_true).size > 1:
        metrics["auroc"] = float(roc_auc_score(y_true, y_score))
        metrics["auprc"] = float(average_precision_score(y_true, y_score))
        metrics["mcc"] = float(matthews_corrcoef(y_true, y_pred))
    metrics["f1"] = float(f1_score(y_true, y_pred, zero_division=0))
    metrics["accuracy"] = float(accuracy_score(y_true, y_pred))
    metrics["specificity"] = compute_specificity(y_true, y_pred)
    return metrics


def build_prediction_table(
    node_manifest,
    label_manifest,
    pred_score,
    pred_label,
    method,
    feature_setting,
):
    predictions = node_manifest[["legacy_gene_id"]].copy()
    predictions["pred_score"] = np.asarray(pred_score, dtype=float)
    predictions["pred_label"] = np.asarray(pred_label).astype(int)
    label_lookup = label_manifest.set_index("legacy_gene_id")
    predictions["split"] = predictions["legacy_gene_id"].map(label_lookup["split"]).fillna("")
    predictions["label"] = predictions["legacy_gene_id"].map(label_lookup["label"]).fillna("")
    predictions["is_labeled"] = predictions["legacy_gene_id"].map(label_lookup["is_labeled"]).fillna("False")
    predictions["method"] = method
    predictions["feature_setting"] = feature_setting
    return predictions
