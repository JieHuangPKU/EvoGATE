from pathlib import Path
from typing import Any, Dict, Tuple, Union

import numpy as np
import pandas as pd
import torch
import yaml

from src.features.epgat_legacy_features import append_degree_block, build_feature_schema, load_feature_table


TABULAR_FEATURE_SETTINGS = {"ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB"}
ESM2_FEATURE_SETTINGS = {
    "ESM2",
    "ORT_ESM2",
    "ORT_EXP_SUB_ESM2",
    "ORT_EXP_SUB_ESM2_GATED",
    "ORT_EXP_SUB_ESM2_OLD_GATED_WBCE",
    "ORT_EXP_SUB_ESM2_GATED_RESIDUAL",
    "ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE",
}
SUPPORTED_FEATURE_SETTINGS = TABULAR_FEATURE_SETTINGS | ESM2_FEATURE_SETTINGS | {"N2V", "NETWORK"}


def strip_feature_setting_variant(feature_setting: str) -> str:
    normalized = str(feature_setting).strip().upper()
    if "_OLD_GATED" in normalized:
        normalized = normalized.replace("_OLD_GATED", "")
    suffixes = ["_WBCE", "_RESIDUAL", "_GATED"]
    stripped = normalized
    changed = True
    while changed:
        changed = False
        for suffix in suffixes:
            if stripped.endswith(suffix):
                stripped = stripped[: -len(suffix)]
                changed = True
                break
    return stripped


def is_gated_feature_setting(feature_setting: str) -> bool:
    normalized = str(feature_setting).strip().upper()
    return "_GATED" in normalized


def normalize_graph_contract(graph_contract):
    normalized = str(graph_contract or "undirected_symmetrized").strip().lower()
    allowed = {"directed_raw", "undirected_symmetrized"}
    if normalized not in allowed:
        raise ValueError(f"Unsupported graph_contract '{graph_contract}'")
    return normalized


def load_config(config_path):
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def protocol_config(config, protocol_name):
    try:
        return dict(config["protocols"][protocol_name])
    except KeyError as exc:
        raise KeyError(f"Unknown frozen protocol '{protocol_name}'") from exc


def parse_feature_setting(feature_setting):
    normalized = str(feature_setting).strip().upper()
    if normalized not in SUPPORTED_FEATURE_SETTINGS:
        raise ValueError(f"Unsupported feature_setting '{feature_setting}'")
    normalized_base = strip_feature_setting_variant(normalized)
    if normalized_base in {"N2V", "NETWORK", "ESM2"}:
        return {"orthologs": False, "expression": False, "sublocalization": False}
    if normalized_base.endswith("_ESM2"):
        normalized_base = normalized_base.replace("_ESM2", "")
    tokens = {token for token in normalized_base.split("_") if token}
    unknown = sorted(tokens.difference({"ORT", "EXP", "SUB"}))
    if unknown:
        raise ValueError(f"Unsupported feature_setting '{feature_setting}' with unknown tokens {unknown}")
    return {
        "orthologs": "ORT" in tokens,
        "expression": "EXP" in tokens,
        "sublocalization": "SUB" in tokens,
    }


def label_manifest_path(config, protocol_name):
    return Path(config["paths"]["labels_dir"]) / str(protocol_config(config, protocol_name)["label_output"])


def split_manifest_path(config, protocol_name):
    return Path(config["paths"]["splits_dir"]) / str(protocol_config(config, protocol_name)["split_output"])


def load_frozen_label_manifest(config, protocol_name):
    path = label_manifest_path(config, protocol_name)
    if not path.exists():
        raise FileNotFoundError(f"Frozen label manifest not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def load_frozen_split_manifest(config, protocol_name):
    path = split_manifest_path(config, protocol_name)
    if not path.exists():
        raise FileNotFoundError(f"Frozen split manifest not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def load_split_manifest_from_path(path: Union[str, Path]) -> pd.DataFrame:
    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Split manifest not found: {resolved}")
    return pd.read_csv(resolved, sep="\t", dtype=str).fillna("")


def _species_feature_path(config, feature_key, data_key, filename):
    return Path(config["feature_roots"][feature_key]) / data_key / filename


def _esm2_cache_path(config, data_key):
    candidates = [
        Path(config["esm2"]["cache_root"]) / str(data_key) / "esm2_pooled.pt",
        Path("data/processed/phase1_esm2_mlp/species") / str(data_key) / "esm2_pooled.pt",
        Path("data/processed/ESM2-bk") / str(data_key) / "esm2_pooled.pt",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _resolve_graph_path(
    config: dict,
    data_key: str,
    graph_source_path: Union[str, Path, None] = None,
) -> Path:
    if graph_source_path not in {None, ""}:
        return Path(graph_source_path)
    return _species_feature_path(config, "ppi", data_key, "string.csv")


def load_string_graph(
    config,
    data_key,
    graph_contract=None,
    graph_source_path: Union[str, Path, None] = None,
    string_threshold: Union[int, float, None] = None,
):
    resolved_graph_contract = normalize_graph_contract(graph_contract or config["runtime"].get("graph_contract", "undirected_symmetrized"))
    ppi_path = _resolve_graph_path(config, data_key, graph_source_path)
    edge_df = pd.read_csv(ppi_path, dtype=str).fillna("")
    if ("A" not in edge_df.columns or "B" not in edge_df.columns) and ppi_path.suffix.lower() in {".tsv", ".txt"}:
        edge_df = pd.read_csv(ppi_path, sep="\t", dtype=str).fillna("")
    source_columns = [str(column) for column in edge_df.columns]
    effective_threshold = float(config["runtime"]["string_threshold"] if string_threshold is None else string_threshold)
    score_column = ""
    if "combined_score" in edge_df.columns:
        score_column = "combined_score"
        scores = pd.to_numeric(edge_df["combined_score"], errors="coerce").fillna(0.0)
        keep_mask = scores >= effective_threshold
        edge_df = edge_df.loc[keep_mask].copy()
        weights = scores.loc[keep_mask].to_numpy(dtype=np.float32) / 1000.0
    elif "edge_weight" in edge_df.columns:
        score_column = "edge_weight"
        weights = pd.to_numeric(edge_df["edge_weight"], errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
        effective_threshold = ""
    else:
        weights = None
        effective_threshold = ""

    edge_df = edge_df[["A", "B"]].copy()
    edge_df["A"] = edge_df["A"].astype(str).str.strip()
    edge_df["B"] = edge_df["B"].astype(str).str.strip()
    edge_df = edge_df[edge_df["A"].ne("") & edge_df["B"].ne("")].copy()
    edge_df = edge_df[edge_df["A"] != edge_df["B"]].drop_duplicates().reset_index(drop=True)
    if resolved_graph_contract == "undirected_symmetrized":
        reverse_df = edge_df.rename(columns={"A": "B", "B": "A"}).copy()
        edge_df = pd.concat([edge_df, reverse_df], ignore_index=True).drop_duplicates().reset_index(drop=True)
    if not bool(config["runtime"].get("use_edge_weights", False)):
        weights = None
    if weights is not None and resolved_graph_contract == "undirected_symmetrized":
        weights = np.concatenate([weights, weights]).astype(np.float32, copy=False)
    graph_metadata = {
        "source_columns": source_columns,
        "score_column": score_column,
        "source_row_count": int(len(edge_df)),
        "has_edge_weights": bool(weights is not None),
    }
    return edge_df, weights, resolved_graph_contract, str(ppi_path), effective_threshold, graph_metadata


def build_node_manifest(labels_df: pd.DataFrame, edge_df: pd.DataFrame) -> pd.DataFrame:
    graph_nodes = set(edge_df["A"].astype(str)) | set(edge_df["B"].astype(str))
    label_nodes = set(labels_df["graph_gene_id"].astype(str))
    node_ids = sorted(graph_nodes | label_nodes)

    label_lookup = labels_df.set_index("graph_gene_id")
    node_manifest = pd.DataFrame({"graph_gene_id": node_ids})
    node_manifest["canonical_gene_id"] = node_manifest["graph_gene_id"].map(label_lookup["canonical_gene_id"]).fillna(
        node_manifest["graph_gene_id"]
    )
    node_manifest["label"] = node_manifest["graph_gene_id"].map(label_lookup["label"]).fillna("")
    node_manifest["split"] = node_manifest["graph_gene_id"].map(label_lookup["split"]).fillna("")
    node_manifest["is_labeled"] = node_manifest["graph_gene_id"].isin(label_nodes)
    node_manifest["in_graph"] = node_manifest["graph_gene_id"].isin(graph_nodes)
    return node_manifest


def _selected_feature_table(config, data_key, flags):
    orth_df = (
        load_feature_table(_species_feature_path(config, "orthologs", data_key, "orthologs.csv"), "ortholog")
        if flags["orthologs"]
        else pd.DataFrame({"legacy_gene_id": []})
    )
    expr_df = (
        load_feature_table(_species_feature_path(config, "expression", data_key, "profile.csv"), "expression")
        if flags["expression"]
        else pd.DataFrame({"legacy_gene_id": []})
    )
    sub_df = (
        load_feature_table(_species_feature_path(config, "sublocalization", data_key, "subloc.csv"), "subloc")
        if flags["sublocalization"]
        else pd.DataFrame({"legacy_gene_id": []})
    )
    return orth_df, expr_df, sub_df


def _normalize_features(features: np.ndarray, train_idx: np.ndarray) -> np.ndarray:
    if features.size == 0:
        return features.astype(np.float32, copy=False)
    if train_idx.size == 0:
        raise ValueError("Frozen protocol train split is empty; cannot normalize features")
    train_block = features[train_idx]
    mean = train_block.mean(axis=0, keepdims=True)
    std = train_block.std(axis=0, keepdims=True)
    std[std < 1e-8] = 1.0
    return ((features - mean) / std).astype(np.float32, copy=False)


def _strip_species_prefix(value: str) -> str:
    text = str(value).strip()
    if "::" in text:
        return text.split("::", 1)[1].strip()
    return text


def _load_esm2_embedding_artifact(path: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Species-level ESM2 cache not found: {path}")
    load_path = str(path)
    try:
        payload = torch.load(load_path, map_location="cpu", weights_only=False)
    except TypeError:
        payload = torch.load(load_path, map_location="cpu")

    if isinstance(payload, dict) and "embeddings" in payload:
        embedding_map = payload["embeddings"]
        metadata = dict(payload.get("metadata", {}))
    elif isinstance(payload, dict):
        embedding_map = payload
        metadata = {}
    else:
        raise ValueError(f"Unsupported ESM2 cache payload structure in {path}")

    vectors = {}
    for key, value in embedding_map.items():
        vector = value.detach().cpu().numpy() if isinstance(value, torch.Tensor) else np.asarray(value)
        vectors[str(key)] = np.ravel(vector).astype(np.float32, copy=False)
    if not vectors:
        raise ValueError(f"No embeddings found in ESM2 cache: {path}")
    dims = {int(vector.shape[0]) for vector in vectors.values()}
    if len(dims) != 1:
        raise ValueError(f"ESM2 cache contains inconsistent embedding dimensions: {sorted(dims)}")
    metadata["embedding_dim"] = int(next(iter(dims)))
    return vectors, metadata


def _resolve_esm2_vectors(config, data_key, node_manifest: pd.DataFrame, train_idx: np.ndarray):
    cache_path = _esm2_cache_path(config, data_key)
    embedding_lookup, metadata = _load_esm2_embedding_artifact(cache_path)
    projection_dim = config.get("runtime", {}).get("esm2_projection_dim")
    projection_dim = int(projection_dim) if projection_dim not in {None, ""} else None

    audit_rows = []
    matched_vectors = []
    missing_rows = []
    for row in node_manifest[["graph_gene_id", "canonical_gene_id"]].itertuples(index=False):
        graph_gene_id = str(row.graph_gene_id).strip()
        canonical_gene_id = str(row.canonical_gene_id).strip()
        candidate_keys = []
        for candidate in [graph_gene_id, canonical_gene_id, _strip_species_prefix(canonical_gene_id)]:
            candidate = str(candidate).strip()
            if candidate and candidate not in candidate_keys:
                candidate_keys.append(candidate)
        resolved_key = next((candidate for candidate in candidate_keys if candidate in embedding_lookup), "")
        status = "matched" if resolved_key else "missing"
        if not resolved_key:
            missing_rows.append(
                {
                    "graph_gene_id": graph_gene_id,
                    "canonical_gene_id": canonical_gene_id,
                    "candidate_keys": "|".join(candidate_keys),
                }
            )
        else:
            matched_vectors.append(embedding_lookup[resolved_key])
        audit_rows.append(
            {
                "graph_gene_id": graph_gene_id,
                "canonical_gene_id": canonical_gene_id,
                "candidate_keys": "|".join(candidate_keys),
                "resolved_embedding_id": resolved_key,
                "mapping_status": status,
                "cache_path": str(cache_path),
            }
        )

    if missing_rows:
        examples = missing_rows[:10]
        raise FileNotFoundError(
            "ESM2 node feature alignment failed for some graph nodes. "
            f"cache={cache_path}; missing_count={len(missing_rows)}; examples={examples}"
        )

    matrix = np.vstack(matched_vectors).astype(np.float32, copy=False)
    original_embedding_dim = int(matrix.shape[1])
    if projection_dim is not None:
        if projection_dim < 1 or projection_dim > original_embedding_dim:
            raise ValueError(
                f"Invalid runtime.esm2_projection_dim={projection_dim}; available embedding_dim={original_embedding_dim}"
            )
        matrix = matrix[:, :projection_dim].astype(np.float32, copy=False)
    matrix = _normalize_features(matrix, train_idx)
    schema = pd.DataFrame(
        [
            {
                "feature_block": "esm2",
                "start_col": 0,
                "end_col": int(matrix.shape[1]) - 1,
                "dimension": int(matrix.shape[1]),
                "data_source": str(cache_path),
                "missing_strategy": "strict_id_alignment_no_silent_drop",
            }
        ]
    )
    metadata.update(
        {
            "cache_path": str(cache_path),
            "matched_nodes": int(len(node_manifest)),
            "missing_nodes": 0,
            "id_strategy": "graph_gene_id_then_canonical_gene_id_then_species_stripped_canonical",
            "original_embedding_dim": original_embedding_dim,
            "projection_dim": projection_dim if projection_dim is not None else original_embedding_dim,
            "embedding_dim": int(matrix.shape[1]),
        }
    )
    return matrix, schema, pd.DataFrame(audit_rows), metadata


def build_feature_bundle(
    config,
    data_key: str,
    feature_setting: str,
    node_manifest: pd.DataFrame,
    edge_df: pd.DataFrame,
    train_idx: np.ndarray,
):
    normalized = str(feature_setting).strip().upper()
    if normalized in {"N2V", "NETWORK"}:
        return np.zeros((len(node_manifest), 0), dtype=np.float32), pd.DataFrame(
            columns=["feature_block", "start_col", "end_col", "dimension", "data_source", "missing_strategy"]
        ), {}

    normalized_base = strip_feature_setting_variant(normalized)

    if normalized_base == "ESM2":
        matrix, schema, audit_df, metadata = _resolve_esm2_vectors(config, data_key, node_manifest, train_idx)
        return matrix, schema, {"esm2_alignment_audit": audit_df, "esm2_metadata": metadata}

    base_feature_setting = normalized_base.replace("_ESM2", "") if normalized_base.endswith("_ESM2") else normalized_base
    flags = parse_feature_setting(base_feature_setting)
    orth_df, expr_df, sub_df = _selected_feature_table(config, data_key, flags)
    feature_df = pd.DataFrame({"legacy_gene_id": node_manifest["graph_gene_id"].astype(str)})
    if flags["orthologs"]:
        feature_df = feature_df.merge(orth_df, on="legacy_gene_id", how="left")
    if flags["expression"]:
        feature_df = feature_df.merge(expr_df, on="legacy_gene_id", how="left")
    if flags["sublocalization"]:
        feature_df = feature_df.merge(sub_df, on="legacy_gene_id", how="left")
    feature_df = feature_df.fillna(0.0)
    include_degree = bool(config["runtime"].get("include_degree", False))
    if include_degree:
        feature_df = append_degree_block(feature_df, edge_df)

    feature_cols = [column for column in feature_df.columns if column != "legacy_gene_id"]
    matrix = feature_df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
    matrix = _normalize_features(matrix, train_idx)
    schema = build_feature_schema(expr_df, orth_df, sub_df, include_degree=include_degree)
    if normalized_base.endswith("_ESM2"):
        esm2_matrix, esm2_schema, audit_df, metadata = _resolve_esm2_vectors(config, data_key, node_manifest, train_idx)
        start_col = int(schema["end_col"].max()) + 1 if not schema.empty else 0
        esm2_schema = esm2_schema.copy()
        esm2_schema["start_col"] = start_col
        esm2_schema["end_col"] = start_col + int(esm2_matrix.shape[1]) - 1
        schema = pd.concat([schema, esm2_schema], ignore_index=True)
        matrix = np.concatenate([matrix, esm2_matrix], axis=1).astype(np.float32, copy=False)
        feature_metadata = {"esm2_alignment_audit": audit_df, "esm2_metadata": metadata}
        if is_gated_feature_setting(normalized):
            omics_rows = schema[schema["feature_block"] != "esm2"].copy()
            esm2_rows = schema[schema["feature_block"] == "esm2"].copy()
            if omics_rows.empty or esm2_rows.empty:
                raise ValueError(f"Gated feature setting '{feature_setting}' requires both omics and esm2 blocks")
            feature_metadata["fusion_partition"] = {
                "omics_start_col": int(omics_rows["start_col"].min()),
                "omics_end_col": int(omics_rows["end_col"].max()),
                "omics_dim": int(omics_rows["dimension"].sum()),
                "esm2_start_col": int(esm2_rows["start_col"].min()),
                "esm2_end_col": int(esm2_rows["end_col"].max()),
                "esm2_dim": int(esm2_rows["dimension"].sum()),
                "partition_source": "feature_schema",
            }
        return matrix, schema, feature_metadata
    return matrix, schema, {}


def load_protocol_dataset(
    config_path,
    protocol_name,
    feature_setting,
    graph_contract=None,
    graph_source_path: Union[str, Path, None] = None,
    string_threshold: Union[int, float, None] = None,
    graph_source_name: Union[str, None] = None,
    split_manifest_override: Union[str, Path, None] = None,
):
    config = load_config(config_path)
    cfg = protocol_config(config, protocol_name)

    labels_df = load_frozen_label_manifest(config, protocol_name)
    split_df = (
        load_split_manifest_from_path(split_manifest_override)
        if split_manifest_override not in {None, ""}
        else load_frozen_split_manifest(config, protocol_name)
    )
    split_df["label"] = pd.to_numeric(split_df["label"], errors="raise").astype(int)

    edge_df, edge_weights, resolved_graph_contract, resolved_graph_source, effective_threshold, graph_metadata = load_string_graph(
        config,
        cfg["data_key"],
        graph_contract=graph_contract,
        graph_source_path=graph_source_path,
        string_threshold=string_threshold,
    )
    node_manifest = build_node_manifest(split_df, edge_df)
    mapping = dict(zip(node_manifest["graph_gene_id"].astype(str), range(len(node_manifest))))

    edge_index = np.column_stack(
        [
            edge_df["A"].astype(str).map(mapping).to_numpy(dtype=np.int64),
            edge_df["B"].astype(str).map(mapping).to_numpy(dtype=np.int64),
        ]
    )

    train_idx = np.array(
        [mapping[gene] for gene in split_df.loc[split_df["split"] == "train", "graph_gene_id"].astype(str)],
        dtype=np.int64,
    )
    val_idx = np.array(
        [mapping[gene] for gene in split_df.loc[split_df["split"] == "val", "graph_gene_id"].astype(str)],
        dtype=np.int64,
    )
    test_idx = np.array(
        [mapping[gene] for gene in split_df.loc[split_df["split"] == "test", "graph_gene_id"].astype(str)],
        dtype=np.int64,
    )

    y_all = np.full(len(node_manifest), np.nan, dtype=np.float32)
    label_map = dict(zip(split_df["graph_gene_id"].astype(str), split_df["label"].astype(int)))
    for gene_id, idx in mapping.items():
        if gene_id in label_map:
            y_all[idx] = float(label_map[gene_id])

    feature_matrix, feature_schema, feature_metadata = build_feature_bundle(
        config=config,
        data_key=cfg["data_key"],
        feature_setting=feature_setting,
        node_manifest=node_manifest,
        edge_df=edge_df,
        train_idx=train_idx,
    )

    return {
        "config": config,
        "protocol_name": protocol_name,
        "species": cfg["species"],
        "regime": cfg["regime"],
        "data_key": cfg["data_key"],
        "feature_setting": str(feature_setting).strip().upper(),
        "label_manifest_path": str(label_manifest_path(config, protocol_name)),
        "split_manifest_path": str(Path(split_manifest_override).resolve()) if split_manifest_override not in {None, ""} else str(split_manifest_path(config, protocol_name)),
        "graph_source": resolved_graph_source,
        "graph_source_name": str(graph_source_name).strip() if graph_source_name not in {None, ""} else Path(resolved_graph_source).stem,
        "string_threshold": effective_threshold,
        "graph_contract": resolved_graph_contract,
        "graph_metadata": graph_metadata,
        "node_manifest": node_manifest,
        "label_manifest": labels_df,
        "split_manifest": split_df,
        "edge_table": edge_df,
        "edge_index": edge_index,
        "edge_weights": edge_weights,
        "mapping": mapping,
        "feature_matrix": feature_matrix,
        "feature_schema": feature_schema,
        "feature_metadata": feature_metadata,
        "train_idx": train_idx,
        "val_idx": val_idx,
        "test_idx": test_idx,
        "y_all": y_all,
        "label_regime": cfg["regime"],
        "split_version": str(split_df["split_version"].iloc[0]),
    }
