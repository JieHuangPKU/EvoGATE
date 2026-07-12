from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import torch


DEFAULT_ID_COLUMNS = ["gene_id", "protein_id", "sequence_id", "id"]
DEFAULT_LABEL_COLUMNS = ["label", "label_binary", "target", "y"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a tensor dataset from pooled embeddings and labels")
    parser.add_argument("--embeddings-pt", required=True, help="Input pooled embedding .pt artifact")
    parser.add_argument("--labels", required=True, help="Input labels TSV/CSV")
    parser.add_argument("--output-pt", required=True, help="Output dataset .pt artifact")
    parser.add_argument("--id-column", default=None, help="Label-table identifier column")
    parser.add_argument("--label-column", default=None, help="Label-table binary label column")
    parser.add_argument("--split-column", default=None, help="Optional split column to carry into the dataset artifact")
    parser.add_argument("--split-manifest", default=None, help="Optional frozen split manifest TSV")
    parser.add_argument("--split-id-column", default=None, help="Identifier column in the split manifest")
    return parser.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def detect_column(columns: list[str], explicit: str | None, candidates: list[str], kind: str) -> str:
    if explicit:
        if explicit not in columns:
            raise ValueError(f"Requested {kind} column '{explicit}' not found. Available columns: {columns}")
        return explicit
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(f"Could not infer {kind} column. Available columns: {columns}")


def normalize_binary_label(value: object) -> int:
    normalized = str(value).strip().lower()
    mapping = {
        "0": 0,
        "1": 1,
        "false": 0,
        "true": 1,
        "negative": 0,
        "positive": 1,
        "non-essential": 0,
        "essential": 1,
        "no": 0,
        "yes": 1,
    }
    if normalized not in mapping:
        raise ValueError(f"Unsupported binary label value: {value}")
    return mapping[normalized]


def load_pooled_embeddings(path: Path) -> tuple[dict[str, torch.Tensor], dict[str, object]]:
    try:
        payload = torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        payload = torch.load(path, map_location="cpu")

    if isinstance(payload, dict) and "embeddings" in payload:
        embeddings = payload["embeddings"]
        metadata = dict(payload.get("metadata", {}))
    elif isinstance(payload, dict):
        embeddings = payload
        metadata = {}
    else:
        raise ValueError(f"Unsupported embedding artifact structure in {path}")

    tensor_map: dict[str, torch.Tensor] = {}
    for key, value in embeddings.items():
        tensor = value if isinstance(value, torch.Tensor) else torch.as_tensor(value)
        tensor = tensor.detach().cpu().to(torch.float32).reshape(-1)
        tensor_map[str(key)] = tensor

    if not tensor_map:
        raise ValueError(f"No pooled embeddings were found in {path}")

    dimensions = {int(tensor.shape[0]) for tensor in tensor_map.values()}
    if len(dimensions) != 1:
        raise ValueError(f"Embedding dimension mismatch detected: {sorted(dimensions)}")

    metadata["feature_dim"] = int(next(iter(dimensions)))
    return tensor_map, metadata


def main() -> None:
    args = parse_args()
    embeddings_path = Path(args.embeddings_pt)
    labels_path = Path(args.labels)
    output_pt = Path(args.output_pt)
    output_pt.parent.mkdir(parents=True, exist_ok=True)

    embedding_map, embedding_metadata = load_pooled_embeddings(embeddings_path)
    label_df = load_table(labels_path).copy()
    label_df.columns = [str(column).strip() for column in label_df.columns]

    id_column = detect_column(label_df.columns.tolist(), args.id_column, DEFAULT_ID_COLUMNS, "ID")
    label_column = detect_column(label_df.columns.tolist(), args.label_column, DEFAULT_LABEL_COLUMNS, "label")
    split_column = args.split_column if args.split_column in label_df.columns else None

    if args.split_manifest:
        split_manifest_path = Path(args.split_manifest)
        split_df = load_table(split_manifest_path).copy()
        split_df.columns = [str(column).strip() for column in split_df.columns]
        split_id_column = detect_column(split_df.columns.tolist(), args.split_id_column, DEFAULT_ID_COLUMNS, "split ID")
        if not args.split_column:
            raise ValueError("--split-column is required when --split-manifest is provided")
        if args.split_column not in split_df.columns:
            raise ValueError(f"Split column '{args.split_column}' not found in split manifest {split_manifest_path}")
        split_df[split_id_column] = split_df[split_id_column].astype(str).str.strip()
        split_df = split_df[[split_id_column, args.split_column]].drop_duplicates(subset=[split_id_column], keep="first")
        label_df = label_df.merge(
            split_df,
            left_on=id_column,
            right_on=split_id_column,
            how="left",
            suffixes=("", "_split_manifest"),
        )
        if split_id_column != id_column and split_id_column in label_df.columns:
            label_df = label_df.drop(columns=[split_id_column])
        split_column = args.split_column

    label_df[id_column] = label_df[id_column].astype(str).str.strip()
    label_df = label_df[label_df[id_column].ne("")].copy()
    label_df["_label_int"] = label_df[label_column].map(normalize_binary_label)

    if split_column:
        missing_splits = int(label_df[split_column].isna().sum())
        if missing_splits > 0:
            raise ValueError(f"Detected {missing_splits} rows without split assignments after label/split alignment")

    embedding_ids = set(embedding_map)
    label_ids = set(label_df[id_column].tolist())
    matched_ids = [gene_id for gene_id in label_df[id_column].tolist() if gene_id in embedding_ids]
    missing_embedding_ids = sorted(label_ids - embedding_ids)
    unused_embedding_ids = sorted(embedding_ids - label_ids)

    if not matched_ids:
        raise ValueError("No overlapping IDs were found between pooled embeddings and label table")

    aligned_df = label_df[label_df[id_column].isin(matched_ids)].copy()
    aligned_df = aligned_df.drop_duplicates(subset=[id_column], keep="first").reset_index(drop=True)
    vectors = [embedding_map[gene_id] for gene_id in aligned_df[id_column]]

    dataset_payload = {
        "X": torch.stack(vectors).to(torch.float32),
        "y": torch.tensor(aligned_df["_label_int"].tolist(), dtype=torch.long),
        "ids": aligned_df[id_column].tolist(),
        "splits": aligned_df[split_column].astype(str).tolist() if split_column else None,
        "table": aligned_df.drop(columns=["_label_int"]).to_dict(orient="records"),
        "metadata": {
            "embeddings_pt": str(embeddings_path),
            "labels_path": str(labels_path),
            "id_column": id_column,
            "label_column": label_column,
            "split_column": split_column,
            "n_labels_total": int(len(label_df)),
            "n_embeddings_total": int(len(embedding_map)),
            "n_matched": int(len(aligned_df)),
            "n_missing_embeddings": int(len(missing_embedding_ids)),
            "n_unused_embeddings": int(len(unused_embedding_ids)),
            "missing_embedding_ids": missing_embedding_ids[:100],
            "unused_embedding_ids": unused_embedding_ids[:100],
            **embedding_metadata,
        },
    }
    torch.save(dataset_payload, output_pt)

    print(
        json.dumps(
            {
                "output_pt": str(output_pt),
                "n_matched": len(aligned_df),
                "n_missing_embeddings": len(missing_embedding_ids),
                "n_unused_embeddings": len(unused_embedding_ids),
                "feature_dim": dataset_payload["metadata"]["feature_dim"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
