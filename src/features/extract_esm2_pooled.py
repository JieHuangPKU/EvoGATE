from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any

import torch


MODEL_DIMENSIONS = {
    "facebook/esm2_t6_8M_UR50D": 320,
    "facebook/esm2_t12_35M_UR50D": 480,
    "facebook/esm2_t30_150M_UR50D": 640,
    "facebook/esm2_t33_650M_UR50D": 1280,
    "facebook/esm2_t36_3B_UR50D": 2560,
    "facebook/esm2_t48_15B_UR50D": 5120,
}

LOCAL_MODEL_REQUIRED_FILES = (
    "config.json",
    "tokenizer_config.json",
    "vocab.txt",
)
LOCAL_MODEL_WEIGHT_FILES = (
    "model.safetensors",
    "pytorch_model.bin",
)
INVALID_LOCAL_PATH_PARTS = {".locks", ".no_exist"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract mean-pooled ESM2 embeddings from a FASTA file")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA path")
    parser.add_argument("--output-pt", required=True, help="Output .pt path")
    parser.add_argument(
        "--model-name-or-path",
        default="facebook/esm2_t6_8M_UR50D",
        help="Hugging Face repo id or local model directory path",
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "mock", "esm2"],
        default="auto",
        help="Embedding backend. 'auto' attempts real ESM2 and falls back to mock.",
    )
    parser.add_argument(
        "--local-files-only",
        type=parse_bool,
        default=False,
        metavar="{true,false}",
        help="Require local Hugging Face files only",
    )
    parser.add_argument("--cache-dir", default=None, help="Explicit Hugging Face cache directory")
    parser.add_argument("--max-length", type=int, default=1024, help="Maximum amino-acid sequence length")
    parser.add_argument("--batch-size", type=int, default=4, help="Batch size for embedding extraction")
    parser.add_argument("--device", default="cpu", help="Torch device, e.g. cpu, cuda, mps")
    parser.add_argument(
        "--pooling",
        choices=["mean"],
        default="mean",
        help="Pooling strategy for per-sequence embeddings",
    )
    parser.add_argument(
        "--mock-embedding-dim",
        type=int,
        default=None,
        help="Override deterministic mock embedding dimension",
    )
    return parser.parse_args()


def parse_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if normalized in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean value, got: {value}")


def read_fasta(path: Path, max_length: int) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_chunks: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequence = "".join(current_chunks).upper()
                    sequences[current_id] = sequence[:max_length]
                current_id = line[1:].split()[0]
                if not current_id:
                    raise ValueError(f"Encountered FASTA header without an identifier in {path}")
                if current_id in sequences:
                    raise ValueError(f"Duplicate FASTA identifier detected: {current_id}")
                current_chunks = []
            else:
                if current_id is None:
                    raise ValueError(f"Invalid FASTA formatting in {path}: sequence before header")
                current_chunks.append(line)

    if current_id is not None:
        sequence = "".join(current_chunks).upper()
        sequences[current_id] = sequence[:max_length]

    if not sequences:
        raise ValueError(f"No sequences found in FASTA: {path}")
    return sequences


def infer_embedding_dim(model_name_or_path: str, mock_override: int | None) -> int:
    if mock_override is not None:
        return int(mock_override)
    candidate = str(model_name_or_path).rstrip("/").split("/")[-1]
    return int(MODEL_DIMENSIONS.get(model_name_or_path, MODEL_DIMENSIONS.get(candidate, 320)))


def looks_like_local_model_path(model_name_or_path: str) -> bool:
    candidate = str(model_name_or_path).strip()
    if not candidate:
        return False
    return candidate.startswith(("/", "./", "../", "~/")) or Path(candidate).exists()


def validate_local_model_dir(model_name_or_path: str) -> Path:
    model_dir = Path(model_name_or_path).expanduser().resolve(strict=False)
    invalid_parts = [part for part in model_dir.parts if part in INVALID_LOCAL_PATH_PARTS]
    if invalid_parts:
        raise ValueError(
            f"Invalid local Hugging Face model directory: {model_name_or_path}. "
            f"Path contains disallowed cache marker(s): {sorted(set(invalid_parts))}"
        )
    if not model_dir.exists():
        raise ValueError(f"Local Hugging Face model directory does not exist: {model_name_or_path}")
    if not model_dir.is_dir():
        raise ValueError(f"Local Hugging Face model path is not a directory: {model_name_or_path}")

    missing_files = [name for name in LOCAL_MODEL_REQUIRED_FILES if not (model_dir / name).is_file()]
    if missing_files:
        raise ValueError(
            f"Invalid local Hugging Face model directory: {model_name_or_path}. "
            f"Missing required file(s): {missing_files}"
        )

    if not any((model_dir / name).is_file() for name in LOCAL_MODEL_WEIGHT_FILES):
        raise ValueError(
            f"Invalid local Hugging Face model directory: {model_name_or_path}. "
            f"Expected one of weight file(s): {list(LOCAL_MODEL_WEIGHT_FILES)}"
        )

    return model_dir


def generate_mock_embedding(sequence_id: str, sequence: str, embedding_dim: int) -> torch.Tensor:
    digest = hashlib.sha256(f"{sequence_id}|{sequence}".encode("utf-8")).hexdigest()
    seed = int(digest[:16], 16) % (2**31)
    generator = torch.Generator(device="cpu")
    generator.manual_seed(seed)
    vector = torch.randn(embedding_dim, generator=generator, dtype=torch.float32)
    if sequence:
        vector = vector * (min(len(sequence), 2048) / max(len(sequence), 1)) ** 0.5
    return vector


def extract_mock_embeddings(
    sequences: dict[str, str],
    model_name_or_path: str,
    mock_embedding_dim: int | None,
) -> tuple[dict[str, torch.Tensor], dict[str, object]]:
    embedding_dim = infer_embedding_dim(model_name_or_path, mock_embedding_dim)
    embeddings = {
        sequence_id: generate_mock_embedding(sequence_id, sequence, embedding_dim)
        for sequence_id, sequence in sequences.items()
    }
    metadata = {"backend_used": "mock", "embedding_dim": embedding_dim}
    return embeddings, metadata


def extract_real_esm2_embeddings(
    sequences: dict[str, str],
    model_name_or_path: str,
    batch_size: int,
    max_length: int,
    device: str,
    local_files_only: bool,
    cache_dir: str | None,
    pooling: str,
) -> tuple[dict[str, torch.Tensor], dict[str, object]]:
    from transformers import AutoModel, AutoTokenizer

    resolved_model_name_or_path = model_name_or_path
    if looks_like_local_model_path(model_name_or_path):
        resolved_model_name_or_path = str(validate_local_model_dir(model_name_or_path))

    resolved_device = resolve_torch_device(device)
    pretrained_kwargs: dict[str, Any] = {"local_files_only": bool(local_files_only)}
    if cache_dir:
        pretrained_kwargs["cache_dir"] = cache_dir

    tokenizer = AutoTokenizer.from_pretrained(resolved_model_name_or_path, **pretrained_kwargs)
    model = AutoModel.from_pretrained(resolved_model_name_or_path, **pretrained_kwargs)
    model.to(resolved_device)
    model.eval()

    cls_token_id = tokenizer.cls_token_id
    eos_token_id = tokenizer.eos_token_id
    pad_token_id = tokenizer.pad_token_id

    items = list(sequences.items())
    n_sequences = len(items)
    n_batches = max(1, math.ceil(n_sequences / batch_size))
    pooled_embeddings: dict[str, torch.Tensor] = {}
    start_time = time.time()

    print(
        (
            f"[extract_esm2_pooled] starting real ESM2 extraction: "
            f"sequences={n_sequences}, batch_size={batch_size}, batches={n_batches}, "
            f"device={resolved_device}, model={resolved_model_name_or_path}"
        ),
        flush=True,
    )

    with torch.inference_mode():
        for start in range(0, len(items), batch_size):
            batch_index = (start // batch_size) + 1
            batch = items[start : start + batch_size]
            batch_ids = [item[0] for item in batch]
            batch_sequences = [item[1] for item in batch]

            encoded = tokenizer(
                batch_sequences,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max_length,
                add_special_tokens=True,
            )
            encoded = {key: value.to(resolved_device) for key, value in encoded.items()}
            outputs = model(**encoded)
            hidden = outputs.last_hidden_state

            attention_mask = encoded["attention_mask"].bool()
            token_mask = attention_mask.clone()
            if cls_token_id is not None:
                token_mask &= encoded["input_ids"] != cls_token_id
            if eos_token_id is not None:
                token_mask &= encoded["input_ids"] != eos_token_id
            if pad_token_id is not None:
                token_mask &= encoded["input_ids"] != pad_token_id

            if pooling != "mean":
                raise ValueError(f"Unsupported pooling strategy: {pooling}")

            token_mask = token_mask.unsqueeze(-1)
            masked_hidden = hidden * token_mask
            pooled = masked_hidden.sum(dim=1) / token_mask.sum(dim=1).clamp(min=1)

            if len(batch_ids) != len(pooled):
                raise RuntimeError(
                    f"Batch output size mismatch: {len(batch_ids)} ids but {len(pooled)} pooled embeddings"
                )

            for sequence_id, vector in zip(batch_ids, pooled):
                pooled_embeddings[sequence_id] = vector.detach().cpu().to(torch.float32)

            processed = len(pooled_embeddings)
            elapsed = time.time() - start_time
            print(
                (
                    f"[extract_esm2_pooled] batch {batch_index}/{n_batches} complete; "
                    f"processed={processed}/{n_sequences}; elapsed_sec={elapsed:.1f}"
                ),
                flush=True,
            )

    embedding_dim = int(next(iter(pooled_embeddings.values())).shape[0])
    return pooled_embeddings, {"backend_used": "esm2", "embedding_dim": embedding_dim}


def resolve_torch_device(device: str) -> torch.device:
    normalized = str(device).strip().lower()
    if normalized == "gpu":
        normalized = "cuda"
    if normalized == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA device requested but torch.cuda.is_available() is False")
    return torch.device(normalized)


def extract_embeddings(
    sequences: dict[str, str],
    model_name_or_path: str,
    backend: str,
    batch_size: int,
    max_length: int,
    device: str,
    mock_embedding_dim: int | None,
    local_files_only: bool,
    cache_dir: str | None,
    pooling: str,
) -> tuple[dict[str, torch.Tensor], dict[str, object]]:
    if backend == "mock":
        return extract_mock_embeddings(sequences, model_name_or_path, mock_embedding_dim)

    if backend in {"auto", "esm2"}:
        try:
            return extract_real_esm2_embeddings(
                sequences=sequences,
                model_name_or_path=model_name_or_path,
                batch_size=batch_size,
                max_length=max_length,
                device=device,
                local_files_only=local_files_only,
                cache_dir=cache_dir,
                pooling=pooling,
            )
        except Exception as exc:
            if backend == "esm2":
                raise
            embeddings, metadata = extract_mock_embeddings(sequences, model_name_or_path, mock_embedding_dim)
            metadata["fallback_reason"] = str(exc)
            return embeddings, metadata

    raise ValueError(f"Unsupported backend: {backend}")


def main() -> None:
    args = parse_args()
    input_fasta = Path(args.input_fasta)
    output_pt = Path(args.output_pt)
    output_pt.parent.mkdir(parents=True, exist_ok=True)

    sequences = read_fasta(input_fasta, max_length=int(args.max_length))
    embeddings, backend_metadata = extract_embeddings(
        sequences=sequences,
        model_name_or_path=args.model_name_or_path,
        backend=args.backend,
        batch_size=int(args.batch_size),
        max_length=int(args.max_length),
        device=args.device,
        mock_embedding_dim=args.mock_embedding_dim,
        local_files_only=bool(args.local_files_only),
        cache_dir=args.cache_dir,
        pooling=args.pooling,
    )

    payload = {
        "metadata": {
            "input_fasta": str(input_fasta),
            "model_name_or_path": args.model_name_or_path,
            "backend_requested": args.backend,
            "device": args.device,
            "local_files_only": bool(args.local_files_only),
            "cache_dir": args.cache_dir,
            "pooling": args.pooling,
            "max_length": int(args.max_length),
            "batch_size": int(args.batch_size),
            "n_sequences": int(len(sequences)),
            **backend_metadata,
        },
        "embeddings": embeddings,
    }
    # Keep the artifact readable by the older PyTorch build in the EPGAT training env.
    torch.save(payload, output_pt, _use_new_zipfile_serialization=False)

    print(f"[extract_esm2_pooled] finished writing embeddings to {output_pt}", flush=True)

    summary = {
        "output_pt": str(output_pt),
        "n_sequences": len(sequences),
        "model_name_or_path": args.model_name_or_path,
        "backend_requested": args.backend,
        "backend_used": payload["metadata"]["backend_used"],
        "embedding_dim": payload["metadata"]["embedding_dim"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
