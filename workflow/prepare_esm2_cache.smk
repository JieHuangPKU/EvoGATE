configfile: "configs/prepare_esm2_cache.yaml"

from pathlib import Path

ESM2 = config["esm2"]
ESM2_PYTHON_BIN = str(ESM2.get("python_bin", "python"))
ESM2_CACHE_ROOT = Path(ESM2["cache_root"])
ESM2_DATA_KEYS = sorted(str(key) for key in ESM2["protein_fastas"])


def esm2_cache_path_for_data_key(data_key):
    return str(ESM2_CACHE_ROOT / str(data_key) / "esm2_pooled.pt")


def esm2_fasta_path_for_data_key(data_key):
    try:
        return str(ESM2["protein_fastas"][data_key])
    except KeyError as exc:
        raise KeyError(f"Missing esm2.protein_fastas entry for data_key '{data_key}'") from exc


wildcard_constraints:
    data_key="|".join(ESM2_DATA_KEYS)


rule all:
    default_target: True
    input:
        expand(str(ESM2_CACHE_ROOT / "{data_key}" / "esm2_pooled.pt"), data_key=ESM2_DATA_KEYS),


rule extract_esm2_pooled_embeddings:
    input:
        fasta=lambda wc: esm2_fasta_path_for_data_key(wc.data_key),
    output:
        pooled=str(ESM2_CACHE_ROOT / "{data_key}" / "esm2_pooled.pt"),
    params:
        output_dir=lambda wc: str(Path(esm2_cache_path_for_data_key(wc.data_key)).parent),
        cache_dir_arg=lambda wc: f'--cache-dir "{ESM2["cache_dir"]}"' if ESM2.get("cache_dir") else "",
        model_name_or_path=str(ESM2["model_name_or_path"]),
        backend=str(ESM2["backend"]),
        local_files_only=str(ESM2["local_files_only"]),
        max_length=int(ESM2["max_length"]),
        batch_size=int(ESM2["batch_size"]),
        device=str(ESM2["device"]),
        pooling=str(ESM2["pooling"]),
        mock_embedding_dim=int(ESM2["mock_embedding_dim"]),
    log:
        str(ESM2_CACHE_ROOT / "{data_key}" / "extract_esm2_pooled.log"),
    shell:
        r"""
        mkdir -p "{params.output_dir}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{ESM2_PYTHON_BIN}" -m src.features.extract_esm2_pooled \
          --input-fasta "{input.fasta}" \
          --output-pt "{output.pooled}" \
          --model-name-or-path "{params.model_name_or_path}" \
          --backend "{params.backend}" \
          --local-files-only "{params.local_files_only}" \
          {params.cache_dir_arg} \
          --max-length {params.max_length} \
          --batch-size {params.batch_size} \
          --device "{params.device}" \
          --pooling "{params.pooling}" \
          --mock-embedding-dim {params.mock_embedding_dim} \
          > "{log}" 2>&1
        """
