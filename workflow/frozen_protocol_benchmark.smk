configfile: "configs/frozen_protocol.yaml"

from pathlib import Path

config_path = workflow.overwrite_configfiles[0] if workflow.overwrite_configfiles else "configs/frozen_protocol.yaml"

PROTOCOLS = config["workflow"]["protocol_order"]
FIGURE1_TRAINABLE_DEFAULTS = [
    "MLP",
    "SVM",
    "NB",
    "N2V_MLP",
    "GAT",
    "GCN",
    "GIN",
    "GraphSAGE_ORT_EXP_SUB",
]
TRAINABLE_MODELS = config["workflow"].get("figure1_trainable_models", FIGURE1_TRAINABLE_DEFAULTS)
DETERMINISTIC_MODELS = config["workflow"]["deterministic_models"]
SEEDS = config["runtime"]["seed_list"]
PYTHON_BIN = config["runtime"]["python_bin"]
EPGAT_PYTHON = config["runtime"]["epgat_python_bin"]
MPLCONFIGDIR = config["runtime"].get("mplconfigdir", ".mplconfig")
XDG_CACHE_HOME = config["runtime"].get("xdg_cache_home", ".cache")
GRAPH_CONTRACT = config["runtime"]["graph_contract"]
ESM2 = config.get("esm2", {})
ESM2_PYTHON_BIN = ESM2.get("python_bin", PYTHON_BIN)

FIGURE1_RESULTS_ROOT = config["paths"].get("figure1_results_root", "results/Figure1")
FIGURE1_LABELS_DIR = config["paths"].get("figure1_labels_dir", f"{FIGURE1_RESULTS_ROOT}/labels")
FIGURE1_SPLITS_DIR = config["paths"].get("figure1_splits_dir", f"{FIGURE1_RESULTS_ROOT}/splits")
SUMMARY_DIR = config["paths"].get("figure1_summary_dir", f"{FIGURE1_RESULTS_ROOT}/summary")
PLOTS_DIR = config["paths"].get("figure1_plots_dir", f"{FIGURE1_RESULTS_ROOT}/plots")
OUTPUT_ROOT = config["paths"].get("figure1_benchmark_output_root", "outputs/Figure1")
OUTPUT_PREFIX = "Figure1"
GROUP_CONFIG_PATH = "workflow/figure1_plot_groups.json"

LABEL_SUMMARY = FIGURE1_LABELS_DIR + "/label_protocol_summary.md"
FUSARIUM_LABEL_SUMMARY = FIGURE1_LABELS_DIR + "/fgraminearum_label_protocol_summary.md"
SPLIT_SUMMARY = FIGURE1_SPLITS_DIR + "/split_protocol_summary.md"
FUSARIUM_SPLIT_SUMMARY = FIGURE1_SPLITS_DIR + "/fgraminearum_split_protocol_summary.md"
LEGACY_LABEL_SUMMARY = config["paths"]["labels_dir"] + "/label_protocol_summary.md"
LEGACY_FUSARIUM_LABEL_SUMMARY = config["paths"]["labels_dir"] + "/fgraminearum_label_protocol_summary.md"
LEGACY_SPLIT_SUMMARY = config["paths"]["splits_dir"] + "/split_protocol_summary.md"
LEGACY_FUSARIUM_SPLIT_SUMMARY = config["paths"]["splits_dir"] + "/fgraminearum_split_protocol_summary.md"
LABEL_SOURCE_FILES = [config["label_sources"][source_key] for source_key in sorted(config["label_sources"])]


def feature_setting_for(model_name):
    return config["models"][model_name]["feature_setting"]


def summary_model_name_for(model_name):
    if model_name == "GraphSAGE_ORT_EXP_SUB":
        return "GraphSAGE"
    return model_name


FEATURE_SETTINGS = sorted({feature_setting_for(model) for model in config["models"]})


def protocol_data_key(protocol_name):
    return str(config["protocols"][protocol_name]["data_key"])


def model_uses_esm2(model_name):
    return str(feature_setting_for(model_name)).strip().upper() in {"ESM2", "ORT_EXP_SUB_ESM2"}


def esm2_cache_path_for_data_key(data_key):
    return str(Path(ESM2["cache_root"]) / str(data_key) / "esm2_pooled.pt")


def esm2_fasta_path_for_data_key(data_key):
    try:
        return str(ESM2["protein_fastas"][data_key])
    except KeyError as exc:
        raise KeyError(f"Missing esm2.protein_fastas entry for data_key '{data_key}'") from exc


def required_esm2_cache_for_run(protocol_name, model_name):
    if not model_uses_esm2(model_name):
        return []
    return [esm2_cache_path_for_data_key(protocol_data_key(protocol_name))]


ESM2_DATA_KEYS = sorted({protocol_data_key(protocol) for protocol in PROTOCOLS})


def run_metrics_path(protocol, model, seed):
    return f"{OUTPUT_ROOT}/{protocol}/{model}/{feature_setting_for(model)}/run_{seed}/metrics.tsv"


def deterministic_metrics_path(protocol, model):
    return f"{OUTPUT_ROOT}/{protocol}/{model}/{feature_setting_for(model)}/deterministic/metrics.tsv"


TRAINABLE_OUTPUTS = [
    run_metrics_path(protocol, model, seed)
    for protocol in PROTOCOLS
    for model in TRAINABLE_MODELS
    for seed in SEEDS
]

DETERMINISTIC_OUTPUTS = [
    deterministic_metrics_path(protocol, model)
    for protocol in PROTOCOLS
    for model in DETERMINISTIC_MODELS
]

REQUIRED_SUMMARY_COMBOS = ",".join(
    sorted(
        {
            f"{summary_model_name_for(model)}:{feature_setting_for(model)}"
            for model in TRAINABLE_MODELS + DETERMINISTIC_MODELS
        }
    )
)

SUMMARY_OUTPUTS = [
    SUMMARY_DIR + f"/{OUTPUT_PREFIX}_per_run_metrics.tsv",
    SUMMARY_DIR + f"/{OUTPUT_PREFIX}_aggregated_metrics.tsv",
    SUMMARY_DIR + f"/{OUTPUT_PREFIX}_publication_summary.tsv",
    SUMMARY_DIR + f"/{OUTPUT_PREFIX}_per_run_metrics.md",
    SUMMARY_DIR + f"/{OUTPUT_PREFIX}_aggregated_metrics.md",
    SUMMARY_DIR + f"/{OUTPUT_PREFIX}_publication_summary.md",
    SUMMARY_DIR + "/per_run_metrics.tsv",
    SUMMARY_DIR + "/aggregated_metrics.tsv",
    SUMMARY_DIR + "/final_summary.tsv",
    SUMMARY_DIR + "/per_run_metrics.md",
    SUMMARY_DIR + "/aggregated_metrics.md",
    SUMMARY_DIR + "/final_summary.md",
]

wildcard_constraints:
    protocol="|".join(PROTOCOLS),
    model="|".join(sorted(set(TRAINABLE_MODELS + DETERMINISTIC_MODELS))),
    feature="|".join(FEATURE_SETTINGS),
    data_key="|".join(ESM2_DATA_KEYS),
    seed=r"\d+"


rule freeze_protocol:
    input:
        LABEL_SOURCE_FILES,
    output:
        LABEL_SUMMARY,
        FUSARIUM_LABEL_SUMMARY,
        SPLIT_SUMMARY,
        FUSARIUM_SPLIT_SUMMARY,
    log:
        SUMMARY_DIR + "/freeze_protocol.log",
    shell:
        r"""
        mkdir -p "{SUMMARY_DIR}" "{FIGURE1_LABELS_DIR}" "{FIGURE1_SPLITS_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.freeze_unified_protocol --config "{config_path}" > "{log}" 2>&1
        cp -f "{LEGACY_LABEL_SUMMARY}" "{LABEL_SUMMARY}"
        cp -f "{LEGACY_FUSARIUM_LABEL_SUMMARY}" "{FUSARIUM_LABEL_SUMMARY}"
        cp -f "{LEGACY_SPLIT_SUMMARY}" "{SPLIT_SUMMARY}"
        cp -f "{LEGACY_FUSARIUM_SPLIT_SUMMARY}" "{FUSARIUM_SPLIT_SUMMARY}"
        """


rule extract_species_esm2_pooled_embeddings:
    input:
        fasta=lambda wc: esm2_fasta_path_for_data_key(wc.data_key),
    output:
        pooled=str(Path(ESM2["cache_root"]) / "{data_key}" / "esm2_pooled.pt"),
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
        str(Path(ESM2["cache_root"]) / "{data_key}" / "extract_esm2_pooled.log"),
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


rule run_trainable_model:
    input:
        LABEL_SUMMARY,
        FUSARIUM_LABEL_SUMMARY,
        SPLIT_SUMMARY,
        FUSARIUM_SPLIT_SUMMARY,
        esm2=lambda wc: required_esm2_cache_for_run(wc.protocol, wc.model),
    output:
        metrics=OUTPUT_ROOT + "/{protocol}/{model}/{feature}/run_{seed}/metrics.tsv",
    log:
        OUTPUT_ROOT + "/{protocol}/{model}/{feature}/run_{seed}/Figure1_materialize_or_run.log",
    params:
        output_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.protocol}/{wc.model}/{wc.feature}/run_{wc.seed}",
    shell:
        r"""
        mkdir -p "{params.output_dir}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}"
        find "{params.output_dir}" -mindepth 1 -delete
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export OMP_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{EPGAT_PYTHON}" -m src.train.run_frozen_protocol_model \
          --config "{config_path}" \
          --protocol "{wildcards.protocol}" \
          --model "{wildcards.model}" \
          --seed {wildcards.seed} \
          --graph-contract "{GRAPH_CONTRACT}" \
          --output-dir "{params.output_dir}" \
          > "{log}" 2>&1
        """


rule run_deterministic_model:
    input:
        LABEL_SUMMARY,
        FUSARIUM_LABEL_SUMMARY,
        SPLIT_SUMMARY,
        FUSARIUM_SPLIT_SUMMARY,
    output:
        metrics=OUTPUT_ROOT + "/{protocol}/{model}/{feature}/deterministic/metrics.tsv",
    log:
        OUTPUT_ROOT + "/{protocol}/{model}/{feature}/deterministic/Figure1_materialize_or_run.log",
    params:
        output_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.protocol}/{wc.model}/{wc.feature}/deterministic",
    shell:
        r"""
        mkdir -p "{params.output_dir}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}"
        find "{params.output_dir}" -mindepth 1 -delete
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export OMP_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{EPGAT_PYTHON}" -m src.train.run_frozen_protocol_model \
          --config "{config_path}" \
          --protocol "{wildcards.protocol}" \
          --model "{wildcards.model}" \
          --graph-contract "{GRAPH_CONTRACT}" \
          --output-dir "{params.output_dir}" \
          > "{log}" 2>&1
        """


rule aggregate_frozen_protocol:
    input:
        trainable=lambda wc: TRAINABLE_OUTPUTS,
        deterministic=lambda wc: DETERMINISTIC_OUTPUTS,
    output:
        SUMMARY_OUTPUTS,
    log:
        SUMMARY_DIR + "/aggregate_frozen_protocol.log",
    shell:
        r"""
        mkdir -p "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.aggregate_frozen_protocol_runs \
          --output-root "{OUTPUT_ROOT}" \
          --summary-dir "{SUMMARY_DIR}" \
          --output-prefix "{OUTPUT_PREFIX}" \
          --include-model-feature-combos "{REQUIRED_SUMMARY_COMBOS}" \
          > "{log}" 2>&1
        """


include: "plots.smk"


rule all:
    default_target: True
    input:
        SUMMARY_OUTPUTS,
        PLOT_OUTPUTS
