configfile: "configs/Figure3b_fusarium_graphsage_esm2_dim_ablation.yaml"

from pathlib import Path

import yaml

BASE_CONFIG_PATH = config["base_config"]
BASE_CONFIG = yaml.safe_load(open(BASE_CONFIG_PATH, "r", encoding="utf-8"))
SUMMARY_PREFIX = config["summary_prefix"]
TARGET_NAME = config["target_name"]
PROTOCOLS = config["protocols"]
MODELS = config["models"]
FEATURE_SETTING = config["feature_setting"]
ESM2_DIMS = [int(value) for value in config["esm2_dimensions"]]
SEEDS = config["seed_list"]
OUTPUT_ROOT = config["output_root"]
SUMMARY_DIR = config["results_root"]
PYTHON_BIN = BASE_CONFIG["runtime"]["python_bin"]
EPGAT_PYTHON = BASE_CONFIG["runtime"]["epgat_python_bin"]
MPLCONFIGDIR = BASE_CONFIG["runtime"].get("mplconfigdir", ".mplconfig")
XDG_CACHE_HOME = BASE_CONFIG["runtime"].get("xdg_cache_home", ".cache")
LABEL_SUMMARY = BASE_CONFIG["paths"]["labels_dir"] + "/label_protocol_summary.md"
FUSARIUM_LABEL_SUMMARY = BASE_CONFIG["paths"]["labels_dir"] + "/fgraminearum_label_protocol_summary.md"
SPLIT_SUMMARY = BASE_CONFIG["paths"]["splits_dir"] + "/split_protocol_summary.md"
FUSARIUM_SPLIT_SUMMARY = BASE_CONFIG["paths"]["splits_dir"] + "/fgraminearum_split_protocol_summary.md"
LABEL_SOURCE_FILES = [
    BASE_CONFIG["label_sources"][source_key]
    for source_key in sorted(BASE_CONFIG["label_sources"])
]


def protocol_data_key(protocol_name):
    return str(BASE_CONFIG["protocols"][protocol_name]["data_key"])


def esm2_cache_path(protocol_name):
    data_key = protocol_data_key(protocol_name)
    return str(Path(BASE_CONFIG["esm2"]["cache_root"]) / data_key / "esm2_pooled.pt")


def dim_label(dim_value):
    return f"ESM2-{int(dim_value)}"


def runtime_config_path(dim_value):
    return f"{SUMMARY_DIR}/runtime/{dim_label(dim_value)}.yaml"


def run_metrics_path(protocol, model, dim_value, seed):
    return f"{OUTPUT_ROOT}/{dim_label(dim_value)}/{protocol}/{model}/{FEATURE_SETTING}/run_{seed}/metrics.tsv"


RUNTIME_CONFIGS = [runtime_config_path(dim_value) for dim_value in ESM2_DIMS]
RUN_OUTPUTS = [
    run_metrics_path(protocol, model, dim_value, seed)
    for protocol in PROTOCOLS
    for model in MODELS
    for dim_value in ESM2_DIMS
    for seed in SEEDS
]


wildcard_constraints:
    protocol="|".join(PROTOCOLS),
    model="|".join(MODELS),
    dim_label="|".join(dim_label(dim_value) for dim_value in ESM2_DIMS),
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
        mkdir -p "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.freeze_unified_protocol --config "{BASE_CONFIG_PATH}" > "{log}" 2>&1
        """


rule render_dim_runtime_config:
    output:
        config=SUMMARY_DIR + "/runtime/{dim_label}.yaml",
    log:
        SUMMARY_DIR + "/runtime/{dim_label}.log",
    params:
        dim_value=lambda wc: int(str(wc.dim_label).split("-")[-1]),
    shell:
        r"""
        mkdir -p "$(dirname "{output.config}")"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.render_figure3_esm2_config \
          --base-config "{BASE_CONFIG_PATH}" \
          --output-config "{output.config}" \
          --esm2-projection-dim {params.dim_value} \
          > "{log}" 2>&1
        """


rule run_figure3b_graphsage_dim:
    input:
        LABEL_SUMMARY,
        FUSARIUM_LABEL_SUMMARY,
        SPLIT_SUMMARY,
        FUSARIUM_SPLIT_SUMMARY,
        config=lambda wc: runtime_config_path(int(str(wc.dim_label).split("-")[-1])),
        esm2=lambda wc: [esm2_cache_path(wc.protocol)],
    output:
        metrics=OUTPUT_ROOT + "/{dim_label}/{protocol}/{model}/" + FEATURE_SETTING + "/run_{seed}/metrics.tsv",
    log:
        OUTPUT_ROOT + "/{dim_label}/{protocol}/{model}/" + FEATURE_SETTING + "/run_{seed}/Figure3b_run.log",
    params:
        output_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.dim_label}/{wc.protocol}/{wc.model}/{FEATURE_SETTING}/run_{wc.seed}",
    shell:
        r"""
        mkdir -p "{params.output_dir}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export OMP_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{EPGAT_PYTHON}" -m src.train.run_frozen_protocol_feature_combo_model \
          --config "{input.config}" \
          --protocol "{wildcards.protocol}" \
          --model "{wildcards.model}" \
          --feature-setting "{FEATURE_SETTING}" \
          --seed {wildcards.seed} \
          --output-dir "{params.output_dir}" \
          > "{log}" 2>&1
        """


rule summarize_figure3b:
    input:
        RUN_OUTPUTS
    output:
        per_run=SUMMARY_DIR + "/Figure3b_per_run_metrics.tsv",
        aggregated=SUMMARY_DIR + "/Figure3b_aggregated_metrics.tsv",
        final=SUMMARY_DIR + "/Figure3b_final_summary.tsv",
        per_run_md=SUMMARY_DIR + "/Figure3b_per_run_metrics.md",
        aggregated_md=SUMMARY_DIR + "/Figure3b_aggregated_metrics.md",
        final_md=SUMMARY_DIR + "/Figure3b_final_summary.md",
    log:
        SUMMARY_DIR + "/Figure3b_generate_summary.log",
    shell:
        r"""
        mkdir -p "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.generate_figure_benchmark_summary \
          --output-root "{OUTPUT_ROOT}" \
          --summary-dir "{SUMMARY_DIR}" \
          --prefix "{SUMMARY_PREFIX}" \
          --target-name "{TARGET_NAME}" \
          > "{log}" 2>&1
        """


rule all:
    default_target: True
    input:
        SUMMARY_DIR + "/Figure3b_final_summary.tsv",
        SUMMARY_DIR + "/Figure3b_final_summary.md",
        SUMMARY_DIR + "/Figure3b_aggregated_metrics.tsv",
        SUMMARY_DIR + "/Figure3b_aggregated_metrics.md",
        SUMMARY_DIR + "/Figure3b_per_run_metrics.tsv",
        SUMMARY_DIR + "/Figure3b_per_run_metrics.md"
