configfile: "configs/Figure2b_fusarium_gnn_feature_ablation.yaml"

import yaml

BASE_CONFIG_PATH = config["base_config"]
BASE_CONFIG = yaml.safe_load(open(BASE_CONFIG_PATH, "r", encoding="utf-8"))
SUMMARY_PREFIX = config["summary_prefix"]
TARGET_NAME = config["target_name"]
PROTOCOLS = config["protocols"]
MODELS = config["models"]
FEATURE_SETTINGS = config["feature_settings"]
SEEDS = config["seed_list"]
OUTPUT_ROOT = config["output_root"]
SUMMARY_DIR = config["results_root"]
REUSE_SOURCE_ROOT = config["reuse_source_root"]
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


def run_metrics_path(protocol, model, feature_setting, seed):
    return f"{OUTPUT_ROOT}/{protocol}/{model}/{feature_setting}/run_{seed}/metrics.tsv"


RUN_OUTPUTS = [
    run_metrics_path(protocol, model, feature_setting, seed)
    for protocol in PROTOCOLS
    for model in MODELS
    for feature_setting in FEATURE_SETTINGS
    for seed in SEEDS
]


wildcard_constraints:
    protocol="|".join(PROTOCOLS),
    model="|".join(MODELS),
    feature_setting="|".join(FEATURE_SETTINGS),
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


rule materialize_or_run_figure2b:
    input:
        LABEL_SUMMARY,
        FUSARIUM_LABEL_SUMMARY,
        SPLIT_SUMMARY,
        FUSARIUM_SPLIT_SUMMARY,
    output:
        metrics=OUTPUT_ROOT + "/{protocol}/{model}/{feature_setting}/run_{seed}/metrics.tsv",
    log:
        OUTPUT_ROOT + "/{protocol}/{model}/{feature_setting}/run_{seed}/Figure2b_materialize_or_run.log",
    params:
        output_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.protocol}/{wc.model}/{wc.feature_setting}/run_{wc.seed}",
        source_dir=lambda wc: f"{REUSE_SOURCE_ROOT}/{wc.protocol}/{wc.model}/{wc.feature_setting}/run_{wc.seed}",
    shell:
        r"""
        mkdir -p "{params.output_dir}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}"
        {{
          echo "summary_prefix={SUMMARY_PREFIX}"
          echo "source_dir={params.source_dir}"
          echo "output_dir={params.output_dir}"
          if [ -f "{params.source_dir}/metrics.tsv" ]; then
            echo "mode=reuse_existing_run_dir"
            cp -a "{params.source_dir}/." "{params.output_dir}/"
          else
            echo "mode=train_missing_run_dir"
            export MPLBACKEND=Agg
            export MPLCONFIGDIR="{MPLCONFIGDIR}"
            export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
            export OMP_NUM_THREADS=1
            export MKL_NUM_THREADS=1
            export OPENBLAS_NUM_THREADS=1
            export NUMEXPR_NUM_THREADS=1
            export PYTHONPATH="${{PYTHONPATH:-.}}:."
            "{EPGAT_PYTHON}" -m src.train.run_frozen_protocol_feature_combo_model \
              --config "{BASE_CONFIG_PATH}" \
              --protocol "{wildcards.protocol}" \
              --model "{wildcards.model}" \
              --feature-setting "{wildcards.feature_setting}" \
              --seed {wildcards.seed} \
              --output-dir "{params.output_dir}"
          fi
          export PYTHONPATH="${{PYTHONPATH:-.}}:."
          "{PYTHON_BIN}" -m src.eval.repair_frozen_protocol_metrics --run-dir "{params.output_dir}"
        }} > "{log}" 2>&1
        """


rule summarize_figure2b:
    input:
        RUN_OUTPUTS
    output:
        per_run=SUMMARY_DIR + "/Figure2b_per_run_metrics.tsv",
        aggregated=SUMMARY_DIR + "/Figure2b_aggregated_metrics.tsv",
        final=SUMMARY_DIR + "/Figure2b_final_summary.tsv",
        per_run_md=SUMMARY_DIR + "/Figure2b_per_run_metrics.md",
        aggregated_md=SUMMARY_DIR + "/Figure2b_aggregated_metrics.md",
        final_md=SUMMARY_DIR + "/Figure2b_final_summary.md",
    log:
        SUMMARY_DIR + "/Figure2b_generate_summary.log",
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
        SUMMARY_DIR + "/Figure2b_final_summary.tsv",
        SUMMARY_DIR + "/Figure2b_final_summary.md",
        SUMMARY_DIR + "/Figure2b_aggregated_metrics.tsv",
        SUMMARY_DIR + "/Figure2b_aggregated_metrics.md",
        SUMMARY_DIR + "/Figure2b_per_run_metrics.tsv",
        SUMMARY_DIR + "/Figure2b_per_run_metrics.md"
