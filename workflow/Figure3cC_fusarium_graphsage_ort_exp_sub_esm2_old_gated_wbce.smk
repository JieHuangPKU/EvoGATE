configfile: "configs/Figure3cC_fusarium_graphsage_ort_exp_sub_esm2_old_gated_wbce.yaml"

from pathlib import Path

import yaml

BASE_CONFIG_PATH = config["base_config"]
BASE_CONFIG = yaml.safe_load(open(BASE_CONFIG_PATH, "r", encoding="utf-8"))
SUMMARY_PREFIX = config["summary_prefix"]
TARGET_NAME = config["target_name"]
PROTOCOLS = config["protocols"]
MODELS = config["models"]
PROTOCOL_FEATURE_SETTINGS = config["protocol_feature_settings"]
FEATURE_SETTINGS = sorted({feature for features in PROTOCOL_FEATURE_SETTINGS.values() for feature in features})
SEEDS = config["seed_list"]
OUTPUT_ROOT = config["output_root"]
SUMMARY_DIR = config["results_root"]
COMPARISON_OUTPUT_ROOT = config["comparison_output_root"]
COMPARISON_RESULTS_ROOT = config["comparison_results_root"]
RUNTIME_CONFIG = SUMMARY_DIR + "/runtime/Figure3cC_runtime_config.yaml"
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


def uses_esm2(feature_setting):
    normalized = str(feature_setting).strip().upper()
    return "ESM2" in normalized


def protocol_data_key(protocol_name):
    return str(BASE_CONFIG["protocols"][protocol_name]["data_key"])


def esm2_cache_path(protocol_name):
    data_key = protocol_data_key(protocol_name)
    return str(Path(BASE_CONFIG["esm2"]["cache_root"]) / data_key / "esm2_pooled.pt")


def run_metrics_path(protocol, model, feature_setting, seed):
    return f"{OUTPUT_ROOT}/{protocol}/{model}/{feature_setting}/run_{seed}/metrics.tsv"


RUN_OUTPUTS = [
    run_metrics_path(protocol, model, feature_setting, seed)
    for protocol in PROTOCOLS
    for model in MODELS
    for feature_setting in PROTOCOL_FEATURE_SETTINGS[protocol]
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


rule render_runtime_config:
    output:
        RUNTIME_CONFIG,
    log:
        SUMMARY_DIR + "/render_runtime_config.log",
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.render_figure3_esm2_config \
          --base-config "{BASE_CONFIG_PATH}" \
          --output-config "{output}" \
          --fusion-mode "{config[fusion_mode]}" \
          --fusion-hidden-dim {config[fusion_hidden_dim]} \
          --fusion-dropout {config[fusion_dropout]} \
          --loss-type "{config[loss_type]}" \
          --pos-weight-mode "{config[pos_weight_mode]}" \
          --pos-weight-scale {config[pos_weight_scale]} \
          > "{log}" 2>&1
        """


rule run_figure3cc_old_gated_wbce:
    input:
        LABEL_SUMMARY,
        FUSARIUM_LABEL_SUMMARY,
        SPLIT_SUMMARY,
        FUSARIUM_SPLIT_SUMMARY,
        config=RUNTIME_CONFIG,
        esm2=lambda wc: [esm2_cache_path(wc.protocol)] if uses_esm2(wc.feature_setting) else [],
    output:
        metrics=OUTPUT_ROOT + "/{protocol}/{model}/{feature_setting}/run_{seed}/metrics.tsv",
    log:
        OUTPUT_ROOT + "/{protocol}/{model}/{feature_setting}/run_{seed}/run_frozen_protocol_model.log",
    params:
        output_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.protocol}/{wc.model}/{wc.feature_setting}/run_{wc.seed}",
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
          --feature-setting "{wildcards.feature_setting}" \
          --seed {wildcards.seed} \
          --output-dir "{params.output_dir}" \
          > "{log}" 2>&1
        """


rule summarize_figure3cc:
    input:
        RUN_OUTPUTS
    output:
        per_run=SUMMARY_DIR + "/Figure3cC_per_run_metrics.tsv",
        aggregated=SUMMARY_DIR + "/Figure3cC_aggregated_metrics.tsv",
        final=SUMMARY_DIR + "/Figure3cC_final_summary.tsv",
        per_run_md=SUMMARY_DIR + "/Figure3cC_per_run_metrics.md",
        aggregated_md=SUMMARY_DIR + "/Figure3cC_aggregated_metrics.md",
        final_md=SUMMARY_DIR + "/Figure3cC_final_summary.md",
    log:
        SUMMARY_DIR + "/Figure3cC_generate_summary.log",
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


rule summarize_figure3cc_gate_statistics:
    input:
        RUN_OUTPUTS
    output:
        SUMMARY_DIR + "/Figure3cC_gate_statistics.tsv",
        SUMMARY_DIR + "/Figure3cC_gate_statistics.md",
    log:
        SUMMARY_DIR + "/Figure3cC_gate_statistics.log",
    shell:
        r"""
        mkdir -p "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.summarize_gate_statistics \
          --output-root "{OUTPUT_ROOT}" \
          --summary-dir "{SUMMARY_DIR}" \
          --prefix "{SUMMARY_PREFIX}" \
          > "{log}" 2>&1
        """


rule summarize_figure3cc_diagnostics:
    input:
        RUN_OUTPUTS,
        SUMMARY_DIR + "/Figure3cC_final_summary.tsv",
        COMPARISON_RESULTS_ROOT + "/Figure3cB_final_summary.tsv",
    output:
        SUMMARY_DIR + "/Figure3cC_comparison_summary.tsv",
        SUMMARY_DIR + "/Figure3cC_comparison_summary.md",
        SUMMARY_DIR + "/Figure3cC_threshold_diagnostics.tsv",
        SUMMARY_DIR + "/Figure3cC_threshold_diagnostics.md",
        SUMMARY_DIR + "/Figure3cC_probability_summary.tsv",
        SUMMARY_DIR + "/Figure3cC_probability_summary.md",
        SUMMARY_DIR + "/Figure3cC_calibration_summary.tsv",
        SUMMARY_DIR + "/Figure3cC_calibration_summary.md",
        SUMMARY_DIR + "/Figure3cC_pr_curve_data.tsv",
        SUMMARY_DIR + "/Figure3cC_pr_curve_data.md",
    log:
        SUMMARY_DIR + "/Figure3cC_diagnostics.log",
    shell:
        r"""
        mkdir -p "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.generate_figure3cc_diagnostics \
          --comparison-root "{COMPARISON_OUTPUT_ROOT}" \
          --new-root "{OUTPUT_ROOT}" \
          --summary-dir "{SUMMARY_DIR}" \
          > "{log}" 2>&1
        """


rule all:
    default_target: True
    input:
        SUMMARY_DIR + "/Figure3cC_final_summary.tsv",
        SUMMARY_DIR + "/Figure3cC_final_summary.md",
        SUMMARY_DIR + "/Figure3cC_aggregated_metrics.tsv",
        SUMMARY_DIR + "/Figure3cC_aggregated_metrics.md",
        SUMMARY_DIR + "/Figure3cC_per_run_metrics.tsv",
        SUMMARY_DIR + "/Figure3cC_per_run_metrics.md",
        SUMMARY_DIR + "/Figure3cC_gate_statistics.tsv",
        SUMMARY_DIR + "/Figure3cC_gate_statistics.md",
        SUMMARY_DIR + "/Figure3cC_comparison_summary.tsv",
        SUMMARY_DIR + "/Figure3cC_comparison_summary.md",
        SUMMARY_DIR + "/Figure3cC_threshold_diagnostics.tsv",
        SUMMARY_DIR + "/Figure3cC_threshold_diagnostics.md",
        SUMMARY_DIR + "/Figure3cC_probability_summary.tsv",
        SUMMARY_DIR + "/Figure3cC_probability_summary.md",
        SUMMARY_DIR + "/Figure3cC_calibration_summary.tsv",
        SUMMARY_DIR + "/Figure3cC_calibration_summary.md",
        SUMMARY_DIR + "/Figure3cC_pr_curve_data.tsv",
        SUMMARY_DIR + "/Figure3cC_pr_curve_data.md"
