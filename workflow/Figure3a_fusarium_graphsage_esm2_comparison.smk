configfile: "configs/Figure3a_fusarium_graphsage_esm2_comparison.yaml"

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
RESULTS_ROOT = config["results_root"]
SUMMARY_DIR = RESULTS_ROOT + "/summary"
DATA_DIR = RESULTS_ROOT + "/data"
PLOTS_DIR = RESULTS_ROOT + "/plots"
RUNTIME_DIR = RESULTS_ROOT + "/runtime"
RUNTIME_CONFIG = RUNTIME_DIR + "/Figure3a_runtime_config.yaml"
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
    return str(feature_setting).strip().upper().endswith("ESM2") or str(feature_setting).strip().upper() == "ESM2"


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
        RESULTS_ROOT + "/freeze_protocol.log",
    shell:
        r"""
        mkdir -p "{RESULTS_ROOT}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.freeze_unified_protocol --config "{BASE_CONFIG_PATH}" > "{log}" 2>&1
        """


rule render_runtime_config:
    output:
        RUNTIME_CONFIG,
    log:
        RESULTS_ROOT + "/render_runtime_config.log",
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.render_figure3_esm2_config \
          --base-config "{BASE_CONFIG_PATH}" \
          --output-config "{output}" \
          > "{log}" 2>&1
        """


rule run_figure3a_graphsage:
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
        OUTPUT_ROOT + "/{protocol}/{model}/{feature_setting}/run_{seed}/Figure3a_run.log",
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


rule generate_figure3a_outputs:
    input:
        RUN_OUTPUTS
    output:
        final_summary=DATA_DIR + "/Figure3a_final_summary.tsv",
        plot_data=DATA_DIR + "/Figure3a_plot_data.tsv",
        panel_a_data=DATA_DIR + "/Figure3a_panelA_fgraminearum_plot_data.tsv",
        panel_b_data=DATA_DIR + "/Figure3a_panelB_scerevisiae_plot_data.tsv",
        main_pdf=PLOTS_DIR + "/Figure3a_main_panels.pdf",
        main_png=PLOTS_DIR + "/Figure3a_main_panels.png",
        panel_a_auprc_pdf=PLOTS_DIR + "/Figure3a_panelA_fgraminearum_auprc.pdf",
        panel_a_auprc_png=PLOTS_DIR + "/Figure3a_panelA_fgraminearum_auprc.png",
        panel_a_mcc_pdf=PLOTS_DIR + "/Figure3a_panelA_fgraminearum_mcc.pdf",
        panel_a_mcc_png=PLOTS_DIR + "/Figure3a_panelA_fgraminearum_mcc.png",
        panel_b_auprc_pdf=PLOTS_DIR + "/Figure3a_panelB_scerevisiae_auprc.pdf",
        panel_b_auprc_png=PLOTS_DIR + "/Figure3a_panelB_scerevisiae_auprc.png",
        panel_b_mcc_pdf=PLOTS_DIR + "/Figure3a_panelB_scerevisiae_mcc.pdf",
        panel_b_mcc_png=PLOTS_DIR + "/Figure3a_panelB_scerevisiae_mcc.png",
        run_manifest=SUMMARY_DIR + "/Figure3a_run_manifest.tsv",
        coverage_audit=SUMMARY_DIR + "/Figure3a_feature_coverage_audit.tsv",
        generation_report=SUMMARY_DIR + "/Figure3a_generation_report.md",
    log:
        RESULTS_ROOT + "/Figure3a_generate_summary.log",
    shell:
        r"""
        mkdir -p "{RESULTS_ROOT}" "{SUMMARY_DIR}" "{DATA_DIR}" "{PLOTS_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.generate_figure3a_outputs \
          --output-root "{OUTPUT_ROOT}" \
          --results-root "{RESULTS_ROOT}" \
          --prefix "{SUMMARY_PREFIX}" \
          > "{log}" 2>&1
        """


rule all:
    default_target: True
    input:
        DATA_DIR + "/Figure3a_final_summary.tsv",
        DATA_DIR + "/Figure3a_plot_data.tsv",
        DATA_DIR + "/Figure3a_panelA_fgraminearum_plot_data.tsv",
        DATA_DIR + "/Figure3a_panelB_scerevisiae_plot_data.tsv",
        PLOTS_DIR + "/Figure3a_main_panels.pdf",
        PLOTS_DIR + "/Figure3a_main_panels.png",
        PLOTS_DIR + "/Figure3a_panelA_fgraminearum_auprc.pdf",
        PLOTS_DIR + "/Figure3a_panelA_fgraminearum_auprc.png",
        PLOTS_DIR + "/Figure3a_panelA_fgraminearum_mcc.pdf",
        PLOTS_DIR + "/Figure3a_panelA_fgraminearum_mcc.png",
        PLOTS_DIR + "/Figure3a_panelB_scerevisiae_auprc.pdf",
        PLOTS_DIR + "/Figure3a_panelB_scerevisiae_auprc.png",
        PLOTS_DIR + "/Figure3a_panelB_scerevisiae_mcc.pdf",
        PLOTS_DIR + "/Figure3a_panelB_scerevisiae_mcc.png",
        SUMMARY_DIR + "/Figure3a_run_manifest.tsv",
        SUMMARY_DIR + "/Figure3a_feature_coverage_audit.tsv",
        SUMMARY_DIR + "/Figure3a_generation_report.md"
