configfile: "configs/frozen_protocol.yaml"

from pathlib import Path


PROTOCOL = str(config.get("label_scarcity_protocol", "fgraminearum_newlabel"))
BENCHMARK_FEATURE = str(config.get("label_scarcity_feature_setting", "ORT_EXP_SUB_ESM2")).strip().upper()
GRAPH_MODEL_VARIANT = str(
    config.get(
        "label_scarcity_graph_model",
        "GraphSAGE_ORT_EXP_SUB_ESM2" if BENCHMARK_FEATURE == "ORT_EXP_SUB_ESM2" else "GraphSAGE_ORT_EXP_SUB",
    )
)
TRAIN_FRACTIONS = [float(value) for value in config.get("label_scarcity_train_fractions", [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90])]
FRACTION_TAGS = [f"{int(round(value * 100)):02d}" for value in TRAIN_FRACTIONS]
TRAIN_FRACTIONS_CSV = ",".join(f"{value:.2f}" for value in TRAIN_FRACTIONS)
SEEDS = [int(value) for value in config["runtime"].get("seed_list", [1029, 1030, 1031, 1032, 1033])]
SEEDS_CSV = ",".join(str(value) for value in SEEDS)
OUTPUT_ROOT = str(config.get("label_scarcity_output_root", "outputs/Figure2_label_scarcity"))
RESULTS_ROOT = str(config.get("label_scarcity_results_root", "results/Figure2_label_scarcity"))
SPLIT_ROOT = RESULTS_ROOT + "/splits/" + PROTOCOL
SUMMARY_DIR = RESULTS_ROOT + "/summary"
PLOT_DIR = RESULTS_ROOT + "/plots"
PYTHON_BIN = config["runtime"]["python_bin"]
EPGAT_PYTHON = config["runtime"]["epgat_python_bin"]
MPLCONFIGDIR = config["runtime"].get("mplconfigdir", ".mplconfig")
XDG_CACHE_HOME = config["runtime"].get("xdg_cache_home", ".cache")
BASE_SPLIT_MANIFEST = config["paths"]["splits_dir"] + f"/{PROTOCOL}_split.tsv"
LABEL_SOURCE_FILES = [config["label_sources"][source_key] for source_key in sorted(config["label_sources"])]
LABEL_SUMMARY = config["paths"]["labels_dir"] + "/label_protocol_summary.md"
FUSARIUM_LABEL_SUMMARY = config["paths"]["labels_dir"] + "/fgraminearum_label_protocol_summary.md"
SPLIT_SUMMARY = config["paths"]["splits_dir"] + "/split_protocol_summary.md"
FUSARIUM_SPLIT_SUMMARY = config["paths"]["splits_dir"] + "/fgraminearum_split_protocol_summary.md"

MODEL_SPECS = [
    {"model_variant": GRAPH_MODEL_VARIANT, "feature_setting": BENCHMARK_FEATURE},
    {"model_variant": "MLP", "feature_setting": BENCHMARK_FEATURE},
    {"model_variant": "SVM", "feature_setting": BENCHMARK_FEATURE},
    {"model_variant": "RF", "feature_setting": BENCHMARK_FEATURE},
    {"model_variant": "N2V_MLP", "feature_setting": "N2V"},
    {"model_variant": "DC", "feature_setting": "NETWORK"},
    {"model_variant": "CC", "feature_setting": "NETWORK"},
]
MODEL_VARIANTS = [spec["model_variant"] for spec in MODEL_SPECS]
MODEL_FEATURE_LOOKUP = {spec["model_variant"]: spec["feature_setting"] for spec in MODEL_SPECS}


def split_manifest_path(fraction_tag, seed):
    return f"{SPLIT_ROOT}/train_fraction_{fraction_tag}/split_seed_{seed}.tsv"


def run_metrics_path(fraction_tag, model_variant, seed):
    feature_setting = MODEL_FEATURE_LOOKUP[model_variant]
    return f"{OUTPUT_ROOT}/{PROTOCOL}/train_fraction_{fraction_tag}/{model_variant}/{feature_setting}/run_{seed}/metrics.tsv"


RUN_OUTPUTS = [
    run_metrics_path(fraction_tag, model_variant, seed)
    for fraction_tag in FRACTION_TAGS
    for model_variant in MODEL_VARIANTS
    for seed in SEEDS
]


wildcard_constraints:
    fraction_tag="|".join(FRACTION_TAGS),
    model_variant="|".join(MODEL_VARIANTS),
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
        "{PYTHON_BIN}" -m src.data.freeze_unified_protocol --config "configs/frozen_protocol.yaml" > "{log}" 2>&1
        """


rule build_label_scarcity_splits:
    input:
        LABEL_SUMMARY,
        FUSARIUM_LABEL_SUMMARY,
        SPLIT_SUMMARY,
        FUSARIUM_SPLIT_SUMMARY,
        BASE_SPLIT_MANIFEST,
    output:
        index=SPLIT_ROOT + "/label_scarcity_split_manifest_index.tsv",
        manifests=expand(
            SPLIT_ROOT + "/train_fraction_{fraction_tag}/split_seed_{seed}.tsv",
            fraction_tag=FRACTION_TAGS,
            seed=SEEDS,
        ),
    log:
        SUMMARY_DIR + "/build_label_scarcity_splits.log",
    shell:
        r"""
        mkdir -p "{SPLIT_ROOT}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.build_label_scarcity_split_manifests \
          --base-split-manifest "{BASE_SPLIT_MANIFEST}" \
          --output-dir "{SPLIT_ROOT}" \
          --train-fractions "{TRAIN_FRACTIONS_CSV}" \
          --seeds "{SEEDS_CSV}" \
          > "{log}" 2>&1
        """


rule run_label_scarcity_model:
    input:
        split_index=SPLIT_ROOT + "/label_scarcity_split_manifest_index.tsv",
        split_manifest=lambda wc: split_manifest_path(wc.fraction_tag, wc.seed),
    output:
        metrics=OUTPUT_ROOT + "/" + PROTOCOL + "/train_fraction_{fraction_tag}/{model_variant}/" + "{feature_setting}" + "/run_{seed}/metrics.tsv",
    log:
        OUTPUT_ROOT + "/" + PROTOCOL + "/train_fraction_{fraction_tag}/{model_variant}/" + "{feature_setting}" + "/run_{seed}/label_scarcity_run.log",
    params:
        feature_setting=lambda wc: MODEL_FEATURE_LOOKUP[wc.model_variant],
        output_dir=lambda wc: f"{OUTPUT_ROOT}/{PROTOCOL}/train_fraction_{wc.fraction_tag}/{wc.model_variant}/{MODEL_FEATURE_LOOKUP[wc.model_variant]}/run_{wc.seed}",
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
        "{EPGAT_PYTHON}" -m src.train.run_frozen_protocol_model \
          --config "configs/frozen_protocol.yaml" \
          --protocol "{PROTOCOL}" \
          --model "{wildcards.model_variant}" \
          --feature-setting "{params.feature_setting}" \
          --seed {wildcards.seed} \
          --split-manifest "{input.split_manifest}" \
          --output-dir "{params.output_dir}" \
          > "{log}" 2>&1
        "{PYTHON_BIN}" -m src.eval.repair_frozen_protocol_metrics --run-dir "{params.output_dir}" >> "{log}" 2>&1
        """


rule summarize_label_scarcity:
    input:
        RUN_OUTPUTS,
        split_index=SPLIT_ROOT + "/label_scarcity_split_manifest_index.tsv",
    output:
        per_run=SUMMARY_DIR + "/label_scarcity_per_run_metrics.tsv",
        summary=SUMMARY_DIR + "/label_scarcity_summary.tsv",
        ranking=SUMMARY_DIR + "/label_scarcity_ranking_table.tsv",
        audit=SUMMARY_DIR + "/label_scarcity_coverage_audit.tsv",
        report=SUMMARY_DIR + "/label_scarcity_report.md",
    log:
        SUMMARY_DIR + "/summarize_label_scarcity.log",
    shell:
        r"""
        mkdir -p "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.summarize_label_scarcity \
          --output-root "{OUTPUT_ROOT}" \
          --split-root "{SPLIT_ROOT}" \
          --summary-dir "{SUMMARY_DIR}" \
          --protocol "{PROTOCOL}" \
          --feature-setting "{BENCHMARK_FEATURE}" \
          --seeds "{SEEDS_CSV}" \
          --train-fractions "{TRAIN_FRACTIONS_CSV}" \
          > "{log}" 2>&1
        """


rule plot_label_scarcity:
    input:
        summary=SUMMARY_DIR + "/label_scarcity_summary.tsv",
        ranking=SUMMARY_DIR + "/label_scarcity_ranking_table.tsv",
    output:
        auprc_pdf=PLOT_DIR + "/label_scarcity_auprc.pdf",
        auprc_png=PLOT_DIR + "/label_scarcity_auprc.png",
        auroc_pdf=PLOT_DIR + "/label_scarcity_auroc.pdf",
        auroc_png=PLOT_DIR + "/label_scarcity_auroc.png",
        retention_pdf=PLOT_DIR + "/label_scarcity_retention.pdf",
        retention_png=PLOT_DIR + "/label_scarcity_retention.png",
    log:
        PLOT_DIR + "/plot_label_scarcity.log",
    shell:
        r"""
        mkdir -p "{PLOT_DIR}"
        ~/anaconda3/bin/Rscript src/plot/plot_label_scarcity.R \
          --summary-dir "{SUMMARY_DIR}" \
          --output-dir "{PLOT_DIR}" \
          > "{log}" 2>&1
        """


rule all:
    default_target: True
    input:
        SUMMARY_DIR + "/label_scarcity_per_run_metrics.tsv",
        SUMMARY_DIR + "/label_scarcity_summary.tsv",
        SUMMARY_DIR + "/label_scarcity_ranking_table.tsv",
        SUMMARY_DIR + "/label_scarcity_coverage_audit.tsv",
        SUMMARY_DIR + "/label_scarcity_report.md",
        PLOT_DIR + "/label_scarcity_auprc.pdf",
        PLOT_DIR + "/label_scarcity_auprc.png",
        PLOT_DIR + "/label_scarcity_auroc.pdf",
        PLOT_DIR + "/label_scarcity_auroc.png",
        PLOT_DIR + "/label_scarcity_retention.pdf",
        PLOT_DIR + "/label_scarcity_retention.png"
