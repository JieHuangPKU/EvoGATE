configfile: "configs/frozen_protocol.yaml"

BASE_CONFIG = config
PYTHON_BIN = BASE_CONFIG["runtime"]["python_bin"]
MPLCONFIGDIR = BASE_CONFIG["runtime"].get("mplconfigdir", ".mplconfig")
XDG_CACHE_HOME = "/tmp/codex_cache"
NUMBA_CACHE_DIR = "/tmp/numba_cache"
TMPDIR = "/tmp"

RESULTS_ROOT = str(BASE_CONFIG.get("figure5_representation", {}).get("results_root", "results/Figure5"))
PLOTS_DIR = RESULTS_ROOT + "/plots"
DATA_DIR = RESULTS_ROOT + "/data"
TABLES_DIR = RESULTS_ROOT + "/tables"
SUMMARY_DIR = RESULTS_ROOT + "/summary"
RUNTIME_CONFIG = "results/Figure3a/runtime/Figure3a_runtime_config.yaml"
UPSTREAM_ROOT = "outputs/Figure3a"
FIGURE5_PROTOCOLS = ["fgraminearum_newlabel", "scerevisiae"]
FIGURE5_SEEDS = [1029, 1030, 1031, 1032, 1033]
FIGURE5_PROTOCOL_ARGS = " ".join(FIGURE5_PROTOCOLS)

REPRESENTATIVE_TABLE = TABLES_DIR + "/Figure5_representative_seed_selection.tsv"
REPRESENTATIVE_SUMMARY = SUMMARY_DIR + "/Figure5_representative_seed_selection.md"

FIG5D_OUTPUTS = [
    DATA_DIR + "/Figure5d_feature_group_manifest.tsv",
    DATA_DIR + "/Figure5d_group_ablation_per_gene.tsv",
    DATA_DIR + "/Figure5d_feature_group_global_importance_plot_data.tsv",
    DATA_DIR + "/Figure5d_feature_group_rescued_vs_other_plot_data.tsv",
    DATA_DIR + "/Figure5d_feature_group_gene_set_heatmap_plot_data.tsv",
    TABLES_DIR + "/Figure5d_group_ablation_global_summary.tsv",
    TABLES_DIR + "/Figure5d_group_ablation_global_metrics.tsv",
    TABLES_DIR + "/Figure5d_group_ablation_by_gene_set.tsv",
    TABLES_DIR + "/Figure5d_group_ablation_stats.tsv",
    TABLES_DIR + "/Figure5d_group_ablation_gene_set_heatmap_values.tsv",
    TABLES_DIR + "/Figure5d_inference_reproducibility_check.tsv",
    TABLES_DIR + "/Figure5d_output_manifest.tsv",
    PLOTS_DIR + "/Figure5d_feature_group_global_importance.pdf",
    PLOTS_DIR + "/Figure5d_feature_group_global_importance.png",
    PLOTS_DIR + "/Figure5d_feature_group_rescued_vs_other.pdf",
    PLOTS_DIR + "/Figure5d_feature_group_rescued_vs_other.png",
    PLOTS_DIR + "/Figure5d_feature_group_gene_set_heatmap.pdf",
    PLOTS_DIR + "/Figure5d_feature_group_gene_set_heatmap.png",
    SUMMARY_DIR + "/Figure5d_feature_group_attribution.md",
]


rule select_figure5_representative_seed:
    input:
        RUNTIME_CONFIG,
        expand("outputs/Figure3a/{protocol}/GraphSAGE/ORT_EXP_SUB/run_{seed}/best_model.pt", protocol=FIGURE5_PROTOCOLS, seed=FIGURE5_SEEDS),
        expand("outputs/Figure3a/{protocol}/GraphSAGE/ORT_EXP_SUB_ESM2/run_{seed}/best_model.pt", protocol=FIGURE5_PROTOCOLS, seed=FIGURE5_SEEDS),
        expand("outputs/Figure3a/{protocol}/GraphSAGE/ORT_EXP_SUB/run_{seed}/predictions.tsv", protocol=FIGURE5_PROTOCOLS, seed=FIGURE5_SEEDS),
        expand("outputs/Figure3a/{protocol}/GraphSAGE/ORT_EXP_SUB_ESM2/run_{seed}/predictions.tsv", protocol=FIGURE5_PROTOCOLS, seed=FIGURE5_SEEDS),
        expand("outputs/Figure3a/{protocol}/GraphSAGE/ORT_EXP_SUB/run_{seed}/metrics.tsv", protocol=FIGURE5_PROTOCOLS, seed=FIGURE5_SEEDS),
        expand("outputs/Figure3a/{protocol}/GraphSAGE/ORT_EXP_SUB_ESM2/run_{seed}/metrics.tsv", protocol=FIGURE5_PROTOCOLS, seed=FIGURE5_SEEDS),
    output:
        REPRESENTATIVE_TABLE,
        REPRESENTATIVE_SUMMARY,
    log:
        SUMMARY_DIR + "/Figure5_representative_seed_selection.log",
    shell:
        r"""
        mkdir -p "{PLOTS_DIR}" "{DATA_DIR}" "{TABLES_DIR}" "{SUMMARY_DIR}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}" "{NUMBA_CACHE_DIR}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export NUMBA_CACHE_DIR="{NUMBA_CACHE_DIR}"
        export TMPDIR="{TMPDIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.analysis.select_figure5_representative_seed \
          --runtime-config "{RUNTIME_CONFIG}" \
          --upstream-root "{UPSTREAM_ROOT}" \
          --protocols {FIGURE5_PROTOCOL_ARGS} \
          --output-table "{output[0]}" \
          --output-summary "{output[1]}" \
          > "{log}" 2>&1
        """


rule build_figure5d:
    input:
        REPRESENTATIVE_TABLE,
        REPRESENTATIVE_SUMMARY,
    output:
        FIG5D_OUTPUTS,
    log:
        SUMMARY_DIR + "/Figure5d_build.log",
    shell:
        r"""
        mkdir -p "{PLOTS_DIR}" "{DATA_DIR}" "{TABLES_DIR}" "{SUMMARY_DIR}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}" "{NUMBA_CACHE_DIR}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export NUMBA_CACHE_DIR="{NUMBA_CACHE_DIR}"
        export TMPDIR="{TMPDIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.analysis.build_figure5d_feature_group_attribution \
          --runtime-config "{RUNTIME_CONFIG}" \
          --upstream-root "{UPSTREAM_ROOT}" \
          --protocols {FIGURE5_PROTOCOL_ARGS} \
          --data-dir "{DATA_DIR}" \
          --table-dir "{TABLES_DIR}" \
          --plot-dir "{PLOTS_DIR}" \
          --summary-dir "{SUMMARY_DIR}" \
          > "{log}" 2>&1
        """


rule figure5d:
    input:
        FIG5D_OUTPUTS


rule all:
    default_target: True
    input:
        REPRESENTATIVE_TABLE,
        REPRESENTATIVE_SUMMARY,
        FIG5D_OUTPUTS
