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

FIG5A_OUTPUTS = [
    DATA_DIR + "/Figure5a_rescue_consensus_per_seed.tsv",
    DATA_DIR + "/Figure5a_rescue_consensus_gene_level.tsv",
    DATA_DIR + "/Figure5a_rescue_frequency_distribution_data.tsv",
    TABLES_DIR + "/Figure5a_rescue_consensus_summary.tsv",
    TABLES_DIR + "/Figure5a_rescue_confidence_tiers.tsv",
    TABLES_DIR + "/Figure5a_rescue_transition_stability_distribution.tsv",
    TABLES_DIR + "/Figure5a_rescue_frequency_distribution.tsv",
    TABLES_DIR + "/Figure5a_transition_counts_all_species.tsv",
    TABLES_DIR + "/Figure5a_output_manifest.tsv",
    PLOTS_DIR + "/Figure5a_rescue_consensus_barplot.pdf",
    PLOTS_DIR + "/Figure5a_rescue_consensus_barplot.png",
    PLOTS_DIR + "/Figure5a_rescue_frequency_distribution.pdf",
    PLOTS_DIR + "/Figure5a_rescue_frequency_distribution.png",
    SUMMARY_DIR + "/Figure5a_representative_visualization_summary.md",
    SUMMARY_DIR + "/Figure5a_rescue_consensus_summary.md",
    SUMMARY_DIR + "/Figure5a_rescue_frequency_distribution.md",
]

FIG5B_OUTPUTS = [
    DATA_DIR + "/Figure5b_hidden_quant_summary_raw.tsv",
    TABLES_DIR + "/Figure5b_hidden_quant_summary_per_seed.tsv",
    TABLES_DIR + "/Figure5b_hidden_quant_summary_stats.tsv",
    TABLES_DIR + "/Figure5b_hidden_quant_summary_all_species_stats.tsv",
    TABLES_DIR + "/Figure5b_hidden_quant_summary_scaling_metadata.tsv",
    TABLES_DIR + "/Figure5b_output_manifest.tsv",
    PLOTS_DIR + "/Figure5b_hidden_quant_summary.pdf",
    PLOTS_DIR + "/Figure5b_hidden_quant_summary.png",
    SUMMARY_DIR + "/Figure5b_hidden_quant_summary.md",
]

FIG5C_OUTPUTS = [
    DATA_DIR + "/Figure5c_input_vs_hidden_compare_coords.tsv",
    DATA_DIR + "/Figure5c_local_neighborhood_analysis.tsv",
    TABLES_DIR + "/Figure5c_input_vs_hidden_compare_metrics.tsv",
    TABLES_DIR + "/Figure5c_input_vs_hidden_compare_all_species_metrics.tsv",
    TABLES_DIR + "/Figure5c_local_neighborhood_summary.tsv",
    TABLES_DIR + "/Figure5c_output_manifest.tsv",
    PLOTS_DIR + "/Figure5c_local_neighborhood_summary.pdf",
    PLOTS_DIR + "/Figure5c_local_neighborhood_summary.png",
    SUMMARY_DIR + "/Figure5c_input_vs_hidden_compare.md",
]

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

FIG5SUMMARY_OUTPUTS = [
    TABLES_DIR + "/Figure5_contract_check.tsv",
    SUMMARY_DIR + "/Figure5_build_summary.md",
    SUMMARY_DIR + "/Figure5_rerun_notes.md",
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


rule build_figure5a:
    input:
        REPRESENTATIVE_TABLE,
    output:
        FIG5A_OUTPUTS,
    log:
        SUMMARY_DIR + "/Figure5a_build.log",
    shell:
        r"""
        mkdir -p "{PLOTS_DIR}" "{DATA_DIR}" "{TABLES_DIR}" "{SUMMARY_DIR}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}" "{NUMBA_CACHE_DIR}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export NUMBA_CACHE_DIR="{NUMBA_CACHE_DIR}"
        export TMPDIR="{TMPDIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.analysis.build_figure5a_final \
          --runtime-config "{RUNTIME_CONFIG}" \
          --upstream-root "{UPSTREAM_ROOT}" \
          --protocols {FIGURE5_PROTOCOL_ARGS} \
          --selection-table "{REPRESENTATIVE_TABLE}" \
          --data-dir "{DATA_DIR}" \
          --table-dir "{TABLES_DIR}" \
          --plot-dir "{PLOTS_DIR}" \
          --summary-dir "{SUMMARY_DIR}" \
          > "{log}" 2>&1
        """


rule build_figure5b:
    input:
        REPRESENTATIVE_TABLE,
    output:
        FIG5B_OUTPUTS,
    log:
        SUMMARY_DIR + "/Figure5b_build.log",
    shell:
        r"""
        mkdir -p "{PLOTS_DIR}" "{DATA_DIR}" "{TABLES_DIR}" "{SUMMARY_DIR}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}" "{NUMBA_CACHE_DIR}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export NUMBA_CACHE_DIR="{NUMBA_CACHE_DIR}"
        export TMPDIR="{TMPDIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.analysis.build_figure5b_final \
          --runtime-config "{RUNTIME_CONFIG}" \
          --upstream-root "{UPSTREAM_ROOT}" \
          --protocols {FIGURE5_PROTOCOL_ARGS} \
          --selection-table "{REPRESENTATIVE_TABLE}" \
          --data-dir "{DATA_DIR}" \
          --table-dir "{TABLES_DIR}" \
          --plot-dir "{PLOTS_DIR}" \
          --summary-dir "{SUMMARY_DIR}" \
          > "{log}" 2>&1
        """


rule build_figure5c:
    input:
        REPRESENTATIVE_TABLE,
        DATA_DIR + "/Figure5a_rescue_consensus_gene_level.tsv",
    output:
        FIG5C_OUTPUTS,
    log:
        SUMMARY_DIR + "/Figure5c_build.log",
    shell:
        r"""
        mkdir -p "{PLOTS_DIR}" "{DATA_DIR}" "{TABLES_DIR}" "{SUMMARY_DIR}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}" "{NUMBA_CACHE_DIR}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export NUMBA_CACHE_DIR="{NUMBA_CACHE_DIR}"
        export TMPDIR="{TMPDIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.analysis.build_figure5c_final \
          --runtime-config "{RUNTIME_CONFIG}" \
          --upstream-root "{UPSTREAM_ROOT}" \
          --protocols {FIGURE5_PROTOCOL_ARGS} \
          --selection-table "{REPRESENTATIVE_TABLE}" \
          --consensus-gene-table "{DATA_DIR}/Figure5a_rescue_consensus_gene_level.tsv" \
          --data-dir "{DATA_DIR}" \
          --table-dir "{TABLES_DIR}" \
          --plot-dir "{PLOTS_DIR}" \
          --summary-dir "{SUMMARY_DIR}" \
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


rule build_figure5_summary:
    input:
        REPRESENTATIVE_TABLE,
        REPRESENTATIVE_SUMMARY,
        FIG5A_OUTPUTS,
        FIG5B_OUTPUTS,
        FIG5C_OUTPUTS,
        FIG5D_OUTPUTS,
    output:
        FIG5SUMMARY_OUTPUTS,
    log:
        SUMMARY_DIR + "/Figure5_build_summary.log",
    shell:
        r"""
        mkdir -p "{PLOTS_DIR}" "{DATA_DIR}" "{TABLES_DIR}" "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.analysis.write_figure5_build_summary \
          --results-root "{RESULTS_ROOT}" \
          --output "{SUMMARY_DIR}/Figure5_build_summary.md" \
          > "{log}" 2>&1
        """


rule figure5a:
    input:
        FIG5A_OUTPUTS


rule figure5b:
    input:
        FIG5B_OUTPUTS


rule figure5c:
    input:
        FIG5C_OUTPUTS


rule figure5d:
    input:
        FIG5D_OUTPUTS


rule all:
    default_target: True
    input:
        REPRESENTATIVE_TABLE,
        REPRESENTATIVE_SUMMARY,
        FIG5A_OUTPUTS,
        FIG5B_OUTPUTS,
        FIG5C_OUTPUTS,
        FIG5D_OUTPUTS,
        FIG5SUMMARY_OUTPUTS
