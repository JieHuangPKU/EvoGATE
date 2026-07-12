configfile: "configs/frozen_protocol.yaml"

FIG4_CFG = config.get("figure4_graph_robustness", {})
PROTOCOL = str(FIG4_CFG.get("protocol", "fgraminearum_newlabel"))
MODEL = str(FIG4_CFG.get("model", "GraphSAGE"))
MODEL_VARIANT = MODEL
FEATURE_SETTING = str(FIG4_CFG.get("feature_setting", "ORT_EXP_SUB_ESM2"))
PYTHON_BIN = config["runtime"]["python_bin"]
EPGAT_PYTHON = config["runtime"]["epgat_python_bin"]
MPLCONFIGDIR = config["runtime"].get("mplconfigdir", ".mplconfig")
XDG_CACHE_HOME = config["runtime"].get("xdg_cache_home", ".cache")
OUTPUT_ROOT = str(FIG4_CFG.get("output_root", "outputs/Figure4"))
RESULTS_ROOT = str(FIG4_CFG.get("results_root", "results/Figure4"))
SUMMARY_DIR = RESULTS_ROOT + "/summary"
DATA_DIR = RESULTS_ROOT + "/data"
SUPP_DIR = RESULTS_ROOT + "/supplementary"
GRAPH_DIR = OUTPUT_ROOT + "/source_graphs"
RUN_OUTPUT_ROOT = OUTPUT_ROOT + "/source_comparison"
THRESHOLD_OUTPUT_ROOT = OUTPUT_ROOT + "/threshold_sweep"
STRING_GRAPH = "data/processed/PPI/fgraminearum/string.csv"
RAW_EFG_GRAPH = "data/processed/PPI/fgraminearum/eFG_ppis.txt"
BRIDGE_PATH = "data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv"
SPLIT_MANIFEST = config["paths"]["splits_dir"] + "/fgraminearum_newlabel_split.tsv"
LABEL_MANIFEST = config["paths"]["labels_dir"] + "/fgraminearum_newlabel.tsv"
SEEDS = [int(seed) for seed in config["runtime"]["seed_list"]]
MAIN_THRESHOLDS = [int(value) for value in FIG4_CFG.get("thresholds_main", [100, 200, 300, 400, 500, 600, 700, 800, 900])]
QUICK_THRESHOLDS = [int(value) for value in FIG4_CFG.get("thresholds_quick", [300, 700])]
THRESHOLD_MODE = str(config.get("figure4_threshold_mode", "main")).lower().strip()
THRESHOLDS = QUICK_THRESHOLDS if THRESHOLD_MODE == "quick" else MAIN_THRESHOLDS
THRESHOLDS_CSV = ",".join(str(value) for value in THRESHOLDS)
SEEDS_CSV = ",".join(str(value) for value in SEEDS)

EFG_GRAPH_PATHS = {
    "eFG": GRAPH_DIR + "/fgraminearum_eFG_canonical.tsv",
    "eFG_HIGH": GRAPH_DIR + "/fgraminearum_eFG_HIGH_canonical.tsv",
    "eFG_HIGH_MEDIUM": GRAPH_DIR + "/fgraminearum_eFG_HIGH_MEDIUM_canonical.tsv",
    "eFG_ALL": GRAPH_DIR + "/fgraminearum_eFG_ALL_canonical.tsv",
}
EFG_CONFIDENCE_FILTERS = {
    "eFG": str(FIG4_CFG.get("efg_main_confidence_filter", "ALL")),
    "eFG_HIGH": "HIGH",
    "eFG_HIGH_MEDIUM": "HIGH_MEDIUM",
    "eFG_ALL": "ALL",
}
EFG_ADAPTER_SUMMARIES = {
    "eFG": SUMMARY_DIR + "/Figure4_eFG_adapter_summary.tsv",
    "eFG_HIGH": SUMMARY_DIR + "/Figure4_eFG_HIGH_adapter_summary.tsv",
    "eFG_HIGH_MEDIUM": SUMMARY_DIR + "/Figure4_eFG_HIGH_MEDIUM_adapter_summary.tsv",
    "eFG_ALL": SUMMARY_DIR + "/Figure4_eFG_ALL_adapter_summary.tsv",
}
EFG_ADAPTER_MAPPINGS = {
    "eFG": SUMMARY_DIR + "/Figure4_eFG_adapter_mapping.tsv",
    "eFG_HIGH": SUMMARY_DIR + "/Figure4_eFG_HIGH_adapter_mapping.tsv",
    "eFG_HIGH_MEDIUM": SUMMARY_DIR + "/Figure4_eFG_HIGH_MEDIUM_adapter_mapping.tsv",
    "eFG_ALL": SUMMARY_DIR + "/Figure4_eFG_ALL_adapter_mapping.tsv",
}
MAIN_EFG_GRAPH = EFG_GRAPH_PATHS["eFG"]
MAIN_EFG_ADAPTER_SUMMARY = EFG_ADAPTER_SUMMARIES["eFG"]
MAIN_EFG_OUTPUT_ROOT = RUN_OUTPUT_ROOT + "/eFG"
SUPPLEMENTARY_SOURCE_ORDER = ["STRING_300", "STRING_700", "eFG_HIGH", "eFG_HIGH_MEDIUM", "eFG_ALL"]
SUPPLEMENTARY_EFG_VARIANTS = ["eFG_HIGH", "eFG_HIGH_MEDIUM", "eFG_ALL"]
SUPPLEMENTARY_OUTPUT_ROOTS_CSV = ",".join(
    [
        f"STRING_300={THRESHOLD_OUTPUT_ROOT}/string_thr_300",
        f"STRING_700={THRESHOLD_OUTPUT_ROOT}/string_thr_700",
        *[f"{variant}={RUN_OUTPUT_ROOT}/{variant}" for variant in SUPPLEMENTARY_EFG_VARIANTS],
    ]
)
SUPPLEMENTARY_EFG_GRAPHS_CSV = ",".join([f"{variant}={EFG_GRAPH_PATHS[variant]}" for variant in SUPPLEMENTARY_EFG_VARIANTS])
SUPPLEMENTARY_EFG_ADAPTER_SUMMARIES_CSV = ",".join([f"{variant}={EFG_ADAPTER_SUMMARIES[variant]}" for variant in SUPPLEMENTARY_EFG_VARIANTS])


def threshold_metrics_path(threshold, seed):
    return f"{THRESHOLD_OUTPUT_ROOT}/string_thr_{threshold}/{PROTOCOL}/{MODEL}/{FEATURE_SETTING}/run_{seed}/metrics.tsv"


def threshold_context_path(threshold, seed):
    return f"{THRESHOLD_OUTPUT_ROOT}/string_thr_{threshold}/{PROTOCOL}/{MODEL}/{FEATURE_SETTING}/run_{seed}/figure4_run_context.tsv"


def source_metrics_path(graph_variant, seed):
    return f"{RUN_OUTPUT_ROOT}/{graph_variant}/{PROTOCOL}/{MODEL}/{FEATURE_SETTING}/run_{seed}/metrics.tsv"


def source_context_path(graph_variant, seed):
    return f"{RUN_OUTPUT_ROOT}/{graph_variant}/{PROTOCOL}/{MODEL}/{FEATURE_SETTING}/run_{seed}/figure4_run_context.tsv"


wildcard_constraints:
    threshold=r"\d+",
    seed=r"\d+",
    graph_variant=r"eFG|eFG_HIGH|eFG_HIGH_MEDIUM|eFG_ALL"


rule freeze_protocol:
    input:
        "configs/frozen_protocol.yaml",
    output:
        LABEL_MANIFEST,
        SPLIT_MANIFEST,
    log:
        RESULTS_ROOT + "/freeze_protocol.log",
    shell:
        r"""
        mkdir -p "{RESULTS_ROOT}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.freeze_unified_protocol --config "configs/frozen_protocol.yaml" > "{log}" 2>&1
        """


rule prepare_efg_graph:
    input:
        SPLIT_MANIFEST,
        BRIDGE_PATH,
        RAW_EFG_GRAPH,
    output:
        graph=GRAPH_DIR + "/fgraminearum_{graph_variant}_canonical.tsv",
        summary=SUMMARY_DIR + "/Figure4_{graph_variant}_adapter_summary.tsv",
        mapping=SUMMARY_DIR + "/Figure4_{graph_variant}_adapter_mapping.tsv",
    log:
        RESULTS_ROOT + "/prepare_{graph_variant}_graph.log",
    params:
        confidence_filter=lambda wc: EFG_CONFIDENCE_FILTERS[wc.graph_variant],
    shell:
        r"""
        mkdir -p "{GRAPH_DIR}" "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.prepare_figure4_efg_graph \
          --raw-efg "{RAW_EFG_GRAPH}" \
          --split-manifest "{SPLIT_MANIFEST}" \
          --bridge "{BRIDGE_PATH}" \
          --confidence-filter "{params.confidence_filter}" \
          --output-graph "{output.graph}" \
          --output-summary "{output.summary}" \
          --output-mapping "{output.mapping}" \
          > "{log}" 2>&1
        """


rule run_threshold_graphsage:
    input:
        LABEL_MANIFEST,
        SPLIT_MANIFEST,
    output:
        metrics=THRESHOLD_OUTPUT_ROOT + "/string_thr_{threshold}/" + PROTOCOL + "/" + MODEL + "/" + FEATURE_SETTING + "/run_{seed}/metrics.tsv",
        context=THRESHOLD_OUTPUT_ROOT + "/string_thr_{threshold}/" + PROTOCOL + "/" + MODEL + "/" + FEATURE_SETTING + "/run_{seed}/figure4_run_context.tsv",
    log:
        THRESHOLD_OUTPUT_ROOT + "/string_thr_{threshold}/" + PROTOCOL + "/" + MODEL + "/" + FEATURE_SETTING + "/run_{seed}/Figure4_run.log",
    params:
        output_dir=lambda wc: f"{THRESHOLD_OUTPUT_ROOT}/string_thr_{wc.threshold}/{PROTOCOL}/{MODEL}/{FEATURE_SETTING}/run_{wc.seed}",
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
          --config "configs/frozen_protocol.yaml" \
          --protocol "{PROTOCOL}" \
          --model "{MODEL_VARIANT}" \
          --feature-setting "{FEATURE_SETTING}" \
          --seed {wildcards.seed} \
          --graph-source "{STRING_GRAPH}" \
          --graph-source-name "STRING" \
          --string-threshold {wildcards.threshold} \
          --output-dir "{params.output_dir}" \
          > "{log}" 2>&1
        "{PYTHON_BIN}" -m src.eval.repair_frozen_protocol_metrics --run-dir "{params.output_dir}" >> "{log}" 2>&1
        printf 'workflow_target\tprotocol\tmodel\tfeature_setting\tgraph_condition\tgraph_source_name\tgraph_threshold\tseed\nFigure4_graph_robustness\t%s\t%s\t%s\tSTRING_%s\tSTRING\t%s\t%s\n' \
          "{PROTOCOL}" "{MODEL}" "{FEATURE_SETTING}" "{wildcards.threshold}" "{wildcards.threshold}" "{wildcards.seed}" \
          > "{output.context}"
        """


rule run_efg_graphsage:
    input:
        LABEL_MANIFEST,
        SPLIT_MANIFEST,
        graph=lambda wc: EFG_GRAPH_PATHS[wc.graph_variant],
    output:
        metrics=RUN_OUTPUT_ROOT + "/{graph_variant}/" + PROTOCOL + "/" + MODEL + "/" + FEATURE_SETTING + "/run_{seed}/metrics.tsv",
        context=RUN_OUTPUT_ROOT + "/{graph_variant}/" + PROTOCOL + "/" + MODEL + "/" + FEATURE_SETTING + "/run_{seed}/figure4_run_context.tsv",
    log:
        RUN_OUTPUT_ROOT + "/{graph_variant}/" + PROTOCOL + "/" + MODEL + "/" + FEATURE_SETTING + "/run_{seed}/Figure4_run.log",
    params:
        output_dir=lambda wc: f"{RUN_OUTPUT_ROOT}/{wc.graph_variant}/{PROTOCOL}/{MODEL}/{FEATURE_SETTING}/run_{wc.seed}",
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
          --config "configs/frozen_protocol.yaml" \
          --protocol "{PROTOCOL}" \
          --model "{MODEL_VARIANT}" \
          --feature-setting "{FEATURE_SETTING}" \
          --seed {wildcards.seed} \
          --graph-source "{input.graph}" \
          --graph-source-name "{wildcards.graph_variant}" \
          --output-dir "{params.output_dir}" \
          > "{log}" 2>&1
        "{PYTHON_BIN}" -m src.eval.repair_frozen_protocol_metrics --run-dir "{params.output_dir}" >> "{log}" 2>&1
        printf 'workflow_target\tprotocol\tmodel\tfeature_setting\tgraph_condition\tgraph_source_name\tgraph_threshold\tseed\nFigure4_graph_robustness\t%s\t%s\t%s\t%s\t%s\t\t%s\n' \
          "{PROTOCOL}" "{MODEL}" "{FEATURE_SETTING}" "{wildcards.graph_variant}" "{wildcards.graph_variant}" "{wildcards.seed}" \
          > "{output.context}"
        """


rule summarize_figure4:
    input:
        threshold_metrics=[threshold_metrics_path(threshold, seed) for threshold in THRESHOLDS for seed in SEEDS],
        threshold_context=[threshold_context_path(threshold, seed) for threshold in THRESHOLDS for seed in SEEDS],
        main_efg_metrics=[source_metrics_path("eFG", seed) for seed in SEEDS],
        main_efg_context=[source_context_path("eFG", seed) for seed in SEEDS],
        supp_efg_metrics=[source_metrics_path(graph_variant, seed) for graph_variant in SUPPLEMENTARY_EFG_VARIANTS for seed in SEEDS],
        supp_efg_context=[source_context_path(graph_variant, seed) for graph_variant in SUPPLEMENTARY_EFG_VARIANTS for seed in SEEDS],
        split_manifest=SPLIT_MANIFEST,
        main_efg_graph=EFG_GRAPH_PATHS["eFG"],
        efg_adapter_summary=EFG_ADAPTER_SUMMARIES["eFG"],
        supp_efg_graphs=[EFG_GRAPH_PATHS[variant] for variant in SUPPLEMENTARY_EFG_VARIANTS],
        supp_efg_adapter_summaries=[EFG_ADAPTER_SUMMARIES[variant] for variant in SUPPLEMENTARY_EFG_VARIANTS],
    output:
        per_run=SUMMARY_DIR + "/Figure4_threshold_per_run_metrics.tsv",
        aggregated=SUMMARY_DIR + "/Figure4_threshold_aggregated_metrics.tsv",
        density=SUMMARY_DIR + "/Figure4_network_density_summary.tsv",
        source_per_run=SUMMARY_DIR + "/Figure4_source_comparison_per_run_metrics.tsv",
        source=SUMMARY_DIR + "/Figure4_source_comparison_metrics.tsv",
        overlap=SUMMARY_DIR + "/Figure4_edge_overlap_summary.tsv",
        supp_per_run=SUMMARY_DIR + "/Figure4_supplementary_source_comparison_per_run_metrics.tsv",
        supp_source=SUMMARY_DIR + "/Figure4_supplementary_source_comparison_metrics.tsv",
    log:
        SUMMARY_DIR + "/Figure4_summary.log",
    shell:
        r"""
        mkdir -p "{SUMMARY_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.generate_figure4_graph_robustness_summary \
          --protocol "{PROTOCOL}" \
          --model "{MODEL}" \
          --feature-setting "{FEATURE_SETTING}" \
          --thresholds "{THRESHOLDS_CSV}" \
          --seeds "{SEEDS_CSV}" \
          --threshold-output-root "{THRESHOLD_OUTPUT_ROOT}" \
          --main-efg-output-root "{MAIN_EFG_OUTPUT_ROOT}" \
          --supplementary-output-roots "{SUPPLEMENTARY_OUTPUT_ROOTS_CSV}" \
          --string-graph "{STRING_GRAPH}" \
          --main-efg-graph "{MAIN_EFG_GRAPH}" \
          --supplementary-efg-graphs "{SUPPLEMENTARY_EFG_GRAPHS_CSV}" \
          --split-manifest "{SPLIT_MANIFEST}" \
          --efg-adapter-summary "{MAIN_EFG_ADAPTER_SUMMARY}" \
          --supplementary-efg-adapter-summaries "{SUPPLEMENTARY_EFG_ADAPTER_SUMMARIES_CSV}" \
          --summary-dir "{SUMMARY_DIR}" \
          > "{log}" 2>&1
        """


rule plot_figure4:
    input:
        SUMMARY_DIR + "/Figure4_threshold_per_run_metrics.tsv",
        SUMMARY_DIR + "/Figure4_threshold_aggregated_metrics.tsv",
        SUMMARY_DIR + "/Figure4_network_density_summary.tsv",
        SUMMARY_DIR + "/Figure4_source_comparison_per_run_metrics.tsv",
        SUMMARY_DIR + "/Figure4_source_comparison_metrics.tsv",
        SUMMARY_DIR + "/Figure4_edge_overlap_summary.tsv",
        SUMMARY_DIR + "/Figure4_supplementary_source_comparison_per_run_metrics.tsv",
        SUMMARY_DIR + "/Figure4_supplementary_source_comparison_metrics.tsv",
    output:
        RESULTS_ROOT + "/Figure4A_threshold_performance_line.pdf",
        RESULTS_ROOT + "/Figure4A_threshold_performance_line.png",
        RESULTS_ROOT + "/Figure4B_network_density_plot.pdf",
        RESULTS_ROOT + "/Figure4B_network_density_plot.png",
        RESULTS_ROOT + "/Figure4C_source_comparison_barplot.pdf",
        RESULTS_ROOT + "/Figure4C_source_comparison_barplot.png",
        RESULTS_ROOT + "/Figure4D_edge_overlap_plot.pdf",
        RESULTS_ROOT + "/Figure4D_edge_overlap_plot.png",
        DATA_DIR + "/Figure4A_threshold_performance_line.tsv",
        DATA_DIR + "/Figure4A_threshold_performance_line_summary.tsv",
        DATA_DIR + "/Figure4B_network_density_plot.tsv",
        DATA_DIR + "/Figure4C_source_comparison_barplot.tsv",
        DATA_DIR + "/Figure4C_source_comparison_barplot_summary.tsv",
        DATA_DIR + "/Figure4D_edge_overlap_plot.tsv",
        SUPP_DIR + "/Figure4S1_eFG_confidence_source_comparison_barplot.pdf",
        SUPP_DIR + "/Figure4S1_eFG_confidence_source_comparison_barplot.png",
        SUPP_DIR + "/Figure4S1_eFG_confidence_source_comparison_barplot.tsv",
        SUPP_DIR + "/Figure4S1_eFG_confidence_source_comparison_barplot_summary.tsv",
    log:
        RESULTS_ROOT + "/Figure4_plot.log",
    shell:
        r"""
        mkdir -p "{RESULTS_ROOT}" "{DATA_DIR}" "{SUPP_DIR}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.plot.plot_figure4_graph_robustness \
          --summary-dir "{SUMMARY_DIR}" \
          --output-dir "{RESULTS_ROOT}" \
          > "{log}" 2>&1
        """


rule all:
    default_target: True
    input:
        SUMMARY_DIR + "/Figure4_threshold_per_run_metrics.tsv",
        SUMMARY_DIR + "/Figure4_threshold_aggregated_metrics.tsv",
        SUMMARY_DIR + "/Figure4_source_comparison_metrics.tsv",
        SUMMARY_DIR + "/Figure4_network_density_summary.tsv",
        SUMMARY_DIR + "/Figure4_edge_overlap_summary.tsv",
        SUMMARY_DIR + "/Figure4_eFG_adapter_summary.tsv",
        SUMMARY_DIR + "/Figure4_supplementary_source_comparison_metrics.tsv",
        RESULTS_ROOT + "/Figure4A_threshold_performance_line.pdf",
        RESULTS_ROOT + "/Figure4B_network_density_plot.pdf",
        RESULTS_ROOT + "/Figure4C_source_comparison_barplot.pdf",
        RESULTS_ROOT + "/Figure4D_edge_overlap_plot.pdf",
        DATA_DIR + "/Figure4A_threshold_performance_line.tsv",
        DATA_DIR + "/Figure4A_threshold_performance_line_summary.tsv",
        DATA_DIR + "/Figure4B_network_density_plot.tsv",
        DATA_DIR + "/Figure4C_source_comparison_barplot.tsv",
        DATA_DIR + "/Figure4C_source_comparison_barplot_summary.tsv",
        DATA_DIR + "/Figure4D_edge_overlap_plot.tsv",
        SUPP_DIR + "/Figure4S1_eFG_confidence_source_comparison_barplot.pdf",
        SUPP_DIR + "/Figure4S1_eFG_confidence_source_comparison_barplot.tsv",
        SUPP_DIR + "/Figure4S1_eFG_confidence_source_comparison_barplot_summary.tsv"
