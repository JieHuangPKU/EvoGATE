import json

with open(GROUP_CONFIG_PATH, "r", encoding="utf-8") as handle:
    FIGURE1_GROUP_CONFIG = json.load(handle)

PLOT_OUTPUTS = []
for group in FIGURE1_GROUP_CONFIG["groups"]:
    group_id = group["id"]
    base_name = f"{OUTPUT_PREFIX}_{group_id}"
    PLOT_OUTPUTS.extend(
        [
            f"{PLOTS_DIR}/{base_name}_radar_2x3.pdf",
            f"{PLOTS_DIR}/{base_name}_radar_2x3.png",
            f"{PLOTS_DIR}/{base_name}_barplot_2x3.pdf",
            f"{PLOTS_DIR}/{base_name}_barplot_2x3.png",
            f"{PLOTS_DIR}/{base_name}_plot_summary.md",
        ]
    )


rule plot_figure1_benchmark_summary:
    input:
        summary=SUMMARY_DIR + f"/{OUTPUT_PREFIX}_publication_summary.tsv",
        group_config=GROUP_CONFIG_PATH,
    output:
        PLOT_OUTPUTS,
    log:
        PLOTS_DIR + "/plot_Figure1.log",
    shell:
        r"""
        mkdir -p "{PLOTS_DIR}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="{MPLCONFIGDIR}"
        export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.eval.plot_classical_baseline_summary \
          --summary-tsv "{input.summary}" \
          --group-config "{input.group_config}" \
          --output-prefix "{OUTPUT_PREFIX}" \
          --output-dir "{PLOTS_DIR}" \
          > "{log}" 2>&1
        """
