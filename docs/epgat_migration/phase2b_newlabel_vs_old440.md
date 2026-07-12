# Phase 2B New-Label vs old440

Update: this document describes a historical Phase 2B comparison path. The current mainline
Fusarium label source is the materialized processed label directory under
`data/processed/essential_gene/fgraminearum/`, not `results/label_rebuild_experiments/`.

## Scope

- old benchmark source: `workflow/fgraminearum_feature_combo_benchmark.smk`
- new-label benchmark source: `workflow/fgraminearum_feature_combo_newlabel_benchmark.smk`
- both workflows run the same Fusarium benchmark grid:
  - 4 graph models: `GAT`, `GCN`, `GIN`, `GraphSAGE`
  - 7 feature combinations: `ORT`, `EXP`, `SUB`, `ORT_EXP`, `ORT_SUB`, `EXP_SUB`, `ORT_EXP_SUB`
  - fixed `string threshold = 300`
  - fixed `include_degree = false`
  - 3 runs per configuration
- the only intentional change in the historical new workflow was the label input:
  - old440 keeps the completed archived benchmark outputs
  - new-label originally used archived `positive_set_P1.tsv` and `negative_set.tsv` from the rebuilt label workflow; current mainline uses materialized processed labels instead

## Why The Two Benchmarks Are Comparable

- feature construction is unchanged
- graph model choices are unchanged
- run count and fixed threshold are unchanged
- `run_epgat_phase2a_single.py` is reused for both
- the companion workflow swaps only the positive/negative label tables, so the label effect is isolated from the feature/model design

## Label Inputs

- old440 labels are the existing archived benchmark result set under `outputs/fgraminearum_feature_combo_benchmark/`
- historical new-label positive snapshot was copied from `results/label_rebuild_experiments/labels/positive_set_P1.tsv`
- historical new-label negative snapshot was copied from `results/label_rebuild_experiments/labels/negative_set.tsv`
- current materialized mainline label assets live under `data/processed/essential_gene/fgraminearum/newlabel/`
- prepared benchmark-ready historical comparison files are emitted to `results/phase2b_new_label/labels/`

## How To Run The New Benchmark

Shell wrapper:

```bash
cd /home/jiehuang/software/fungi/ProGATE_v2
bash scripts/run_fgraminearum_feature_combo_newlabel_benchmark.sh 4
```

Manual Snakemake:

```bash
cd /home/jiehuang/software/fungi/ProGATE_v2
conda activate EPGAT
export PYTHONPATH="${PYTHONPATH:-.}:."
export MPLBACKEND=Agg
export MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig
export XDG_CACHE_HOME=/home/jiehuang/software/fungi/ProGATE_v2/.cache
/home/jiehuang/anaconda3/bin/snakemake \
  --snakefile workflow/fgraminearum_feature_combo_newlabel_benchmark.smk \
  --cores 4 \
  --rerun-incomplete \
  --keep-going
```

## Outputs

New-label benchmark runs:

- `outputs/fgraminearum_feature_combo_newlabel_benchmark/{model}/{feature_combo}/thr_300/run_{run_id}/`
- `outputs/fgraminearum_feature_combo_newlabel_benchmark/{model}/{feature_combo}/thr_300/aggregated_metrics.tsv`

Comparison outputs:

- `results/phase2b_new_label/final_summary.tsv`
- `results/phase2b_new_label/old440_vs_newlabel_comparison.tsv`
- `results/phase2b_new_label/phase2b_summary.md`
- `results/phase2b_new_label/figures/heatmap_old440_auprc.pdf`
- `results/phase2b_new_label/figures/heatmap_newlabel_auprc.pdf`
- `results/phase2b_new_label/figures/delta_heatmap_auprc.pdf`
- `results/phase2b_new_label/figures/delta_heatmap_mcc.pdf`
- `results/phase2b_new_label/figures/best_feature_old_vs_new_barplot.pdf`

## How To Read The Comparison

- use `final_summary.tsv` for the new-label aggregated benchmark table
- use `old440_vs_newlabel_comparison.tsv` for per-model x per-feature-combo metric deltas
- use `phase2b_summary.md` for the main biological/model-selection interpretation
- use the figure PDFs for publication-style old-vs-new visual comparison
