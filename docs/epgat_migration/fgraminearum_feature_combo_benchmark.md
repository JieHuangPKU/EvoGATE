# F. graminearum Feature Combo Benchmark

## Scope

- species: `fgraminearum`
- graph models: `GAT`, `GCN`, `GIN`, `GraphSAGE`
- feature combos: `ORT`, `EXP`, `SUB`, `ORT_EXP`, `ORT_SUB`, `EXP_SUB`, `ORT_EXP_SUB`
- fixed config only: `string threshold = 300`, `include_degree = false`
- execution grid: `7 combos x 4 models x 3 runs`
- analysis includes summary tables, AUROC/AUPRC/MCC heatmaps, and combination-level feature contribution analysis

## Environment

- conda environment: `EPGAT`
- Snakemake file: `workflow/fgraminearum_feature_combo_benchmark.smk`
- shell launcher: `scripts/run_fgraminearum_feature_combo_benchmark.sh`

## How To Run

Shell wrapper:

```bash
bash scripts/run_fgraminearum_feature_combo_benchmark.sh 4
```

Manual Snakemake:

```bash
cd /home/jiehuang/software/fungi/ProGATE_v2
conda activate EPGAT
export PYTHONPATH="${PYTHONPATH:-.}:."
/home/jiehuang/anaconda3/bin/snakemake \
  --snakefile workflow/fgraminearum_feature_combo_benchmark.smk \
  --cores 4 \
  --rerun-incomplete \
  --keep-going
```

## Output Locations

Benchmark outputs:

- `outputs/fgraminearum_feature_combo_benchmark/{model}/{feature_combo}/thr_300/run_{run_id}/metrics.tsv`
- `outputs/fgraminearum_feature_combo_benchmark/{model}/{feature_combo}/thr_300/run_{run_id}/resolved_config.yaml`
- `outputs/fgraminearum_feature_combo_benchmark/{model}/{feature_combo}/thr_300/run_{run_id}/feature_schema.tsv`
- `outputs/fgraminearum_feature_combo_benchmark/{model}/{feature_combo}/thr_300/run_{run_id}/feature_summary.tsv`
- `outputs/fgraminearum_feature_combo_benchmark/{model}/{feature_combo}/thr_300/aggregated_metrics.tsv`
- `outputs/fgraminearum_feature_combo_benchmark/final_summary.tsv`

Analysis outputs:

- `results/fgraminearum_feature_combo_analysis/model_feature_combo_summary.tsv`
- `results/fgraminearum_feature_combo_analysis/fgraminearum_auroc_heatmap.png`
- `results/fgraminearum_feature_combo_analysis/fgraminearum_auprc_heatmap.png`
- `results/fgraminearum_feature_combo_analysis/fgraminearum_mcc_heatmap.png`
- `results/fgraminearum_feature_combo_analysis/feature_main_effects.tsv`
- `results/fgraminearum_feature_combo_analysis/fgraminearum_feature_combo_summary.md`
