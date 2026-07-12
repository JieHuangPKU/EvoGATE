# Phase 2A Fixed300 Multirun

## Scope

- fixed config reproduction only
- 4 species x 4 models x `EXP_SUB_ORT` x `thr_300` x 3 runs
- no threshold sweep
- no alternate feature combinations
- no archived replay benchmark outputs in the formal summary

## Environment

- training runtime: `conda activate EPGAT`
- workflow launcher: `scripts/run_phase2a_fixed300_multirun.sh`
- Snakemake file: `workflow/phase2a_fixed300_multirun.smk`

## Output Root

- `outputs/epgat_phase2a_fixed300_multirun`

Per-run outputs:

- `outputs/epgat_phase2a_fixed300_multirun/{species}/{model}/EXP_SUB_ORT/thr_300/run_{run_id}/metrics.tsv`
- `outputs/epgat_phase2a_fixed300_multirun/{species}/{model}/EXP_SUB_ORT/thr_300/run_{run_id}/resolved_config.yaml`
- `outputs/epgat_phase2a_fixed300_multirun/{species}/{model}/EXP_SUB_ORT/thr_300/run_{run_id}/feature_schema.tsv`
- `outputs/epgat_phase2a_fixed300_multirun/{species}/{model}/EXP_SUB_ORT/thr_300/run_{run_id}/feature_summary.tsv`

Per-species x model aggregates:

- `outputs/epgat_phase2a_fixed300_multirun/{species}/{model}/EXP_SUB_ORT/thr_300/aggregated_metrics.tsv`

Final summary:

- `outputs/epgat_phase2a_fixed300_multirun/final_summary.tsv`

## How To Run

Shell wrapper:

```bash
bash scripts/run_phase2a_fixed300_multirun.sh 4
```

Manual Snakemake:

```bash
conda activate EPGAT
export PYTHONPATH="${PYTHONPATH:-.}:."
/home/jiehuang/anaconda3/bin/snakemake \
  --snakefile workflow/phase2a_fixed300_multirun.smk \
  --cores 4 \
  --rerun-incomplete \
  --keep-going
```

## How To Read Results

- Inspect one run under its `run_{id}` directory for raw metrics and resolved config.
- Inspect `aggregated_metrics.tsv` for each species x model mean/std across 3 runs.
- Inspect `final_summary.tsv` for the cross-species final table.
