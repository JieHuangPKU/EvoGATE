# Label Rebuild Archive Status

## Current Status

`results/label_rebuild_experiments/` is now treated as a historical archive / provenance directory.

It preserves:

- the old label rebuild comparison outputs
- the archived `positive_set_P1.tsv` / `negative_set.tsv` snapshots
- the old440 replay diagnostics

It is no longer the mainline source of Fusarium label assets.

## Current Mainline

The current mainline Fusarium label pipeline is the label materialization workflow:

- `scripts/run_fgraminearum_label_materialization.sh`
- `workflow/fgraminearum_label_materialization.smk`

That workflow materializes the current label assets under:

- `data/processed/essential_gene/fgraminearum/oldlabel/`
- `data/processed/essential_gene/fgraminearum/newlabel/`
- `data/processed/essential_gene/fgraminearum/bridge/`

Those processed assets are then consumed by:

- `configs/frozen_protocol.yaml`
- `src/data/freeze_unified_protocol.py`
- `workflow/frozen_protocol_benchmark.smk`
- `workflow/Figure2a_fusarium_label_compare_graphsage.smk`
- `workflow/Figure4_graph_robustness.smk`

## Remaining Legacy References

The following legacy scripts still reference `results/label_rebuild_experiments/` directly:

- `src/train/run_label_rebuild_compare.py`
- `src/train/run_old440_diagnostic.py`

These scripts are retained for historical experiment replay and provenance only.

## Why It Is Not Deleted Yet

The archive is being kept in place because:

- it still documents the historical label rebuild outputs used in earlier comparison work
- legacy diagnostic scripts still read and write this directory
- several refactor / migration notes still refer to these archived outputs for provenance

Current policy: keep the directory as archive, do not use it as a new workflow input, and do not delete it until the remaining legacy references are retired.
