# Legacy Consolidation

## Goal

Keep historical protocol code available without allowing it to contaminate the new frozen benchmark path.

## Mainline Boundary

The new mainline boundary is:

- `configs/frozen_protocol.yaml`
- `workflow/frozen_protocol_benchmark.smk`
- `scripts/run_frozen_protocol_benchmark.sh`
- `src/data/freeze_unified_protocol.py`
- `src/data/frozen_protocol_loader.py`
- `src/train/run_frozen_protocol_model.py`
- `src/eval/aggregate_frozen_protocol_runs.py`

Everything outside this boundary is non-authoritative unless explicitly documented otherwise.

## Legacy Replay Retained

Retained for replay / comparison:

- `workflow/classical_baseline_benchmark.smk`
- `workflow/phase2a_fixed300_multirun.smk`
- `src/train/run_epgat_phase2a_single.py`
- `src/eval/aggregate_phase2a_multirun.py`
- `src/train/train_epgat_legacy.py`
- `src/train/train_epgat_graph_models.py`

## Deprecated Fusarium-Specific Comparison Assets

Retained only for archival comparison:

- `workflow/fgraminearum_feature_combo_benchmark.smk`
- `workflow/fgraminearum_feature_combo_newlabel_benchmark.smk`
- `scripts/run_fgraminearum_feature_combo_benchmark.sh`
- `scripts/run_fgraminearum_feature_combo_newlabel_benchmark.sh`
- `src/train/run_old440_diagnostic.py`
- `src/train/run_label_rebuild_compare.py`
- `src/eval/compare_old440_vs_newlabel.py`
- `configs/label_rebuild_compare.yaml`

## Reference Updates Made

The new mainline no longer depends on:

- `data_registry`
- ranking-only subset files
- implicit old440 defaulting
- ambiguous Fusarium filenames

The new mainline now depends on:

- repo-local processed tables
- explicit frozen manifests
- explicit protocol names

## Cleanup Action in This Refactor

The ambiguous residual split file:

- `results/frozen_protocol/splits/fgraminearum_split.tsv`

is treated as a deprecated residual from the earlier frozen-protocol trial and should not remain in the mainline output surface.

## Future Optional Move

If desired after full validation, the deprecated Fusarium comparison assets can be moved under a dedicated subtree such as:

- `legacy/`
- `archive/`

This refactor does not mass-move or mass-delete them yet.
