# Protocol State Audit

This audit records the repository state after the dependency-aware frozen-protocol refactor and before any large-scale legacy deletion.

## Current Mainline Files

The authoritative new benchmark path is now:

- `configs/frozen_protocol.yaml`
- `workflow/frozen_protocol_benchmark.smk`
- `scripts/run_frozen_protocol_benchmark.sh`
- `src/data/freeze_unified_protocol.py`
- `src/data/frozen_protocol_loader.py`
- `src/train/run_frozen_protocol_model.py`
- `src/eval/aggregate_frozen_protocol_runs.py`

These files define the two-layer architecture now used by mainline:

- Layer 1: frozen protocol inputs under `results/frozen_protocol/labels/` and `results/frozen_protocol/splits/`
- Layer 2: model families consuming the frozen contract through `src/data/frozen_protocol_loader.py`

Mainline benchmark protocols are now:

- `human`
- `celegans`
- `scerevisiae`
- `dmelanogaster`
- `fgraminearum_newlabel`

`fgraminearum_oldlabel` remains available only for explicit legacy replay and controlled comparison.

## Mainline Data Roots

Mainline now resolves from repo-local paths only:

- `data/processed/essential_gene/...`
- `data/processed/PPI/...`
- `data/processed/OR/...`
- `data/processed/EXP/...`
- `data/processed/LC/...`
- `results/frozen_protocol/...`

`data_registry` is no longer part of the required mainline loader contract.

## Legacy Replay Files

These remain useful for replay, comparison, or migration provenance, but they are not the forward mainline benchmark entry points:

- `workflow/phase2a_fixed300_multirun.smk`
- `workflow/classical_baseline_benchmark.smk`
- `src/train/run_epgat_phase2a_single.py`
- `src/eval/aggregate_phase2a_multirun.py`
- `src/eval/analyze_phase2a_fixed300_multirun.py`
- `src/train/run_classical_baseline_single.py`
- `src/eval/aggregate_classical_baselines.py`
- `src/train/train_epgat_legacy.py`
- `src/train/train_epgat_graph_models.py`
- `src/eval/evaluate_epgat_legacy.py`
- `src/eval/evaluate_epgat_graph_models.py`
- `scripts/run_epgat_legacy.sh`

Reason:

- they preserve historical EPGAT or transition-stage benchmark behavior
- they may still rebuild labels or splits internally
- they do not define the frozen publication protocol

## Old Ranking-Only Fusarium Files

These files still encode deprecated ranking-only Fusarium assumptions and must not be used by the new mainline:

- `workflow/fgraminearum_feature_combo_benchmark.smk`
- `workflow/fgraminearum_feature_combo_newlabel_benchmark.smk`
- `scripts/run_fgraminearum_feature_combo_benchmark.sh`
- `scripts/run_fgraminearum_feature_combo_newlabel_benchmark.sh`
- `src/data/build_baseline_dataset.py`
- `src/data/gene_graph_adapter.py`
- `src/data/build_graph_manifest.py`
- `src/data/build_fusarium_graph_inference_inputs.py`
- `src/train/run_old440_diagnostic.py`
- `src/train/run_label_rebuild_compare.py`
- `src/train/run_fusarium_orthology_patch_compare.py`
- `src/eval/evaluate_graph_model.py`
- `src/eval/compare_old440_vs_newlabel.py`
- `configs/baseline.yaml`
- `configs/graph_prototype.yaml`
- `configs/label_rebuild_compare.yaml`

Deprecated ranking-only artifacts referenced by these paths:

- `broad79`
- `strict29`
- `conflict8`
- `fgraminearum_gold_positive*`
- implicit `old440` defaulting

## Safe Archive Candidates

These are strong archive candidates once the frozen protocol benchmark is fully rerun and validated:

- `workflow/fgraminearum_feature_combo_benchmark.smk`
- `workflow/fgraminearum_feature_combo_newlabel_benchmark.smk`
- `scripts/run_fgraminearum_feature_combo_benchmark.sh`
- `scripts/run_fgraminearum_feature_combo_newlabel_benchmark.sh`
- `src/eval/compare_old440_vs_newlabel.py`
- `src/train/run_old440_diagnostic.py`
- `src/train/run_label_rebuild_compare.py`
- `configs/label_rebuild_compare.yaml`

## Data-Registry Status

`data_registry` usage is now classified as follows:

- legacy-only:
  `src/data/build_fusarium_graph_inference_inputs.py`
  `src/data/build_support_feature_matrices.py`
  `src/data/build_support_graph_admission.py`
  `src/train/run_old440_diagnostic.py`
  `src/train/run_label_rebuild_compare.py`
- needs migration but not required by mainline:
  any replay path still reading `master_evidence_table.preliminary.tsv`
- archive-compatible:
  bridge and historical mapping tables retained for provenance

The only mirrored table created during this refactor is:

- `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`

This mirror supports provenance documentation but is not required by the new mainline loader.

## dmelanogaster Wiring Status

`dmelanogaster` is now wired in the new mainline protocol.

Verified available repo-local assets:

- `data/processed/essential_gene/melanogaster/labels.standard.tsv`
- `data/processed/PPI/melanogaster/string.csv`
- `data/processed/OR/melanogaster/orthologs.csv`
- `data/processed/EXP/melanogaster/profile.csv`
- `data/processed/LC/melanogaster/subloc.csv`

Verified new outputs:

- `results/frozen_protocol/labels/dmelanogaster_labels.tsv`
- `results/frozen_protocol/splits/dmelanogaster_split.tsv`
- smoke-tested model outputs under `outputs/frozen_protocol_benchmark_v2/dmelanogaster/...`

## Fusarium Role Separation

The refactor now makes the Fusarium role split explicit:

- `fgraminearum_newlabel` = intended mainline regime
- `fgraminearum_oldlabel` = legacy replay / historical comparison only

The new mainline no longer uses:

- ambiguous `fgraminearum_split.tsv`
- ambiguous `fgraminearum_labels.tsv`
- any implicit fallback from new labels to old440

## Recommendations

- Keep the new frozen protocol files as the only documented main benchmark path.
- Keep old Fusarium ranking and reconstruction code available only under explicit legacy documentation.
- Do not delete replay code until the full `outputs/frozen_protocol_benchmark_v2` benchmark has been rerun and archived summaries have been validated.
