# Data Root Migration

## Goal

Mainline protocol inputs must resolve from the ProGATE_v2 repository itself and must not require `data_registry` or hard-coded absolute paths.

## New Mainline Roots

The refactored mainline now resolves from:

- `data/processed/essential_gene/...`
- `data/processed/PPI/...`
- `data/processed/OR/...`
- `data/processed/EXP/...`
- `data/processed/LC/...`
- `results/frozen_protocol/...`

## Mainline Migration Result

The following mainline files no longer depend on `data_registry`:

- `configs/frozen_protocol.yaml`
- `src/data/freeze_unified_protocol.py`
- `src/data/frozen_protocol_loader.py`
- `src/train/run_frozen_protocol_model.py`
- `workflow/frozen_protocol_benchmark.smk`

Mainline now uses repo-local processed tables and frozen manifests only.

## Legacy-Only `data_registry` Usage

The following code paths still read `data_registry` and are therefore legacy-only or provenance-only:

- `src/data/build_fusarium_graph_inference_inputs.py`
- `src/data/build_support_feature_matrices.py`
- `src/data/build_support_graph_admission.py`
- `src/train/run_old440_diagnostic.py`
- `src/train/run_label_rebuild_compare.py`
- `configs/label_rebuild_compare.yaml`

Classification:

- legacy-only:
  ranking-era graph inference and support-graph utilities
- needs migration:
  provenance rebuild code for old label reconstruction
- archive-compatible:
  bridge tables and historical registry snapshots

## Mirrored Files

To keep provenance auditable inside the repo, this refactor created:

- `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`

Role:

- mirror of `data_registry/master_evidence_table.preliminary.tsv`
- supports provenance documentation of the Fusarium new-label negative filtering logic
- not required by the new mainline loader

## Remaining Absolute-Path Issue

One legacy config still contains a hard-coded absolute path:

- `configs/label_rebuild_compare.yaml`
- key: `paths.yeast_labels`

Current value points outside the repository and should not be used in mainline.

Recommended migration target:

- `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`

## Mainline Policy After Refactor

- No new mainline code should read from `data_registry`.
- No new mainline code should depend on absolute external paths.
- Any remaining registry-based or absolute-path workflow is legacy replay only.
