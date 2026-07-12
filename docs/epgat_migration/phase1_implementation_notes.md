# Phase 1 Implementation Notes

## New Files

- `src/data/epgat_legacy_id_adapter.py`
- `src/data/epgat_legacy_label_adapter.py`
- `src/features/epgat_legacy_features.py`
- `src/data/build_epgat_legacy_dataset.py`
- `src/models/epgat_original.py`
- `src/train/train_epgat_legacy.py`
- `src/eval/evaluate_epgat_legacy.py`
- `configs/epgat_legacy.yaml`
- `scripts/run_epgat_legacy.sh`

## Reused Legacy Sources

- `/home/jiehuang/software/fungi/EPGAT/utils/utils.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/tools.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/run_gat.py`
- `/home/jiehuang/software/fungi/EPGAT/models/gat/gat_pytorch.py`
- `/home/jiehuang/software/fungi/EPGAT/models/gat/params.py`

## Faithfully Migrated Logic

- original feature contract
- original feature ordering observed in `utils.data`
- original GAT layer behavior
- original split pattern: train/test then train/val stratified split
- original result semantics: AUC, AUPR, F1, Accuracy, MCC

## Adapterized Logic

- explicit node manifest
- explicit feature schema
- explicit label manifest
- explicit alignment audit
- unified output directory under `outputs/epgat_legacy/<run_name>/`

## Deliberately Not Done

- no ESM / ESMC / ProtT5 / ProtBERT migration
- no extended layer
- no unified full runner migration
- no Fusarium full canonical integration replay
- no batch matrix orchestration

## Current Implementation Status

- `build_epgat_legacy_dataset.py`: implemented
- `epgat_original.py`: implemented
- `train_epgat_legacy.py`: implemented
- `evaluate_epgat_legacy.py`: implemented
- `export_epgat_legacy_results.py`: implemented
- `configs/epgat_legacy.yaml`: implemented
- `scripts/run_epgat_legacy.sh`: implemented

## Real Smoke Test

- smoke species: `scerevisiae`
- output root: `outputs/epgat_legacy/scerevisiae_original_smoke`
- legacy-compatible summary export: `outputs/epgat_legacy/GAT_results.csv`
