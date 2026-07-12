# Phase 1 Smoke Test

## Selected Species

- `scerevisiae`

## Why This Species

- all four original feature blocks are present
- smaller than human and Fusarium
- original `gat_yeast` parameter preset exists
- suitable for a fast original-compatible smoke test

## Inputs Used

- `EPGAT/data/essential_genes/yeast/EssentialGenes/ogee.csv`
- `EPGAT/data/essential_genes/yeast/PPI/STRING/string.csv`
- `EPGAT/data/essential_genes/yeast/Expression/profile.csv`
- `EPGAT/data/essential_genes/yeast/Orthologs/orthologs.csv`
- `EPGAT/data/essential_genes/yeast/SubLocalizations/subloc.csv`

## Smoke Test Targets

- dataset builder success
- original-compatible GAT training success
- metrics / predictions export success

## Actual Smoke Test Status

- species/cohort: `scerevisiae_original_smoke`
- dataset build: success
- training: success
- evaluation export: success
- legacy-compatible exporter: success

## Real Output Snapshot

- output_dir: `outputs/epgat_legacy/scerevisiae_original_smoke`
- feature dimension: `110`
- test_count: `1091`
- AUROC: `0.8079`
- AUPRC: `0.5331`
- Accuracy: `0.6819`
- F1: `0.4949`
- MCC: `0.3664`

## Current Known Limits

- smoke test uses only one support species
- epochs reduced for fast validation
- legacy-compatible exporter is minimal, not a full old results tree clone
- no ESM / ESMC / extended layer in Phase 1
