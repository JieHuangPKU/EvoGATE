# Phase 2A Official Definition

## Status

- Phase 2A no longer means the 4-species x 4-model replay benchmark.
- Effective 2026-04-07, the old replay benchmark is archive-only and cannot be cited as the formal Phase 2A result line.
- Phase 2A now means: reproduce the old EPGAT original result contract, starting with human.

## Formal Goal

- Primary species: human
- Primary target family: legacy GAT threshold-sweep results in `EPGAT/results/GAT_results-string.csv`
- Immediate target point: `GAT_SUB_ORT` at `string_thr=300`
- Immediate numerical objective: approach the old human AUROC `~0.91`, concretely `0.9141008214`

## What Counts As Phase 2A

- Matching the old result-table condition, not the old replay runtime matrix by itself
- Respecting the legacy feature subset encoded by the target row:
  - `EXP`
  - `SUB`
  - `ORT`
  - and their combinations such as `SUB_ORT`, `EXP_ORT`, `EXP_SUB`, `EXP_SUB_ORT`
- Respecting legacy threshold sweep semantics from `GAT_results-string.csv`
- Using human as the first formal reproduction target

## What No Longer Counts As Phase 2A

- `outputs/epgat_graph_benchmark/*`
- the previous 4-species x 4-model comparison
- fixed full-contract replay inputs used as a proxy for old final results
- any benchmark that assumes `degree` is on by default without old-result evidence

## Degree Policy

- Degree is no longer a default Phase 2A feature.
- It must be enabled only when the target old result definition explicitly requires it.
- The current formal human Phase 2A configs keep `include_degree: false`.

## Official Inputs And Entry Points

- Single target reproduction config:
  - [`configs/phase2a_human_old_result.yaml`](/home/jiehuang/software/fungi/ProGATE_v2/configs/phase2a_human_old_result.yaml)
- Human threshold-sweep config:
  - [`configs/phase2a_human_threshold_sweep.yaml`](/home/jiehuang/software/fungi/ProGATE_v2/configs/phase2a_human_threshold_sweep.yaml)
- Single target runner:
  - [`src/train/run_epgat_human_exact.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_epgat_human_exact.py)
- Threshold-sweep runner:
  - [`src/train/run_epgat_human_threshold_sweep.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_epgat_human_threshold_sweep.py)
- Feature builder with explicit subset flags:
  - [`src/data/build_epgat_legacy_dataset.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/data/build_epgat_legacy_dataset.py)

## Current Formal Target

- Source table: `EPGAT/results/GAT_results-string.csv`
- Species: `human`
- Row name: `GAT_SUB_ORT`
- Feature subset: `SUB_ORT`
- STRING threshold: `300`
- Aggregation level: `mean_over_2_runs_from_old_string_threshold_table`
- Degree: disabled unless separately proven necessary for this result line

## Archive Notice

- The earlier replay benchmark remains useful as an engineering archive only:
  - feature-contract stress test
  - cross-model migration smoke benchmark
  - implementation validation for GAT/GCN/GIN/GraphSAGE loaders
- It is not an acceptable proxy for old EPGAT final-result reproduction.
