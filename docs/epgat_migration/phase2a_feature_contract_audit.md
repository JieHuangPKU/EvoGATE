# Phase 2A Feature Contract Audit

## Status Note

- This document audits the archived replay benchmark feature contract.
- It is retained as provenance only.
- It no longer defines the formal meaning of Phase 2A after the 2026-04-07 scope reset.

## Scope

- Project: `ProGATE_v2`
- Date: 2026-04-07
- Objective: audit the legacy EPGAT feature contract only
- Explicitly excluded: any new training or benchmark reruns

## Executive Verdict

- The archived replay benchmark reuses the legacy replay node matrix directly from `outputs/epgat_legacy/*_original_replay/feature_matrix.npy` for human, celegans, and fgraminearum, plus `outputs/epgat_legacy/scerevisiae_original_smoke/feature_matrix.npy` for yeast.
- That replay matrix is a fixed four-block contract: `orthologs -> expression -> sublocalization -> degree`, then z-scored column-wise.
- Therefore, the archived graph benchmark is consistent with the raw legacy non-PLM runtime feature matrix.
- It is not fully consistent with the remembered old `final.csv` / `human_final.csv` interpretation, because the old published figure path is a post-hoc processed summary and the human `~0.91` point aligns more closely to `SUB_ORT` threshold-sweep rows than to the fixed full-contract `EXP_SUB_ORT + degree` replay condition.

## Short Answers

### 1. In old `final.csv`, what do `EXP / SUB / ORT` correspond to?

- `EXP` = `EPGAT/data/essential_genes/<species>/Expression/profile.csv`
  - Source columns: all columns except `Gene`
  - In ProGATE_v2 replay builder these are renamed to `expression_0 ... expression_n`
- `SUB` = `EPGAT/data/essential_genes/<species>/SubLocalizations/subloc.csv`
  - Source columns: all columns except `Gene`
  - In ProGATE_v2 replay builder these are renamed to `subloc_0 ... subloc_n`
- `ORT` = `EPGAT/data/essential_genes/<species>/Orthologs/orthologs.csv`
  - Source columns: all columns except `Gene`
  - Important file-format note: in the ortholog CSVs, `Gene` is the last column, not the first
  - In ProGATE_v2 replay builder these are renamed to `ortholog_0 ... ortholog_n`

### 2. What are the columns in current `outputs/epgat_legacy/*_original_replay/feature_matrix.npy`?

- Source assembly is implemented in [`src/data/build_epgat_legacy_dataset.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/data/build_epgat_legacy_dataset.py) and [`src/features/epgat_legacy_features.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/features/epgat_legacy_features.py).
- Column construction order is:
  1. all ortholog columns from `Orthologs/orthologs.csv`
  2. all expression columns from `Expression/profile.csv`
  3. all sublocalization columns from `SubLocalizations/subloc.csv`
  4. one explicit `degree_0` column computed from thresholded STRING PPI edges in `PPI/STRING/string.csv`
- Merge semantics:
  - left join onto the PPI/label gene universe
  - missing values filled with `0.0`
  - final numeric matrix z-scored per column

### 3. Is degree really part of the old final benchmark’s explicit feature set?

- Yes at runtime, with one caveat.
- It is not part of the `EXP / SUB / ORT` shorthand itself.
- But it is part of the actual node feature matrix consumed by the legacy runners, because legacy `runners/tools.py` appended degree after the omics blocks, and the replay builder preserves that behavior.
- So the most accurate statement is:
  - old shorthand names describe omics subset selection
  - old runtime feature matrix for those runs still included an appended degree column when using the shared legacy graph loader

### 4. Is the current migration benchmark definitionally identical to old `final.csv`?

- No.
- At the raw feature-matrix level, it is close to the legacy full-contract runtime path.
- At the result-provenance level, it is not the same object as old `final.csv`:
  - old figures were read from `outputs/results/human_final.csv`
  - that file was built by legacy `utils/prepare_data.py`
  - `prepare_data.py` aggregated multiple result cohorts and applied post-processing
  - the remembered human `~0.91` point aligns better with `GAT_SUB_ORT` threshold-sweep rows at STRING threshold 300 than with the current fixed-threshold full replay

## Evidence Chain

### A. Legacy runtime feature contract

- [`docs/epgat_migration/01_epgat_repo_audit.md`](/home/jiehuang/software/fungi/ProGATE_v2/docs/epgat_migration/01_epgat_repo_audit.md) records the legacy loader flow:
  - load PPI
  - load labels
  - optionally join `Orthologs/orthologs.csv`, `Expression/profile.csv`, `SubLocalizations/subloc.csv`
  - fill missing with zero
  - append degree in `runners/tools.py`
- Older audit work established that the archived replay benchmark used `orthologs, expression, sublocalization, then degree`.

### B. Current replay builder

- [`src/data/build_epgat_legacy_dataset.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/data/build_epgat_legacy_dataset.py) loads:
  - `Orthologs/orthologs.csv`
  - `Expression/profile.csv`
  - `SubLocalizations/subloc.csv`
  - `PPI/STRING/string.csv`
- It merges orthologs, expression, sublocalization, fills missing with zero, appends degree via `append_degree_block()`, then saves `feature_matrix.npy`.
- [`src/features/epgat_legacy_features.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/features/epgat_legacy_features.py) confirms the preserved order:
  - orthologs
  - expression
  - sublocalization
  - degree

### C. Current graph benchmark input path

- [`configs/epgat_graph_benchmark_human.yaml`](/home/jiehuang/software/fungi/ProGATE_v2/configs/epgat_graph_benchmark_human.yaml)
- [`configs/epgat_graph_benchmark_celegans.yaml`](/home/jiehuang/software/fungi/ProGATE_v2/configs/epgat_graph_benchmark_celegans.yaml)
- [`configs/epgat_graph_benchmark_fgraminearum.yaml`](/home/jiehuang/software/fungi/ProGATE_v2/configs/epgat_graph_benchmark_fgraminearum.yaml)
- [`configs/epgat_graph_benchmark_scerevisiae.yaml`](/home/jiehuang/software/fungi/ProGATE_v2/configs/epgat_graph_benchmark_scerevisiae.yaml)

These configs point the benchmark to:

- `outputs/epgat_legacy/human_original_replay`
- `outputs/epgat_legacy/celegans_original_replay`
- `outputs/epgat_legacy/fgraminearum_original_replay`
- `outputs/epgat_legacy/scerevisiae_original_smoke`

[`src/train/train_epgat_graph_models.py`](/home/jiehuang/software/fungi/ProGATE_v2/src/train/train_epgat_graph_models.py) then loads `feature_matrix.npy` from those directories directly, with no feature-contract rewrite.

### D. Why old `final.csv` is not the same evaluation target

- [`docs/epgat_migration/why_replay_diff_report.md`](/home/jiehuang/software/fungi/ProGATE_v2/docs/epgat_migration/why_replay_diff_report.md) documents that:
  - old plots read `outputs/results/human_final.csv`
  - `human_final.csv` came from `utils/prepare_data.py`
  - that script aggregated multiple result directories and applied post-processing
- [`outputs/epgat_legacy/human_original_replay/old_vs_replay_condition_diff.tsv`](/home/jiehuang/software/fungi/ProGATE_v2/outputs/epgat_legacy/human_original_replay/old_vs_replay_condition_diff.tsv) further shows:
  - the human `~0.91` neighborhood matches `GAT_SUB_ORT`
  - current replay is `EXP_SUB_ORT` with fixed threshold 500
  - degree itself is not the main mismatch

## Feature Source Summary

| species | benchmark source dir | ORT source | ORT dim | EXP source | EXP dim | SUB source | SUB dim | degree source | degree dim | replay matrix shape | audit conclusion |
| --- | --- | --- | ---: | --- | ---: | --- | ---: | --- | ---: | --- | --- |
| human | `outputs/epgat_legacy/human_original_replay` | `EPGAT/data/essential_genes/human/Orthologs/orthologs.csv` except trailing `Gene` | 162 | `EPGAT/data/essential_genes/human/Expression/profile.csv` except `Gene` | 64 | `EPGAT/data/essential_genes/human/SubLocalizations/subloc.csv` except `Gene` | 11 | degree count from thresholded `EPGAT/data/essential_genes/human/PPI/STRING/string.csv` | 1 | `(18822, 238)` | replay matrix matches raw legacy full-contract loader; not equivalent to processed `human_final.csv` |
| celegans | `outputs/epgat_legacy/celegans_original_replay` | `EPGAT/data/essential_genes/celegans/Orthologs/orthologs.csv` except trailing `Gene` | 273 | `EPGAT/data/essential_genes/celegans/Expression/profile.csv` except `Gene` | 72 | `EPGAT/data/essential_genes/celegans/SubLocalizations/subloc.csv` except `Gene` | 11 | degree count from thresholded `EPGAT/data/essential_genes/celegans/PPI/STRING/string.csv` | 1 | `(5766, 357)` | replay matrix matches raw legacy full-contract loader |
| fgraminearum | `outputs/epgat_legacy/fgraminearum_original_replay` | `EPGAT/data/essential_genes/fgraminearum/Orthologs/orthologs.csv` except trailing `Gene` | 274 | `EPGAT/data/essential_genes/fgraminearum/Expression/profile.csv` except `Gene` | 344 | `EPGAT/data/essential_genes/fgraminearum/SubLocalizations/subloc.csv` except `Gene` | 11 | degree count from thresholded `EPGAT/data/essential_genes/fgraminearum/PPI/STRING/string.csv` | 1 | `(6996, 630)` | replay matrix matches raw legacy full-contract loader |
| scerevisiae | `outputs/epgat_legacy/scerevisiae_original_smoke` | `EPGAT/data/essential_genes/yeast/Orthologs/orthologs.csv` except trailing `Gene` | 62 | `EPGAT/data/essential_genes/yeast/Expression/profile.csv` except `Gene` | 36 | `EPGAT/data/essential_genes/yeast/SubLocalizations/subloc.csv` except `Gene` | 11 | degree count from thresholded `EPGAT/data/essential_genes/yeast/PPI/STRING/string.csv` | 1 | `(6049, 110)` | graph benchmark uses the same raw-contract matrix, but from `original_smoke` rather than `original_replay` |

## Detailed Conclusion

### What is reproduced correctly

- The current ProGATE_v2 legacy replay builder faithfully reconstructs the raw legacy non-PLM node feature matrix contract:
  - STRING PPI defines graph topology
  - `ORT`, `EXP`, `SUB` come from the three legacy flat CSV tables
  - degree is appended explicitly as one numeric feature
  - feature order is preserved

### What is not reproduced by the current graph benchmark

- The benchmark does not reproduce the exact provenance of old `final.csv` / `human_final.csv`.
- In particular, for human:
  - the remembered old `~0.91` figure is not the same condition as the current fixed full replay
  - it is tied to threshold-sweep and/or processed summary outputs
  - the most likely old comparison point is `SUB_ORT`, not fixed `EXP_SUB_ORT`

## Bottom Line

- If the question is "did the archived replay benchmark use the same raw feature blocks as old legacy runtime?", the answer is yes.
- If the question is "did that archived replay benchmark reproduce the exact old `final.csv` result definition?", the answer is no.
- The contract mismatch is not mainly about degree.
- The main mismatch is result provenance and condition selection: processed `final.csv` and `SUB_ORT` threshold-sweep behavior versus current fixed full-contract replay inputs.
