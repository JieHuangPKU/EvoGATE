# Classical vs Frozen Parameter Audit

## Scope

This compares effective runtime behavior between:

- legacy `workflow/classical_baseline_benchmark.smk`
- `workflow/frozen_protocol_benchmark.smk`

and records the post-refactor resolution:

- `workflow/classical_baseline_benchmark.smk` now delegates to the frozen runtime contract.

Compared code:

- `configs/classical_baseline_benchmark.yaml`
- `configs/frozen_protocol.yaml`
- `src/train/run_classical_baseline_single.py`
- `src/network/run_network_heuristics.py`
- `src/train/run_frozen_protocol_model.py`
- `src/data/frozen_protocol_loader.py`
- `src/classical_baselines/common.py`

## Authority Decision

Authoritative mainline settings for overlapping models are the frozen settings, because:

- the frozen workflow is the explicit publication-grade protocol workflow
- it owns the current protocol list
- it owns the frozen label and split manifests
- it owns the current output and aggregation schema

Therefore:

- legacy classical parameter differences should be treated as legacy drift unless there is strong evidence they were intentional benchmark design
- the refactor updates classical benchmarking to the frozen settings by delegation

## General Parity

### Python Entry Point

- `classical_pre_refactor`
  - trainable: `src.train.run_classical_baseline_single`
  - deterministic: `src.network.run_network_heuristics`
- `frozen_protocol`
  - all overlapping models: `src.train.run_frozen_protocol_model`

Classification:

- `inconsistent_benchmark_breaking`

Reason:

- The entrypoint difference also implies different dataset loading and output schemas.

### Config Source

- classical: `configs/classical_baseline_benchmark.yaml`
- frozen: `configs/frozen_protocol.yaml`

Classification:

- `inconsistent_benchmark_breaking`

### Random Seed Handling

- classical:
  - `base_seed: 1029`
  - run IDs: `0,1,2`
  - effective seeds: `1029,1030,1031`
- frozen:
  - `seed_list: [1029,1030,1031,1032,1033]`

Classification:

- `inconsistent_benchmark_breaking`

### Train / Val / Test Split Source

- classical:
  - split is rebuilt per run through `build_epgat_legacy_dataset`
  - controlled by `seed`, `test_fraction=0.20`, `val_fraction=0.04`
- frozen:
  - split is loaded from `results/frozen_protocol/splits/*.tsv`
  - manifests were frozen by `src.data.freeze_unified_protocol`
  - config records `test_fraction=0.20`, `val_fraction=0.10`

Classification:

- `inconsistent_benchmark_breaking`

### Class Imbalance Handling

- MLP and N2V_MLP MLP head:
  - both use `BalancedBCE`
- RF:
  - both use `class_weight=balanced`
- SVM:
  - both use `class_weight=balanced`
- NB:
  - no explicit class weighting in either

Classification:

- `identical`

### Output Schema

- classical metrics row:
  - `species`
  - `method`
  - `feature_setting`
  - `label_regime`
  - `run_id`
  - `seed`
  - `test_count`
  - `val_auroc`
  - `val_auprc`
  - `val_mcc`
  - test metrics `auroc/auprc/mcc/f1/accuracy/specificity`
- frozen metrics row:
  - protocol metadata plus
  - `val_auroc`
  - `val_auprc`
  - `val_mcc`
  - `val_f1`
  - `val_accuracy`
  - `test_auroc`
  - `test_auprc`
  - `test_mcc`
  - `test_f1`
  - `test_accuracy`
  - counts
  - artifact references

Classification:

- `inconsistent_benchmark_breaking`

### Evaluation Metrics Written

- classical wrote `specificity` but did not write `val_f1`, `val_accuracy`, or protocol metadata.
- frozen writes explicit `test_*` and `val_*` metrics and richer metadata.

Classification:

- `inconsistent_benchmark_breaking`

## Classical Models

### MLP

Parameters:

- feature_setting:
  - classical: many feature sweeps
  - frozen: `ORT_EXP_SUB`
- hidden_dim:
  - `32` vs `32`
- dropout:
  - `0.2` vs `0.2`
- learning_rate:
  - `0.001` vs `0.001`
- weight_decay:
  - `0.0` vs `0.0`
- epochs:
  - `500` vs `200`
- patience:
  - `20` vs `20`
- early_stopping:
  - classical: implicit through patience break
  - frozen: explicit `early_stopping: true`
- preprocessing / scaling path:
  - classical uses legacy dataset builder output
  - frozen uses frozen protocol loader normalization over train split

Classification:

- same hyperparameters except epochs: `different_but_acceptable`
- effective benchmark behavior overall: `inconsistent_benchmark_breaking`

Reason:

- feature contract, split contract, and loader contract differ.

### RF

Parameters:

- feature_setting:
  - `ORT_EXP_SUB` vs `ORT_EXP_SUB`
- n_estimators:
  - `500` vs `500`
- class_weight:
  - `balanced` vs `balanced`
- preprocessing:
  - different loader and feature construction

Classification:

- model hyperparameters: `identical`
- effective benchmark behavior overall: `inconsistent_benchmark_breaking`

### SVM

Parameters:

- feature_setting:
  - `ORT_EXP_SUB` vs `ORT_EXP_SUB`
- kernel:
  - `rbf` vs `rbf`
- gamma:
  - `scale` vs `scale`
- C:
  - `1.0` vs `1.0`
- class_weight:
  - `balanced` vs `balanced`
- preprocessing:
  - different loader and feature construction

Classification:

- model hyperparameters: `identical`
- effective benchmark behavior overall: `inconsistent_benchmark_breaking`

### NB

Parameters:

- feature_setting:
  - `ORT_EXP_SUB` vs `ORT_EXP_SUB`
- preprocessing:
  - different loader and feature construction

Classification:

- model hyperparameters: `identical`
- effective benchmark behavior overall: `inconsistent_benchmark_breaking`

## N2V_MLP

Parameters:

- node2vec backend:
  - classical: no config-level backend passed
  - frozen: `runtime.node2vec_backend=auto`
- embedding_dim:
  - `64` vs `64`
- walk_length:
  - `32` vs `32`
- context_size:
  - `32` vs `32`
- walks_per_node:
  - `16` vs `16`
- num_negative_samples:
  - `1` vs `1`
- node2vec training epochs:
  - `50` vs `50`
- batch_size:
  - `256` vs `256`
- node2vec learning rate:
  - `0.005` vs `0.005`
- MLP hidden_dim:
  - `32` vs `32`
- MLP dropout:
  - `0.2` vs `0.2`
- MLP learning rate:
  - `0.001` vs `0.001`
- MLP weight_decay:
  - `0.0` vs `0.0`
- MLP epochs:
  - `500` vs `200`
- MLP patience:
  - `20` vs `20`
- early_stopping:
  - implicit vs explicit true
- input mix:
  - pure N2V in both

Classification:

- node2vec hyperparameters: mostly `identical`
- backend handling: `different_but_acceptable`
- MLP epochs: `different_but_acceptable`
- effective benchmark behavior overall: `inconsistent_benchmark_breaking`

Reason:

- split source, protocol scope, output schema, and loader contract still differed.

## Deterministic Models

### DC

Compared behavior:

- graph source:
  - classical: graph from legacy dataset bundle
  - frozen: graph from frozen protocol loader
- edge filtering threshold:
  - `300` vs `300`
- weighted/unweighted:
  - unweighted in both
- ranking-to-classification:
  - top K equals total labeled positives in both
- evaluation path:
  - different metrics schema and protocol metadata

Classification:

- centrality logic: `identical`
- effective benchmark behavior overall: `inconsistent_benchmark_breaking`

### CC

Same result as `DC`, replacing degree centrality with clustering coefficient.

Classification:

- centrality logic: `identical`
- effective benchmark behavior overall: `inconsistent_benchmark_breaking`

## Graph Models

### GraphSAGE

- classical workflow did not invoke it.
- frozen uses `src.train.run_frozen_protocol_model` and `EPGATOriginalSAGE` with:
  - `feature_setting: ORT_EXP_SUB`
  - `n_hidden: 64`
  - `n_layers: 2`
  - `dropout: 0.4`
  - `lr: 0.005`
  - `weight_decay: 0.0001`
  - `aggregator_type: mean`
  - `epochs: 200`
  - `patience: 20`
  - `early_stopping: true`

Classification:

- `not_applicable` for cross-workflow parity because classical did not include graph models.

### GCN

- classical workflow did not invoke it.
- frozen uses `EPGATOriginalGCN` with:
  - `hidden_dim: 64`
  - `dropout: 0.5`
  - `lr: 0.001`
  - `weight_decay: 0.0001`
  - `epochs: 200`
  - `patience: 20`
  - `early_stopping: true`

Classification:

- `not_applicable`

### GAT

- classical workflow did not invoke it.
- frozen uses `EPGATOriginalGAT` with:
  - species-conditioned hidden sizes in code:
    - human: `[16, 1]`
    - others: `[12, 1]`
  - heads `[8,1]`
  - `dropout: 0.3`
  - `lr: 0.005`
  - `weight_decay: 0.0002`
  - `epochs: 200`
  - `patience: 20`
  - `early_stopping: true`

Classification:

- `not_applicable`

### GIN

- classical workflow did not invoke it.
- frozen uses `EPGATOriginalGIN` with:
  - `dim_h: 64`
  - `dropout: 0.4`
  - `lr: 0.005`
  - `weight_decay: 0.0005`
  - `epochs: 200`
  - `patience: 20`
  - `early_stopping: true`

Classification:

- `not_applicable`

## Legacy Drift vs Intentional Differences

### Clear legacy drift

- separate config source
- only 3 runs instead of 5
- missing `dmelanogaster`
- implicit single Fusarium regime
- non-frozen split generation
- non-frozen summary schema
- MLP and N2V_MLP using `500` epochs while frozen uses `200`
- `include_degree: false` in classical config vs `true` in frozen config
- `val_fraction: 0.04` in classical vs `0.10` in frozen

These should be considered drift, not intentional benchmark design, because they move the workflow away from the stated frozen benchmark protocol.

### Acceptable implementation differences after convergence

- classical-only workflow retains a distinct output root:
  - `outputs/classical_baseline_benchmark_v2`
- classical-only workflow retains a distinct summary dir:
  - `results/classical_baseline_benchmark`
- graph models remain excluded from classical workflow by design

These are acceptable because they preserve workflow responsibility without changing benchmark semantics for overlapping models.

## Refactor Resolution

Post-refactor, `workflow/classical_baseline_benchmark.smk` now:

- uses `configs/frozen_protocol.yaml`
- uses all 6 frozen protocol targets
- uses the explicit 5-seed list
- uses `src.data.freeze_unified_protocol`
- uses `src.train.run_frozen_protocol_model`
- uses `src.eval.aggregate_frozen_protocol_runs`
- excludes graph models intentionally

Result for overlapping models:

- `MLP`, `RF`, `SVM`, `NB`, `N2V_MLP`, `DC`, and `CC` are now aligned with frozen runtime behavior
- the authoritative workflow for graph models remains `workflow/frozen_protocol_benchmark.smk`

## Recommendation

Going forward:

- all-model benchmark: use `workflow/frozen_protocol_benchmark.smk`
- classical-only benchmark: use `workflow/classical_baseline_benchmark.smk`

Rationale:

- the frozen workflow is the authoritative benchmark implementation
- the classical workflow should remain only a subset delegate, not an independent implementation branch
