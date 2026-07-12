# Classical vs Frozen Workflow Audit

## Scope

Compared code paths:

- `workflow/classical_baseline_benchmark.smk`
- `workflow/frozen_protocol_benchmark.smk`
- `configs/classical_baseline_benchmark.yaml`
- `configs/frozen_protocol.yaml`
- `src/train/run_classical_baseline_single.py`
- `src/network/run_network_heuristics.py`
- `src/train/run_frozen_protocol_model.py`
- `src/data/frozen_protocol_loader.py`
- `src/classical_baselines/common.py`
- `src/eval/aggregate_classical_baselines.py`
- `src/eval/aggregate_frozen_protocol_runs.py`
- `scripts/run_classical_baseline_benchmark.sh`
- `scripts/run_frozen_protocol_benchmark.sh`

This audit distinguishes:

- `classical_pre_refactor`: the original legacy classical workflow.
- `classical_post_refactor`: the updated classical workflow implemented in this change.
- `frozen_protocol`: the existing authoritative all-model workflow.

## A. Workflow Comparison

### Configfile

- `classical_pre_refactor` used `configs/classical_baseline_benchmark.yaml`.
- `frozen_protocol` uses `configs/frozen_protocol.yaml`.
- `classical_post_refactor` now also uses `configs/frozen_protocol.yaml`.

Implementation consequence:

- Before refactor, classical ran from a separate parameter space and split contract.
- After refactor, classical inherits the same protocol definitions, model settings, runtime settings, and frozen manifests as the mainline workflow.

### Protocols / Species Coverage

- `classical_pre_refactor` hard-coded `["human", "celegans", "scerevisiae", "fgraminearum"]`.
- It omitted `dmelanogaster`.
- It represented Fusarium as one implicit species entry with `label_regime: new_label`.
- `frozen_protocol` reads `workflow.protocol_order` from `configs/frozen_protocol.yaml`:
  - `human`
  - `celegans`
  - `scerevisiae`
  - `dmelanogaster`
  - `fgraminearum_oldlabel`
  - `fgraminearum_newlabel`
- `classical_post_refactor` now uses the same 6 frozen protocol targets.

Conclusion:

- Pre-refactor classical protocol scope was benchmark-breaking relative to frozen.
- Post-refactor classical has full required 6-target coverage with explicit Fusarium old/new separation.

### Seed Expansion

- `classical_pre_refactor` used `RUN_IDS = [0, 1, 2]` and derived actual seed as `base_seed + run_id`.
- Effective seeds were `1029`, `1030`, `1031`.
- `frozen_protocol` uses `runtime.seed_list = [1029, 1030, 1031, 1032, 1033]`.
- `classical_post_refactor` now uses the same explicit 5-seed list.

Conclusion:

- Pre-refactor classical had only 3 seeds.
- Post-refactor classical supports the required 5 seeds.

### Rules

`classical_pre_refactor` rules:

- `prepare_new_labels`
- `run_trainable_baseline`
- `run_network_heuristic`
- `aggregate_trainable_feature_setting`
- `aggregate_network_feature_setting`
- `final_summary`

`frozen_protocol` rules:

- `freeze_protocol`
- `run_trainable_model`
- `run_deterministic_model`
- `aggregate_frozen_protocol`
- `all`

`classical_post_refactor` rules:

- `freeze_protocol`
- `run_classical_trainable_model`
- `run_classical_deterministic_model`
- `aggregate_classical_protocol`
- `all`

Conclusion:

- Pre-refactor classical had independent orchestration and independent manifest generation logic.
- Post-refactor classical is a thin subset wrapper over the frozen protocol contract.

### Python Modules Called

`classical_pre_refactor`:

- `src.eval.prepare_phase2b_new_labels`
- `src.train.run_classical_baseline_single`
- `src.network.run_network_heuristics`
- `src.eval.aggregate_classical_baselines`

`frozen_protocol`:

- `src.data.freeze_unified_protocol`
- `src.train.run_frozen_protocol_model`
- `src.eval.aggregate_frozen_protocol_runs`

`classical_post_refactor`:

- `src.data.freeze_unified_protocol`
- `src.train.run_frozen_protocol_model`
- `src.eval.aggregate_frozen_protocol_runs`

Conclusion:

- The refactor removes legacy classical-only training and aggregation entrypoints from the benchmark workflow path.

### Outputs

`classical_pre_refactor` output root:

- `outputs/classical_baseline_benchmark/{species}/{method}/{feature}/run_{run_id}/...`
- heuristics used `{species}/{method}/network/...`

`frozen_protocol` output root:

- `outputs/frozen_protocol_benchmark_v2/{protocol}/{model}/{feature}/run_{seed}/...`
- deterministic models use `{protocol}/{model}/{feature}/deterministic/...`

`classical_post_refactor` output root:

- `outputs/classical_baseline_benchmark_v2/{protocol}/{model}/{feature}/run_{seed}/...`
- deterministic models use `{protocol}/{model}/{feature}/deterministic/...`

Conclusion:

- Pre-refactor classical used a different directory contract and lower-case `network`.
- Post-refactor classical uses the frozen-style output contract with its own output root.

### Aggregation

`classical_pre_refactor`:

- Aggregated per feature group with `src.eval.aggregate_classical_baselines`.
- Final summary grouped by `species`, `method`, `feature_setting`, `label_regime`.
- Metrics schema used `auroc`, `auprc`, `mcc`, `f1`, `accuracy`, `specificity`.

`frozen_protocol`:

- Aggregates all per-run metrics with `src.eval.aggregate_frozen_protocol_runs`.
- Groups by protocol metadata including:
  - `protocol`
  - `species`
  - `regime`
  - `model`
  - `feature_setting`
  - `label_regime`
  - `split_version`
  - `graph_source`
  - `label_manifest`
  - `split_manifest`
  - `config_used`
  - `is_deterministic`

`classical_post_refactor`:

- Uses `src.eval.aggregate_frozen_protocol_runs`.
- Produces frozen-compatible `per_run_metrics.tsv`, `aggregated_metrics.tsv`, and `final_summary.tsv`.

Conclusion:

- Pre-refactor classical summary schema was not frozen-compatible.
- Post-refactor classical summary output is frozen-compatible.

### Logging

`classical_pre_refactor`:

- Trainable log: `run_classical_baseline_single.log`
- Heuristic log: `run_network_heuristics.log`
- Aggregation steps did not persist rich metadata per run beyond local resolved YAML and feature summary files.

`frozen_protocol`:

- Uses `run_frozen_protocol_model.log`.
- Persists richer run outputs:
  - `resolved_config.yaml`
  - `predictions.tsv`
  - `feature_schema.tsv`
  - `edge_table.tsv`
  - `split_manifest.tsv`
  - `training_log.tsv` when applicable
  - `best_model.pt` when applicable

`classical_post_refactor`:

- Uses the same `run_frozen_protocol_model.log` and same run artifact contract as frozen.

Conclusion:

- Pre-refactor diagnosability was materially weaker and not schema-aligned.
- Post-refactor classical has the same diagnostics as frozen for overlapping models.

### Rerun / Keep-going / Incomplete Behavior

`scripts/run_classical_baseline_benchmark.sh` pre-refactor:

- did not use `--keep-going`
- manually enumerated legacy matrix targets from `configs/classical_baseline_benchmark.yaml`
- then ran a second Snakemake invocation for `results/classical_baseline_benchmark/final_summary.tsv`

`scripts/run_frozen_protocol_benchmark.sh`:

- uses `--rerun-incomplete`
- uses `--keep-going`
- runs a single workflow invocation against the Snakefile default target

`scripts/run_classical_baseline_benchmark.sh` post-refactor:

- now mirrors the frozen invocation style
- uses `--rerun-incomplete`
- uses `--keep-going`
- runs the default target of `workflow/classical_baseline_benchmark.smk`

Conclusion:

- Pre-refactor classical execution behavior was less robust and more brittle.
- Post-refactor classical execution behavior matches frozen.

### Frozen Manifest Dependency

- `classical_pre_refactor` did not use frozen manifests.
- For non-Fusarium species, it rebuilt splits inside `build_epgat_legacy_dataset`.
- For Fusarium, it called `prepare_phase2b_new_labels` and then also used legacy dataset-building logic.
- `frozen_protocol` uses `src.data.freeze_unified_protocol` plus `src.data.frozen_protocol_loader`.
- `classical_post_refactor` now depends on the same frozen label and split manifests as frozen.

Conclusion:

- Pre-refactor classical depended on independent split/label generation.
- Post-refactor classical reuses frozen labels and splits as required.

### Same Output Root Style

- `classical_pre_refactor`: no.
- `frozen_protocol` and `classical_post_refactor`: yes.

The only intended difference is root path:

- all-model benchmark: `outputs/frozen_protocol_benchmark_v2`
- classical-only benchmark: `outputs/classical_baseline_benchmark_v2`

## B. Model-Family Parity Audit

### Overlapping baseline / embedding / heuristic models

#### MLP

- Pre-refactor classical and frozen did not match.
- Different workflow entrypoints:
  - classical: `src.train.run_classical_baseline_single`
  - frozen: `src.train.run_frozen_protocol_model`
- Different dataset loaders:
  - classical: `src.classical_baselines.common.build_dataset_for_benchmark` -> `src.data.build_epgat_legacy_dataset.build_dataset`
  - frozen: `src.data.frozen_protocol_loader.load_protocol_dataset`
- Different protocol scope, split source, and output schema.
- Different feature coverage:
  - classical benchmark mixed feature ablations
  - frozen benchmark uses `ORT_EXP_SUB`
- Same hidden dimension, dropout, LR, weight decay, patience.
- Different epochs:
  - classical: `500`
  - frozen: `200`

Parity verdict:

- `classical_pre_refactor` vs frozen: `inconsistent_benchmark_breaking`
- `classical_post_refactor` vs frozen: `identical`

#### RF

- Same sklearn family and same model hyperparameters.
- Different entrypoints, dataset loader, split source, protocol scope, output schema.

Parity verdict:

- pre-refactor: `inconsistent_benchmark_breaking`
- post-refactor: `identical`

#### SVM

- Same `C`, `kernel`, `gamma`, `class_weight`.
- Different entrypoints, dataset loader, split source, protocol scope, output schema.

Parity verdict:

- pre-refactor: `inconsistent_benchmark_breaking`
- post-refactor: `identical`

#### NB

- Same GaussianNB implementation.
- Different entrypoints, dataset loader, split source, protocol scope, output schema.

Parity verdict:

- pre-refactor: `inconsistent_benchmark_breaking`
- post-refactor: `identical`

#### N2V_MLP

- Both use `src.graph.run_node2vec_embedding.train_node2vec_embeddings`.
- Same embedding and MLP hyperparameters except MLP epochs:
  - classical: `mlp_epochs=500`
  - frozen: `mlp_epochs=200`
- Frozen additionally passes `backend=config["runtime"]["node2vec_backend"]`.
- Different dataset loader, split source, protocol scope, output schema.

Parity verdict:

- pre-refactor: `inconsistent_benchmark_breaking`
- post-refactor: `identical`

#### DC

- Pre-refactor classical used `src.network.run_network_heuristics`.
- Frozen uses deterministic branch inside `src.train.run_frozen_protocol_model`.
- Core logic is the same:
  - build unweighted NetworkX graph from edge list
  - compute `nx.degree_centrality`
  - predict top-K positives where K equals total labeled positives
- But graph source and labeled universe differed because the dataset loader differed.

Parity verdict:

- pre-refactor: logic similar, benchmark contract inconsistent
- post-refactor: identical

#### CC

- Same conclusion as `DC`, using `nx.clustering`.

Parity verdict:

- pre-refactor: logic similar, benchmark contract inconsistent
- post-refactor: identical

### Graph models

#### GraphSAGE, GCN, GAT, GIN

- `classical_pre_refactor` did not invoke them at all.
- `frozen_protocol` invokes all four through `src.train.run_frozen_protocol_model`.
- Therefore there is no same-way invocation across the two workflows because classical never covered graph models.
- `classical_post_refactor` intentionally still excludes graph models.

Parity verdict:

- between original classical and frozen: `not the same`, because classical does not invoke them
- going forward: graph models remain authoritative only in `frozen_protocol_benchmark.smk`

## C. Exact Answers to the Main Questions

### 1. How do the two workflows differ?

Protocol scope:

- Classical pre-refactor covered 4 species and fused Fusarium into one implicit regime.
- Frozen covers all 6 required protocol targets.

Seed handling:

- Classical pre-refactor used 3 run IDs and derived seeds.
- Frozen uses 5 explicit seeds.

Model coverage:

- Classical pre-refactor covered only `MLP`, `RF`, `SVM`, `NB`, `N2V_MLP`, `DC`, `CC`.
- Frozen covers those plus `GraphSAGE`, `GCN`, `GAT`, `GIN`.

Output structure:

- Classical pre-refactor used `{species}/{method}/{feature}/run_{run_id}` and `{species}/{method}/network`.
- Frozen uses `{protocol}/{model}/{feature}/run_{seed}` and deterministic subdirectories.

Aggregation logic:

- Classical pre-refactor used `aggregate_classical_baselines`.
- Frozen uses `aggregate_frozen_protocol_runs`.

Dependency on frozen labels/splits:

- Classical pre-refactor: no.
- Frozen: yes.

Environment setup:

- Classical pre-refactor hard-coded `EPGAT_PYTHON`, `MPLCONFIGDIR`, and `XDG_CACHE_HOME` in Snakefile.
- Frozen reads runtime environment paths from config.

Logging / diagnosability:

- Frozen persists richer metadata and a normalized run artifact set.
- Classical pre-refactor was thinner.

Rerun / keep-going / incomplete behavior:

- Frozen runner used `--keep-going --rerun-incomplete`.
- Classical runner pre-refactor only used `--rerun-incomplete` and split execution into two invocations.

### 2. Are GraphSAGE / GCN / GAT / GIN invoked in the same way?

- No.
- `classical_pre_refactor` did not invoke any of them.
- `frozen_protocol` invokes all of them through `src.train.run_frozen_protocol_model`.
- `classical_post_refactor` intentionally continues to exclude them.

### 3. Are MLP / RF / SVM / NB / N2V_MLP / DC / CC invoked in the same way?

- Pre-refactor: no.
- Post-refactor: yes.

Pre-refactor exact differences:

- different Python entry points
- different config source
- different split sources
- different protocol scope
- different output directory conventions
- different summary schemas
- different Fusarium handling

### 4. Where exactly are the differences?

Different Python entry points:

- `src.train.run_classical_baseline_single` vs `src.train.run_frozen_protocol_model`
- `src.network.run_network_heuristics` vs deterministic branch in `src.train.run_frozen_protocol_model`
- `src.eval.aggregate_classical_baselines` vs `src.eval.aggregate_frozen_protocol_runs`

Different config keys:

- classical had `runtime.base_seed`, `species.*`, `methods.*`
- frozen has `runtime.seed_list`, `protocols.*`, `workflow.*`, `models.*`

Different feature settings:

- classical benchmark path let MLP sweep multiple feature combinations
- frozen benchmark path pins models to one benchmark feature setting per model

Different split sources:

- classical built splits dynamically through legacy dataset builder
- frozen consumes frozen split manifests

Different output directory conventions:

- classical: species-oriented
- frozen: protocol-oriented

Different summary schemas:

- classical: species/method/feature with classical metrics
- frozen: protocol-aware schema with frozen metadata

Different Fusarium handling:

- classical: one implicit `fgraminearum` species with new_label
- frozen: explicit `fgraminearum_oldlabel` and `fgraminearum_newlabel`

Different logging behavior:

- frozen records richer metadata artifacts and standardized logs

### 5. Cleanest way to make classical support the required behavior

The cleanest implementation is:

- keep `classical_baseline_benchmark.smk` as a classical-only workflow
- make it a thin wrapper over frozen protocol machinery
- reuse:
  - `configs/frozen_protocol.yaml`
  - `src.data.freeze_unified_protocol`
  - `src.train.run_frozen_protocol_model`
  - `src.eval.aggregate_frozen_protocol_runs`
- exclude graph models intentionally
- give classical its own output root and summary dir so it does not collide with the all-model benchmark

That is exactly what was implemented.

### 6. Should classical stay separate, be upgraded, or delegate?

Recommendation:

- keep it separate as a classical-only workflow
- but make it a thin delegate to the frozen machinery

Why:

- classical-only benchmarking remains useful operationally
- the frozen workflow should remain authoritative for all-model benchmarking
- duplicating split, label, model-runner, and aggregation logic would reintroduce drift

## Recommended Refactor

Recommendation category:

- `wrapper/delegate approach`

Implementation reason:

- The frozen runner already owns the authoritative implementations for:
  - `MLP`
  - `RF`
  - `SVM`
  - `NB`
  - `N2V_MLP`
  - `DC`
  - `CC`
  - graph models
- The frozen loader already owns:
  - explicit protocol names
  - frozen manifests
  - current graph source
  - current node universe
  - current feature construction
- The frozen aggregator already owns:
  - frozen-compatible summary schema

Therefore the maintainable change is not to modernize the old classical stack in place, but to retire it from benchmark orchestration and delegate to the frozen stack.

## Implemented Refactor

Modified files:

- `workflow/classical_baseline_benchmark.smk`
- `configs/frozen_protocol.yaml`
- `scripts/run_classical_baseline_benchmark.sh`
- `docs/protocol_refactor/classical_vs_frozen_workflow_audit.md`
- `docs/protocol_refactor/classical_vs_frozen_workflow_audit.tsv`
- `docs/protocol_refactor/classical_vs_frozen_parameter_audit.md`
- `docs/protocol_refactor/classical_vs_frozen_parameter_audit.tsv`

Behavior after refactor:

- supports 5 seeds
- supports all 6 protocol targets
- keeps `fgraminearum_oldlabel` and `fgraminearum_newlabel` explicit
- uses frozen labels and splits
- uses frozen-compatible output schema and aggregation
- intentionally excludes graph models from the classical workflow

## Validation

Dry-run commands used:

```bash
XDG_CACHE_HOME=/tmp/codex_cache TMPDIR=/tmp /home/jiehuang/anaconda3/bin/snakemake -s workflow/classical_baseline_benchmark.smk --configfile configs/frozen_protocol.yaml --cores 1 --dry-run --rerun-incomplete --keep-going
XDG_CACHE_HOME=/tmp/codex_cache TMPDIR=/tmp /home/jiehuang/anaconda3/bin/snakemake -s workflow/frozen_protocol_benchmark.smk --configfile configs/frozen_protocol.yaml --cores 1 --dry-run --rerun-incomplete --keep-going
```

Results:

- classical dry-run succeeded and expanded:
  - `150` trainable jobs
  - `12` deterministic jobs
  - `1` aggregate job
  - `1` final target
- frozen dry-run succeeded and reported the DAG up to date.

Interpretation of classical counts:

- `6 protocols * 5 trainable classical models * 5 seeds = 150`
- `6 protocols * 2 deterministic models = 12`

## Going Forward

Use:

- all-model benchmark: `workflow/frozen_protocol_benchmark.smk`
- classical-only benchmark: `workflow/classical_baseline_benchmark.smk`

Authoritative implementations:

- graph models: frozen workflow only
- overlapping classical models: frozen runner and frozen config
