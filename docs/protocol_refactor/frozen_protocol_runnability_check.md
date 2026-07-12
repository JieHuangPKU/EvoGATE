# Frozen Protocol Runnability Check

Date: 2026-04-10
Repo root: `/home/jiehuang/software/fungi/ProGATE_v2`

## Scope

This check validates mainline runnability of the current frozen protocol workflow without redesigning the workflow.

Validated areas:

- authoritative entry point
- current config expansion scope
- fresh workflow dry-run
- fresh workflow-driven representative target reruns
- `N2V_MLP` path and output correctness
- aggregation module output generation
- per-run log file creation

## Entry Point

Authoritative full-run entry point remains:

```bash
./scripts/run_frozen_protocol_benchmark.sh 32
```

Verified in [scripts/run_frozen_protocol_benchmark.sh](/home/jiehuang/software/fungi/ProGATE_v2/scripts/run_frozen_protocol_benchmark.sh):

- changes directory to repo root
- uses `workflow/frozen_protocol_benchmark.smk`
- uses `configs/frozen_protocol.yaml`
- creates local `.mplconfig` and `.cache`
- exports `MPLBACKEND=Agg`
- exports `MPLCONFIGDIR=$PROJECT_ROOT/.mplconfig`
- exports `XDG_CACHE_HOME=$PROJECT_ROOT/.cache`
- passes `--rerun-incomplete`
- passes `--keep-going`
- activates `conda activate EPGAT` if `/home/jiehuang/anaconda3/etc/profile.d/conda.sh` is available
- invokes `/home/jiehuang/anaconda3/bin/snakemake`

Availability check:

- `/home/jiehuang/anaconda3/etc/profile.d/conda.sh`: present
- `/home/jiehuang/anaconda3/bin/snakemake`: executable
- `/home/jiehuang/anaconda3/envs/EPGAT/bin/python`: executable

## Config Scope

Verified in [configs/frozen_protocol.yaml](/home/jiehuang/software/fungi/ProGATE_v2/configs/frozen_protocol.yaml):

- protocols:
  - `human`
  - `celegans`
  - `scerevisiae`
  - `dmelanogaster`
  - `fgraminearum_oldlabel`
  - `fgraminearum_newlabel`
- trainable models:
  - `MLP`
  - `RF`
  - `SVM`
  - `NB`
  - `N2V_MLP`
  - `GAT`
  - `GCN`
  - `GIN`
  - `GraphSAGE`
- deterministic models:
  - `DC`
  - `CC`
- seeds:
  - `1029`
  - `1030`
  - `1031`
  - `1032`
  - `1033`

Training budget remains intact:

- `MLP`: `epochs=200`, `patience=20`, `early_stopping=true`
- `N2V_MLP` MLP head: `mlp_epochs=200`, `mlp_patience=20`, `early_stopping=true`
- `GAT`: `epochs=200`, `patience=20`, `early_stopping=true`
- `GCN`: `epochs=200`, `patience=20`, `early_stopping=true`
- `GIN`: `epochs=200`, `patience=20`, `early_stopping=true`
- `GraphSAGE`: `epochs=200`, `patience=20`, `early_stopping=true`

## Dry-Run

Fresh dry-run executed via the main entry script:

```bash
./scripts/run_frozen_protocol_benchmark.sh 32 -n
```

Result:

- exit code `0`
- full DAG built successfully
- DAG size at validation time: `285` jobs
  - `freeze_protocol`: `1`
  - `run_trainable_model`: `270`
  - `run_deterministic_model`: `12`
  - `aggregate_frozen_protocol`: `1`
  - `all`: `1`

## Fresh Workflow-Driven Validation Targets

To keep validation scoped while staying workflow-driven, representative targets were rerun directly through Snakemake with explicit target files and restricted rule sets.

### Deterministic

Command:

```bash
export MPLBACKEND=Agg MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig XDG_CACHE_HOME=/home/jiehuang/software/fungi/ProGATE_v2/.cache
/home/jiehuang/anaconda3/bin/snakemake -s workflow/frozen_protocol_benchmark.smk --configfile configs/frozen_protocol.yaml --cores 1 --nolock --allowed-rules freeze_protocol run_deterministic_model --forcerun outputs/frozen_protocol_benchmark_v2/human/DC/NETWORK/deterministic/metrics.tsv outputs/frozen_protocol_benchmark_v2/human/DC/NETWORK/deterministic/metrics.tsv
```

Result:

- success
- output: [metrics.tsv](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/human/DC/NETWORK/deterministic/metrics.tsv)
- log exists: [run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/human/DC/NETWORK/deterministic/run_frozen_protocol_model.log)

### Classical Trainable

Command:

```bash
export MPLBACKEND=Agg MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig XDG_CACHE_HOME=/home/jiehuang/software/fungi/ProGATE_v2/.cache
/home/jiehuang/anaconda3/bin/snakemake -s workflow/frozen_protocol_benchmark.smk --configfile configs/frozen_protocol.yaml --cores 1 --nolock --allowed-rules freeze_protocol run_trainable_model --forcerun outputs/frozen_protocol_benchmark_v2/celegans/MLP/ORT_EXP_SUB/run_1031/metrics.tsv outputs/frozen_protocol_benchmark_v2/celegans/MLP/ORT_EXP_SUB/run_1031/metrics.tsv
```

Result:

- success
- output: [metrics.tsv](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/celegans/MLP/ORT_EXP_SUB/run_1031/metrics.tsv)
- log exists: [run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/celegans/MLP/ORT_EXP_SUB/run_1031/run_frozen_protocol_model.log)

### N2V_MLP

Command:

```bash
export MPLBACKEND=Agg MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig XDG_CACHE_HOME=/home/jiehuang/software/fungi/ProGATE_v2/.cache
/home/jiehuang/anaconda3/bin/snakemake -s workflow/frozen_protocol_benchmark.smk --configfile configs/frozen_protocol.yaml --cores 1 --nolock --allowed-rules freeze_protocol run_trainable_model --forcerun outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/metrics.tsv outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/metrics.tsv
```

Result:

- success
- output dir: [run_1032](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032)
- log exists: [run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/run_frozen_protocol_model.log)

### GNN

Command:

```bash
export MPLBACKEND=Agg MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig XDG_CACHE_HOME=/home/jiehuang/software/fungi/ProGATE_v2/.cache
/home/jiehuang/anaconda3/bin/snakemake -s workflow/frozen_protocol_benchmark.smk --configfile configs/frozen_protocol.yaml --cores 1 --nolock --allowed-rules freeze_protocol run_trainable_model --forcerun outputs/frozen_protocol_benchmark_v2/scerevisiae/GCN/ORT_EXP_SUB/run_1029/metrics.tsv outputs/frozen_protocol_benchmark_v2/scerevisiae/GCN/ORT_EXP_SUB/run_1029/metrics.tsv
```

Result:

- success
- output: [metrics.tsv](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/scerevisiae/GCN/ORT_EXP_SUB/run_1029/metrics.tsv)
- log exists: [run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/scerevisiae/GCN/ORT_EXP_SUB/run_1029/run_frozen_protocol_model.log)

## N2V_MLP Revalidation

Relevant implementation:

- [src/graph/run_node2vec_embedding.py](/home/jiehuang/software/fungi/ProGATE_v2/src/graph/run_node2vec_embedding.py)
- [src/train/run_frozen_protocol_model.py](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_frozen_protocol_model.py)

Current config:

- `runtime.node2vec_backend=auto`

Observed runtime behavior for `dmelanogaster / N2V_MLP / seed 1032`:

- backend resolved to `svd_fallback`
- generated summary: [node2vec_summary.tsv](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/node2vec_summary.tsv)
- generated metrics: [metrics.tsv](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/metrics.tsv)
- generated feature schema: [feature_schema.tsv](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/feature_schema.tsv)

Verified values:

- `node2vec_summary.tsv` backend: `svd_fallback`
- `metrics.tsv` `feature_dim`: `64`
- `feature_schema.tsv`: non-empty

Observed log content:

- [run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/run_frozen_protocol_model.log) contains only a SciPy sparse efficiency warning

## Aggregation

Fresh aggregation command executed:

```bash
mkdir -p results/frozen_protocol/summary
export PYTHONPATH="${PYTHONPATH:-.}:."
/home/jiehuang/anaconda3/bin/python -m src.eval.aggregate_frozen_protocol_runs --output-root outputs/frozen_protocol_benchmark_v2 --summary-dir results/frozen_protocol/summary > results/frozen_protocol/summary/aggregate_frozen_protocol.log 2>&1
```

Result:

- success
- [per_run_metrics.tsv](/home/jiehuang/software/fungi/ProGATE_v2/results/frozen_protocol/summary/per_run_metrics.tsv)
- [aggregated_metrics.tsv](/home/jiehuang/software/fungi/ProGATE_v2/results/frozen_protocol/summary/aggregated_metrics.tsv)
- [final_summary.tsv](/home/jiehuang/software/fungi/ProGATE_v2/results/frozen_protocol/summary/final_summary.tsv)
- [aggregate_frozen_protocol.log](/home/jiehuang/software/fungi/ProGATE_v2/results/frozen_protocol/summary/aggregate_frozen_protocol.log)

Validation-time file counts:

- `per_run_metrics.tsv`: `79` lines
- `aggregated_metrics.tsv`: `35` lines
- `final_summary.tsv`: `35` lines

Note:

- the aggregation module works on the currently available contents of `outputs/frozen_protocol_benchmark_v2`
- the Snakemake aggregation rule itself still expects the full benchmark output set as declared inputs, which is normal for a full run and was not exercised here because this check intentionally reran representative targets rather than all 282 benchmark outputs

## Per-Run Logs

Verified these files exist:

- [outputs/frozen_protocol_benchmark_v2/human/DC/NETWORK/deterministic/run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/human/DC/NETWORK/deterministic/run_frozen_protocol_model.log)
- [outputs/frozen_protocol_benchmark_v2/celegans/MLP/ORT_EXP_SUB/run_1031/run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/celegans/MLP/ORT_EXP_SUB/run_1031/run_frozen_protocol_model.log)
- [outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/dmelanogaster/N2V_MLP/N2V/run_1032/run_frozen_protocol_model.log)
- [outputs/frozen_protocol_benchmark_v2/scerevisiae/GCN/ORT_EXP_SUB/run_1029/run_frozen_protocol_model.log](/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2/scerevisiae/GCN/ORT_EXP_SUB/run_1029/run_frozen_protocol_model.log)

Observed sizes:

- deterministic DC log: `0` bytes
- classical MLP log: `0` bytes
- N2V_MLP log: `242` bytes
- GCN log: `0` bytes

This means the log-path wiring is working. Silent-success runs simply do not emit stdout/stderr content.

## Blockers

No current functional blockers were reproduced in this check.

| blocker_id | severity | affected_protocols | affected_models | exact failing command | exact failing file/module | root cause | minimal fix | whether partial benchmark can proceed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| none | none | none | none | none | none | none | none | yes |

## Verdict

Verdict: `READY_FOR_FULL_RUN`

Rationale:

- main entry script is still the correct full-run entry point
- full DAG dry-run succeeds
- deterministic, classical, `N2V_MLP`, and GNN representative targets all reran successfully through the workflow
- `N2V_MLP` is fixed inside the main workflow and now produces correct `backend`, `feature_dim`, and non-empty schema outputs
- aggregation module runs successfully and rewrites the expected summary files
