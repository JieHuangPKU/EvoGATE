# EvoGATE runbook

_Read-only inspection, lightweight validation, canonical entry points, blocked paths, and execution safety._

---

## Operating boundary

This repository is not a portable software release. Inspect before running. Do not run training, Snakemake workflows, ESM2 extraction, Figure regeneration, or large data processing without explicit user approval. Do not write to `data/` or existing `results/` unless explicitly authorized.

## Read-only inspection

```bash
find . -maxdepth 2 -type d
find configs docs scripts src workflow -maxdepth 2 -type f
du -sh data results data/* results/*
rg -n "ProGATE_v2|EPGAT|Bingo" configs docs scripts src workflow
head -n 5 data/processed/essential_gene/fgraminearum/newlabel/summary.tsv
head -n 5 results/Figure3a/data/Figure3a_final_summary.tsv
```

These commands do not validate scientific correctness; they only inspect repository state and small text artifacts.

## Lightweight checks

The following checks parse source or configuration without launching the pipeline:

```bash
python -m src.data.freeze_unified_protocol --help
python -m src.train.run_frozen_protocol_model --help
python -m src.features.extract_esm2_pooled --help
python -m src.eval.build_figure5_candidate_prioritization --help
```

Import checks may fail if the active environment lacks graph or machine-learning dependencies. Such failure indicates environment incompleteness and must not be repaired by installing packages without approval.

## Canonical entry points

The entries below are canonical code paths, not blanket authorization to execute them.

| Task | Entry | Expected output | Status |
|---|---|---|---|
| Materialize Fusarium labels | `workflow/fgraminearum_label_materialization.smk` | `data/processed/essential_gene/fgraminearum/{oldlabel,newlabel}/` | Partially reproducible |
| Freeze labels and splits | `python -m src.data.freeze_unified_protocol --config configs/frozen_protocol.yaml` | `results/frozen_protocol/labels/`, `results/frozen_protocol/splits/` | Implemented; writes frozen results |
| Run one model task | `python -m src.train.run_frozen_protocol_model ...` | Per-run output directory | Implemented; requires full environment |
| Full frozen benchmark | `workflow/frozen_protocol_benchmark.smk` | `outputs/Figure1/`, `results/Figure1/` | Large job; approval required |
| Prepare ESM2 cache | `workflow/prepare_esm2_cache.smk` | `data/processed/ESM2/<species>/esm2_pooled.pt` | Large job; approval required |
| Candidate prioritization | `python -m src.eval.build_figure5_candidate_prioritization` | `results/Figure5_new_candidate_prioritization/` | Currently Blocked |

## Dry-run guidance

A Snakemake dry-run is the preferred preflight for a workflow, but it may still evaluate configured paths and imports. It requires explicit approval under the current project rules because workflow execution was excluded from the documentation phase.

When approved, use the workflow file directly from the repository root rather than a historical wrapper. Confirm all planned outputs are in a new directory before proceeding.

## Historical and non-portable entry points

Shell wrappers including `scripts/run_Figure1_frozen_protocol_benchmark.sh`, `scripts/run_fgraminearum_label_materialization.sh`, `scripts/run_label_scarcity_benchmark.sh`, and several Figure2-Figure4 wrappers hard-code `/home/jiehuang/software/fungi/ProGATE_v2` or machine-specific environments.

Status: **Historical / non-portable**. Do not present these as recommended EvoGATE commands. Report the path first; do not batch-replace it.

Some Figure5 wrappers derive the repository root relative to their own path and are structurally more portable, but their upstream `outputs/` dependencies are absent.

## Currently blocked operations

| Operation | Blocker |
|---|---|
| Rebuild all published summaries | Repository-local `outputs/` is absent |
| Recreate exact software environment | No environment lock or package manifest |
| Audit all evaluation implementations | Some evaluation source files are missing while `.pyc` remains |
| Rebuild yeast-transfer confidence | Upstream confidence generator is missing |
| Use Git history for provenance | `.git/` is empty |
| Run historical wrappers from current root | Hard-coded ProGATE_v2/macOS paths |

## Large jobs requiring confirmation

- full frozen benchmark
- any Figure benchmark or ablation workflow
- label materialization that writes under `data/processed/`
- ESM2 extraction
- graph reconstruction or threshold sweeps
- candidate regeneration that overwrites existing results
- any job using GPUs, many CPU cores, or large matrices

Before approval, report the exact command, inputs, outputs, overwrite behavior, environment, estimated scale, and recovery plan.

## Preflight checklist

1. Confirm the current working directory is the EvoGATE root.
2. Read `AGENTS.md`, `docs/INCONSISTENCIES.md`, and the relevant config.
3. Verify every input path exists.
4. Verify the output directory is new and isolated.
5. Record protocol version, split version, split seed, and training seed.
6. Confirm no test-split optimization is performed.
7. Confirm the command does not delete or overwrite prior artifacts.
8. Obtain explicit approval for a large job.

## Common failure modes

| Symptom | Likely cause | Response |
|---|---|---|
| Shell wrapper changes to a missing directory | Historical `PROJECT_ROOT` | Stop and report; use no batch replacement |
| ESM2 cache not found | Missing local model or cache path | Report configured path; do not download automatically |
| Graph import failure | Missing DGL/PyG/torch-scatter environment | Record environment failure; do not install automatically |
| Candidate builder cannot find predictions | Missing `outputs/Figure3a/` | Treat as Blocked |
| Metrics differ between summaries | Conflicting historical artifacts | Record in `INCONSISTENCIES.md`; do not select a preferred number |
| Snakemake proposes deletion or overwrite | Existing workflow cleanup logic | Stop before execution and request confirmation |

## Run recording

For any approved future run, retain the command, timestamp, working directory, environment description, resolved config, input checksums, output directory, protocol version, split version, seed, exit status, and logs. Results without these fields should not become authoritative manuscript evidence.
