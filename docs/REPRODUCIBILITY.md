# EvoGATE reproducibility

_Current reproducibility contract, known gaps, and requirements for a release-grade reconstruction._

---

## Current reproducibility level

EvoGATE is **Partially reproducible**. Frozen labels, splits, configurations, processed modalities, aggregated metrics, and many Figure artifacts are present. Full end-to-end reproduction is **Blocked** because the software environment, some upstream producers, some source modules, Git history, and the main per-run output tree are missing.

## Frozen protocol

| Item | Authoritative value | Evidence |
|---|---|---|
| Protocol version | `frozen_protocol_v1` | `configs/frozen_protocol.yaml` |
| Split strategy | Stratified fixed split | `src/data/freeze_unified_protocol.py` |
| Train fraction | 70% | Derived from test `0.20` and validation `0.10` |
| Validation fraction | 10% | Frozen config and split artifacts |
| Test fraction | 20% | Frozen config and split artifacts |
| Split seed | `20260409` | Frozen config and split version |
| Training seeds | `1029`-`1033` | Frozen config and result summaries |
| Main Fusarium regime | `fgraminearum_newlabel` | Frozen protocol definitions |
| Main newlabel counts | 1,097 positive; 10,868 negative | Newlabel summary |

All formal comparisons must use the same frozen split. A model may not generate a private split for a formal comparison.

## Randomness and model selection

`src/train/run_frozen_protocol_model.py` seeds Python, NumPy, PyTorch, and available CUDA devices. Neural model checkpoints are selected by validation AUPRC when defined, with AUROC or accuracy only as fallback behavior in the implementation.

The five training seeds quantify optimization variability on one frozen split; they do not replace independent data splits or biological replication. Existing summaries generally report mean and standard deviation. Not all primary comparisons include paired confidence intervals or hypothesis tests.

## Threshold policy

Standard trainable-model outputs use a fixed threshold of `0.5`. Network heuristics use a top-k rule tied to the labeled positive count. Threshold-tuned Figure3c analyses derive F1- or MCC-optimal thresholds from validation predictions and then apply them to test predictions.

The test split must never select a threshold, model, feature combination, epoch, or hyperparameter. Descriptive test-split threshold curves must be labeled as descriptive and cannot define a final operating point.

## Data and feature determinism

Frozen label and split TSVs are materialized. Numeric ORT/EXP/SUB/degree features and ESM2 vectors are normalized using training-node statistics. ESM2 alignment fails on missing node embeddings instead of silently dropping nodes.

The current repository contains ESM2 `.pt` artifacts and extraction logs but should not load these large binary files for routine documentation or inspection. Exact ESM2 reproduction depends on the configured local model checkpoint and software environment.

## Environment status

| Area | Current evidence | Status |
|---|---|---|
| Python interpreter | Absolute Python paths in config/workflows | Historical / non-portable |
| Conda environment | `EPGAT` and `progate` names in configs/scripts | Historical |
| Python dependencies | Imports in source | Partially documented |
| R dependencies | R scripts present, no lock | Blocked |
| Snakemake version | Not locked | Blocked |
| PyTorch/graph backend | Multiple legacy environments and implementations | Blocked |
| ESM2 model | Local `esm2_t33_650M_UR50D` path evidenced | Partially documented |
| Hardware | ESM2 config requests CUDA; frozen benchmark config uses CPU | Partially documented |

No environment lock, package manifest, container, or reproducible installation procedure exists. Package installation must not be inferred from imports.

## Missing artifacts

- repository-local `outputs/` with predictions, checkpoints, and resolved configs
- source for several bytecode-only evaluation modules
- upstream producer for `weak_positive_confidence`
- non-empty Git metadata
- fully portable raw-data build inputs and commands
- environment lock files

These gaps prevent a full claim of reproducibility.

## Existing reproducibility assets

- frozen configuration and model hyperparameters
- materialized labels and splits
- source manifests and selected checksums for label construction
- processed modality summaries and mapping audits
- per-Figure aggregated metrics and reports
- ESM2 extraction logs and embedding metadata
- runtime config snapshots for several Figure3 experiments

## Release-grade requirements

1. Recover or establish version control with a documented baseline.
2. Lock Python, R, Snakemake, PyTorch, graph, and plotting dependencies.
3. Replace machine-specific runtime roots with reviewed portable configuration.
4. Recover all source modules used for authoritative results.
5. Recover or regenerate per-run outputs without overwriting historical artifacts.
6. Restore the upstream evolutionary-confidence producer and exact rules.
7. Create an explicit exclusion manifest.
8. Record input and output checksums for a named release.
9. Run contract-matched, paired validation for manuscript claims.
10. Rebuild all manuscript tables and Figures from the named release.

## Reproduction claim policy

Until all release requirements are met, documentation must say **Partially reproducible** or **Blocked for full reproduction**. It must not say fully reproducible, production-ready, or independently reproduced.

