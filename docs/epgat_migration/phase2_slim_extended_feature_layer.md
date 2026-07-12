# Phase 2 Slim Extended Feature Layer

## 1. PLMs integrated in this round

- ESM2: HDF5 embeddings, 1280 dimensions.
- ProtT5: HDF5 embeddings, 1024 dimensions.
- No ESMC, ESM1b, or ProtBERT were included in this slim round.

## 2. Integration in ProGATE_v2

- Manifest: `src/features/plm_manifest.py` records PLM source path, format, dimension, key type, and sample keys.
- Loader: `src/features/plm_loaders.py` loads top-level HDF5 datasets into a key->vector lookup.
- Extended feature assembly: `src/features/epgat_extended_features.py` appends one PLM block and re-zscores the concatenated matrix.
- Dataset builder: `src/data/build_epgat_extended_dataset.py` starts from an existing Phase 1/1.5 replay dataset, appends `esm2` or `prott5`, and writes `plm_coverage.tsv` plus `plm_mapping_audit.tsv`.
- Train/eval: `src/train/train_epgat_extended.py` and `src/eval/evaluate_epgat_extended.py` reuse the existing legacy-compatible GAT path without changing the model family.

## 3. Coverage

- human: PLM coverage is `0.9496` for both ESM2 and ProtT5.
- celegans: PLM coverage is `0.9405` for both ESM2 and ProtT5.
- fgraminearum_canonical: PLM coverage is `0.9999` for both ESM2 and ProtT5 after canonical->raw `FGRAMPH1_*` projection.

## 4. Initial experiment results

- human baseline vs +ESM2 vs +ProtT5: baseline remained best; both PLMs underperformed baseline on AUROC/AUPRC, with ProtT5 less harmful than ESM2.
- celegans baseline vs +ESM2 vs +ProtT5: both PLMs strongly improved over baseline; ESM2 gave the largest AUROC/AUPRC gain, while ProtT5 yielded the highest MCC.
- fgraminearum_canonical baseline vs +ESM2 vs +ProtT5: ESM2 improved AUROC, AUPRC, and MCC over baseline; ProtT5 slightly improved AUPRC but reduced AUROC and MCC.
- Best PLM overall by mean AUROC gain across the three species: `esm2`.
- Best PLM for Fusarium by AUROC on `fgraminearum_canonical`: `baseline_plus_esm2`.

## 5. Current limits

- No ESMC integration in this round.
- No baseline zoo migration in this round.
- No multi-task upgrade in this round.
- No deeper fusion module; PLMs are appended as explicit flat feature blocks only.