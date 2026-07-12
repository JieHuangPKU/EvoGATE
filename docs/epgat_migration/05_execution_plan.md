# Execution Plan

## Phase 1. Original-Compatible Migration

Objective:
- create a strict EPGAT legacy replay path inside ProGATE_v2

Deliverables:
- `src/data/build_epgat_legacy_dataset.py`
- `src/data/epgat_legacy_id_adapter.py`
- `src/data/epgat_legacy_label_adapter.py`
- `src/features/epgat_legacy_features.py`
- `src/models/epgat_original_gat.py`
- `src/train/train_epgat_legacy.py`
- `src/eval/export_epgat_legacy_results.py`
- `configs/epgat_legacy.yaml`

Tasks:
1. Freeze the exact legacy input file list per species.
2. Encode original feature ordering and split behavior.
3. Port original GAT model behavior with minimal change.
4. Export a `results/GAT_results.csv`-compatible table from ProGATE_v2 outputs.
5. Validate a small replay on one support species and one Fusarium legacy setting.

Exit criteria:
- one legacy GAT run can be executed entirely from ProGATE_v2
- result table schema matches legacy expectations

## Phase 2. Extended Feature Migration

Objective:
- migrate post-paper feature and model expansions without contaminating the original-compatible layer

Deliverables:
- `src/features/epgat_plm_features.py`
- `src/features/epgat_feature_schema.py`
- `src/models/epgat_gcn.py`
- `src/models/epgat_gin.py`
- `src/models/epgat_graphsage.py`
- `src/models/epgat_baselines.py`
- `src/train/train_epgat_extended.py`
- `configs/epgat_extended_embeddings.yaml`
- `configs/epgat_extended_omics.yaml`

Tasks:
1. Convert PLM assets into a manifest-driven feature block.
2. Define one explicit extended feature schema.
3. Re-express embedding and omics matrices as config grids.
4. Move model selection from runner scripts into config.
5. Export legacy-compatible summary CSVs only as a downstream view.

Exit criteria:
- ESM/ESMC/ProtT5/ProtBERT runs are config-driven
- omics and PLM features share the same feature-schema contract

## Phase 3. Fusarium Integration

Objective:
- make Fusarium the canonical long-term target without breaking legacy replay

Deliverables:
- canonicalized Fusarium legacy feature adapters
- Fusarium-specific config such as `configs/epgat_fusarium_legacy_replay.yaml`
- feature mapping audit tables

Tasks:
1. Map legacy Fusarium IDs through `fgraminearum_canonical_id_resolved.tsv`.
2. Convert legacy PPI/expression/orthology/subloc assets into canonical feature blocks.
3. Preserve unresolved or ambiguous mappings as audit rows.
4. Align benchmark-facing evaluation with `broad79`, `strict29`, and `conflict8`.
5. Separate training labels from evidence-only upstream files.

Exit criteria:
- all Fusarium legacy feature blocks can be joined on canonical IDs
- benchmark evaluation uses ProGATE_v2 subset tables, not ad hoc legacy label files

## Phase 4. Unified Training / Eval / Result Protocol

Objective:
- make EPGAT legacy and EPGAT extended first-class citizens in the ProGATE_v2 runtime

Deliverables:
- shared evaluation exporters
- shared output directory policy
- shared run manifest policy

Tasks:
1. Standardize training output directories.
2. Standardize prediction table columns.
3. Standardize evaluation table schema.
4. Provide optional exporters for legacy CSV compatibility.
5. Deprecate reliance on duplicated legacy artifact trees.

Exit criteria:
- no new experiments need to write directly into `results/*.csv`
- legacy CSVs are generated from structured outputs as compatibility artifacts

## Recommended Order

1. Phase 1
2. Phase 3
3. Phase 2
4. Phase 4

Reason:
- strict replay is needed before safe abstraction
- Fusarium ID and label correctness is higher-risk than model variety
