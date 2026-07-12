# EPGAT To ProGATE_v2 Migration Blueprint

## Goal

Design a migration that simultaneously:
- reproduces legacy EPGAT behavior
- preserves paper-original EPGAT logic
- absorbs post-paper expansions such as ESM/ESMC/omics
- aligns with ProGATE_v2 canonical registry and evaluation contracts
- supports Fusarium as the long-term priority

## C1. Migration Principles

### Layer 1. EPGAT-Original-Compatible Layer

Purpose:
- preserve the closest possible behavior of original EPGAT
- support strict legacy reruns

Contents:
- original-style data assembly
- original GAT model defaults
- original feature flags: PPI, expression, orthologs, sublocalization
- legacy result export compatibility

Design rule:
- behavior-preserving adapter layer, not a speculative rewrite

### Layer 2. EPGAT-Extended Feature Layer

Purpose:
- preserve the repository's real post-paper capability set

Contents:
- ESM2
- ESM1b
- ESMC
- ProtT5
- ProtBERT
- GraphSAGE / GCN / GIN / classical baselines if retained
- omics combination schema

Design rule:
- separate feature preparation from model runtime
- turn dynamic CLI feature toggles into explicit config/schema objects

### Layer 3. ProGATE Unified Training/Eval Layer

Purpose:
- make legacy and extended EPGAT variants executable under a shared ProGATE_v2 protocol

Contents:
- YAML configs
- structured train/eval outputs
- canonical ID adapters
- evaluation exporters

Design rule:
- new work happens here
- old EPGAT behavior is exposed through wrappers, not by preserving old repository chaos

## C2. Migration Strategy

### Directly Migrate

Direct candidates:
- core GAT architecture logic
- active GCN/GIN/GraphSAGE model implementations if still needed
- PLM loading semantics
- species feature tables as upstream reference assets

### Rewrite As Adapters

Adapter candidates:
- legacy data assembly from `data/essential_genes/<species>/...`
- ID alignment logic
- label ingestion from `ogee.csv` and related evidence files
- legacy results CSV export
- batch experiment matrix generation

### Keep As Reference Only

Reference-only candidates:
- `runners/backup/`
- `runners/20250714/`
- `utils/prepare_data.py`
- dated model variants
- duplicated result cohorts such as `results/Third/`

## C3. Proposed New Directory Structure

Recommended additions inside ProGATE_v2:

- `src/data/build_epgat_legacy_dataset.py`
- `src/data/build_epgat_extended_dataset.py`
- `src/data/epgat_legacy_id_adapter.py`
- `src/data/epgat_legacy_label_adapter.py`
- `src/features/epgat_legacy_features.py`
- `src/features/epgat_plm_features.py`
- `src/features/epgat_feature_schema.py`
- `src/models/epgat_original_gat.py`
- `src/models/epgat_gcn.py`
- `src/models/epgat_gin.py`
- `src/models/epgat_graphsage.py`
- `src/models/epgat_baselines.py`
- `src/train/train_epgat_legacy.py`
- `src/train/train_epgat_extended.py`
- `src/eval/export_epgat_legacy_results.py`
- `src/eval/evaluate_epgat_extended.py`
- `configs/epgat_legacy.yaml`
- `configs/epgat_extended_embeddings.yaml`
- `configs/epgat_extended_omics.yaml`
- `configs/epgat_fusarium_legacy_replay.yaml`

## C4. Compatibility Strategy

### Reproduce Old Original EPGAT Results

Method:
1. freeze legacy input file list
2. freeze legacy feature ordering
3. freeze legacy train/val/test split behavior
4. freeze legacy model defaults from `models/gat/params.py`
5. export legacy-format metrics CSVs

Compatibility object:
- `epgat_legacy.yaml` should explicitly pin:
  - source data paths
  - feature flags
  - string threshold
  - split seed
  - model params
  - output export format

### Replace Old Runner Functions With New Entrypoints

Method:
- map each legacy runner to a ProGATE_v2 train config + exporter
- keep thin wrappers only if needed for command-line parity

Example:
- legacy `python -m runners.run_gat --organism fgraminearum --ppi string --expression --orthologs --sublocs`
- new `python -m src.train.train_epgat_legacy --config configs/epgat_legacy.yaml --override feature_scope.expression=true ...`

### Preserve Old Core CSV Formats

Method:
- implement export adapters producing:
  - `GAT_results.csv`
  - `GCN_results.csv`
  - `GIN_results.csv`
  - `SAGE_results.csv`
  - `MLP_results.csv`
  - `SVM_results.csv`
  - `RandomForest_results.csv`
  - `NaiveBayes_results.csv`
  - `N2V_MLP_results.csv`
  - `NetworkCC_results.csv`
  - `NetworkDC_results.csv`

Important boundary:
- these should be exported from structured ProGATE_v2 outputs, not used as the primary native storage format

### Support Old Species Data Paths Through An Adapter Layer

Method:
- do not relocate legacy EPGAT data immediately
- implement path adapters that ingest `EPGAT/data/essential_genes/<species>/...`
- normalize to ProGATE_v2 canonical sample/feature tables

## C5. ESM / ESMC Migration Strategy

### How ESM / ESMC Are Currently Integrated In EPGAT

Current pattern:
- generated externally into HDF5 files
- loaded via `runners/tools.py`
- matched directly to gene IDs
- appended to the node feature matrix

This means:
- they are preprocessing-time assets
- not train-time neural encoders

### Where They Should Land In ProGATE_v2

Recommended landing:
- preparation and indexing: `src/features/plm_prepare.py`
- feature lookup / manifests: `src/features/load_embeddings.py`
- EPGAT-specific feature assembly: `src/features/epgat_plm_features.py`

### Preprocessing Vs Runtime Encoder

Recommendation:
- keep them as preprocessing / feature-loading assets
- do not introduce runtime PLM encoding in the migration phase

Reason:
- matches current EPGAT behavior
- preserves reproducibility
- avoids unnecessary compute and dependency drift

### Fusion Location

Recommended fusion point:
- `src/features/epgat_feature_schema.py`

Fusion contract should explicitly define:
- canonical gene ID
- degree feature
- expression block
- orthology block
- sublocalization block
- PLM block
- missingness masks per block

## C6. Fusarium-Specific Migration Strategy

### Source Legacy Layout

Legacy location:
- `EPGAT/data/essential_genes/fgraminearum/`

Important legacy sources:
- `EssentialGenes/ogee.csv`
- `EssentialGenes/gene_list.txt`
- `PPI/STRING/string.csv`
- `PPI/mapping.tab`
- `Expression/profile.csv`
- `Expression/profile_fc.csv`
- `Orthologs/orthologs.csv`
- `SubLocalizations/subloc.csv`

### Target Contract In ProGATE_v2

Target should be canonical-ID-centered:
- registry tables remain source of truth for canonical IDs and benchmark subsets
- legacy feature assets are mapped into canonical sample space by adapters

### ID Alignment Plan

1. use `data_registry/fgraminearum_canonical_id_resolved.tsv` as the canonical bridge
2. build `epgat_legacy_id_adapter.py` to resolve:
   - raw gene IDs
   - raw protein IDs
   - old PLM embedding keys
   - old PPI node IDs
3. attach mapping confidence and unresolved flags
4. keep unresolved rows auditable instead of silently dropping them

### Avoiding Legacy ID Chaos

Rules:
- no new runtime path should match Fusarium features by naked string equality alone
- every feature table must pass through a canonicalization step
- each adapted feature table should emit an audit table:
  - raw_id
  - canonical_gene_id
  - mapping_rule
  - mapping_confidence
  - unresolved_flag

### Unified Fusarium Feature Schema

Recommended single feature schema:
- `canonical_gene_id`
- `legacy_gene_id`
- `has_ppi_node`
- `degree`
- expression features
- orthology features
- subcellular localization features
- PLM embedding vector
- missingness masks for each feature block

Optional but recommended:
- provenance columns for each block

## Recommended End State

### Preserve

- a strict EPGAT legacy replay path
- original GAT defaults
- legacy summary CSV export compatibility

### Modernize

- data loading
- ID mapping
- feature manifests
- training/evaluation orchestration

### Do Not Migrate As-Is

- old backup directories
- ad hoc post-hoc summary-merging scripts
- duplicate result-directory structures

## Direct Blueprint Conclusion

The migration should be adapter-first:
- preserve original behavior through explicit legacy dataset/model/export wrappers
- absorb extended features into a structured feature layer
- promote all new training/evaluation work into ProGATE_v2's config- and registry-driven runtime
