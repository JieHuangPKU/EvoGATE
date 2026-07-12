# ProGATE_v2 Gap Analysis

## Scope

This report audits `/home/jiehuang/software/fungi/ProGATE_v2` against the legacy EPGAT runtime and data contracts.

Inspected areas:
- `configs/`
- `data_registry/`
- `outputs/`
- `scripts/`
- `src/`

Key inspected files:
- `src/data/build_baseline_dataset.py`
- `src/data/graph_dataset_adapter.py`
- `src/data/gene_graph_adapter.py`
- `src/features/load_embeddings.py`
- `src/features/plm_prepare.py`
- `src/registry/load_registry.py`
- `src/registry/support_graph_registry.py`
- `src/train/train_baseline.py`
- `src/train/train_graph_model.py`
- `src/train/train_support_graph_baseline.py`
- `src/eval/evaluate_baseline.py`
- `src/eval/evaluate_graph_model.py`

## B1. Current Modularity Of ProGATE_v2

### Configuration

Location:
- `configs/`

Observed role:
- YAML-driven execution for baseline, graph-ready, support-graph, and experiment variants.

Assessment:
- Much more explicit than EPGAT.
- Better suited as migration landing zone for legacy configs.

### Data Construction

Primary locations:
- `src/data/build_baseline_dataset.py`
- `src/data/gene_graph_adapter.py`
- `src/data/graph_dataset_adapter.py`
- other graph/input builders under `src/data/`

Observed role:
- build dataset tables from registry assets
- convert raw graph assets into canonical graph-ready tables
- summarize graph compatibility

Assessment:
- This is the correct landing layer for EPGAT data/feature adapters.

### Registry

Primary locations:
- `data_registry/`
- `src/registry/load_registry.py`
- `src/registry/support_graph_registry.py`

Observed role:
- explicit registry bundle loading
- canonical label tables
- evidence tables
- PPI bridge tables
- species-aware graph registry

Assessment:
- This is the biggest structural upgrade over EPGAT.
- It should absorb legacy EPGAT identifiers and evidence provenance rather than bypass them.

### Feature Loading

Primary locations:
- `src/features/load_embeddings.py`
- `src/features/plm_prepare.py`

Observed role:
- embedding manifest discovery
- embedding loading and pooling
- embedding availability tracking

Assessment:
- This is already aligned with migrating ESM/ESMC/ProtT5/ProtBERT into a first-class feature layer.

### Models

Primary locations:
- `src/models/baseline_models.py`
- `src/models/gnn_gat.py`
- `src/models/gnn_gcn.py`
- `src/models/gnn_graphsage.py`
- `src/models/support_graph_baseline.py`

Assessment:
- Models are already separated from runners and evaluation.
- This is the correct target for EPGAT model-class migration.

### Training / Evaluation

Primary locations:
- `src/train/train_baseline.py`
- `src/train/train_graph_model.py`
- `src/train/train_support_graph_baseline.py`
- `src/eval/evaluate_baseline.py`
- `src/eval/evaluate_graph_model.py`

Assessment:
- ProGATE_v2 already has a cleaner train/eval protocol than EPGAT.
- The missing piece is a legacy-compatible EPGAT adapter/training contract.

## B2. Existing Data Assets Reusable For Migration

### `data_registry/fgraminearum_canonical_id_resolved.tsv`

Observed fields:
- `raw_gene_id`
- `raw_protein_id`
- `raw_transcript_id`
- `final_canonical_gene_id`
- `mapping_confidence`
- `mapping_rule`
- `conflict_resolved`
- `needs_manual_review`

Migration value:
- critical for replacing EPGAT's implicit ID matching with explicit canonical mapping
- should become the bridge between old Fusarium features and new canonical runtime

### `data_registry/fgraminearum_gold_positive.broad79.tsv`

Observed role:
- canonical Fusarium gold positive benchmark
- includes evidence provenance and notes linking back to legacy Fusarium assets

Migration value:
- should replace legacy `ogee.csv` as the benchmark-facing positive set in ProGATE_v2
- useful for validating whether legacy labels are evidence-compatible or only legacy-training artifacts

### `data_registry/fgraminearum_gold_positive.strict29.tsv`

Observed role:
- stricter canonical Fusarium positive subset

Migration value:
- evaluation standard for ranking-compatible migration

### `data_registry/master_label_table.preliminary.tsv`

Observed role:
- canonical multi-species label table
- contains source file, source rule, confidence, notes

Migration value:
- should absorb the legacy EPGAT label universe into a traceable registry contract

### `data_registry/master_evidence_table.preliminary.tsv`

Observed role:
- evidence-level provenance table

Migration value:
- ideal place to preserve legacy `gene_list.txt`, Ogee-derived files, PHI-base linkage, and other upstream sources without overloading runtime feature tables

### `*_ppi_canonical_bridge.tsv`

Observed:
- `human_ppi_canonical_bridge.tsv`
- `scerevisiae_ppi_canonical_bridge.tsv`
- `celegans_ppi_canonical_bridge.tsv`

Migration value:
- these bridge tables are exactly the sort of explicit adapter that legacy EPGAT lacks
- same pattern should be extended to Fusarium legacy PPI and possibly other legacy feature assets

## B3. Structural Differences Between EPGAT And ProGATE_v2

### EPGAT

Dominant pattern:
- script-driven
- runner-driven
- experiment-flag-driven
- artifact-directory-driven

Key characteristics:
- data assumptions embedded in code
- feature schema assembled on the fly
- IDs assumed compatible after ad hoc preprocessing
- result summaries appended in-place

### ProGATE_v2

Dominant pattern:
- registry-driven
- module-driven
- dataset-builder-driven
- config-driven

Key characteristics:
- canonical IDs explicit
- evidence provenance explicit
- dataset and evaluation contracts explicit
- outputs separated by purpose

### Biggest Interface Incompatibilities

1. ID contract mismatch
- EPGAT uses per-file processed IDs with mostly implicit agreement
- ProGATE_v2 expects canonical gene IDs

2. Label contract mismatch
- EPGAT uses `EssentialGenes/ogee.csv` directly as training labels
- ProGATE_v2 separates raw evidence, label tables, and benchmark subsets

3. Feature contract mismatch
- EPGAT builds features dynamically from file presence and CLI flags
- ProGATE_v2 expects explicit dataset-build steps and feature manifests

4. Output contract mismatch
- EPGAT writes summary CSVs directly by model family in `results/`
- ProGATE_v2 writes structured training/evaluation outputs by config/run directory

5. Runtime orchestration mismatch
- EPGAT uses Python batch launchers and many runners
- ProGATE_v2 uses YAML configs plus modular builders/trainers/evaluators

## Where EPGAT Should Land Inside ProGATE_v2

Best landing zones:
- legacy dataset assembly: `src/data/`
- legacy feature schema normalization: `src/features/`
- legacy/original model implementations: `src/models/`
- legacy-compatible train wrappers: `src/train/`
- legacy result exporters and converters: `src/eval/`
- config snapshots: `configs/`
- migration/audit docs: `docs/epgat_migration/`

## Gap Summary

### Already Present In ProGATE_v2

- canonical registry layer
- evidence-aware label tables
- embedding manifest logic
- modular train/eval pipeline
- graph-ready adapter pattern
- ranking-aware Fusarium evaluation protocol

### Missing For A Seamless EPGAT Migration

- explicit EPGAT legacy dataset builder
- explicit EPGAT legacy feature schema definition
- explicit legacy label adapter from `ogee.csv` and companion tables
- explicit export layer that can emit legacy `results/*.csv`-style summaries
- explicit separation between original EPGAT core and extended EPGAT add-ons

## Direct Conclusions

1. ProGATE_v2 is already structurally better than EPGAT for long-term maintenance.
2. The migration should not transplant EPGAT's directory layout. It should translate EPGAT behavior into ProGATE_v2 contracts.
3. The main engineering task is adapter design, not model porting alone.
