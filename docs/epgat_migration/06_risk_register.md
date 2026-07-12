# Risk Register

## R1. Original And Extended Behavior Are Mixed In One Legacy Tree

- Risk description: EPGAT no longer cleanly separates paper-original behavior from later extensions.
- Impact scope: model migration, benchmark interpretation, reproducibility claims.
- Priority: high
- Mitigation: create explicit `epgat_legacy` and `epgat_extended` layers before any refactor of shared code.

## R2. Legacy Feature Schema Is Not Frozen

- Risk description: base omics, PLM embeddings, and runner-specific options form a dynamic feature contract.
- Impact scope: reproducibility, apples-to-apples comparisons, Fusarium migration.
- Priority: high
- Mitigation: define explicit feature-schema objects and persist schema manifests per run.

## R3. Legacy IDs Are Implicitly Aligned

- Risk description: EPGAT mostly assumes processed node IDs already match labels, features, and embeddings.
- Impact scope: Fusarium, PLM joins, PPI bridges, graph migration.
- Priority: high
- Mitigation: route all legacy feature ingestion through canonical ID adapters and produce mapping audit tables.

## R4. Results Are Split Across Multiple Output Systems

- Risk description: `results/`, `outputs/results/`, `outputs/evaluation/`, `results/time_logs/`, and `training_logs/` all carry overlapping provenance.
- Impact scope: reproducibility, migration validation, auditability.
- Priority: high
- Mitigation: make structured ProGATE_v2 outputs primary and generate legacy CSVs only as exporters.

## R5. Historical Residue Can Be Mistaken For Live Runtime

- Risk description: `runners/backup/`, `runners/20250714/`, and dated model variants sit inside the live code tree.
- Impact scope: wrong file migration, regression risk, maintenance burden.
- Priority: high
- Mitigation: classify these as reference-only and exclude them from the executable migration path.

## R6. Fusarium Uses Mixed Upstream Evidence Types

- Risk description: legacy Fusarium runtime labels and evidence files are not equivalent and should not be treated as interchangeable.
- Impact scope: benchmark construction, label correctness, scientific claims.
- Priority: high
- Mitigation: keep evidence tables and training labels separate; map `gene_list.txt` into evidence provenance, not direct runtime labels unless explicitly justified.

## R7. Batch Experiment Logic Is Encoded As Imperative Scripts

- Risk description: experiment matrices live in Python constants and filename conventions.
- Impact scope: reproducibility, rerun automation, future extensibility.
- Priority: medium
- Mitigation: encode experiment spaces in YAML configs and config grids.

## R8. `utils/prepare_data.py` Can Corrupt Interpretation

- Risk description: the post-hoc merger uses hard-coded local paths and manual metric swapping logic.
- Impact scope: result trustworthiness, comparison reports.
- Priority: high
- Mitigation: do not migrate it as runtime logic; treat it as legacy analysis residue only.

## R9. Visualization Side Effects Are Mixed Into Training Runners

- Risk description: t-SNE, ROC/PR export, fold summary export, and training happen inside the same runner path.
- Impact scope: runtime complexity, failure modes, reproducibility.
- Priority: medium
- Mitigation: move evaluation and visualization into dedicated ProGATE_v2 eval/export steps.

## R10. Species-Specific Exceptions Are Hidden In Code

- Risk description: examples include `coli` subloc suppression and Fusarium-specific file variants.
- Impact scope: portability, adapter correctness.
- Priority: medium
- Mitigation: move species exceptions into data contracts or config, not runner code branches.
