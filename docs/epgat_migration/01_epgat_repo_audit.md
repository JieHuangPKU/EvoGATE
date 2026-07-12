# EPGAT Repository Audit

## Scope And Method

This audit is based on direct inspection of `/home/jiehuang/software/fungi/EPGAT`.

Inspected top-level targets:
- `README.md`
- `config/`
- `data/`
- `logs/`
- `models/`
- `outputs/`
- `results/`
- `run_all_experiments.py`
- `run_all_experiments_omics.py`
- `runners/`
- `scripts/`
- `training_logs/`
- `utils/`

Inspected code paths in detail:
- `runners/run_gat.py`
- `runners/run_gat_string.py`
- `runners/run_gcn.py`
- `runners/run_gin.py`
- `runners/run_sage.py`
- `runners/run_mlp.py`
- `runners/run_rf.py`
- `runners/run_svm.py`
- `runners/run_n2v_mlp.py`
- `runners/run_bayes.py`
- `runners/run_cc.py`
- `runners/run_dc.py`
- `runners/run_bc.py`
- `runners/run_gb.py`
- `runners/tools.py`
- `utils/utils.py`
- `utils/prepare_data.py`
- representative Fusarium files under `data/essential_genes/fgraminearum/`

Sampling note:
- `data/essential_genes/*/Orthologs/processed/` contains many FASTA assets; those were sampled by path pattern and not line-read exhaustively because they are bulk data, not active runtime logic.
- `outputs/evaluation/` and `logs/` contain many repeated artifacts; audit focused on naming contracts and representative files rather than full-file reading of every artifact.

## A1. Top-Level Directory Summary

### `README.md`

Represents the paper-original positioning of the repository. It describes EPGAT as code for the published GAT-based essential gene prediction paper and points users to `python runners/run_gat.py ...`.

Assessment:
- This is the closest thing to the original-paper contract.
- It does not document later additions such as ESM, ESMC, ProtBERT, ProtT5, GraphSAGE, GIN, RF, SVM, Naive Bayes, node2vec+MLP, centrality baselines, batch experiment drivers, or Fusarium-specific data preparation.

### `config/`

Contains `default.yaml` plus config package bootstrap.

Observed role:
- `utils/utils.py` imports `config.get_config_manager()` and uses it to resolve `DATA_ROOT` and STRING threshold defaults.

Assessment:
- This is not a full experiment configuration system.
- Core runtime settings still live inside runner code and batch launcher constants.

### `data/`

Contains the active data lake for legacy EPGAT experiments, centered under `data/essential_genes/<species>/`.

Observed schema for each species:
- `EssentialGenes/ogee.csv`
- `PPI/STRING/string.csv` and sometimes `PPI/BIOGRID/biogrid.csv`, `PPI/DIP/dip.csv`
- `Expression/profile.csv` or species variants
- `Orthologs/orthologs.csv`
- `SubLocalizations/subloc.csv`

Fusarium-specific observations:
- `data/essential_genes/fgraminearum/EssentialGenes/ogee.csv`
- `data/essential_genes/fgraminearum/EssentialGenes/gene_list.txt`
- `data/essential_genes/fgraminearum/PPI/229533.protein.links.detailed.v12.0.txt`
- `data/essential_genes/fgraminearum/PPI/mapping.tab`
- `data/essential_genes/fgraminearum/PPI/process.py`
- `data/essential_genes/fgraminearum/PPI/STRING/string.csv`
- `data/essential_genes/fgraminearum/Expression/profile.csv`
- `data/essential_genes/fgraminearum/Expression/profile_fc.csv`
- `data/essential_genes/fgraminearum/Orthologs/orthologs.csv`
- `data/essential_genes/fgraminearum/SubLocalizations/subloc.csv`

Assessment:
- This directory is both source-of-truth and staging area.
- Raw assets, processed assets, and runtime-ready assets are mixed together.

### `logs/`

Contains per-run log files with naming pattern like:
- `fgraminearum_run_gcn_esm1b.log`
- `human_run_gat_exp_ort_sub.log`
- `coli_run_n2v_mlp_esmc.log`

Assessment:
- This is a batch-run operational log directory.
- Naming reflects expanded feature space beyond original EPGAT.
- It is useful for audit provenance but not a stable machine-readable experiment registry.

### `models/`

Contains model implementations:
- `models/gat/`
- `models/gcn/`
- `models/gin/`
- `models/graphsage/`
- `models/node2vec/`

Assessment:
- The original-paper core is `models/gat/`.
- Everything else is a post-paper expansion or baseline extension.
- There is no separate package boundary for “original” vs “extended”.

### `outputs/`

Contains runtime artifacts:
- `outputs/evaluation/`
- `outputs/preds/`
- `outputs/results/`
- `outputs/weights/`

Observed artifacts:
- fold-level ROC/PR TSVs
- t-SNE PDFs and TSVs
- per-run prediction CSVs
- model checkpoints
- model-specific results CSVs such as `outputs/results/SAGE_<organism>_<ppi>_<embedding>_results.csv`

Assessment:
- This is a second result system in parallel to `results/`.
- It contains richer per-run artifacts than `results/` but without a single unified manifest.

### `results/`

Contains aggregated summary tables:
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
- `GAT_results-string.csv`
- `Third/`
- `time_logs/`

Assessment:
- This is the main legacy summary layer used to compare experiments across organisms and feature sets.
- It coexists with `outputs/results/`, causing duplication.
- `Third/` strongly suggests repeated manual rerun cohorts instead of a versioned experiment protocol.

### `run_all_experiments.py`

Batch scheduler for embedding-centric experiments.

Observed defaults:
- organisms: `celegans`, `coli`, `fgraminearum`, `yeast`, `melanogaster`, `human`
- methods currently narrowed to `run_svm`, `run_sage`
- embeddings: `esm2_embeddings`, `esm1b_embeddings`, `esmc_embeddings`, `prott5_embeddings`, `protbert_embeddings`
- writes time summary to `results/time_logs/experiment_summary_*.csv`
- writes logs to `logs/`

Assessment:
- Not paper-original.
- This is a post-expansion orchestration layer.

### `run_all_experiments_omics.py`

Batch scheduler for omics feature combinations.

Observed features:
- `sublocs`
- `expression`
- `orthologs`

Assessment:
- Also an expanded orchestration layer.
- Handles combinatorial feature toggles and species-specific exceptions like removing `sublocs` for `coli`.

### `runners/`

Contains the executable experiment entry points and helpers.

Observed categories:
- primary runners: `run_gat.py`, `run_gcn.py`, `run_gin.py`, `run_sage.py`, `run_mlp.py`, `run_rf.py`, `run_svm.py`, `run_n2v_mlp.py`, `run_bayes.py`
- network heuristics: `run_cc.py`, `run_dc.py`, `run_bc.py`
- extra baseline: `run_gb.py`
- helpers: `tools.py`, `utils.py`, `feature_selection.py`
- historical residue: `backup/`, `20250714/`

Assessment:
- This is the operational heart of the repository.
- It mixes live code, old code, and dated snapshots.

### `scripts/`

Observed important files:
- `scripts/shared_visualizer.py`
- `scripts/generate_embeddings.py`
- `scripts/data.py`
- subdirs `embeddings/`, `essential_genes/`, `sequence/`

Assessment:
- This houses post-paper support tooling, especially PLM embedding generation and evaluation visualization.
- It is not just auxiliary; parts are runtime dependencies for expanded experiments.

### `training_logs/`

Contains:
- training log text files
- `experiment_summary_20250719_103817.csv`
- prior variants like `parallel_training.log`, `run_all_experiments-v1.py`

Assessment:
- Another historical artifact area.
- Some files duplicate `results/time_logs/`.

### `utils/`

Observed important files:
- `utils/utils.py`
- `utils/prepare_data.py`
- `utils/training.py`
- stats collectors and caches

Assessment:
- `utils/utils.py` still defines the original-style data contract.
- `utils/prepare_data.py` is not a runtime data builder; it is a legacy post-hoc summary merger, and it includes hard-coded macOS absolute paths.

## A2. Entry Scripts And Experiment Scheduling

### True Main Entrypoints

There is no single global main. The runtime breaks into three layers:

1. Single-model runners:
- `runners/run_gat.py`
- `runners/run_gcn.py`
- `runners/run_gin.py`
- `runners/run_sage.py`
- `runners/run_mlp.py`
- `runners/run_rf.py`
- `runners/run_svm.py`
- `runners/run_n2v_mlp.py`
- `runners/run_bayes.py`
- `runners/run_cc.py`
- `runners/run_dc.py`
- `runners/run_bc.py`
- `runners/run_gb.py`

2. Batch experiment schedulers:
- `run_all_experiments.py`
- `run_all_experiments_omics.py`

3. Shared data/CLI helper:
- `runners/tools.py`

### Which Is The Paper-Original Entry

Most likely:
- `runners/run_gat.py`
- with `README.md` explicitly advertising it
- backed by `models/gat/gat_pytorch.py`
- backed by `models/gat/params.py`
- backed by `utils/utils.py` and `runners/tools.py`

### Which Are Single-Model Entrypoints

Single-model execution entrypoints:
- `run_gat.py`
- `run_gat_string.py`
- `run_gcn.py`
- `run_gin.py`
- `run_sage.py`
- `run_mlp.py`
- `run_rf.py`
- `run_svm.py`
- `run_n2v_mlp.py`
- `run_bayes.py`
- `run_cc.py`
- `run_dc.py`
- `run_bc.py`
- `run_gb.py`

### Which Are Batch Entrypoints

Batch drivers:
- `run_all_experiments.py`: embedding combinations
- `run_all_experiments_omics.py`: omics combinations

### Which Are Legacy/Backup/Temporary Residues

Clear residues:
- `runners/backup/*`
- `runners/20250714/*`
- `training_logs/run_all_experiments-v1.py`
- `training_logs/1.parallel_training-v3.py`
- `training_logs/1.parallel_training_string.py`
- `utils/neighborhood_weights-bk.py`
- `utils/utils-20250710.py`

These should not be migrated as first-class runtime modules.

### Result Write Locations

Observed write surfaces:
- `results/*.csv`: aggregated benchmark tables
- `results/time_logs/experiment_summary_*.csv`: batch scheduler timing summaries
- `outputs/preds/*.csv`: per-run predictions
- `outputs/results/*.csv`: model/run-specific result tables
- `outputs/evaluation/*.tsv` and PDFs: fold metrics and t-SNE artifacts
- `outputs/weights/<model>/...`: checkpoints
- `logs/*.log`: execution logs

### Metrics Aggregation Logic

Aggregation is decentralized:
- each runner appends or updates its own summary CSV in `results/`
- batch launchers aggregate timing only
- `scripts/shared_visualizer.py` aggregates fold ROC/PR summaries for visualization outputs
- `utils/prepare_data.py` performs a separate post-hoc merge across “First / Second / Third” result cohorts using hard-coded absolute paths and even metric-swapping logic for GAT/GCN

This means there is no single authoritative metrics protocol.

### Result Overwrite / Naming / Artifact Drift Problems

Observed issues:
- `results/` and `outputs/results/` both store metrics
- `results/GAT_results.csv` and `results/GAT_results-string.csv` reflect overlapping but different protocols
- `results/Third/` duplicates summary files by cohort rather than by explicit config identity
- `results/time_logs/` and `training_logs/experiment_summary_*.csv` duplicate timing-summary responsibilities
- `outputs/evaluation/` file names use a mix of `Base`, embedding names, and feature abbreviations
- `run_all_experiments.py` currently restricts `METHODS` to `run_svm` and `run_sage`, which means the file no longer truthfully represents the full intended batch space

## A3. Model Layer Audit

### Paper-Original EPGAT Core

Core original components:
- `models/gat/gat_pytorch.py`
- `models/gat/params.py`
- `runners/run_gat.py`
- `utils/utils.py`
- `runners/tools.py`

Original feature contract:
- PPI graph
- optional expression
- optional orthologs
- optional subcellular localization

### GAT Variants

Observed:
- `run_gat.py`: current expanded main GAT runner
- `run_gat_string.py`: older or alternate GAT path
- `models/gat/gat_pytorch0617.py`
- `models/gat/gat_pytorch20250710.py`

Assessment:
- `gat_pytorch.py` is the live implementation.
- dated variant files are historical residues and should be treated as reference only.

### GCN / GIN / GraphSAGE

Observed modules:
- `models/gcn/gcn_pytorch.py`
- `models/gin/gin_pytorch.py`
- `models/graphsage/graphsage_adapter.py`

Assessment:
- These are post-paper extensions and baselines.
- `models/gcn/gcn.py` is a DGL version and looks like older/reference code; `gcn_pytorch.py` is the active path used by `run_gcn.py`.
- `graphsage_adapter.py` exists to adapt edge lists to DGL GraphSAGE, which is clearly not paper-original.

### MLP / RF / SVM / Node2Vec+MLP / Naive Bayes / Gradient Boost

Observed runtime families:
- `run_mlp.py`
- `run_rf.py`
- `run_svm.py`
- `run_n2v_mlp.py`
- `run_bayes.py`
- `run_gb.py`

Assessment:
- All are expansion-era baselines.
- They consume the same or similar assembled feature matrix from `runners/tools.py`.

### Centrality Baselines

Observed:
- `run_cc.py`
- `run_dc.py`
- `run_bc.py`

Assessment:
- These are heuristic/network baselines, not original EPGAT core.
- `run_dc.py` and `run_cc.py` include explicit embedding-type-dependent feature weighting logic, which is strongly post-paper and hacky.

### ESM / ESMC / PLM-Related Modules

Observed integration points:
- `scripts/generate_embeddings.py`
- `scripts/data.py`
- runner flags: `--esm2_embeddings`, `--esm1b_embeddings`, `--esmc_embeddings`, `--prott5_embeddings`, `--protbert_embeddings`
- embedding loader in `runners/tools.py` using HDF5

Assessment:
- ESM/ESMC are not separate model encoders inside `models/`.
- They are precomputed feature assets loaded at runtime and appended into node features.

### Multimodal Fusion Logic

Observed fusion style:
- concatenate base omics features from `utils/utils.data()`
- append PLM embeddings inside runner/data assembly path
- allow `--only_embeddings` in some runners

Assessment:
- Fusion is feature-level early concatenation.
- There is no explicit multimodal fusion module or schema object.

## A4. Feature Engineering And Data Flow Audit

### Core Data Loader Contract

The legacy base contract is defined by `utils/utils.py:data()` and wrapped by `runners/tools.py:get_data()`.

Base flow:
1. load PPI edges from `data/essential_genes/<organism>/PPI/...`
2. if STRING, threshold by `combined_score`
3. load labels from `EssentialGenes/ogee.csv`
4. intersect labels with PPI gene universe
5. construct gene union
6. start empty feature matrix indexed by gene
7. optionally join:
   - `Orthologs/orthologs.csv`
   - `Expression/profile.csv`
   - `SubLocalizations/subloc.csv`
8. fill missing values with zero
9. split labels into train/test in `utils.utils.data()`
10. in `runners/tools.py`, convert gene IDs to integer node indices, remove self-loops/duplicates, add degree feature, and append embeddings depending on runner flags

### Where Essential Gene Labels Come From

Legacy runtime source:
- `data/essential_genes/<organism>/EssentialGenes/ogee.csv`

Fusarium example:
- `data/essential_genes/fgraminearum/EssentialGenes/ogee.csv`

Also present but not part of the main loader contract:
- `data/essential_genes/fgraminearum/EssentialGenes/gene_list.txt`

Assessment:
- `ogee.csv` is the active runtime label file.
- `gene_list.txt` is an upstream evidence artifact, not the direct runtime label input in current EPGAT code.

### Where PPI Comes From

Active runtime paths:
- `PPI/STRING/string.csv`
- optionally `PPI/BIOGRID/biogrid.csv`
- optionally `PPI/DIP/dip.csv`

For Fusarium:
- raw source `229533.protein.links.detailed.v12.0.txt`
- bridge file `mapping.tab`
- processor `PPI/process.py`
- runtime-ready file `PPI/STRING/string.csv`

### Where Expression Comes From

Active runtime path:
- `Expression/profile.csv`

Fusarium also contains:
- `Expression/profile_fc.csv`

Assessment:
- The base loader expects `profile.csv`.
- Presence of `profile_fc.csv` indicates alternate or newer expression preprocessing not formalized in the loader contract.

### Where Orthologs Come From

Active runtime path:
- `Orthologs/orthologs.csv`

Assessment:
- This is treated as a flat feature table, not a graph relation in legacy EPGAT.

### Where Subcellular Localization Comes From

Active runtime path:
- `SubLocalizations/subloc.csv`

### Where ESM / ESMC Features Come From

Observed contract:
- precomputed `.h5` files under `scripts/embeddings/`
- file names inferred by runner code, e.g. `<organism>_esm2_embeddings.h5`, `<organism>_esmc_embeddings.h5`
- loader function `runners/tools.py:load_protein_embeddings()`

Assessment:
- These are file-backed feature matrices, not registry-managed assets.

### Feature Concatenation Contract

Original-paper contract:
- node degree
- optional orthologs
- optional expression
- optional sublocalizations
- graph topology from PPI

Expanded contract:
- all of the above
- plus one PLM embedding family among ESM2 / ESM1b / ESMC / ProtT5 / ProtBERT
- plus runner-specific `only_embeddings` modes in extended paths
- plus centrality baselines using embedding-derived weighting heuristics

### Different Species Data Path Organization

Uniform legacy pattern:
- `data/essential_genes/<species>/...`

But with species-specific deviations:
- `coli` skips sublocalization
- Fusarium has extra raw processing files and variant expression file names

### Fusarium Current Placement

Current legacy location:
- `data/essential_genes/fgraminearum/`

Observed runtime-relevant files:
- `EssentialGenes/ogee.csv`
- `PPI/STRING/string.csv`
- `Expression/profile.csv`
- `Orthologs/orthologs.csv`
- `SubLocalizations/subloc.csv`

Observed additional upstream or alternative files:
- `EssentialGenes/gene_list.txt`
- `Expression/profile_fc.csv`
- raw STRING and mapping assets
- many ortholog FASTA intermediates

### Is The Feature Schema Frozen

No.

Evidence:
- `utils/utils.data()` defines one base schema
- runner-specific embedding flags extend it dynamically
- `run_sage.py` and peers add model-specific output logic and `only_embeddings`
- `profile.csv` vs `profile_fc.csv` shows uncontrolled feature variant drift
- centrality heuristics consume embeddings differently than classifiers

### Are Different Runners Reading Different Feature Schemas

Yes.

Concrete examples:
- GNN and classical ML runners all share `tools.get_data()` but may append different flags and embed handling
- some runners support `only_embeddings`
- centrality runners alter behavior based on active embedding type
- `run_gat_string.py` is older and does not align cleanly with the current expanded argument/result contract

## A5. Result Artifact Audit

### Core Products

Most important legacy aggregate products:
- `results/GAT_results.csv`
- `results/GCN_results.csv`
- `results/GIN_results.csv`
- `results/SAGE_results.csv`
- `results/MLP_results.csv`
- `results/SVM_results.csv`
- `results/RandomForest_results.csv`
- `results/NaiveBayes_results.csv`
- `results/N2V_MLP_results.csv`
- `results/NetworkCC_results.csv`
- `results/NetworkDC_results.csv`

These are the closest thing to a legacy benchmark table.

### Intermediate / Run-Level Products

Observed intermediates:
- `outputs/preds/*.csv`
- `outputs/evaluation/*`
- `outputs/weights/*`
- `logs/*.log`
- `results/time_logs/experiment_summary_*.csv`

### Duplicate Or Repeated Result Cohorts

Observed duplicates:
- `results/Third/`
- `results/time_logs/`
- `training_logs/experiment_summary_*.csv`
- `outputs/results/`

Assessment:
- Repeated “cohort” directories such as `Third/` indicate manual rerun layers rather than controlled immutable run IDs.

### What `GAT_results.csv` Represents

It stores aggregated metrics across:
- organism
- PPI
- boolean feature flags
- embedding flags
- `n_runs`

Columns observed:
- `mean`, `std`
- `auc_pr`, `auc_pr_std`
- `precision`, `precision_std`
- `mcc`, `mcc_std`
- `accuracy`, `accuracy_std`

### What `GAT_results-string.csv` Represents

It is a separate STRING-threshold-specific summary table containing:
- feature flags
- metrics
- explicit `string_thr`

Assessment:
- This is another protocol variant, not merely another slice of `GAT_results.csv`.

### Which Files Are Needed For Strict Legacy Reproduction

Minimum critical set:
- runtime code:
  - `runners/*.py`
  - `runners/tools.py`
  - `utils/utils.py`
  - active `models/*`
- data:
  - `data/essential_genes/<species>/EssentialGenes/ogee.csv`
  - `data/essential_genes/<species>/PPI/...`
  - `data/essential_genes/<species>/Expression/profile.csv`
  - `data/essential_genes/<species>/Orthologs/orthologs.csv`
  - `data/essential_genes/<species>/SubLocalizations/subloc.csv`
  - PLM `.h5` files from `scripts/embeddings/`
- result expectations:
  - `results/*.csv` for legacy summary-format compatibility

## A6. Technical Debt And Migration Risk

### Risk 1. Original And Extended Logic Are Not Separated

Impact:
- hard to prove which results are “paper-original”
- hard to preserve exact legacy behavior during migration

### Risk 2. Configuration Is Fragmented

Evidence:
- `config/default.yaml`
- hard-coded constants in batch launchers
- hard-coded file naming in runners
- hard-coded metric and path logic in `utils/prepare_data.py`

Impact:
- migration can easily change behavior accidentally

### Risk 3. Feature Schema Is Dynamic And Unfrozen

Evidence:
- base omics in `utils/utils.py`
- dynamic embedding append in runners
- `profile.csv` vs `profile_fc.csv`
- `only_embeddings` path

Impact:
- old runs are not reproducible unless schema snapshots are explicit

### Risk 4. Results Are Split Across Multiple Output Systems

Evidence:
- `results/`
- `outputs/results/`
- `outputs/evaluation/`
- `results/time_logs/`
- `training_logs/`

Impact:
- difficult to reconstruct authoritative result provenance

### Risk 5. Historical Residue Is Mixed Into Live Runtime Tree

Evidence:
- `runners/backup/`
- `runners/20250714/`
- dated model files

Impact:
- migration may copy stale code or wrong variants

### Risk 6. Species-Specific Hacks Are Embedded In Orchestration

Evidence:
- `coli` subloc suppression in `run_all_experiments_omics.py`
- Fusarium-specific raw processing files in active data tree

Impact:
- generalized adapters in ProGATE_v2 can break hidden assumptions

### Risk 7. ID Alignment Is Implicit In Legacy Runtime

Evidence:
- legacy EPGAT mostly assumes processed files already share node IDs
- Fusarium `PPI/process.py` depends on `mapping.tab`
- PLM loader matches directly on gene IDs

Impact:
- migration to canonical IDs needs an explicit adapter layer or legacy equivalence will be lost

### Risk 8. `utils/prepare_data.py` Is Not Trustworthy As A Reproducibility Asset

Evidence:
- hard-coded local macOS absolute paths
- post-hoc merging
- manual metric swapping between GAT and GCN rows

Impact:
- this file should not be treated as a migration target for runtime

## Original vs Extended Classification

### Best-Effort “Original EPGAT” Core

- `README.md`
- `models/gat/gat_pytorch.py`
- `models/gat/params.py`
- `runners/run_gat.py`
- `utils/utils.py`
- `runners/tools.py`
- base data contract:
  - PPI
  - expression
  - orthologs
  - sublocalizations

### Clear Post-Paper Extensions

- all PLM embedding support:
  - ESM2
  - ESM1b
  - ESMC
  - ProtT5
  - ProtBERT
- `scripts/generate_embeddings.py`
- `scripts/shared_visualizer.py`
- `run_all_experiments.py`
- `run_all_experiments_omics.py`
- `run_gcn.py`
- `run_gin.py`
- `run_sage.py`
- `run_mlp.py`
- `run_rf.py`
- `run_svm.py`
- `run_n2v_mlp.py`
- `run_bayes.py`
- `run_cc.py`
- `run_dc.py`
- `run_bc.py`
- `run_gb.py`
- Fusarium-specific data processing scripts and files
- backup/date-stamped runtime copies

## Direct Audit Conclusions

1. EPGAT is no longer a paper-only repository. It is a hybrid of original EPGAT, classical baseline zoo, PLM feature experiments, visualization tooling, Fusarium-specific curation, and batch automation.
2. The most important migration task is not copying models. It is freezing and externalizing the legacy feature/data/result contract.
3. The strongest migration boundary is:
   - preserve original EPGAT-compatible runtime behavior as a legacy layer
   - migrate data/ID/feature assembly into explicit ProGATE_v2 adapters
   - do not import legacy result-directory chaos into the new repo as-is
