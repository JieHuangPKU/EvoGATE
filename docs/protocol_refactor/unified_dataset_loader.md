# Unified Dataset Loader

## Authoritative Module

- `src/data/frozen_protocol_loader.py`

## Purpose

This module is the mandatory contract between:

- Layer 1: frozen protocol inputs
- Layer 2: model families

No mainline model is allowed to bypass this loader and read protocol inputs ad hoc.

## Required Inputs

`load_protocol_dataset(config_path, protocol_name, feature_setting)`

Inputs:

- `config_path`
- `protocol_name`
- `feature_setting`

Mainline protocol names:

- `human`
- `celegans`
- `scerevisiae`
- `dmelanogaster`
- `fgraminearum_newlabel`

Legacy replay protocol name:

- `fgraminearum_oldlabel`

## Frozen Inputs Consumed

Labels:

- `results/frozen_protocol/labels/*.tsv`

Splits:

- `results/frozen_protocol/splits/*_split.tsv`

Graphs:

- `data/processed/PPI/<data_key>/string.csv`

Features:

- `data/processed/OR/<data_key>/orthologs.csv`
- `data/processed/EXP/<data_key>/profile.csv`
- `data/processed/LC/<data_key>/subloc.csv`

## Expected Columns

Frozen label manifests must provide:

- `protocol_name`
- `species`
- `regime`
- `data_key`
- `canonical_gene_id`
- `graph_gene_id`
- `label`

Frozen split manifests must provide:

- all frozen label columns
- `split`
- `split_seed`
- `split_strategy`
- `split_version`

PPI graphs must provide:

- `A`
- `B`
- optional `combined_score`

Feature tables must provide:

- `Gene` as the gene identifier column

## Species and Regime Naming

Species names used in protocol manifests:

- `human`
- `celegans`
- `scerevisiae`
- `dmelanogaster`
- `fgraminearum`

Fusarium regimes are explicit:

- `newlabel`
- `oldlabel`

The loader does not infer the Fusarium regime from filename fragments or file presence.

## Feature Resolution Rules

Feature settings supported:

- `ORT`
- `EXP`
- `SUB`
- `ORT_EXP`
- `ORT_SUB`
- `EXP_SUB`
- `ORT_EXP_SUB`
- `N2V`
- `NETWORK`

Rules:

- `N2V` and `NETWORK` load no tabular feature blocks from OR/EXP/LC
- all other settings resolve directly from repo-local `data/processed/...`
- optional degree is computed from the loaded PPI graph
- normalization is fit on the frozen training subset only

## Graph Resolution Rules

Graph nodes are taken from the union of:

- all PPI nodes above the configured STRING threshold
- all labeled nodes in the frozen split manifest

This prevents dropping labeled nodes silently when a label is present but graph coverage is incomplete.

## Returned Bundle

The loader returns:

- protocol metadata
- label and split manifest paths
- graph source path
- node manifest
- split manifest
- edge table
- edge index
- feature matrix
- feature schema
- `train_idx`
- `val_idx`
- `test_idx`
- `y_all`
- `split_version`

## Mainline Consumers

The new mainline runner is:

- `src/train/run_frozen_protocol_model.py`

Mainline model families using this contract:

- Classical: `MLP`, `RF`, `SVM`, `NB`
- Network / embedding: `N2V_MLP`, `DC`, `CC`
- GNN: `GAT`, `GCN`, `GIN`, `GraphSAGE`

## Failure Modes

The loader fails fast when:

- the protocol name is unknown
- a frozen label or split manifest is missing
- required graph or feature files are missing
- split manifests are missing `split` or `label`
- feature normalization would operate on an empty train split
- unsupported feature settings are requested

## Contract Guarantee

If a model entry point still loads labels, splits, graphs, or features outside this loader, the protocol is not unified. The current mainline benchmark runner is wired through this loader and should remain the only supported path.
