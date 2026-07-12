# Fgraminearum Newlabel Bridge Gap Audit

## Scope

This document records the last provenance gap that remained after the mainline benchmark stopped reading `results/label_rebuild_experiments/...` directly.

That gap was specific to the high-confidence yeast-transfer component of `fgraminearum_newlabel`.

## What The Gap Was

The missing artifact was an explicit bridge from:

- PH-1 `XP_*` protein accessions in `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`

to:

- the final canonical Fusarium gene space
- `fgraminearum::FGRAMPH1_*`

Because that bridge was not protocolized, the materialization workflow previously had to recover the high-confidence transfer-supported canonical positives from:

- `results/label_rebuild_experiments/labels/positive_set_P1.tsv`

In effect, `positive_set_P1.tsv` was hiding the missing mapping layer.

## Why `positive_set_P1.tsv` Was Historically Needed

The historical rebuild logic defined:

- `positive_set_P1.tsv = lethal_set | high_set`

The `high_set` part originated from the PH-1 yeast-transfer table in `XP_*` accession space.
However, the transformation from `XP_*` accessions to canonical `FGRAMPH1_*` genes had never been written down as a processed artifact.

As a result:

- the final canonical positive set existed
- but the bridge that produced it did not

That is why `positive_set_P1.tsv` was still needed upstream even after mainline benchmark execution no longer depended on it directly.

## What The Missing Bridge Actually Is

The bridge turns out to be the following chain:

1. `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
   This table defines PH-1 yeast-transfer candidates in `XP_*` protein accession space.

2. `data/derived_labels/proteome_manifest.tsv`
   This repo-local manifest records the PH-1 anchor proteome source used to derive the yeast-transfer table.

3. PH-1 NCBI proteome:
   `/data276/jiehuang/fungi/Fusarium/fusarium-protein/ncbi_dataset/data/GCF_000240135.3/protein.faa`

4. PH-1 legacy `229533` proteome:
   `/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/Orthologs/229533.fasta`

5. legacy accession-to-canonical mapping:
   `/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/PPI/mapping.tab`

6. supplemental repo-local FGSG support:
   - `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`
   - `data/processed/PPI/fgraminearum/string_id_mapping.tsv`
   - `data/processed/EXP/fgraminearum/exp_id_mapping.tsv`
   - `data/processed/LC/fgraminearum/subloc_id_mapping.tsv`

In other words, the missing bridge was:

- `XP_*` protein accession
- to PH-1 legacy protein accession
- to final canonical `FGRAMPH1_*`

with `FGSG_*` locus tags acting as an auxiliary support channel when they are present in the NCBI header.

## What Canonical ID Space Is Expected

The final expected output space is:

- `fgraminearum::FGRAMPH1_*`

This remains the authoritative ID space used by the frozen protocol manifests and mainline benchmark.

## What Was Protocolized

The bridge is now materialized under:

- `data/processed/essential_gene/fgraminearum/bridge/`

Key files:

- `protein_to_canonical_bridge.tsv`
- `source_to_canonical_mapping.tsv`
- `high_confidence_yeast_transfer_candidates.tsv`
- `unresolved_high_confidence_ids.tsv`
- `bridge_summary.tsv`
- `bridge_summary.md`
- `bridge_source_manifest.tsv`

## New Protocolized Bridge Logic

The bridge uses two explicit evidence channels:

### 1. Exact Sequence Bridge

The primary route is exact protein-sequence identity:

- `XP_*` protein in the PH-1 NCBI proteome
- exact same amino-acid sequence in the PH-1 legacy `229533` proteome
- `229533.<protein accession>` translated to `FGRAMPH1_*` via `mapping.tab`

This is the strongest mapping path and is marked as exact.

### 2. FGSG Header Support

When the NCBI protein header itself contains an `FGSG_*` locus tag, that `FGSG_*` identifier is mapped to canonical `FGRAMPH1_*` IDs using repo-local processed/mirrored mapping tables.

This route is retained as explicit supporting evidence and is marked as inferred when it is the sole resolver.

## What Plan Eliminated The Dependency

The implemented refactor is:

1. protocolize the missing XP-to-canonical bridge under `data/processed`
2. protocolize the rebuilt high-confidence yeast-transfer-supported candidate table
3. rebuild `fgraminearum_newlabel` from:
   - lethal PHI-supported positives
   - protocolized bridge-derived high-confidence transfer positives
   - repo-local evidence-based negative reconstruction
4. remove the direct dependency of the materialization workflow on:
   - `results/label_rebuild_experiments/labels/positive_set_P1.tsv`

## Current Limitation

The hidden bridge is now explicit and reproducible, but not every historical `XP_*` row resolves through currently protocolized local evidence.

Therefore:

- the provenance gap is closed
- the materialization workflow no longer directly depends on `positive_set_P1.tsv`
- but rebuilt newlabel counts can differ from the historical snapshot when unresolved `XP_*` rows are excluded by the explicit bridge

Those unresolved rows are now written to audit outputs instead of remaining silently hidden inside a legacy materialized positive set.
