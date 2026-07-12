# Phase 1.5 Fusarium Canonical Replay

## Goal

Preserve the legacy-compatible EPGAT feature contract while projecting
Fusarium features, labels, and graph IDs into current ProGATE_v2
canonical gene space.

## Registry Files Used

- `data_registry/master_label_table.preliminary.tsv`
- `data_registry/master_evidence_table.preliminary.tsv`
- `data_registry/fgraminearum_gold_positive.broad79.tsv`
- `data_registry/fgraminearum_gold_positive.strict29.tsv`
- `data_registry/fgraminearum_gold_reconciliation.v3.tsv`

## Mapping Strategy

- legacy `FGRAMPH1_*` gene IDs are resolved through exact alias lookup
- canonical adapter uses `raw_gene_id`, `raw_protein_id`, and
  `raw_transcript_id` from the preliminary master registry
- replay dataset stores `canonical_gene_id`, but keeps a compatible
  `legacy_gene_id` column so the Phase 1 trainer can be reused

## Label Strategy

- minimal replay-compatible label source: legacy
  `EPGAT/data/essential_genes/fgraminearum/EssentialGenes/ogee.csv`
- this was chosen to keep the replay close to legacy EPGAT behavior
- broader benchmark-aware label upgrades remain outside Phase 1.5

## Mapping Audit Snapshot

- canonical mapping audit written to:
  `outputs/epgat_legacy/fgraminearum_canonical_replay/canonical_mapping_audit.tsv`
- mapping status counts:
  - exact: `21141`
  - unresolved: `1`

## Current Metrics

- AUROC: `0.5567`
- AUPRC: `0.0377`
- Accuracy: `0.7529`
- F1: `0.0617`
- MCC: `0.0124`

## Current Limits

- canonical replay still uses the legacy upstream label contract
- benchmark-facing Fusarium label redesign is deliberately deferred
- this is a minimal replay-compatible canonicalization layer, not a
  full Phase 2 Fusarium integration

## Practical Interpretation

- canonical replay is now runnable end-to-end
- compared with `fgraminearum_original_replay`, canonical replay slightly improves AUROC
- it therefore establishes a valid canonical gene-space bridge for later Fusarium-specific upgrades
