# Phase 1.5 Multispecies Replay

## Scope

Phase 1.5 extends the original-compatible EPGAT layer from the initial
`scerevisiae` smoke test to:

- `celegans_original_replay`
- `human_original_replay`
- `fgraminearum_original_replay`
- `fgraminearum_canonical_replay`

## Input Completeness Summary

### `celegans`

- labels: present
- PPI: present
- expression: present
- orthologs: present
- sublocalization: present
- build status: success
- train/eval/export status: success

### `human`

- labels: present
- PPI: present
- expression: present
- orthologs: present
- sublocalization: present
- build status: success
- train/eval/export status: success

### `fgraminearum`

- labels: present
- PPI: present
- expression: present
- orthologs: present
- sublocalization: present
- build status: success
- train/eval/export status: success

## Key Metrics

- `celegans_original_replay`: AUROC `0.2418`, AUPRC `0.0589`
- `human_original_replay`: AUROC `0.7965`, AUPRC `0.3610`
- `fgraminearum_original_replay`: AUROC `0.5307`, AUPRC `0.0577`
- `scerevisiae_original_smoke`: AUROC `0.8079`, AUPRC `0.5331`

## Replay Outcome

- multispecies legacy replay is now operational inside ProGATE_v2
- all three requested original replays completed dataset build, training, evaluation, and compatibility export
- the compatibility layer remains intentionally limited to the original feature contract and original GAT family

## Observed Limits By Species

- `celegans`: pipeline is runnable, but current replay quality is weak
- `human`: replay is strongest among the new Phase 1.5 additions, but also the heaviest
- `fgraminearum original`: runnable, but much weaker than support-species replay quality

## Known Limits

- no ESM / ESMC / ProtT5 / ProtBERT migration
- no extended feature layer
- no baseline zoo migration
- no runner unification
