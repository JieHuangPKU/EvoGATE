# Fusarium Newlabel Processed Label Summary

## Regime Definition
The `newlabel` regime is the current mainline Fusarium essential-gene protocol.
It corresponds to the newer lethal plus evolution regime used for the publication-grade frozen benchmark.

## Exact Source Files
- Lethal-positive provenance table: `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`
- Repo-local yeast-transfer table: `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
- Protocolized protein-to-canonical bridge: `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv`
- Protocolized high-confidence transfer candidate table: `data/processed/essential_gene/fgraminearum/bridge/high_confidence_yeast_transfer_candidates.tsv`
- Mirrored evidence table used to document virulence exclusion logic: `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`

## Construction Logic
The final positive set is defined as the union of two components.
First, the lethal component is taken from the preserved lethal-positive table, which records the PHI-supported direct lethal, failed-deletion, and non-viable evidence that survived canonical mapping review.
Second, the transfer component is rebuilt from the repo-local PH-1 yeast-transfer table after explicit XP-to-canonical bridge reconstruction.
Each PH-1 `XP_*` protein accession is mapped into the final `fgraminearum::FGRAMPH1_*` space through the protocolized bridge under `data/processed/essential_gene/fgraminearum/bridge/`.
The rebuilt high-confidence transfer-positive component then retains only those PH-1 proteins whose bridge result resolves to exactly one canonical Fusarium gene.

The negative set follows the same protocol definition as the mainline rebuild.
It is the weak-confidence `none` pool after canonical bridging, followed by removal of genes flagged by virulence/pathogenicity evidence and removal of all genes retained in the positive set.
The mirrored master evidence table is recorded here so the biological exclusion rule remains explicit in the processed artifact lineage.

## How PHI Evidence Contributes
PHI evidence contributes through the lethal component.
Genes carrying direct lethal, failed-deletion, or non-viable mutant evidence and surviving canonical mapping review are retained in the lethal-positive provenance table and therefore enter the final positive set directly.
Those genes anchor the experimentally supported essential core of the regime.

## How Derived Labels Contribute
Derived labels contribute through the high-confidence yeast-transfer component.
Genes with strong orthology-based support from yeast essentiality transfer are first bridged from PH-1 `XP_*` protein IDs to canonical Fusarium genes, then deduplicated in canonical gene space, and finally combined with the lethal PHI-supported set.
This removes the previous dependence on the historical `positive_set_P1.tsv` snapshot and makes the bridge logic explicit.

## Final Counts
- positives: 1097
- negatives: 10868
- total labeled genes: 11965
- lethal PHI-supported positives: 77
- high-confidence yeast-transfer-supported positives retained after subtracting lethal genes: 1020
- unresolved high-confidence transfer rows excluded by the protocolized bridge: 8
- historical snapshot positive count for comparison: 1096
- historical snapshot negative count for comparison: 10270

## Split Policy
- split seed: 20260409
- test fraction: 0.20
- val fraction: 0.10
- split version: `frozen_protocol_v1_seed20260409_test20_val10`
- train/val/test counts: 8375/1197/2393

## Caveats and Limitations
This workflow is now rebuilt from the protocolized bridge rather than from `positive_set_P1.tsv`.
Because the protocolized bridge intentionally keeps only exact or auditable inferred mappings, rebuilt counts can differ from the historical snapshot when some `XP_*` rows still lack enough local evidence to resolve to a unique canonical Fusarium gene.
Any unresolved transfer rows are written to the bridge audit outputs and can be revisited later when additional annotation evidence is protocolized.
