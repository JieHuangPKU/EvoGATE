# Fusarium Newlabel Count Gap Audit

## Scope

This audit compares:

- historical expected positives from `results/label_rebuild_experiments/labels/positive_set_P1.tsv`
- the broken pre-repair processed newlabel positive set
- the rebuilt post-repair processed newlabel positive set

Pre-repair baseline snapshot used for this audit:

- `/tmp/fgraminearum_newlabel_positive_genes.pre_reaudit.tsv`

## Historical Reference

Historical expected positive count:

- `1096`

Historical composition:

- `55` lethal-only
- `1019` weak-positive-only
- `22` lethal-plus-weak overlap

## Broken State Before Repair

Pre-repair processed newlabel positive count:

- `950`

Gap versus historical expected:

- `158` expected genes absent from the processed newlabel positive set

Breakdown of the `158` missing expected genes:

- `145` weak-positive-only
- `13` lethal-only
- `0` overlap genes

This means the main failure was the yeast-transfer bridge, not the final union rule itself.

## Where The 145 Yeast-Transfer Genes Were Lost

Before the repair, the bridge showed:

- `1056` high-confidence `XP_*` rows in the yeast-transfer table
- only `899` protocolized unique high-confidence canonical genes
- `157` unresolved high-confidence rows

Unresolved reason breakdown before repair:

- `67` `no_sequence_or_fgsg_bridge`
- `47` `header_fgsg_not_mapped_to_canonical`
- `43` `unclassified`

Additional issue:

- `2` expected genes were not unresolved; they were mis-mapped to neighboring canonical genes by lower-quality bridge precedence

So the historical weak-positive gap was:

- `143` truly unresolved high-confidence genes
- `2` wrong-gene assignments

## Repair Outcome

Recovered expected genes after rebuild:

- `145`

Recovered composition:

- `144` weak-positive-only
- `1` overlap gene (`lethal;weak_positive`)

Recovered bridge-method breakdown:

- `106` `ncbi_protein_id_via_unified_map`
- `38` `unified_plus_fgsg_support`
- `1` `unified_map_preferred_over_legacy_sequence_bridge`

Remaining expected genes still absent after rebuild:

- `13`

All `13` remaining unresolved expected genes are lethal-only historical positives.

These are listed in:

- `data/processed/essential_gene/fgraminearum/newlabel/still_unresolved_genes.tsv`

## Final Counts

Post-repair processed newlabel positive count:

- `1097`

Post-repair processed newlabel negative count:

- `10868`

Post-repair total labeled count:

- `11965`

Relative to the historical expected `1096` positives:

- `13` expected historical genes remain absent
- `14` protocolized positives are now present that were not in the historical `positive_set_P1.tsv`
- net positive count difference = `+1`

## Interpretation

The historical `145`-gene yeast-transfer gap has been repaired.

The remaining mismatch to the historical expectation is not due to the yeast-transfer upstream anymore. It is due to lethal-side provenance differences:

- the current protocolized lethal source still does not exactly reproduce the historical 77-gene lethal table
- the remaining unresolved expected genes are all historical lethal-only genes

## Files Produced

- `data/processed/essential_gene/fgraminearum/newlabel/missing_from_expected_1096.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/recovered_genes.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/still_unresolved_genes.tsv`
