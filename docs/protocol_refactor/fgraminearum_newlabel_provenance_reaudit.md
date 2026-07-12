# Fusarium Newlabel Provenance Re-Audit

## Conclusion

`data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` is the true upstream source table for the yeast-transfer component.

The repo-local table is byte-identical to the external historical table:

- repo-local: `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
- external historical location: `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`

The original derivation logic lives in:

- `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/workflow/scripts/derive_ph1_yeast_essential_labels.py`

That script does not produce the final newlabel positives directly. It produces a PH-1 `XP_*` candidate table with columns such as:

- `ph1_gene_id`
- `orthogroup_id`
- `yeast_essential_support_class`
- `weak_positive_confidence`

## Relationship To `positive_set_P1.tsv`

`results/label_rebuild_experiments/labels/positive_set_P1.tsv` is a downstream materialized set, not the raw yeast-transfer source.

Historical construction logic was confirmed in `src/train/run_label_rebuild_compare.py`:

- read `ph1_yeast_essential_ortholog_labels.tsv`
- merge it with the older `fusarium_orthology_id_bridge.tsv`
- keep `weak_positive_confidence == high` and `mapping_status == exact`
- define `high_set` from the resulting canonical genes
- define `lethal_set` separately from the PHI-derived lethal table
- define `positive_set_P1 = lethal_set ∪ high_set`

So the provenance chain is:

1. `derive_ph1_yeast_essential_labels.py`
2. `ph1_yeast_essential_ortholog_labels.tsv`
3. bridge from `XP_*` to canonical `FGRAMPH1_*`
4. `high_set`
5. `positive_set_P1.tsv`
6. `phase2b_new_label` and later frozen/newlabel materializations

## What Defines The Historical Expected 1096 Count

The expected `1096` does not come directly from `ph1_yeast_essential_ortholog_labels.tsv`.

It is defined by the historical `positive_set_P1.tsv` composition:

- `55` lethal-only positives
- `1019` weak-positive-only positives
- `22` overlap positives

Equivalent summary:

- historical lethal set size: `77`
- historical yeast-transfer-supported set size: `1041`
- historical union: `1096`

## What Was Wrong In The Current Protocolized Build

The protocolized build was using the repo-local yeast-transfer table, but the new bridge was missing one of the historical mapping paths:

- the direct `XP_* -> FGRAMPH1_*` route preserved in `/data276/jiehuang/fungi/Fusarium/Evidence/00_idmap/FG_gene_id_unified_map.tsv`

As a result, many high-confidence `XP_*` rows that were historically mapped into `positive_set_P1.tsv` were either:

- left unresolved
- or mapped through weaker local heuristics to the wrong neighboring canonical gene

## Repair Applied

The bridge was updated so that:

- the unified ID map is an explicit workflow input
- `XP_* -> FGRAMPH1_*` unique mappings from the unified map are protocolized
- unified-map FGSG aliases are added into the source-to-canonical mapping table
- the unified path is preferred over weaker `FGSG`-header-only mappings
- the unified path is also preferred over legacy sequence-only mappings when the direct `XP_*` mapping is unique and there is no higher-confidence sequence-plus-FGSG agreement

Files changed for this repair:

- `configs/fgraminearum_label_materialization.yaml`
- `workflow/fgraminearum_label_materialization.smk`
- `src/data/build_fgraminearum_newlabel_bridge.py`

## Result

After the repair:

- bridge high-confidence resolved unique genes: `1045`
- bridge unresolved high-confidence rows: `8`
- processed `newlabel` positive count: `1097`

The remaining gap versus the historical `1096` expectation is no longer on the yeast-transfer side. It is now a PHI lethal provenance mismatch:

- `145` historically missing expected genes were recovered
- all recovered expected genes came from the yeast-transfer side except one `lethal;weak_positive` overlap gene
- the remaining `13` unresolved expected genes are all historical lethal-only genes

## Output Audits

Generated audit outputs:

- `data/processed/essential_gene/fgraminearum/newlabel/missing_from_expected_1096.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/recovered_genes.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/still_unresolved_genes.tsv`
