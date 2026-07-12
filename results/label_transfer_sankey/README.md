# Label Transfer Sankey Data

## Actual Inputs Used
- Canonical yeast-transfer source: `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
- Canonical bridge source: `data/processed/essential_gene/fgraminearum/bridge/high_confidence_yeast_transfer_candidates.tsv` and `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv`
- Canonical final label source: `data/processed/essential_gene/fgraminearum/newlabel/label_construction_audit.tsv`
- Canonical final positive subset: `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv`
- Canonical final label audit source: `data/processed/essential_gene/fgraminearum/newlabel/label_construction_audit.tsv`
- Lethal positive provenance source: `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`
- Orthogroup membership source: `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.tsv`
- Orthogroup gene-count source: `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.GeneCount.tsv`
- OrthoFinder single-copy list: `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups_SingleCopyOrthologues.txt`
- OrthoFinder source provenance summary used for path validation: `data/derived_labels/ph1_yeast_essential_ortholog_labels.summary.md`

### Full File List
- `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/high_confidence_yeast_transfer_candidates.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/bridge_summary.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/bridge_source_manifest.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/unresolved_high_confidence_ids.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/labels.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/label_construction_audit.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/source_manifest.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/summary.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`
- `data/derived_labels/proteome_manifest.tsv`
- `data/derived_labels/ph1_yeast_essential_ortholog_labels.summary.md`
- `data/derived_labels/yeast_essential/Scerevisiae.essential_genes.tsv`
- `data/derived_labels/yeast_essential/Spombe.essential_genes.tsv`
- `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.tsv`
- `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.GeneCount.tsv`
- `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups_SingleCopyOrthologues.txt`

## Canonical Source Decisions
- Canonical final label source was assigned to `data/processed/essential_gene/fgraminearum/newlabel/label_construction_audit.tsv` because it is the full materialized `newlabel` construction audit spanning all retained and excluded genes in the final label universe, while `data/processed/essential_gene/fgraminearum/newlabel/labels.tsv` is a narrower downstream label manifest.
- Canonical orthogroup membership source was assigned to `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.tsv` because the repo-local summary `data/derived_labels/ph1_yeast_essential_ortholog_labels.summary.md` records that exact OrthoFinder results directory as the upstream source for `ph1_yeast_essential_ortholog_labels.tsv`.
- Canonical XP to `fgraminearum::FGRAMPH1_*` mapping source was assigned to `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv` plus the protocolized positive subset `data/processed/essential_gene/fgraminearum/bridge/high_confidence_yeast_transfer_candidates.tsv`.

## Layer Definitions And Counting Unit
- Counting unit for `sankey_gene_level_long.tsv`: one aggregated real evidence record per `(yeast_support_source, orthogroup_id, canonical_fusarium_gene_id)` after joining the upstream yeast-transfer table to the protocolized bridge. Multiple `XP_*` rows collapsing to the same orthogroup and canonical gene are aggregated, with `source_row_count` preserving the original row multiplicity.
- Layer 1 `Yeast support source`: taken directly from `yeast_essential_support_class` in `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` and restricted to `scer_only`, `spom_only`, `both`.
- Layer 2 `Orthogroups with essential support`: `orthogroup_id` from `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`; occupancy and copy metrics come from the same table (`fg_presence_count_18`, `fg_occupancy_18`, `fg_single_copy_fraction_18`, `is_strict_core_18`, `is_exact_single_copy_core_18`) and were cross-linked to `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.tsv` / `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups.GeneCount.tsv` / `/data276/jiehuang/fungi/Fusarium/orthofinder_essential_workflow/results/orthofinder_results/run_20260405T213342_139369/Results_Apr05/Orthogroups/Orthogroups_SingleCopyOrthologues.txt`.
- Layer 3 `Mapped Fusarium genes`: `canonical_fusarium_gene_id` from the protocolized bridge join to `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv`; only bridge-resolved rows can enter the drawn Sankey edges because layer 3 is constrained to canonical `fgraminearum::FGRAMPH1_*` IDs.
- Layer 4 `Final label class`: assigned by joining canonical genes to the materialized newlabel outputs.

## Final Class Rules
- `Positive_High`: bridge-resolved canonical gene is present in `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv` and also present in `data/processed/essential_gene/fgraminearum/bridge/high_confidence_yeast_transfer_candidates.tsv`.
- `Positive`: bridge-resolved canonical gene is present in `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv` but not present in the protocolized high-confidence transfer set. In practice these are yeast-supported mapped genes retained in final positives through non-high-confidence logic, typically lethal-backed positives.
- `Excluded`: everything else in the yeast-support pool. This includes bridge-resolved genes that ended up in final negatives or outside the final label table, plus unresolved bridge rows in the long table.

## Fusarium Filtering Definition
- `passes_fusarium_filtering = yes` means `weak_positive_confidence != none` in `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`.
- `passes_fusarium_filtering = no` means `weak_positive_confidence == none`.
- This script does not re-infer the filter. It reuses the upstream materialized confidence field exactly.

## Copy-Status Definition
- `exact_single_copy_core`: `is_exact_single_copy_core_18 == 1`
- `near_single_copy_core`: strict core plus `fg_max_copy_18 <= 2` and `fg_single_copy_fraction_18 >= 0.90`
- `strict_core_multicopy`: strict core but not exact/near single-copy
- `high_occupancy_mostly_single_copy`: `fg_occupancy_18 >= 0.80` and `fg_single_copy_fraction_18 >= 0.75`
- `noncore_or_multicopy`: remaining rows

## Output Summary
- Upstream yeast-supported rows (`has_any_yeast_essential == 1`): 1378
- Aggregated gene-level evidence rows: 1373
- Bridge-resolved aggregated rows: 1363
- Bridge-unresolved aggregated rows: 10
- `Positive_High` aggregated rows: 1045
- `Positive` aggregated rows: 7
- `Excluded` aggregated rows: 321
- Unique canonical genes in final positive set within this Sankey pool: 1052
- Unique canonical genes in protocolized high-confidence set within this Sankey pool: 1045
- Sankey edge rows written: 3977
- Stage-count rows written: 2624

## Caveats
- The repo does not carry a repo-local copy of `Orthogroups.tsv`; the canonical membership table was read directly from the upstream OrthoFinder path documented in `data/derived_labels/ph1_yeast_essential_ortholog_labels.summary.md`.
- Unresolved bridge rows are retained in `sankey_gene_level_long.tsv` with `bridge_resolved = no` and `final_label_class = Excluded`, but they are not drawn into `sankey_aggregated_edges.tsv` because layer 3 is restricted to canonical `fgraminearum::FGRAMPH1_*` gene IDs.
- `Orthogroups_SingleCopyOrthologues.txt` is used as supplemental context only. The primary occupancy and copy metrics still come from the materialized yeast-transfer table so the Sankey remains traceable to the exact project pipeline outputs.
