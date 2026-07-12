# Figure 1C Source-Resolved Transfer Sankey

Figure 1C is a source-resolved transfer Sankey for the `fgraminearum` new label positive set. It is designed to answer both source composition and transfer path questions: where the final positives came from, and how yeast essential evidence passed through supported orthogroups and canonical Fusarium mapping into final retained positives.

## Actual Input Files
- `data/processed/essential_gene/fgraminearum/newlabel/labels.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/label_construction_audit.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/summary.tsv`
- `data/processed/essential_gene/fgraminearum/oldlabel/labels.tsv`
- `data/processed/essential_gene/fgraminearum/oldlabel/summary.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/high_confidence_yeast_transfer_candidates.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/unresolved_high_confidence_ids.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/bridge_summary.tsv`
- `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`

## Canonical Final Label Source
`data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv` is used as the canonical final positive source because it is the materialized newlabel positive subset with final retained canonical IDs, construction buckets, supporting XP IDs, orthogroup IDs, and bridge methods. `data/processed/essential_gene/fgraminearum/newlabel/labels.tsv` is used for the canonical newlabel total, and `data/processed/essential_gene/fgraminearum/newlabel/label_construction_audit.tsv` is retained as the construction audit table that records the final label-construction decision fields. The summary box totals come from `data/processed/essential_gene/fgraminearum/newlabel/summary.tsv` and `data/processed/essential_gene/fgraminearum/oldlabel/summary.tsv`.

## Source Assignment
The Sankey uses a primary-source unique assignment, so path counts close and no positive gene is silently double counted.

- `S. cerevisiae only`: final positive has `support_from_protocolized_bridge=true` and its high-confidence bridge row has `yeast_essential_support_class=scer_only`.
- `S. pombe only`: final positive has `support_from_protocolized_bridge=true` and `yeast_essential_support_class=spom_only`.
- `Shared by both yeasts`: final positive has `support_from_protocolized_bridge=true` and `yeast_essential_support_class=both`.
- `PHI-base essential/lethal`: final positive has no protocolized yeast-transfer primary support but is present in `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`.
- Yeast + PHI overlap genes are assigned to the yeast primary source, while `figure1c_source_audit.tsv` records `has_phi_essential_lethal_support=true`.

## Transfer Path Definition
- Layer 1 `Support source`: primary source category above.
- Layer 2 `Supported orthogroups`: yeast sources use the real `orthogroup_id` support summarized into source-specific orthogroup pools for plotting; PHI-only positives use `Direct PHI evidence (no yeast orthogroup)`.
- Layer 3 `Mapped Fusarium genes`: canonical `fgraminearum::FGRAMPH1_*` genes retained in `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv`; unresolved high-confidence yeast rows are displayed as `Unresolved bridge IDs (not retained)` and terminate at `Non-essential / excluded`.
- Layer 4 `Final label class`: `Essential (positive)` for retained positives and `Non-essential / excluded` for unresolved high-confidence transfer candidates that did not enter the final positive set.

## Positive Source Composition
- PHI-base essential/lethal: 52
- S. cerevisiae only: 117
- S. pombe only: 399
- Shared by both yeasts: 529

## Summary Box Values
- old positives: 439
- new positives: 1097
- high-confidence transferred positives: 1045
- old total: 11189
- new total: 11965
- PHI-supported positives: 77
- yeast-transfer-supported positives: 1045

## Stage Counts
- Final label class | Essential (positive): 1097
- Final label class | Non-essential / excluded: 8
- Mapped Fusarium genes | PHI-mapped Fusarium genes (n=52): 52
- Mapped Fusarium genes | S. cerevisiae only
mapped Fusarium genes (n=117): 117
- Mapped Fusarium genes | S. pombe only
mapped Fusarium genes (n=399): 399
- Mapped Fusarium genes | Shared by both yeasts
mapped Fusarium genes (n=529): 529
- Mapped Fusarium genes | Unresolved bridge IDs
(not retained): 8
- Support source | PHI-base essential/lethal: 52
- Support source | S. cerevisiae only: 118
- Support source | S. pombe only: 402
- Support source | Shared by both yeasts: 533
- Supported orthogroups | Direct PHI evidence
(no yeast orthogroup) (n=52): 52
- Supported orthogroups | S. cerevisiae only
supported orthogroups (n=1): 1
- Supported orthogroups | S. cerevisiae only
supported orthogroups (n=117): 117
- Supported orthogroups | S. pombe only
supported orthogroups (n=3): 3
- Supported orthogroups | S. pombe only
supported orthogroups (n=398): 399
- Supported orthogroups | Shared by both yeasts
supported orthogroups (n=4): 4
- Supported orthogroups | Shared by both yeasts
supported orthogroups (n=529): 529

## Outputs
- `figure1c_sankey_long.tsv`: ggsankey-ready transition table with `x`, `node`, `next_x`, `next_node`, `value`, and source metadata.
- `figure1c_stage_counts.tsv`: node counts for each plotted layer.
- `figure1c_summary_box.tsv`: right-side summary box statistics.
- `figure1c_source_audit.tsv`: one row per final positive gene, preserving yeast support class, PHI support, old positive overlap, orthogroup IDs, XP IDs, and bridge method.
