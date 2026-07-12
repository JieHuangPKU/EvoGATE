# Fgraminearum Label Audit

Update: `results/label_rebuild_experiments/` is now a historical archive / provenance directory.
Current mainline Fusarium label assets are materialized under
`data/processed/essential_gene/fgraminearum/{oldlabel,newlabel}/`.

## Scope

This audit separates the historical Fusarium label regime from the intended new mainline regime.

## Old Label Regime

Suggested explicit name:

- `fgraminearum_oldlabel`

Role:

- legacy replay
- historical comparison
- manuscript back-comparison against older work

Current frozen outputs:

- `results/frozen_protocol/labels/fgraminearum_oldlabel.tsv`
- `results/frozen_protocol/splits/fgraminearum_oldlabel_split.tsv`

Historical archive sources retained for provenance:

- `results/label_rebuild_experiments/old440/labels/positive_old440.tsv`
- `results/label_rebuild_experiments/old440/labels/negative_old440.tsv`
- `results/label_rebuild_experiments/old440/labels/old440_label_summary.tsv`
- `results/label_rebuild_experiments/old440/labels/old440_mapping_audit.tsv`

Generation logic:

- positive genes were reconstructed from the historical `gene_list.txt` mapping audit
- mapping rule recorded in the materialized audit: `raw_gene_exact_match`
- negatives were taken from the historical negative pool after removing overlap with the positive set

Positive / negative definition:

- positives: historical old440 positive set, conceptually the preserved old lethal + virulence regime
- negatives: historical old440 negative set after overlap removal

Counts:

- positives: 439
- negatives: 10154
- total: 10593

Still used in code/workflows:

- consumed only through the new explicit protocol name `fgraminearum_oldlabel`
- replay provenance remains in `src/train/run_old440_diagnostic.py`

## New Label Regime

Suggested explicit name:

- `fgraminearum_newlabel`

Role:

- intended new mainline Fusarium regime

Current frozen outputs:

- `results/frozen_protocol/labels/fgraminearum_newlabel.tsv`
- `results/frozen_protocol/splits/fgraminearum_newlabel_split.tsv`

Historical archive files that recorded the older rebuild snapshot:

- `results/label_rebuild_experiments/labels/positive_set_P1.tsv`
- `results/label_rebuild_experiments/labels/negative_set.tsv`
- `results/label_rebuild_experiments/labels/positive_set_summary.md`
- `results/label_rebuild_experiments/labels/negative_set_summary.md`
- `results/label_rebuild_experiments/final_recommendation.md`
- `src/eval/prepare_phase2b_new_labels.py`

Is `phase2b_new_label` already materialized?

- yes
- the current mainline new-label regime is materialized under `data/processed/essential_gene/fgraminearum/newlabel/`
- `src/eval/prepare_phase2b_new_labels.py` now copies from those processed materialized tables into the historical `results/phase2b_new_label/labels/` location

Was it reconstructed during this refactor?

- the current mainline no longer uses `results/label_rebuild_experiments/` as the source of frozen manifests
- frozen manifests are built from the processed materialized label assets

Generation logic recorded in code:

- `src/train/run_label_rebuild_compare.py` recorded the historical rebuild logic
- lethal positives are built from the union of `strict29`, `broad79`, and `conflict8`, then filtered to `primary_evidence_term == lethal`
- weak positives are derived from yeast-transfer labels with `weak_positive_confidence == high`
- `positive_set_P1.tsv` = `lethal_set | high_set`
- `negative_set.tsv` = `none_set - virulence_genes - lethal_set - high_set`

Does it correspond to `phase2b_new_label`?

- yes in biological intent
- but the current mainline source of truth is the materialized processed label directory, while `phase2b_new_label` is retained as a historical comparison output location

Is it based on `derived_labels + PHI 94`?

- partially confirmed
- confirmed from current in-repo materialized provenance:
  `positive_set_P1.tsv` combines lethal positives with high-confidence yeast-transfer positives
- not confirmed exactly as “derived_labels + PHI 94” in the currently materialized files
- the current materialized summary reports:
  - lethal-only count: 55
  - yeast-transfer-supported positives: 1041
  - overlap: 22
- therefore the present materialized regime is better described as:
  lethal union high-confidence yeast-transfer positives
- it is not safe to restate this as exactly “derived_labels + PHI 94” without stronger upstream provenance than the repo currently records

Exact positive / negative definition:

- positives:
  all genes in `results/label_rebuild_experiments/labels/positive_set_P1.tsv`
  with `positive_sources` showing `lethal`, `weak_positive`, or both
- negatives:
  all genes in `results/label_rebuild_experiments/labels/negative_set.tsv`
  defined from weak-confidence-none genes after excluding virulence genes and all positive genes

Counts:

- positives: 1096
- negatives: 10270
- total: 11366

## Missing Inputs

No additional file is currently required from the user to materialize the frozen new-label manifest because the processed materialized source tables are already present in-repo.

However, exact end-to-end raw reconstruction provenance is incomplete because:

- `configs/label_rebuild_compare.yaml` still points to an external absolute-path yeast label source
- the negative filtering provenance historically referenced `data_registry/master_evidence_table.preliminary.tsv`

These are migration issues, not blockers for freezing the current publication protocol.
