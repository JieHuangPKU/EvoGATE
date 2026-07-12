# Fgraminearum Missing Inputs Checklist

Current status:

- no additional user-supplied file is required to materialize the current frozen `fgraminearum_newlabel` manifest
- the materialized source tables already exist in-repo

## Non-Blocking Migration Gaps

| suggested_path | expected_filename | role | expected_columns | source_or_derived | why_required | downstream_dependency | can_pipeline_proceed_without_it |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` | `ph1_yeast_essential_ortholog_labels.tsv` | repo-local replacement for the legacy external yeast-transfer label table | `ph1_gene_id`, `weak_positive_confidence`, `yeast_essential_support_class` | source table | needed only to rerun raw provenance reconstruction without the old absolute path | `src/train/run_label_rebuild_compare.py` | yes |
| `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv` | `master_evidence_table.preliminary.tsv` | mirrored negative-filter provenance table | `species`, `canonical_gene_id`, `evidence_class`, `evidence_term_raw` | derived mirror | needed only to document and replay virulence exclusion logic inside the repo | `docs/protocol_refactor/fgraminearum_label_audit.md`; legacy rebuild code | yes |

## Conclusion

The frozen mainline benchmark can proceed without any additional user file.
