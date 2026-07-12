# Yeast Essential Gene List Interface

This directory contains one generated input table per yeast proteome used by the workflow:

- `Scerevisiae.essential_genes.tsv`
- `Spombe.essential_genes.tsv`

## Source

The current workflow builds these files from:

- `/data276/jiehuang/fungi/Fusarium/Evidence/OGEE/gene_essentiality.txt`

using only:

- `Ref_db == SGD` and `essentiality == E` for `Scerevisiae`
- `Ref_db == PomBase` and `essentiality == E` for `Spombe`

## Output format

Each generated file is tab-delimited with:

- `gene_id`
- `source_ref_db`
- `source_locus`
- `source_gene`
- `match_type`

## Matching rule

`gene_id` is the FASTA-header-compatible ID that will be used downstream.

- For `Scerevisiae`, OGEE `locus` values are SGD IDs such as `S000000032`, which are matched against `SGDID:` tokens in the local FASTA headers and then converted to the local primary sequence IDs such as `YAL034W-A`.
- For `Spombe`, OGEE `locus` values such as `SPAC1002.04c` are matched against the local FASTA primary IDs after removing the transcript suffix from headers like `SPAC1002.04c.1:pep`, and the final `gene_id` is written as the real FASTA primary ID `SPAC1002.04c.1:pep`.

Detailed audit outputs are written to:

- `results/derived_labels/yeast_essential_mapping_audit.scerevisiae.tsv`
- `results/derived_labels/yeast_essential_mapping_audit.spombe.tsv`
- `results/derived_labels/yeast_essential_mapping_summary.tsv`
- `results/derived_labels/yeast_essential_mapping_summary.md`
