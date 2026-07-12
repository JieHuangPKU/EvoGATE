# Fusarium Oldlabel Processed Label Summary

## Regime Definition
The `oldlabel` regime is the historical lethal plus virulence replay used for manuscript back-comparison.
It is retained as an explicit legacy branch and is not the intended mainline Fusarium benchmark.

## Exact Source Files
- Historical gene-list replay audit: `data/interim/protocol_refactor/fgraminearum_label_materialization/old440_mapping_audit.tsv`
- Historical old440 summary table: `data/interim/protocol_refactor/fgraminearum_label_materialization/old440_label_summary.tsv`
- Preserved historical source gene list path recorded for provenance: `/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/EssentialGenes/gene_list.txt`
- Base negative pool reused for the replay: `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` via the rebuilt newlabel none-derived negatives

## Construction Logic
The oldlabel positive set is reconstructed from the preserved old440 replay audit.
Each retained positive gene corresponds to a source entry from the historical `gene_list.txt` whose `Target` label was 1 and whose canonical mapping was resolved successfully.
The replay audit preserves the mapping rule and mapping source, allowing the processed artifact to point back to the exact historical replay lineage.

The oldlabel negative set is defined as the preserved negative pool after removal of any gene that appears in the old440 positive replay set.
This reproduces the old lethal plus virulence comparison regime without leaving the mainline benchmark directly dependent on the historical results directory.

## Positive and Negative Definitions
- Positives: genes with `target_label == 1` in the old440 replay audit after canonical mapping.
- Negatives: genes from the preserved negative pool after overlap removal against the old440 positive set.

## Final Counts
- positives: 439
- negatives: 10750
- total labeled genes: 11189
- overlap removed from the negative pool during replay: 118

## Split Policy
- split seed: 20260409
- test fraction: 0.20
- val fraction: 0.10
- split version: `frozen_protocol_v1_seed20260409_test20_val10`
- train/val/test counts: 7832/1119/2238

## Caveats and Limitations
This regime is a controlled historical replay, not a newly curated essentiality protocol.
Its biological meaning should be interpreted as the preserved old lethal plus virulence comparison branch rather than a current mainline recommendation.
