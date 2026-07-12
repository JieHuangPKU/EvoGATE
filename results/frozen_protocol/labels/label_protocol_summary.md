# Frozen Label Protocol Summary

- protocol_version: `frozen_protocol_v1`
- benchmark species set: human, celegans, scerevisiae, dmelanogaster, fgraminearum
- Fusarium mainline regime: `fgraminearum_newlabel`
- Fusarium legacy replay regime: `fgraminearum_oldlabel`
- Deprecated mainline inputs: `broad79`, `strict29`, `conflict8`, `fgraminearum_gold_positive*`, implicit `old440` defaults

| protocol | species | regime | positives | negatives | total | source_manifest | output | mainline |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| human | human | standard | 1827 | 16631 | 18458 | `data/processed/essential_gene/human/labels.standard.tsv` | `human_labels.tsv` | true |
| celegans | celegans | standard | 578 | 13104 | 13682 | `data/processed/essential_gene/celegans/labels.standard.tsv` | `celegans_labels.tsv` | true |
| scerevisiae | scerevisiae | standard | 1050 | 4578 | 5628 | `data/processed/essential_gene/scerevisiae/labels.standard.tsv` | `scerevisiae_labels.tsv` | true |
| dmelanogaster | dmelanogaster | standard | 402 | 6138 | 6540 | `data/processed/essential_gene/melanogaster/labels.standard.tsv` | `dmelanogaster_labels.tsv` | true |
| fgraminearum_oldlabel | fgraminearum | oldlabel | 439 | 10750 | 11189 | `data/processed/essential_gene/fgraminearum/oldlabel/positive_genes.tsv` | `fgraminearum_oldlabel.tsv` | false |
| fgraminearum_newlabel | fgraminearum | newlabel | 1097 | 10868 | 11965 | `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv` | `fgraminearum_newlabel.tsv` | true |
