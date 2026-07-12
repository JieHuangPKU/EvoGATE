# Supplemental Table 1 species feature statistics

This file follows the legacy four-species feature-table organization and expands it to the six protocol objects used in the current manuscript.

## Supplemental Table 1A. Feature dimension summary for the six protocol objects used in this study

Feature Type | F. graminearum (new label) | F. graminearum (old label) | S. cerevisiae | H. sapiens | C. elegans | D. melanogaster
--- | --- | --- | --- | --- | --- | ---
Ortholog | 167 | 167 | 62 | 162 | 233 | 203
Expression | 60 | 60 | 36 | 64 | 30 | 66
Localization | 12 | 12 | 13 | 13 | 13 | 13

## Supplemental Table 1B. Network and label statistics for the six protocol objects used in this study

Statistic | F. graminearum (new label) | F. graminearum (old label) | S. cerevisiae | H. sapiens | C. elegans | D. melanogaster
--- | --- | --- | --- | --- | --- | ---
N.nodes | 14146 | 14146 | 5628 | 18458 | 13682 | 6540
N.edges | 1439790 | 1439790 | 1281599 | 6321712 | 1937760 | 517640
N.labeled genes | 11965 | 11189 | 5628 | 18458 | 13682 | 6540
N.essential genes | 1097 | 439 | 1050 | 1827 | 578 | 402
N.train genes | 8375 | 7832 | 3939 | 12920 | 9576 | 4578
N.val genes | 1197 | 1119 | 563 | 1846 | 1369 | 654
N.test genes | 2393 | 2238 | 1126 | 3692 | 2737 | 1308

## Source audit

### F. graminearum (new label)

- protocol id: `fgraminearum_newlabel`
- orthologs: `data/processed/OR/fgraminearum/orthologs.csv`
- expression: `data/processed/EXP/fgraminearum/profile.csv`
- sublocalization: `data/processed/LC/fgraminearum/subloc.csv`
- ppi: `data/processed/PPI/fgraminearum/string.csv`
- labels: `results/frozen_protocol/labels/fgraminearum_newlabel.tsv`
- splits: `results/frozen_protocol/splits/fgraminearum_newlabel_split.tsv`

### F. graminearum (old label)

- protocol id: `fgraminearum_oldlabel`
- orthologs: `data/processed/OR/fgraminearum/orthologs.csv`
- expression: `data/processed/EXP/fgraminearum/profile.csv`
- sublocalization: `data/processed/LC/fgraminearum/subloc.csv`
- ppi: `data/processed/PPI/fgraminearum/string.csv`
- labels: `results/frozen_protocol/labels/fgraminearum_oldlabel.tsv`
- splits: `results/frozen_protocol/splits/fgraminearum_oldlabel_split.tsv`

### S. cerevisiae

- protocol id: `scerevisiae`
- orthologs: `data/processed/OR/scerevisiae/orthologs.csv`
- expression: `data/processed/EXP/scerevisiae/profile.csv`
- sublocalization: `data/processed/LC/scerevisiae/subloc.csv`
- ppi: `data/processed/PPI/scerevisiae/string.csv`
- labels: `results/frozen_protocol/labels/scerevisiae_labels.tsv`
- splits: `results/frozen_protocol/splits/scerevisiae_split.tsv`

### H. sapiens

- protocol id: `human`
- orthologs: `data/processed/OR/human/orthologs.csv`
- expression: `data/processed/EXP/human/profile.csv`
- sublocalization: `data/processed/LC/human/subloc.csv`
- ppi: `data/processed/PPI/human/string.csv`
- labels: `results/frozen_protocol/labels/human_labels.tsv`
- splits: `results/frozen_protocol/splits/human_split.tsv`

### C. elegans

- protocol id: `celegans`
- orthologs: `data/processed/OR/celegans/orthologs.csv`
- expression: `data/processed/EXP/celegans/profile.csv`
- sublocalization: `data/processed/LC/celegans/subloc.csv`
- ppi: `data/processed/PPI/celegans/string.csv`
- labels: `results/frozen_protocol/labels/celegans_labels.tsv`
- splits: `results/frozen_protocol/splits/celegans_split.tsv`

### D. melanogaster

- protocol id: `dmelanogaster`
- orthologs: `data/processed/OR/melanogaster/orthologs.csv`
- expression: `data/processed/EXP/melanogaster/profile.csv`
- sublocalization: `data/processed/LC/melanogaster/subloc.csv`
- ppi: `data/processed/PPI/melanogaster/string.csv`
- labels: `results/frozen_protocol/labels/dmelanogaster_labels.tsv`
- splits: `results/frozen_protocol/splits/dmelanogaster_split.tsv`

