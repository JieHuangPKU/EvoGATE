# Frozen Split Protocol Summary

- split type: stratified fixed split
- split seed: 20260409
- test fraction: 0.20
- val fraction: 0.10
- no model may generate its own split internally after this refactor

| protocol | species | regime | train | val | test | train_pos | val_pos | test_pos | output |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| human | human | standard | 12920 | 1846 | 3692 | 1279 | 183 | 365 | `human_split.tsv` |
| celegans | celegans | standard | 9576 | 1369 | 2737 | 404 | 58 | 116 | `celegans_split.tsv` |
| scerevisiae | scerevisiae | standard | 3939 | 563 | 1126 | 735 | 105 | 210 | `scerevisiae_split.tsv` |
| dmelanogaster | dmelanogaster | standard | 4578 | 654 | 1308 | 282 | 40 | 80 | `dmelanogaster_split.tsv` |
| fgraminearum_oldlabel | fgraminearum | oldlabel | 7832 | 1119 | 2238 | 307 | 44 | 88 | `fgraminearum_oldlabel_split.tsv` |
| fgraminearum_newlabel | fgraminearum | newlabel | 8375 | 1197 | 2393 | 768 | 110 | 219 | `fgraminearum_newlabel_split.tsv` |
