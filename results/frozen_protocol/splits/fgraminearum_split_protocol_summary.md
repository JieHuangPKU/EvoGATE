# Fusarium Split Protocol Summary

## Shared Rules
- split type: stratified fixed split
- split seed: 20260409
- test fraction: 0.20
- val fraction: 0.10

## Mainline
- mainline split: `fgraminearum_newlabel_split.tsv`
- counts: train=8375, val=1197, test=2393
- class balance: train_pos=768, val_pos=110, test_pos=219

## Legacy Comparison
- legacy split: `fgraminearum_oldlabel_split.tsv`
- counts: train=7832, val=1119, test=2238
- class balance: train_pos=307, val_pos=44, test_pos=88

## Caveats
- The old regime exists only for replay and manuscript comparison.
- The new regime is the only Fusarium split used by the new mainline benchmark protocol.
