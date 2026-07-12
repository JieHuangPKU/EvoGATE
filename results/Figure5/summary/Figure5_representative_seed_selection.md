# Figure5 Representative Seed Selection

- Seeds considered: `1029, 1030, 1031, 1032, 1033`.
- Protocols considered: `fgraminearum_newlabel, scerevisiae`.
- Base metric families: `rescued_gene_count, rescued_gene_fraction, rescued_mean_delta_probability, all_gene_mean_delta_probability, test_auprc_esm2, test_mcc_esm2, test_f1_esm2, test_auroc_esm2`.
- Species-specific metric columns are scored separately, then summed into a single cross-species deviation score.
- Selected representative seed: `1029`.
- Selection principle: closest overall behavior to the five-seed median, not the visually best-looking seed.
- Exact rule: for every protocol-specific selection metric, compute the scaled absolute deviation from the five-seed median and choose the seed with the smallest aggregate cross-species deviation score.
- Representative seed score: `7.159868`.
- Full per-seed table: `results/Figure5/tables/Figure5_representative_seed_selection.tsv`.

## Representative seed metrics

- `fgraminearum` rescued gene count: `22`.
- `fgraminearum` rescued gene fraction: `0.009193`.
- `fgraminearum` mean rescued delta probability: `0.489279`.
- `fgraminearum` Bio+ESM2 test AUPRC: `0.538844`.
- `fgraminearum` Bio+ESM2 test MCC: `0.499461`.
- `scerevisiae` rescued gene count: `11`.
- `scerevisiae` rescued gene fraction: `0.009769`.
- `scerevisiae` mean rescued delta probability: `0.426175`.
- `scerevisiae` Bio+ESM2 test AUPRC: `0.664879`.
- `scerevisiae` Bio+ESM2 test MCC: `0.500094`.
