# Label Scarcity Benchmark Report

- protocol: `fgraminearum_newlabel`
- feature_setting: `ORT_EXP_SUB_ESM2`
- main robustness answer at 10% labels: `MLP`
- strongest retention answer: `DC`
- models with retention_AUPRC < 0.50 at 10%: `none`

## How To Run
- Run `scripts/run_label_scarcity_benchmark.sh` to build the benchmark.
- Run `Rscript src/plot/plot_label_scarcity.R --summary-dir results/Figure2_label_scarcity/summary --output-dir results/Figure2_label_scarcity/plots` to regenerate the figures.

## Interpretation
- Judge robustness primarily from `performance_retention_AUPRC` and secondarily from `AUPRC@10%`.
- If `GraphSAGE` is top-ranked in both low-label AUPRC and retention, the benchmark supports the claim that graph-based learning is the most robust under label scarcity.
- Classical ML degradation is reflected by lower `AUPRC@10%` and larger `performance_drop_AUPRC` relative to `GraphSAGE`.
- `N2V_MLP` should be interpreted as a shallow topology-aware buffer model between GraphSAGE and classical/tabular baselines.
- `DC` and `CC` provide topology-only reference floors rather than competitive predictors.
