# Fusarium Label Regime Comparison

## Scope
This file contrasts the preserved historical `oldlabel` regime with the current mainline `newlabel` regime after protocolized materialization into `data/processed`.

## Count Summary
| regime   |   positive_count |   negative_count |   total_count |   train_count |   val_count |   test_count | definition_summary                                                | positive_definition                                                                     | negative_definition                                                 |
|:---------|-----------------:|-----------------:|--------------:|--------------:|------------:|-------------:|:------------------------------------------------------------------|:----------------------------------------------------------------------------------------|:--------------------------------------------------------------------|
| oldlabel |              439 |            10750 |         11189 |          7832 |        1119 |         2238 | historical lethal plus virulence replay from old440 mapping audit | mapped historical gene_list Target=1 positives                                          | preserved negative pool after overlap removal with old440 positives |
| newlabel |             1097 |            10868 |         11965 |          8375 |        1197 |         2393 | current lethal plus evolution regime                              | lethal PHI-supported positives union high-confidence yeast-transfer-supported positives | weak-none pool after virulence/pathogenicity and positive exclusion |

## Key Definitional Differences
- `oldlabel` is a historical replay. Its positives come from the old440 gene-list mapping audit and its negatives are the preserved negative pool after overlap removal.
- `newlabel` is the current mainline regime. Its positives combine lethal PHI-supported genes with the high-confidence yeast-transfer-supported component, and its negatives follow the weak-none minus virulence/pathogenicity and positive exclusion rule preserved in the canonical materialized set.
- The two regimes should remain separate because they answer different scientific questions: legacy back-comparison versus current publication-grade benchmark evaluation.
