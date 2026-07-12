# Figure5 new design and results summary

Theme: "ESM2 expands the essentiality prediction manifold and identifies ESM2-specific essential gene candidates."

Protocol: `fgraminearum_newlabel`; model: `GraphSAGE`; seeds: 1029, 1030, 1031, 1032, 1033.
Representative Figure5A seed: `1029`.

## Panel inputs and outputs
### Panel 5A
Inputs:
- `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB/run_1029/predictions.tsv`
- `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1029/predictions.tsv`
Outputs:
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_essential_vs_nonessential.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_essential_vs_nonessential.png`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_esm2_specific.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_esm2_specific.png`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_space_umap_two_panel.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_space_umap_two_panel.png`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_shift_summary.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5A_ESM2_specific_essential_genes.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5A_shared_essential_genes.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_space_audit.md`

### Panel 5B
Inputs:
- `results/Figure5/data/Figure5c_input_vs_hidden_compare_fgraminearum_seed1029_coords.tsv`
- `results/Figure5/data/Figure5a_hidden_umap_error_transition_fgraminearum_seed1029_coords.tsv`
Outputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_B_input_vs_hidden_umap.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5_new_B_input_vs_hidden_umap.png`
- `results/Figure5_new_candidate_prioritization/Figure5_new_B_input_vs_hidden_umap_plot_data.tsv`

### Panel 5C
Inputs:
- `results/Figure5/tables/Figure5d_group_ablation_global_summary.tsv`
- `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB/run_1029/feature_schema.tsv`
- `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1029/feature_schema.tsv`
Outputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_C_grouped_feature_importance.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_C_grouped_feature_importance.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5_new_C_grouped_feature_importance.png`

### Panel 5D
Inputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_D_per_seed_rank_percentiles.tsv`
Outputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_candidate_rank_table.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_D_topk_overlap_summary.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_D_topk_candidate_overlap.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5_new_D_topk_candidate_overlap.png`

### Panel 5E
Inputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_candidate_rank_table.tsv`
Outputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_E_esm2_rescued_candidates.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_E_rank_shift_summary.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_E_rank_shift_rescue_plot.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5_new_E_rank_shift_rescue_plot.png`

### Panel 5F
Inputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_candidate_rank_table.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_E_esm2_rescued_candidates.tsv`
Outputs:
- `results/Figure5_new_candidate_prioritization/Figure5_new_F_biological_profile_gene_level.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_F_biological_profile_summary.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_F_enrichment_results.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5_new_F_biological_profile_heatmap.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5_new_F_biological_profile_heatmap.png`

## Figure5A Top-K summary

| K   |   n_top_noESM2 |   n_top_ESM2 |   n_shared |   n_ESM2_specific |   n_lost_after_ESM2 |   jaccard |   median_delta_score_ESM2_specific |   mean_delta_score_ESM2_specific |
|:----|---------------:|-------------:|-----------:|------------------:|--------------------:|----------:|-----------------------------------:|---------------------------------:|
| 100 |            100 |          100 |         23 |                77 |                  77 |  0.129944 |                          0.0910028 |                         0.154747 |
| 200 |            200 |          200 |         63 |               137 |                 137 |  0.186944 |                          0.106276  |                         0.180848 |
| 500 |            500 |          500 |        222 |               278 |                 278 |  0.285347 |                          0.154395  |                         0.224935 |
| 5%  |            652 |          652 |        332 |               320 |                 320 |  0.341564 |                          0.179755  |                         0.25653  |

## Figure5A top 20 ESM2-specific genes by delta_score

| gene_id           |   score_ORT_EXP_SUB |   score_ORT_EXP_SUB_ESM2 |   rank_ORT_EXP_SUB |   rank_ORT_EXP_SUB_ESM2 |   percentile_ORT_EXP_SUB |   percentile_ORT_EXP_SUB_ESM2 |   delta_score |   delta_rank |   true_label | in_top100   | in_top200   | in_top500   | in_top5percent   |
|:------------------|--------------------:|-------------------------:|-------------------:|------------------------:|-------------------------:|------------------------------:|--------------:|-------------:|-------------:|:------------|:------------|:------------|:-----------------|
| FGRAMPH1_01G18469 |          0.00139232 |                 0.951532 |              11386 |                     496 |                 0.125643 |                      0.961984 |      0.950139 |        10890 |          nan | False       | False       | True        | True             |
| FGRAMPH1_01G18017 |          0.0529729  |                 0.987998 |               7494 |                      81 |                 0.424545 |                      0.993856 |      0.935025 |         7413 |            0 | True        | True        | True        | True             |
| FGRAMPH1_01G24479 |          0.105725   |                 0.953154 |               6011 |                     469 |                 0.538438 |                      0.964058 |      0.847429 |         5542 |            1 | False       | False       | True        | True             |
| FGRAMPH1_01G17083 |          0.21707    |                 0.995509 |               4378 |                      14 |                 0.663851 |                      0.999002 |      0.778438 |         4364 |            0 | True        | True        | True        | True             |
| FGRAMPH1_01G01237 |          0.217452   |                 0.984958 |               4374 |                     109 |                 0.664158 |                      0.991706 |      0.767506 |         4265 |            0 | False       | True        | True        | True             |
| FGRAMPH1_01G26671 |          0.195932   |                 0.952338 |               4586 |                     481 |                 0.647877 |                      0.963136 |      0.756406 |         4105 |            0 | False       | False       | True        | True             |
| FGRAMPH1_01G06067 |          0.213955   |                 0.96266  |               4409 |                     368 |                 0.66147  |                      0.971815 |      0.748704 |         4041 |            0 | False       | False       | True        | True             |
| FGRAMPH1_01G26141 |          0.251464   |                 0.982153 |               4040 |                     145 |                 0.689809 |                      0.988941 |      0.730688 |         3895 |            0 | False       | True        | True        | True             |
| FGRAMPH1_01G22917 |          0.270628   |                 0.980571 |               3894 |                     165 |                 0.701021 |                      0.987405 |      0.709943 |         3729 |            1 | False       | True        | True        | True             |
| FGRAMPH1_01G07249 |          0.291725   |                 0.985656 |               3741 |                     105 |                 0.712772 |                      0.992013 |      0.693931 |         3636 |            1 | False       | True        | True        | True             |
| FGRAMPH1_01G17373 |          0.31436    |                 0.972041 |               3595 |                     265 |                 0.723984 |                      0.979725 |      0.65768  |         3330 |            1 | False       | False       | True        | True             |
| FGRAMPH1_01G23327 |          0.319492   |                 0.976199 |               3578 |                     222 |                 0.72529  |                      0.983027 |      0.656707 |         3356 |            1 | False       | False       | True        | True             |
| FGRAMPH1_01G19487 |          0.310447   |                 0.953469 |               3619 |                     461 |                 0.722141 |                      0.964672 |      0.643022 |         3158 |            1 | False       | False       | True        | True             |
| FGRAMPH1_01G11323 |          0.331441   |                 0.958797 |               3497 |                     412 |                 0.731511 |                      0.968436 |      0.627357 |         3085 |            0 | False       | False       | True        | True             |
| FGRAMPH1_01G20303 |          0.364066   |                 0.983059 |               3302 |                     137 |                 0.746486 |                      0.989555 |      0.618993 |         3165 |            1 | False       | True        | True        | True             |
| FGRAMPH1_01G18641 |          0.352238   |                 0.960521 |               3370 |                     392 |                 0.741264 |                      0.969972 |      0.608283 |         2978 |            0 | False       | False       | True        | True             |
| FGRAMPH1_01G02051 |          0.37337    |                 0.969846 |               3247 |                     290 |                 0.75071  |                      0.977805 |      0.596477 |         2957 |            1 | False       | False       | True        | True             |
| FGRAMPH1_01G04705 |          0.396962   |                 0.987109 |               3121 |                      92 |                 0.760387 |                      0.993011 |      0.590147 |         3029 |            0 | True        | True        | True        | True             |
| FGRAMPH1_01G03115 |          0.379779   |                 0.960869 |               3214 |                     389 |                 0.753245 |                      0.970202 |      0.58109  |         2825 |            1 | False       | False       | True        | True             |
| FGRAMPH1_01G06775 |          0.392782   |                 0.9683   |               3144 |                     310 |                 0.758621 |                      0.976269 |      0.575518 |         2834 |          nan | False       | False       | True        | True             |

## Top-K overlap

|   top_k |   shared_candidates |   ESM2_unique_candidates |   baseline_only_candidates |   jaccard |
|--------:|--------------------:|-------------------------:|---------------------------:|----------:|
|     100 |                  22 |                       78 |                         78 |  0.123596 |
|     200 |                  75 |                      125 |                        125 |  0.230769 |
|     500 |                 237 |                      263 |                        263 |  0.310616 |

## ESM2-rescued candidates

Definition: not in Top500 baseline, in Top500 ESM2, and rank_gain >= 500 or percentile_gain >= 0.2.
Count: 189

## Biological profile

| comparison_group       | metric                     |           value | status    |
|:-----------------------|:---------------------------|----------------:|:----------|
| Shared_top500          | n_genes                    |   237           | available |
| Shared_top500          | known_essential_fraction   |     0.921659    | available |
| Shared_top500          | mean_ppi_degree            |   664.143       | available |
| Shared_top500          | median_ppi_degree          |   628           | available |
| Shared_top500          | mean_orthology_support     |     0.847437    | available |
| Shared_top500          | mean_expression_percentile |     0.499792    | available |
| Shared_top500          | mean_subcellular_signal    |     0.752385    | available |
| Shared_top500          | mean_rank_gain             |    -6.173       | available |
| Shared_top500          | mean_delta_score           |    -0.00160999  | available |
| ESM2_unique_top500     | n_genes                    |   263           | available |
| ESM2_unique_top500     | known_essential_fraction   |     0.764706    | available |
| ESM2_unique_top500     | mean_ppi_degree            |   281.665       | available |
| ESM2_unique_top500     | median_ppi_degree          |   258           | available |
| ESM2_unique_top500     | mean_orthology_support     |     0.31946     | available |
| ESM2_unique_top500     | mean_expression_percentile |     0.352424    | available |
| ESM2_unique_top500     | mean_subcellular_signal    |     0.705146    | available |
| ESM2_unique_top500     | mean_rank_gain             |  1175.64        | available |
| ESM2_unique_top500     | mean_delta_score           |     0.0890496   | available |
| Baseline_unique_top500 | n_genes                    |   263           | available |
| Baseline_unique_top500 | known_essential_fraction   |     0.513514    | available |
| Baseline_unique_top500 | mean_ppi_degree            |   619.24        | available |
| Baseline_unique_top500 | median_ppi_degree          |   534           | available |
| Baseline_unique_top500 | mean_orthology_support     |     0.45103     | available |
| Baseline_unique_top500 | mean_expression_percentile |     0.457713    | available |
| Baseline_unique_top500 | mean_subcellular_signal    |     0.799672    | available |
| Baseline_unique_top500 | mean_rank_gain             | -1240.46        | available |
| Baseline_unique_top500 | mean_delta_score           |    -0.0957561   | available |
| Neither                | n_genes                    | 12259           | available |
| Neither                | known_essential_fraction   |     0.0532424   | available |
| Neither                | mean_ppi_degree            |    85.9574      | available |
| Neither                | median_ppi_degree          |    24           | available |
| Neither                | mean_orthology_support     |     0.223678    | available |
| Neither                | mean_expression_percentile |     0.504118    | available |
| Neither                | mean_subcellular_signal    |     0.73804     | available |
| Neither                | mean_rank_gain             |     1.50999     | available |
| Neither                | mean_delta_score           |     0.000175005 | available |

## Incomplete analyses

- GO_enrichment: not_completed; No local F. graminearum GO annotation/background table was found in the inspected project paths.
- KEGG_enrichment: not_completed; No local F. graminearum KEGG annotation/background table was found in the inspected project paths.

## Draft manuscript figure legend

Figure 5. ESM2 expands the essentiality prediction manifold and identifies ESM2-specific essential gene candidates. (A) GraphSAGE penultimate hidden-embedding UMAP from the ORT+EXP+SUB+ESM2 model highlights predicted essential genes and separates shared versus ESM2-specific Top500 essential predictions by final score-based ranking. (B) Matched input- and hidden-representation UMAPs show how rescued genes reorganize across representation levels. (C) Grouped feature ablation shows that ESM2 provides an independent contribution beyond ORT, EXP, and SUB feature groups. (D) Genome-wide Top-K overlap compares candidates prioritized by baseline GraphSAGE and GraphSAGE+ESM2 across seeds. (E) Rank-shift analysis identifies ESM2-rescued candidates that enter Top500 only after ESM2 and show large rank or percentile gains. (F) Biological profiling summarizes graph, orthology, expression, localization, and label-support characteristics of shared, ESM2-unique, baseline-unique, and non-top500 genes.
