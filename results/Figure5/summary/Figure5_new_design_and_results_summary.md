# Figure5 new design and results summary

Theme: "ESM2 reshapes genome-wide candidate prioritization and rescues ORT/EXP/SUB blind-spot genes".

Protocol: `fgraminearum_newlabel`; model: `GraphSAGE`; seeds: 1029, 1030, 1031, 1032, 1033.

## Panel inputs and outputs
### Panel 5A
Inputs:
- `results/Figure5/data/Figure5a_hidden_umap_error_transition_fgraminearum_seed1029_coords.tsv`
Outputs:
- `results/Figure5/plots/Figure5_new_A_hidden_space_rescue_umap.pdf`
- `results/Figure5/plots/Figure5_new_A_hidden_space_rescue_umap.png`
- `results/Figure5/data/Figure5_new_A_hidden_space_rescue_umap_plot_data.tsv`

### Panel 5B
Inputs:
- `results/Figure5/data/Figure5c_input_vs_hidden_compare_fgraminearum_seed1029_coords.tsv`
- `results/Figure5/data/Figure5a_hidden_umap_error_transition_fgraminearum_seed1029_coords.tsv`
Outputs:
- `results/Figure5/plots/Figure5_new_B_input_vs_hidden_umap.pdf`
- `results/Figure5/plots/Figure5_new_B_input_vs_hidden_umap.png`
- `results/Figure5/data/Figure5_new_B_input_vs_hidden_umap_plot_data.tsv`

### Panel 5C
Inputs:
- `results/Figure5/tables/Figure5d_group_ablation_global_summary.tsv`
- `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB/run_1029/feature_schema.tsv`
- `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1029/feature_schema.tsv`
Outputs:
- `results/Figure5/tables/Figure5_new_C_grouped_feature_importance.tsv`
- `results/Figure5/plots/Figure5_new_C_grouped_feature_importance.pdf`
- `results/Figure5/plots/Figure5_new_C_grouped_feature_importance.png`

### Panel 5D
Inputs:
- `results/Figure5/data/Figure5_new_D_per_seed_rank_percentiles.tsv`
Outputs:
- `results/Figure5/tables/Figure5_new_candidate_rank_table.tsv`
- `results/Figure5/tables/Figure5_new_D_topk_overlap_summary.tsv`
- `results/Figure5/plots/Figure5_new_D_topk_candidate_overlap.pdf`
- `results/Figure5/plots/Figure5_new_D_topk_candidate_overlap.png`

### Panel 5E
Inputs:
- `results/Figure5/tables/Figure5_new_candidate_rank_table.tsv`
Outputs:
- `results/Figure5/tables/Figure5_new_E_esm2_rescued_candidates.tsv`
- `results/Figure5/tables/Figure5_new_E_rank_shift_summary.tsv`
- `results/Figure5/plots/Figure5_new_E_rank_shift_rescue_plot.pdf`
- `results/Figure5/plots/Figure5_new_E_rank_shift_rescue_plot.png`

### Panel 5F
Inputs:
- `results/Figure5/tables/Figure5_new_candidate_rank_table.tsv`
- `results/Figure5/tables/Figure5_new_E_esm2_rescued_candidates.tsv`
Outputs:
- `results/Figure5/tables/Figure5_new_F_biological_profile_summary.tsv`
- `results/Figure5/tables/Figure5_new_F_enrichment_results.tsv`
- `results/Figure5/plots/Figure5_new_F_biological_profile_heatmap.pdf`
- `results/Figure5/plots/Figure5_new_F_biological_profile_heatmap.png`
- `results/Figure5/data/Figure5_new_F_biological_profile_gene_level.tsv`

## Top-K overlap

|   top_k |   shared_candidates |   ESM2_unique_candidates |   baseline_only_candidates |   jaccard |
|--------:|--------------------:|-------------------------:|---------------------------:|----------:|
|     100 |                  22 |                       78 |                         78 |  0.123596 |
|     200 |                  75 |                      125 |                        125 |  0.230769 |
|     500 |                 237 |                      263 |                        263 |  0.310616 |

## ESM2-rescued candidates

Definition: not in Top500 baseline, in Top500 ESM2, and rank_gain >= 500 or percentile_gain >= 0.2.
Count: 189

## Top 20 genes by rank_gain

| gene_id           |   score_baseline |   score_ESM2 |   rank_baseline |   rank_ESM2 |   delta_score |   rank_gain |
|:------------------|-----------------:|-------------:|----------------:|------------:|--------------:|------------:|
| FGRAMPH1_01G28139 |        0.086208  |     0.871141 |           11945 |        1670 |      0.784933 |       10275 |
| FGRAMPH1_01G26189 |        0.0816925 |     0.865474 |           11997 |        1747 |      0.783781 |       10250 |
| FGRAMPH1_01G19309 |        0.0498387 |     0.812993 |           12416 |        2433 |      0.763155 |        9983 |
| FGRAMPH1_01G03751 |        0.0742436 |     0.830779 |           12094 |        2211 |      0.756535 |        9883 |
| FGRAMPH1_01G26169 |        0.0476117 |     0.800476 |           12446 |        2590 |      0.752864 |        9856 |
| FGRAMPH1_01G17351 |        0.0400706 |     0.788834 |           12531 |        2728 |      0.748764 |        9803 |
| FGRAMPH1_01G06931 |        0.0479496 |     0.792659 |           12442 |        2690 |      0.744709 |        9752 |
| FGRAMPH1_01G01707 |        0.11522   |     0.858209 |           11561 |        1837 |      0.742989 |        9724 |
| FGRAMPH1_01G01479 |        0.0777914 |     0.818476 |           12048 |        2366 |      0.740685 |        9682 |
| FGRAMPH1_01G17629 |        0.0248195 |     0.76595  |           12704 |        3030 |      0.74113  |        9674 |
| FGRAMPH1_01G11317 |        0.0516357 |     0.775073 |           12387 |        2889 |      0.723437 |        9498 |
| FGRAMPH1_01G00605 |        0.0436339 |     0.766764 |           12494 |        3021 |      0.72313  |        9473 |
| FGRAMPH1_01G18551 |        0.098372  |     0.822347 |           11777 |        2319 |      0.723975 |        9458 |
| FGRAMPH1_01G23731 |        0.123806  |     0.843373 |           11457 |        2020 |      0.719567 |        9437 |
| FGRAMPH1_01G19041 |        0.0151282 |     0.738719 |           12819 |        3387 |      0.723591 |        9432 |
| FGRAMPH1_01G14359 |        0.0522961 |     0.765351 |           12373 |        3037 |      0.713055 |        9336 |
| FGRAMPH1_01G19743 |        0.0400706 |     0.752634 |           12532 |        3206 |      0.712563 |        9326 |
| FGRAMPH1_01G14285 |        0.144678  |     0.854585 |           11208 |        1882 |      0.709906 |        9326 |
| FGRAMPH1_01G27429 |        0.143281  |     0.847689 |           11236 |        1964 |      0.704408 |        9272 |
| FGRAMPH1_01G20539 |        0.0431885 |     0.746997 |           12498 |        3275 |      0.703809 |        9223 |

## Biological profile

| comparison_group       | metric                     |         value | status    |
|:-----------------------|:---------------------------|--------------:|:----------|
| Shared_top500          | n_genes                    |   237         | available |
| Shared_top500          | known_essential_fraction   |     0.921659  | available |
| Shared_top500          | mean_ppi_degree            |   664.143     | available |
| Shared_top500          | median_ppi_degree          |   628         | available |
| Shared_top500          | mean_orthology_support     |     0.847437  | available |
| Shared_top500          | mean_expression_percentile |     0.499792  | available |
| Shared_top500          | mean_subcellular_signal    |     0.752385  | available |
| Shared_top500          | mean_rank_gain             |    -6.16034   | available |
| Shared_top500          | mean_delta_score           |    -0.0016099 | available |
| ESM2_unique_top500     | n_genes                    |   263         | available |
| ESM2_unique_top500     | known_essential_fraction   |     0.764706  | available |
| ESM2_unique_top500     | mean_ppi_degree            |   281.665     | available |
| ESM2_unique_top500     | median_ppi_degree          |   258         | available |
| ESM2_unique_top500     | mean_orthology_support     |     0.31946   | available |
| ESM2_unique_top500     | mean_expression_percentile |     0.352424  | available |
| ESM2_unique_top500     | mean_subcellular_signal    |     0.705146  | available |
| ESM2_unique_top500     | mean_rank_gain             |  1175.6       | available |
| ESM2_unique_top500     | mean_delta_score           |     0.0890387 | available |
| Baseline_unique_top500 | n_genes                    |   263         | available |
| Baseline_unique_top500 | known_essential_fraction   |     0.513514  | available |
| Baseline_unique_top500 | mean_ppi_degree            |   619.24      | available |
| Baseline_unique_top500 | median_ppi_degree          |   534         | available |
| Baseline_unique_top500 | mean_orthology_support     |     0.45103   | available |
| Baseline_unique_top500 | mean_expression_percentile |     0.457713  | available |
| Baseline_unique_top500 | mean_subcellular_signal    |     0.799672  | available |
| Baseline_unique_top500 | mean_rank_gain             | -1240.45      | available |
| Baseline_unique_top500 | mean_delta_score           |    -0.0957488 | available |

## Incomplete analyses

- GO_enrichment: not_completed; No local F. graminearum GO annotation/background table was found in the inspected project paths.
- KEGG_enrichment: not_completed; No local F. graminearum KEGG annotation/background table was found in the inspected project paths.
- Baseline grouped SHAP/feature-ablation values: not completed because no baseline ORT_EXP_SUB grouped attribution file was found; the table records these rows as unavailable rather than imputing values.

## Draft manuscript figure legend

Figure 5. ESM2 reshapes genome-wide candidate prioritization and rescues genes missed by orthology, expression, and subcellular localization features. (A) Hidden-space UMAP of F. graminearum test genes highlights FN-to-TP rescued genes after adding ESM2. (B) Matched input- and hidden-representation UMAPs show that rescued genes are reorganized in the learned hidden space. (C) Grouped feature ablation shows that ESM2 provides an independent contribution beyond ORT, EXP, and SUB feature groups. (D) Genome-wide Top-K overlap compares candidates prioritized by baseline GraphSAGE and GraphSAGE+ESM2 across seeds. (E) Rank-shift analysis identifies ESM2-rescued candidates that enter Top500 only after ESM2 and show large rank or percentile gains. (F) Biological profiling summarizes graph, orthology, expression, localization, and label-support characteristics of shared, ESM2-unique, baseline-unique, and non-top500 genes.
