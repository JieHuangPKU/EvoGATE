# Figure5A prediction UMAP audit

- Input file path ORT_EXP_SUB: `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB/run_1029/predictions.tsv`.
- Input file path ORT_EXP_SUB_ESM2: `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1029/predictions.tsv`.
- Using seed/run: `run_1029`.
- UMAP input file: baseline `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB/run_1029/best_model.pt`; ESM2 `outputs/Figure3a/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1029/best_model.pt`.
- UMAP input dimensionality: baseline `13022 x 64`; ESM2 `13022 x 64`.
- Score column name used only for coloring/highlight and Top-K definition: `pred_score`.
- Aligned gene count: `13022`.
- UMAP feature columns: `checkpoint-derived GraphSAGE penultimate hidden embedding`.
- UMAP metric: baseline `cosine`; ESM2 `cosine`.
- K definition for main figure: `Top500` by descending prediction score within each model.
- Additional summaries also include `Top100`, `Top200`, and `Top5%`.
- This figure uses hidden embedding for UMAP coordinates and prediction score only for coloring/highlight.
- Hidden embedding source: `checkpoint_penultimate_hidden` for both panels.
- ESM2-specific essential candidates colored green (`#2CA25F`).
- Left ORT+EXP+SUB panel keeps original y-axis limits `(-8.703361511230469, 16.868297576904297)`.
- Right ORT+EXP+SUB+ESM2 panel display is restricted to `UMAP2 >= 4.0` with y-axis limits `(4.0, 14.316974639892578)`.
- Underlying UMAP coordinates are unchanged unless a filtered plot version is explicitly used.
- No score/rank/delta fallback was used.

## Output files

- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_essential_vs_nonessential.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_essential_vs_nonessential.png`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_esm2_specific.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_esm2_specific.png`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_space_umap_two_panel.pdf`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_space_umap_two_panel.png`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_shift_summary.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5A_ESM2_specific_essential_genes.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5A_shared_essential_genes.tsv`
- `results/Figure5_new_candidate_prioritization/Figure5A_prediction_umap_plot_data.tsv`

## Top 20 ESM2-specific genes by delta_score

| gene_id           |   score_ORT_EXP_SUB |   score_ORT_EXP_SUB_ESM2 |   delta_score |   delta_rank |   rank_ORT_EXP_SUB |   rank_ORT_EXP_SUB_ESM2 |
|:------------------|--------------------:|-------------------------:|--------------:|-------------:|-------------------:|------------------------:|
| FGRAMPH1_01G18469 |          0.00139232 |                 0.951532 |      0.950139 |        10890 |              11386 |                     496 |
| FGRAMPH1_01G18017 |          0.0529729  |                 0.987998 |      0.935025 |         7413 |               7494 |                      81 |
| FGRAMPH1_01G24479 |          0.105725   |                 0.953154 |      0.847429 |         5542 |               6011 |                     469 |
| FGRAMPH1_01G17083 |          0.21707    |                 0.995509 |      0.778438 |         4364 |               4378 |                      14 |
| FGRAMPH1_01G01237 |          0.217452   |                 0.984958 |      0.767506 |         4265 |               4374 |                     109 |
| FGRAMPH1_01G26671 |          0.195932   |                 0.952338 |      0.756406 |         4105 |               4586 |                     481 |
| FGRAMPH1_01G06067 |          0.213955   |                 0.96266  |      0.748704 |         4041 |               4409 |                     368 |
| FGRAMPH1_01G26141 |          0.251464   |                 0.982153 |      0.730688 |         3895 |               4040 |                     145 |
| FGRAMPH1_01G22917 |          0.270628   |                 0.980571 |      0.709943 |         3729 |               3894 |                     165 |
| FGRAMPH1_01G07249 |          0.291725   |                 0.985656 |      0.693931 |         3636 |               3741 |                     105 |
| FGRAMPH1_01G17373 |          0.31436    |                 0.972041 |      0.65768  |         3330 |               3595 |                     265 |
| FGRAMPH1_01G23327 |          0.319492   |                 0.976199 |      0.656707 |         3356 |               3578 |                     222 |
| FGRAMPH1_01G19487 |          0.310447   |                 0.953469 |      0.643022 |         3158 |               3619 |                     461 |
| FGRAMPH1_01G11323 |          0.331441   |                 0.958797 |      0.627357 |         3085 |               3497 |                     412 |
| FGRAMPH1_01G20303 |          0.364066   |                 0.983059 |      0.618993 |         3165 |               3302 |                     137 |
| FGRAMPH1_01G18641 |          0.352238   |                 0.960521 |      0.608283 |         2978 |               3370 |                     392 |
| FGRAMPH1_01G02051 |          0.37337    |                 0.969846 |      0.596477 |         2957 |               3247 |                     290 |
| FGRAMPH1_01G04705 |          0.396962   |                 0.987109 |      0.590147 |         3029 |               3121 |                      92 |
| FGRAMPH1_01G03115 |          0.379779   |                 0.960869 |      0.58109  |         2825 |               3214 |                     389 |
| FGRAMPH1_01G06775 |          0.392782   |                 0.9683   |      0.575518 |         2834 |               3144 |                     310 |
