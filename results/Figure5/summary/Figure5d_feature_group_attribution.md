# Figure5d Feature-Group Attribution

## Method

- Primary attribution method: fixed-model leave-one-group-out inference ablation on the saved `ORT_EXP_SUB_ESM2` GraphSAGE checkpoint.
- Masking rule: `zero_out_selected_standardized_columns`.
- Masking rule detail: The saved frozen-protocol bundle is already z-scored by the training split; masking sets the selected standardized columns to 0.0.
- Attribution score definition: `contribution(group) = probability_full - probability_drop_group`.
- Stored prediction quantity: probabilities from `sigmoid(logit)`.
- No model retraining was performed for the ablated groups.
- The additional graph-degree scalar present in the frozen bundle was detected and held constant; Figure5d attributes only the four requested groups `ORT / EXP / SUB / ESM2`.

## Gene subsets used

- `stable_rescued`: essential genes with `stable_rescued_ge2 = True` (Figure5a consensus threshold >= 2 seeds).
- `always_correct_TP`: essential genes with `TP_stable == n_seeds_observed`.
- `persistent_FN`: essential genes with `FN_persistent == n_seeds_observed`.
- `corrected_FP`: non-essential genes with `FP_to_TN_corrected == n_seeds_observed`.

## Main findings

- Fusarium global importance: the largest mean `ΔAUPRC` drop was observed for `ESM2` (`0.1928`), with mean `ΔMCC = 0.2002`.
- Fusarium ESM2 dependence: stable rescued genes had `not higher` mean ESM2 contribution than always-correct TPs (`stable_rescued = 0.2357`, `always_correct_TP = 0.4124`).
- Fusarium `stable_rescued` heatmap pattern: highest mean contribution came from `ESM2` (`0.2357`).
- Fusarium `always_correct_TP` heatmap pattern: highest mean contribution came from `ESM2` (`0.4124`).
- Fusarium `persistent_FN` heatmap pattern: highest mean contribution came from `ORT` (`-0.0128`).
- Fusarium `corrected_FP` heatmap pattern: highest mean contribution came from `SUB` (`-0.0044`).

## Limitations

- Figure5d is group-level attribution only; it does not identify individual raw dimensions, residues, or causal sequence motifs.
- Zeroing standardized columns is an inference-time perturbation on the frozen normalized feature matrix, not a retrained drop-group model.
- Contribution values can be negative when masking a group increases predicted probability or improves a metric for a specific gene/run.

## Outputs

- Feature-group manifest: `results/Figure5/data/Figure5d_feature_group_manifest.tsv`.
- Per-gene attribution table: `results/Figure5/data/Figure5d_group_ablation_per_gene.tsv`.
- Global summary table: `results/Figure5/tables/Figure5d_group_ablation_global_summary.tsv`.
- Gene-set summary table: `results/Figure5/tables/Figure5d_group_ablation_by_gene_set.tsv`.
- Statistical comparison table: `results/Figure5/tables/Figure5d_group_ablation_stats.tsv`.
- Global importance plot: `results/Figure5/plots/Figure5d_feature_group_global_importance.pdf` / `results/Figure5/plots/Figure5d_feature_group_global_importance.png`.
- Rescued-vs-other plot: `results/Figure5/plots/Figure5d_feature_group_rescued_vs_other.pdf` / `results/Figure5/plots/Figure5d_feature_group_rescued_vs_other.png`.
- Gene-set heatmap: `results/Figure5/plots/Figure5d_feature_group_gene_set_heatmap.pdf` / `results/Figure5/plots/Figure5d_feature_group_gene_set_heatmap.png`.
