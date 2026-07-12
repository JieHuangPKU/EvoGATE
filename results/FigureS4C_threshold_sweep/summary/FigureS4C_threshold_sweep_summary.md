# Figure S4C threshold sweep summary

Decision metrics such as MCC and F1 vary with the classification threshold applied to the model output score, whereas AUROC and AUPRC remain threshold-invariant ranking metrics.

## Prediction files used

- `Concat`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1029/predictions.tsv`
- `Concat`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1030/predictions.tsv`
- `Concat`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1031/predictions.tsv`
- `Concat`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1032/predictions.tsv`
- `Concat`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2/run_1033/predictions.tsv`
- `Gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED/run_1029/predictions.tsv`
- `Gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED/run_1030/predictions.tsv`
- `Gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED/run_1031/predictions.tsv`
- `Gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED/run_1032/predictions.tsv`
- `Gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED/run_1033/predictions.tsv`
- `Gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cC/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_OLD_GATED_WBCE/run_1029/predictions.tsv`
- `Gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cC/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_OLD_GATED_WBCE/run_1030/predictions.tsv`
- `Gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cC/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_OLD_GATED_WBCE/run_1031/predictions.tsv`
- `Gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cC/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_OLD_GATED_WBCE/run_1032/predictions.tsv`
- `Gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cC/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_OLD_GATED_WBCE/run_1033/predictions.tsv`
- `Residual gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL/run_1029/predictions.tsv`
- `Residual gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL/run_1030/predictions.tsv`
- `Residual gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL/run_1031/predictions.tsv`
- `Residual gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL/run_1032/predictions.tsv`
- `Residual gated`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL/run_1033/predictions.tsv`
- `Residual gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE/run_1029/predictions.tsv`
- `Residual gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE/run_1030/predictions.tsv`
- `Residual gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE/run_1031/predictions.tsv`
- `Residual gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE/run_1032/predictions.tsv`
- `Residual gated+WBCE`: `/home/jiehuang/software/fungi/ProGATE_v2/outputs/Figure3cB/fgraminearum_newlabel/GraphSAGE/ORT_EXP_SUB_ESM2_GATED_RESIDUAL_WBCE/run_1033/predictions.tsv`

## Split used for threshold optimization

- `validation`

## Split used for threshold sweep curves

- `test`

## Fusion variants included

- `Concat`
- `Gated`
- `Gated+WBCE`
- `Residual gated`
- `Residual gated+WBCE`

## Number of genes

- `validation`: 1197
- `test`: 2393

## Number of seeds/runs

- `Concat`: 5
- `Gated`: 5
- `Gated+WBCE`: 5
- `Residual gated`: 5
- `Residual gated+WBCE`: 5

## Optimal F1 thresholds per fusion method

| fusion_method | threshold | metric_value |
|:--|--:|--:|
| Concat | 0.71 | 0.593128 |
| Gated | 0.83 | 0.586487 |
| Gated+WBCE | 0.61 | 0.586508 |
| Residual gated | 0.76 | 0.572257 |
| Residual gated+WBCE | 0.66 | 0.565271 |

## Optimal MCC thresholds per fusion method

| fusion_method | threshold | metric_value |
|:--|--:|--:|
| Concat | 0.58 | 0.553651 |
| Gated | 0.83 | 0.544142 |
| Gated+WBCE | 0.61 | 0.543205 |
| Residual gated | 0.76 | 0.529547 |
| Residual gated+WBCE | 0.48 | 0.524914 |
