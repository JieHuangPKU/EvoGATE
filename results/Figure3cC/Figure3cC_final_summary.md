# Figure3cC Final Summary

Publication-facing benchmark summary with fixed split, fixed seeds, and mean ± std test metrics.

| Target                | Model     | Feature_Setting                 |   Runs | AUROC         | AUPRC         | MCC           | F1            | Precision     | Recall        | Specificity   | Species      | Regime   |   ESM2_Dim | Label_Regime   | Split_Version                                | Seed_List                |
|:----------------------|:----------|:--------------------------------|-------:|:--------------|:--------------|:--------------|:--------------|:--------------|:--------------|:--------------|:-------------|:---------|-----------:|:---------------|:---------------------------------------------|:-------------------------|
| fgraminearum_newlabel | GraphSAGE | ORT_EXP_SUB_ESM2_OLD_GATED_WBCE |      5 | 0.916 ± 0.001 | 0.580 ± 0.006 | 0.505 ± 0.014 | 0.549 ± 0.012 | 0.469 ± 0.015 | 0.665 ± 0.037 | 0.924 ± 0.008 | fgraminearum | newlabel |       1280 | newlabel       | frozen_protocol_v1_seed20260409_test20_val10 | 1029,1030,1031,1032,1033 |
