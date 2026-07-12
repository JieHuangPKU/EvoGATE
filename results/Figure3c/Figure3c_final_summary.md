# Figure3c Final Summary

Publication-facing benchmark summary with fixed split, fixed seeds, and mean ± std test metrics.

| Target                | Model     | Feature_Setting        |   Runs | AUROC         | AUPRC         | MCC           | F1            | Precision     | Recall        | Specificity   | Species      | Regime   |   ESM2_Dim | Label_Regime   | Split_Version                                | Seed_List                |
|:----------------------|:----------|:-----------------------|-------:|:--------------|:--------------|:--------------|:--------------|:--------------|:--------------|:--------------|:-------------|:---------|-----------:|:---------------|:---------------------------------------------|:-------------------------|
| fgraminearum_newlabel | GraphSAGE | ORT_EXP_SUB_ESM2       |      5 | 0.918 ± 0.002 | 0.550 ± 0.023 | 0.502 ± 0.005 | 0.540 ± 0.010 | 0.437 ± 0.040 | 0.720 ± 0.063 | 0.904 ± 0.022 | fgraminearum | newlabel |       1280 | newlabel       | frozen_protocol_v1_seed20260409_test20_val10 | 1029,1030,1031,1032,1033 |
| fgraminearum_newlabel | GraphSAGE | ORT_EXP_SUB_ESM2_GATED |      5 | 0.913 ± 0.002 | 0.563 ± 0.013 | 0.490 ± 0.008 | 0.512 ± 0.011 | 0.371 ± 0.014 | 0.827 ± 0.015 | 0.858 ± 0.011 | fgraminearum | newlabel |       1280 | newlabel       | frozen_protocol_v1_seed20260409_test20_val10 | 1029,1030,1031,1032,1033 |
