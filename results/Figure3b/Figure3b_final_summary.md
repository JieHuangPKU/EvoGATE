# Figure3b Final Summary

Publication-facing benchmark summary with fixed split, fixed seeds, and mean ± std test metrics.

| Target       | Model     | Feature_Setting   |   Runs | AUROC         | AUPRC         | MCC           | F1            | Precision   | Recall   | Specificity   | Species      | Regime   |   ESM2_Dim | Label_Regime   | Split_Version                                | Seed_List                |
|:-------------|:----------|:------------------|-------:|:--------------|:--------------|:--------------|:--------------|:------------|:---------|:--------------|:-------------|:---------|-----------:|:---------------|:---------------------------------------------|:-------------------------|
| fgraminearum | GraphSAGE | ORT_EXP_SUB_ESM2  |      5 | 0.918 ± 0.002 | 0.550 ± 0.023 | 0.502 ± 0.005 | 0.540 ± 0.010 | NA ± NA     | NA ± NA  | NA ± NA       | fgraminearum | newlabel |       1280 | newlabel       | frozen_protocol_v1_seed20260409_test20_val10 | 1029,1030,1031,1032,1033 |
| fgraminearum | GraphSAGE | ORT_EXP_SUB_ESM2  |      5 | 0.910 ± 0.001 | 0.552 ± 0.003 | 0.486 ± 0.015 | 0.524 ± 0.020 | NA ± NA     | NA ± NA  | NA ± NA       | fgraminearum | newlabel |        160 | newlabel       | frozen_protocol_v1_seed20260409_test20_val10 | 1029,1030,1031,1032,1033 |
| fgraminearum | GraphSAGE | ORT_EXP_SUB_ESM2  |      5 | 0.912 ± 0.001 | 0.549 ± 0.017 | 0.490 ± 0.012 | 0.526 ± 0.011 | NA ± NA     | NA ± NA  | NA ± NA       | fgraminearum | newlabel |        320 | newlabel       | frozen_protocol_v1_seed20260409_test20_val10 | 1029,1030,1031,1032,1033 |
| fgraminearum | GraphSAGE | ORT_EXP_SUB_ESM2  |      5 | 0.914 ± 0.003 | 0.534 ± 0.028 | 0.497 ± 0.017 | 0.526 ± 0.024 | NA ± NA     | NA ± NA  | NA ± NA       | fgraminearum | newlabel |        640 | newlabel       | frozen_protocol_v1_seed20260409_test20_val10 | 1029,1030,1031,1032,1033 |
