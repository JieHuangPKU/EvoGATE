# EvoGATE 论文映射

_从科学 claim 到 evidence、code、config、data、output 和 Figure 的追溯。_

---

## Claim-to-artifact 矩阵

| Claim | Evidence | Code | Config | Input data | Output artifact | Figure | Status | Caveat |
|---|---|---|---|---|---|---|---|---|
| Evolution-aware label reconstruction 产生当前 Fusarium supervision regime | Materialized count 与 source manifest | `build_fgraminearum_newlabel_bridge.py`；`prepare_fgraminearum_label_materialization_sources.py`；`materialize_fgraminearum_label_regimes.py` | `fgraminearum_label_materialization.yaml` | Transfer table、evidence mirror、unified/legacy map | `data/processed/essential_gene/fgraminearum/newlabel/` | Figure1B/1C | Validated | 上游 confidence producer 缺失 |
| Newlabel 含 1,097 positives 和 10,868 negatives | `summary.tsv`、frozen label summary | 同一 materializer 与 freeze module | Label/frozen configs | Materialized positive/negative table | Newlabel/frozen summary | Figure1C | Validated | Exclusion set 不是 standalone complete manifest |
| Newlabel 与 oldlabel 的组成和 model behavior 不同 | Old/new comparison table | Figure2a aggregation 与 feature-combo runner | `Figure2a_fusarium_label_compare_graphsage.yaml` | Oldlabel/newlabel frozen manifest | `results/Figure2a/` | Figure2a | Partially validated | Figure2 artifact version 冲突；见 `INCONSISTENCIES.md` |
| 已 benchmark 多种 classical、topology 和 GNN family | Frozen model list 与 Figure1 summary | `run_frozen_protocol_model.py`、aggregator | `frozen_protocol.yaml` | Frozen label/split 与 processed feature | `results/Figure1/summary/` | Figure1 benchmark | Partially validated | Per-run `outputs/` 缺失；缺少完整 paired inference |
| Biological feature group 的贡献不同 | Feature-combination result | Feature-combo runner 与 Figure2b workflow | `Figure2b_fusarium_gnn_feature_ablation.yaml` | ORT、EXP、SUB、graph | `results/Figure2b/` | Figure2b | Partially validated | 与 Figure2a 名义相似 setting 冲突 |
| ESM2 为 Fusarium GraphSAGE 增加信息 | Four-setting comparison summary | Figure3a workflow/aggregation | `Figure3a_fusarium_graphsage_esm2_comparison.yaml` | ESM2 cache、ORT、EXP、SUB、graph | `results/Figure3a/` | Figure3a | Partially validated | 五 seed descriptive evidence；per-run output 缺失 |
| ESM2 dimension truncation 改变性能 | Dimension summary | Figure3b workflow | `Figure3b_fusarium_graphsage_esm2_dim_ablation.yaml` | Full ESM2 cache 与 frozen bundle | `results/Figure3b/` | Figure3b | Partially validated | 截断 embedding dimension；未比较 ESM2 layer/model |
| Gated fusion 相对 concatenation 改变性能 | Figure3c summary | Fusion model 与 Figure3c workflow | Figure3c YAML files | ORT/EXP/SUB/ESM2 frozen bundle | `results/Figure3c*/` | Figure3c | Partially validated | AUPRC 与 threshold metric 方向不同；不可声称全面优越 |
| 已测量 reduced-label 条件下的性能 | Ranking 与 per-run scarcity table | Split builder、workflow、summarizer | Scarcity workflow 使用 base frozen config | Frozen newlabel split 与 feature | `results/Figure2_label_scarcity/` | Label-scarcity Figure | Partially validated | 10% 时 MLP 第一；旧 GraphSAGE-superiority 叙述被否定 |
| PPI threshold/source 影响性能 | Threshold/source comparison table | Figure4 graph adapter、runner、summarizer | `frozen_protocol.yaml` Figure4 section | STRING/eFG graph 与 frozen bundle | `results/Figure4/` | Figure4 | Partially validated | 主要限于 Fusarium；inference scope 有限 |
| ESM2 改变 hidden representation 与 feature-group sensitivity | UMAP、hidden summary、zero-out table | `src/analysis/` Figure5 modules | Figure5 workflow 与 Figure3 runtime config | Prediction、checkpoint、feature schema | `results/Figure5/` | Figure5a-5d | Partially validated | Zero-out 不是 SHAP/GNNExplainer 或 causal attribution |
| EvoGATE 优先排序 candidate essential genes | Candidate rank、top-k overlap、rank-shift table | `build_figure5_candidate_prioritization.py` | CLI/runtime argument | 跨 seed Figure3a prediction | `results/Figure5_new_candidate_prioritization/` | Figure5 | Partially implemented | 缺少 upstream `outputs/`；候选无 wet-lab validation |
| EvoGATE 发现有效 RNA target | None | None | None | None | None | None | Planned | 禁止提出该 claim |

## 建议论文结构

| 章节 | 主要证据 |
|---|---|
| Background and problem | Figure1A 与 label-source summary |
| Scientific hypothesis | Frozen protocol 与预先规定的 claim hierarchy |
| Evolution-aware label reconstruction | Figure1B/1C 与 label materialization artifact |
| Framework | Loader、representation、graph 和 model architecture |
| Benchmark | Figure1、Figure2、Figure4，并披露 conflict |
| ESM2 | Figure3 与 Figure5 group perturbation |
| Interpretability | Figure5 representation 与 perturbation limitation |
| Candidate prioritization | Figure5 rank，作为 computational hypothesis |
| Discussion | Label provenance、imbalance、graph scope、portability 与 missing validation |
| Future RNA target discovery | 仅 Planned |

## Claim 语言规范

代码路径使用 “implemented”；只有已有证据支持明确 contract/result 时使用 “validated”；候选使用 “predicted” 或 “prioritized”；RNA application 使用 “planned”。Candidate ranking 不得写为 “experimentally validated”。受 conflicting artifact 影响的 claim 必须链接 [INCONSISTENCIES.md](INCONSISTENCIES.md)，不得选择有利版本。

