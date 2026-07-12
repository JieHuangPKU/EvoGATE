# EvoGATE 项目状态

_基于证据的实现与验证矩阵。_

---

## 状态词

本文只使用：**Implemented**、**Validated**、**Partially implemented**、**Partially validated**、**Planned**、**Unknown**、**Blocked** 和 **Historical**。

## 模块状态矩阵

| 模块 | 科学目的 | 实现状态 | 验证状态 | 规范入口 | 预期输出 | 当前 blocker | 证据路径 |
|---|---|---|---|---|---|---|---|
| Yeast essentiality transfer artifact | 提供进化监督 | Partially implemented | Partially validated | Unknown upstream generator | `ph1_yeast_essential_ortholog_labels.tsv` | Confidence generator 缺失 | `data/derived_labels/` |
| PH-1 canonical bridge | 解析 protein 与 legacy ID | Implemented | Validated | `src.data.build_fgraminearum_newlabel_bridge` | Bridge 与 unresolved audit | External configured source 影响 rerun | `data/processed/essential_gene/fgraminearum/bridge/` |
| Lethal evidence preparation | 选择直接 lethal support | Implemented | Validated | `src.data.prepare_fgraminearum_label_materialization_sources` | `lethal_positive_gene_list.tsv` | 部分 configured input 是 historical path | `data/interim/protocol_refactor/fgraminearum_label_materialization/` |
| Newlabel materialization | 定义当前监督信号 | Implemented | Validated | `src.data.materialize_fgraminearum_label_regimes` | 1,097 positives；10,868 negatives | 完整上游重建不完整 | `data/processed/essential_gene/fgraminearum/newlabel/` |
| Oldlabel replay | 保留历史比较 | Implemented | Historical | 同一 materialization module | Oldlabel table | 非 mainline regime | `data/processed/essential_gene/fgraminearum/oldlabel/` |
| Frozen labels and splits | 防止 model-specific resampling | Implemented | Validated | `src.data.freeze_unified_protocol` | Frozen label/split manifest | 会写 existing frozen result | `results/frozen_protocol/` |
| Orthology features | 表示进化背景 | Implemented | Partially validated | `src.data.build_inparanoid_ortholog_matrix` | `orthologs.csv` | Builder 有历史 machine path | `data/processed/OR/` |
| Expression features | 表示 expression profile | Implemented | Partially validated | `src.data.build_expression_profile_csv` | `profile.csv` | Raw/source build portability 不完整 | `data/processed/EXP/` |
| Localization features | 表示 cellular context | Implemented | Partially validated | `src.data.build_subloc_csv_from_compartments` | `subloc.csv` | Raw/source build portability 不完整 | `data/processed/LC/` |
| ESM2 embeddings | 表示 protein sequence | Implemented | Validated | `src.features.extract_esm2_pooled` | `esm2_pooled.pt` | Local model/environment path | `data/processed/ESM2/` |
| STRING graph | 定义默认 PPI topology | Implemented | Validated | `src.data.frozen_protocol_loader` | Filtered edge table/index | 完整 raw rebuild 不可移植 | `data/processed/PPI/` |
| eFG graph comparison | 检验 graph-source robustness | Implemented | Partially validated | `workflow/Figure4_graph_robustness.smk` | Figure4 summary | 大型 rerun 需批准 | `results/Figure4/` |
| Unified frozen loader | 对齐 label、graph 与 feature | Implemented | Validated | `src.data.frozen_protocol_loader` | In-memory benchmark bundle | Environment 未锁定 | `src/data/frozen_protocol_loader.py` |
| Classical baselines | 提供 non-graph comparison | Implemented | Partially validated | `src.train.run_frozen_protocol_model` | Per-run metric/prediction | `outputs/` 缺失 | `results/Figure1/summary/` |
| Topology baselines | 测量 network-only information | Implemented | Partially validated | 同一 runner | node2vec/DC/CC metrics | Backend/environment reproducibility | `results/Figure1/summary/` |
| GNN families | 比较 message-passing model | Implemented | Partially validated | 同一 runner | GAT/GCN/GIN/GraphSAGE runs | Figure2 artifact 冲突 | `results/Figure2a/`、`results/Figure2b/` |
| ESM2 comparison | 量化 sequence representation contribution | Implemented | Partially validated | `workflow/Figure3a_fusarium_graphsage_esm2_comparison.smk` | Figure3a summary | Per-run output 缺失 | `results/Figure3a/` |
| ESM2 dimension ablation | 比较截断 embedding dimension | Implemented | Partially validated | `workflow/Figure3b_fusarium_graphsage_esm2_dim_ablation.smk` | Figure3b summary | Wrapper 不可移植 | `results/Figure3b/` |
| Fusion ablation | 比较 concatenation 与 gate | Implemented | Partially validated | Figure3c workflows | Figure3c summary | Metric 方向混合；claim 未解决 | `results/Figure3c*/` |
| Label-scarcity benchmark | 测试减少标签后的行为 | Implemented | Partially validated | `workflow/label_scarcity_benchmark.smk` | Scarcity metric/plot | 旧叙述与排名冲突 | `results/Figure2_label_scarcity/` |
| Graph robustness | 测试 PPI threshold/source | Implemented | Partially validated | `workflow/Figure4_graph_robustness.smk` | Figure4 table/plot | Statistical claim scope 不完整 | `results/Figure4/` |
| Representation analysis | 检查 hidden/input space | Implemented | Partially validated | `workflow/Figure5_representation_mechanism.smk` | UMAP/summary artifact | 重建依赖缺失 run output | `results/Figure5/` |
| Group zero-out analysis | 估计 feature-group sensitivity | Implemented | Partially validated | `workflow/Figure5d_feature_group_attribution.smk` | Group perturbation summary | 不是 causal attribution 或 SHAP | `results/Figure5/` |
| Candidate prioritization | 排序实验候选 | Partially implemented | Partially validated | `src.eval.build_figure5_candidate_prioritization` | Candidate rank table | 缺少 `outputs/`；无 wet-lab validation | `results/Figure5_new_candidate_prioritization/` |
| RNA target discovery | 评价 RNA target suitability | Planned | Unknown | None | None | 无实现 | None |
| Off-target filtering | 移除 host/non-target match | Planned | Unknown | None | None | 无实现 | None |
| dsRNA design | 设计 silencing sequence | Planned | Unknown | None | None | 无实现 | None |

## 仓库成熟度

EvoGATE 是具有 validated supervision core 和大量实验 artifact 的科研仓库。发布级 portability 与 reproducibility 为 **Blocked**，不得描述为 production-ready software release 或 validated RNA-target platform。

Frozen contract 保存在 `results/frozen_protocol/`，已验证的 ESM2 artifact 保存在 `data/processed/ESM2/`。
