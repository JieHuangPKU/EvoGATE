# EvoGATE

_面向植物病原真菌的进化感知必需基因预测与优先级排序框架。_

---

> EvoGATE 当前是科研仓库，尚不是生产就绪的软件发布。

## 项目摘要

EvoGATE 从表型与进化证据重建可审计的必需基因监督标签，将其与蛋白互作网络、生物学特征和 ESM2 embedding 结合，并在禾谷镰刀菌（*Fusarium graminearum*）中优先排序候选必需基因。

## 当前科学范围

当前范围是必需基因预测和计算候选优选。主要目标为 *F. graminearum* PH-1；人、秀丽隐杆线虫、酿酒酵母和黑腹果蝇作为 benchmark 物种保留。

### EvoGATE 当前能做什么

- **Validated**：物化 `fgraminearum_newlabel`，包含 1,097 个正样本和 10,868 个负样本
- **Validated**：冻结分层的 70/10/20 train/validation/test 划分
- **Implemented**：整合 PPI、orthology（ORT）、expression（EXP）、subcellular localization（SUB）和 ESM2 表征
- **Implemented**：比较 MLP、RF、SVM、NB、node2vec、GAT、GCN、GIN、GraphSAGE 和网络启发式方法
- **Partially implemented**：进行全基因组候选排序并识别与 ESM2 相关的排名变化

### EvoGATE 尚不能做什么

- **Planned**：RNA 靶点发现
- **Planned**：宿主与非靶标脱靶过滤
- **Planned**：dsRNA 或 siRNA 设计
- **Blocked**：在当前工作区进行发布级端到端复现

候选预测是计算假设，不是经实验验证的必需基因或有效 RNA 靶点。

## 核心科学贡献

首要贡献是 evolution-aware label reconstruction。项目在明确的 PH-1 ID bridge 后，将 PHI 支持的 lethal evidence 与高置信酵母必需直系同源迁移结合。GraphSAGE 是主要建模载体，不是唯一创新。

## 仓库状态警告

仓库由历史名称 ProGATE 和 ProGATE_v2 迁移而来。EPGAT 是方法学前身和遗留实现来源；Bingo 是外部比较项目及部分标准标签的历史来源。二者都不是 EvoGATE 的别名或内部组件。

若干 wrapper 仍包含历史绝对路径，仓库内缺少 `outputs/`，`.git/` 为空，缺少环境锁文件，部分源码或上游生成步骤也缺失。详见[已知不一致](docs/INCONSISTENCIES.md)与[可复现性](docs/REPRODUCIBILITY.zh-CN.md)。

这些迁移残留的状态为 **Historical**，不代表当前项目身份。

## 规范工作流入口

| 目的 | 规范入口 | 状态 |
|---|---|---|
| 标签物化 | `workflow/fgraminearum_label_materialization.smk` | Partially reproducible |
| Frozen benchmark | `workflow/frozen_protocol_benchmark.smk` | Implemented；大型任务 |
| 单个 benchmark 任务 | `python -m src.train.run_frozen_protocol_model` | Implemented |
| ESM2 cache 准备 | `workflow/prepare_esm2_cache.smk` | Implemented；大型任务 |
| 候选优选 | `python -m src.eval.build_figure5_candidate_prioritization` | 因缺少 `outputs/` 输入而 Blocked |

硬编码 `/home/jiehuang/software/fungi/ProGATE_v2` 的历史 Shell wrapper 不可移植，不是推荐入口。

## Frozen protocol

| 设置 | 值 |
|---|---|
| Protocol | `frozen_protocol_v1` |
| Split | 70% train / 10% validation / 20% test |
| Split seed | `20260409` |
| Training seeds | `1029`、`1030`、`1031`、`1032`、`1033` |
| 主要指标 | AUPRC、MCC、AUROC |
| 阈值策略 | 标准运行固定 `0.5`；任何调优阈值都必须由 validation 数据决定 |

## 仓库地图

| 路径 | 职责 |
|---|---|
| `configs/` | Frozen protocol 和实验配置 |
| `data/` | Manifest、派生标签、中间 provenance 和处理后模态 |
| `docs/` | 项目知识库、迁移审计和科学文档 |
| `scripts/` | 薄执行 wrapper、历史 wrapper 和 Figure builder |
| `src/` | 数据、特征、图、模型、训练、评价和分析模块 |
| `workflow/` | Snakemake workflow |
| `results/` | Frozen manifest、摘要、Figure 和历史实验 artifact |

## 文档索引

- [项目概览](docs/PROJECT_OVERVIEW.zh-CN.md)
- [科学设计](docs/SCIENTIFIC_DESIGN.zh-CN.md)
- [架构](docs/ARCHITECTURE.zh-CN.md)
- [数据流](docs/DATA_FLOW.zh-CN.md)
- [运行手册](docs/RUNBOOK.zh-CN.md)
- [项目状态](docs/PROJECT_STATUS.zh-CN.md)
- [可复现性](docs/REPRODUCIBILITY.zh-CN.md)
- [数据字典](docs/DATA_DICTIONARY.zh-CN.md)
- [项目历史](docs/PROJECT_HISTORY.zh-CN.md)
- [迁移指南](docs/MIGRATION_GUIDE.zh-CN.md)
- [论文映射](docs/MANUSCRIPT_MAPPING.zh-CN.md)
- [术语表](docs/GLOSSARY.zh-CN.md)
- [已知不一致](docs/INCONSISTENCIES.md)
- [英文 README](README.md)

## 可复现性、引用和许可状态

- **可复现性**：Partially reproducible；尚未达到发布级
- **引用**：尚未建立规范 citation file
- **许可**：尚未建立仓库 license
- **项目成熟度**：具有已验证核心 artifact，但仍存在移植与 provenance 阻断项的科研仓库
