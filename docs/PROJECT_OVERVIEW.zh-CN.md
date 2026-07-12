# EvoGATE 项目概览

_科学范围、贡献、应用边界与当前成熟度。_

---

## 生物学问题

植物病原真菌的实验必需基因证据稀缺。Lethality、virulence、pathogenicity 和 reduced fitness 相互关联但并不等价。蛋白、转录本、旧注释与 PH-1 canonical gene 之间的 ID 变化进一步增加了标签构建难度。

EvoGATE 先处理监督信号问题，再应用机器学习。当前目标为 *Fusarium graminearum* PH-1，并保留四个模式生物用于比较 benchmark。

## 标签稀缺

主线 Fusarium regime 包含 1,097 个正样本和 10,868 个负样本。尽管规模大于历史 oldlabel，正类仍为少数，因此 accuracy 可能掩盖较差的正类性能；AUPRC 和 MCC 是主要指标。

现有 label-scarcity 实验作为 artifact 已 **Validated**，但旧叙述与实际排名不一致。在保留 10% 标签时，MLP 的 AUPRC 排名第一，GraphSAGE 排名第三，因此不能据此声称 GraphSAGE 全面优越。

## 进化感知监督

当前正样本集合由两类证据并集构成：

1. 77 个满足 protocolized evidence 与 canonical-ID 要求的 PHI-supported lethal genes
2. 1,045 个经过 PH-1 protein-to-gene bridge 的高置信 yeast-essentiality-transfer genes

两部分重叠 25 个基因，得到 1,097 个唯一正样本。负样本是 `weak_positive_confidence == none` 且成功解析的基因，再排除具有 virulence/pathogenicity 证据的基因和全部正样本，最终为 10,868 个。

低/中置信迁移、未解析 mapping 及被生物学规则排除的基因不进入监督标签。当前没有单一、完整的 exclusion manifest，因此 exclusion layer 为 **Partially implemented**。

## 多模态图学习

EvoGATE 组合互补输入：

| 模态 | 生物学作用 | 当前状态 |
|---|---|---|
| PPI | 用于邻域聚合的图结构 | Implemented |
| ORT | 跨物种 orthology presence 与进化背景 | Implemented |
| EXP | Expression profile | Implemented |
| SUB | Subcellular localization indicator | Implemented |
| ESM2 | 蛋白序列表征 | Implemented |

GraphSAGE 是主要模型载体。GAT、GCN、GIN、MLP、RF、SVM、NB、node2vec、degree centrality 和 closeness centrality 提供比较。科学主张关注监督质量和互补证据，而不是 GraphSAGE 的唯一性。

## 候选优选

项目已计算全基因组预测、跨 seed 排名摘要、ESM2 相关排名变化和 feature-group perturbation。因此 candidate prioritization 为 **Partially implemented**。

候选尚未经过 gene deletion、conditional knockdown、infection assay 或 RNA interference 实验验证，必须表述为 predicted 或 prioritized candidates。

## 预期应用

近期应用是在 *F. graminearum* 中优先排序供实验跟进的必需基因假设。RNA target discovery、off-target filtering 和 dsRNA design 是 **Planned** 的下游应用。当前没有模块能够证明 RNA 靶点适用性或 dsRNA 有效性。

## 贡献

### 科学贡献

- 面向真菌 essentiality 的 evolution-aware label reconstruction
- 明确区分 lethal evidence、transfer evidence、negative 和 exclusion
- 在 frozen protocol 下评价网络、组学、进化和 protein-language 表征
- 对候选必需基因进行计算优选

### 工程贡献

- 可审计的 PH-1 ID bridge 和 source manifest
- 多模型共享的 frozen label 与 split contract
- 统一的 multimodal loader 和结构化 per-run output
- 可复用的 Figure 聚合与解释 workflow

## 当前状态

| 领域 | 状态 | 解释 |
|---|---|---|
| 标签重构 | Validated | Count、source 和 materialized table 均存在 |
| Frozen evaluation contract | Validated | Split 与 seeds 明确并已物化 |
| 模型与特征 workflow | Implemented | 代码和结果摘要存在 |
| 主要性能 claim | Partially validated | 仍有结果版本冲突且统计检验不完整 |
| 候选优选 | Partially implemented | 已计算但未实验验证 |
| 发布级复现 | Blocked | 缺少 outputs、lock、部分源码和可移植入口 |
| RNA target discovery | Planned | 没有实现证据 |

## 局限性

负责赋值 `weak_positive_confidence` 的上游 generator 缺失；部分 evaluation module 仅有 bytecode；当前工作区缺少 `outputs/`；若干 workflow 引用历史绝对路径；名义上相似的 Figure2 设置存在多个不一致结果。详见 [INCONSISTENCIES.md](INCONSISTENCIES.md)。

