# EvoGATE 科学设计

_科学假设、监督设计、表征、评价、解释与局限。_

---

## 科学问题

在实验标签稀缺的植物病原真菌中，能否结合 lethal phenotype evidence 与保守 essential orthology 重建监督信号，再整合互作网络、生物学特征和 protein-language 表征，预测并优先排序必需基因？

## 核心假设

EvoGATE 检验三个相互关联的假设：

1. 标签重构对可学习生物学问题的影响大于单纯更换 GNN family。
2. PPI topology、evolutionary feature、expression、localization 和 ESM2 含有互补信息。
3. 跨 seed 稳定预测可用于候选实验优选，但在独立验证前仍只是计算假设。

## 进化感知标签重构

### 正样本定义

`fgraminearum_newlabel` 的正样本集合为：

```text
PHI-supported lethal canonical genes
UNION
resolved high-confidence yeast-essentiality-transfer canonical genes
```

物化 count 为 77 个 lethal genes、1,045 个 transfer-supported genes、25 个重叠 genes 和 1,097 个唯一 positives。`src/data/materialize_fgraminearum_label_regimes.py` 在 `configs/fgraminearum_label_materialization.yaml` 驱动下消费这些证据。

### 负样本定义

负样本从 transfer confidence 为 `none` 且成功 bridge 的基因开始，移除具有 virulence/pathogenicity evidence 或出现在正样本集合中的基因。最终物化 10,868 个负样本。

### 排除逻辑

下列基因不会成为监督负样本：

- 高置信正样本
- 具有 virulence 或 pathogenicity evidence 的基因
- 低/中置信 transfer genes
- 未解析或 ambiguous ID mappings
- 其他未被正负样本规则接纳的基因

相关逻辑分布在 source preparation 和 materialization code 中，但缺少完整的 standalone exclusion table。状态：**Partially implemented**。

### Canonical ID bridge

`src/data/build_fgraminearum_newlabel_bridge.py` 将 PH-1 蛋白 ID 映射至 canonical `fgraminearum::FGRAMPH1_*` gene space。它使用配置中的 unified map、legacy mapping、sequence/header evidence 及各模态 mapping table。未解析的高置信 row 会保存在 audit output 中，而不是静默加入标签。

### 进化证据

`data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` 记录 *S. cerevisiae* 和 *S. pombe* essential-ortholog support、fungal occupancy、copy statistics、single-copy indicator、support class 和 confidence。最初赋值 `weak_positive_confidence` 的脚本缺失，因此生成过程为 **Unknown**；下游 artifact 及其使用方式已有记录。

## 特征模态

| 模态 | 输入 | 表征 | 作用 |
|---|---|---|---|
| ORT | `data/processed/OR/<species>/orthologs.csv` | Binary orthology matrix | Evolutionary context |
| EXP | `data/processed/EXP/<species>/profile.csv` | Numeric expression profile | Condition-associated activity |
| SUB | `data/processed/LC/<species>/subloc.csv` | Binary localization matrix | Cellular context |
| ESM2 | `data/processed/ESM2/<species>/esm2_pooled.pt` | Mean-pooled protein embedding | Sequence-derived representation |
| Degree | PPI edge table | Scalar graph feature | Local topology |

`src/data/frozen_protocol_loader.py` 使用 training node 计算 feature normalization。ESM2 alignment 是严格的：缺少 graph-node embedding 会报错，不会静默丢弃 row。

## 图构建

默认 graph source 是 `data/processed/PPI/<species>/string.csv`。Frozen protocol 按 `combined_score >= 300` 过滤 STRING edge，移除 self-loop 和 duplicate，默认 graph contract 为 `undirected_symmetrized`，且关闭 edge weight。Figure4 另外评价 STRING threshold 和 eFG graph source。

## 模型家族

| 家族 | 模型 | 科学作用 |
|---|---|---|
| Tabular | MLP、RF、SVM、NB | Non-graph feature baseline |
| Topology | node2vec+MLP、degree、closeness | Network-only 或 embedding baseline |
| Graph neural network | GAT、GCN、GIN、GraphSAGE | Message-passing comparison |
| Fusion variants | Concatenation、gated、residual gated、weighted BCE variants | Multimodal design ablation |

GraphSAGE 是主要实现载体。现有 artifact 不支持将其描述为唯一科学贡献。

## Frozen evaluation

| 项目 | Frozen value |
|---|---|
| Split | 70% train / 10% validation / 20% test |
| Split seed | `20260409` |
| Training seeds | `1029`-`1033` |
| 标准 decision threshold | `0.5` |
| Tuned-threshold source | 仅 validation split |
| 主要指标 | AUPRC、MCC、AUROC |

由于正类远少于负类，项目优先报告 AUPRC。MCC 同时汇总 confusion matrix 的四类结果，在不平衡任务中仍具有解释力。AUROC 作为补充 ranking metric。Test data 不得用于选择 model、threshold 或 hyperparameter。

## 消融与解释

已实现的分析包括 oldlabel vs newlabel、feature combination、ESM2 inclusion、ESM2 dimension truncation、fusion variant、label scarcity、graph threshold/source、hidden representation 和 feature-group zero-out perturbation。

Feature-group attribution 不是 SHAP 或 GNNExplainer。它测量 frozen model 中指定 standardized column 置零后的变化，只支持 group-level sensitivity，不支持 causal residue-level interpretation。

## 候选排序

`src/eval/build_figure5_candidate_prioritization.py` 组合 baseline 与 ESM2 prediction，计算 gene rank、top-k overlap、rank change 和 candidate-group profile。当前 candidate table 在 label field 指示的位置包含 labeled 与 unlabeled graph nodes；下游用户必须按实验问题显式过滤。

状态：**Partially implemented**，且未通过 wet-lab validation。

Candidate-ranking code 本身为 **Implemented**，端到端应用仍因输入与验证缺口而不完整；规范入口为 `src/eval/build_figure5_candidate_prioritization.py`，ID bridge 的规范实现为 `src/data/build_fgraminearum_newlabel_bridge.py`。

## 贡献与边界

### 科学贡献

- 面向 fungal essentiality 的 evolution-aware supervision
- 可审计的 positive、negative 和 exclusion rule
- 对互补 biological representation 的 frozen comparison
- 基于 prediction stability 与 evidence profile 的 candidate prioritization

### 工程贡献

- Canonical identifier bridge
- 共享 frozen loader 与 model contract
- Per-run provenance、prediction、metric 和 feature-schema output
- 结构化 result aggregation 和 Figure workflow

### 当前局限

- 缺少上游 confidence generator
- 缺少显式完整 exclusion manifest
- 部分 comparison 存在结果版本冲突
- 主要 claim 尚无完整的发布级统计推断
- 当前工作区缺少 run-level `outputs/`
- 优选候选没有实验验证

### 未来工作

RNA target discovery、off-target filtering、conserved target-region selection、dsRNA design、cross-species transfer evaluation 和 experimental validation 均为 **Planned**，不是当前 EvoGATE 能力。
