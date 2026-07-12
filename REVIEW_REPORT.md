# ProGATE_v2 项目系统性审核报告

> 生成日期: 2026-07-07
> 审核方法: 6阶段系统性审核

---

## Phase 0: 项目现状摘要

### 0.1 项目概况

ProGATE_v2 是一个基于图神经网络（GraphSAGE）和蛋白语言模型（ESM2）的真菌植物病原体必需基因预测框架。项目从EPGAT（Schapke et al., 2021, GAT架构）演化而来，经历了从脚本驱动到注册表/配置驱动的架构迁移。

| 维度 | 详情 |
|------|------|
| **目标物种** | *Fusarium graminearum* (小麦赤霉病菌, PH-1菌株) |
| **支持物种** | Human, C. elegans, S. cerevisiae, D. melanogaster (4种模式生物) |
| **主模型** | GraphSAGE (DGL/PyG双后端) + Gated Residual Multimodal Fusion |
| **特征模态** | 直系同源(ORT) + 基因表达(EXP) + 亚细胞定位(SUB) + ESM2蛋白embedding |
| **标签来源** | PHI-base lethal基因 + 酵母必需基因直系同源迁移 (newlabel: 1097正/10868负) |
| **ESM2模型** | esm2_t33_650M_UR50D (650M参数, embedding dim 1280) |
| **PPI网络** | STRING v12 (threshold=300) + eFG (Fusarium专用) |
| **评估协议** | 5种子(1029-1033), 80/10/10 train/val/test分层划分 |

### 0.2 当前5组Figure验证的假设

| Figure | 验证假设 | 核心配置 |
|--------|---------|---------|
| **Figure1** | 标签稀缺是领域核心挑战；标签重建流程可追溯；跨物种标签迁移可行 | 标签来源Sankey、物种特征表、经典基线benchmark |
| **Figure2** (label_scarcity) | GraphSAGE+ESM2在低标签率(10%-90%)下显著优于传统ML方法和纯网络中心性方法 | 7种模型 x 9个标签率 x 5种子 |
| **Figure3** (3a/3b/3c) | ESM2特征对GraphSAGE有额外增益; gated residual融合优于简单拼接; 不同融合策略和损失函数有显著差异 | 早期拼接 vs gated vs residual_gated_concat x BCE vs wBCE |
| **Figure4** | 图结构质量(STRING阈值、eFG来源)影响模型鲁棒性,但GraphSAGE+ESM2在合理范围内稳定 | 多阈值扫描(100-900), 多PPI来源对比 |
| **Figure5** (5a-5d) | ESM2引入后,隐层表征可解释: UMAP错误迁移、特征组归因、新候选致病基因有生物学意义 | UMAP投影、输入vs隐层对比、特征组消融(zero-out masking)、候选基因优先级排序 |

### 0.3 架构迁移历史 (EPGAT -> ProGATE_v2)

基于 `docs/epgat_migration/` 文档还原的关键事实:

1. **迁移动机**: EPGAT仓库已成为"原始GAT论文代码 + 后论文扩展"的混合体(GCN/GIN/GraphSAGE/ESM2/ESMC/ProtT5/ProtBERT/RF/SVM/NB/centerality baselines)，脚本和结果目录混乱，不可复现
2. **迁移内容**: 将EPGAT的功能移植到ProGATE_v2的注册表驱动、配置驱动架构中(非简单模型替换)
3. **关键发现**: ProGATE_v2 **同时保留了GAT、GCN、GIN作为baseline**，GraphSAGE是经过消融实验筛选后的主模型选择
4. **已知问题**: `evaluate_graph_model.py`等核心评估文件源代码丢失(仅保留.pyc缓存); EPGAT原版结果表的`prepare_data.py`存在人工指标交换等不可靠后处理

---

## Phase 1: 文献定位与研究空白

### 1.1 7篇已知相关工作的代码可用性

| 论文 | 期刊(IF) | 代码可用 | 备注 |
|------|---------|---------|------|
| **EPGAT** (Schapke 2021) | IEEE/ACM TCBB (~3.0) | [github.com/JSchapke/essential-gene-detection](https://github.com/JSchapke/essential-gene-detection) | 本项目直接前身 |
| **Bingo** (Ma 2024) | Briefings in Bioinformatics (~9.5) | **未找到公开仓库** | 头号对标论文 |
| **AraPathogen2.0** (Lei 2024) | J. Proteome Research (~4.4) | 未确认 | 序列+图+XGBoost同类方法 |
| **Fungtion** (Li 2024) | J. Molecular Biology (~5.6) | 未确认 | 真菌效应蛋白预测竞品 |
| **MFGAC-PPI** (2024) | Frontiers in Plant Science (~5.7) | 未确认 | 含AlphaFold结构特征 |
| **CLEF** (2025) | Nature Communications (~16.6) | 可能可用 | 对比学习融合 |
| **Meta-TGLink** (Yu 2025) | Genome Biology (~12.3) | [github.com/Yoyiming/Meta-TGLink](https://github.com/Yoyiming/Meta-TGLink) | 图元学习+小样本 |

### 1.2 文献定位: 创新点对比

| 维度 | EPGAT | Bingo | Meta-TGLink | **ProGATE_v2** |
|------|-------|-------|-------------|----------------|
| **GNN架构** | GAT | GNN(未详述) | GNN+Transformer | **GraphSAGE (含GAT/GCN/GIN对比)** |
| **蛋白语言模型** | 无 | ESM-2 | 无 | **ESM2 (esm2_t33_650M)** |
| **多模态融合** | 早期拼接 | 注意力融合 | 位置编码+结构增强 | **Gated Residual (含3种融合消融)** |
| **标签迁移** | 无 | Zero-shot跨物种 | 跨域meta-learning | **直系同源标签迁移 + 标签稀缺梯度** |
| **可解释性** | 无 | GNNExplainer | 注意力权重 | **特征组消融 + UMAP + 误差迁移** |
| **目标领域** | 模式生物 | 后生动物模式生物 | 人类细胞系GRN | **真菌植物病原体(非模式生物)** |

### 1.3 方法演化链

```
EPGAT (2021, GAT, 无PLM)
  -> Bingo (2024, ESM2+GNN, zero-shot, GNNExplainer)
    -> ProGATE_v2 (2026, GraphSAGE+ESM2+gated residual+标签稀缺+真菌病原体)
```

### 1.4 关键空白与差异点分析

#### 相比Bingo (头号对标)

| Bingo特征 | ProGATE_v2状态 | 差距 |
|-----------|---------------|------|
| ESM-2+GNN融合 | 已实现 | 对齐 |
| Zero-shot跨物种 | Figure2 label_scarcity部分覆盖 | 概念不同(Bingo跨物种不训练目标物种,本项目在目标物种上减少训练标签) |
| GNNExplainer可解释性 | Figure5d特征组消融 | 特征组消融不等于GNNExplainer,需说明选择理由 |
| 对抗训练 | 未实现 | 缺失 |
| 应用于真菌 | Bingo未涉及 | **本项目的核心差异化** |

**最关键问题**: 缺少与Bingo的直接方法论对比或实验对比。`src/train/run_bingo_440_replay.py` 是MLP+旧440标签回放,**并非Bingo模型实现**。

#### 相比EPGAT (直接前身)

- ProGATE_v2已保留GAT作为baseline模型,但**未做GAT vs GraphSAGE在同一特征配置下的严格消融对比**
- `docs/epgat_migration/`集中于数据/架构迁移,不包含GAT vs GraphSAGE的消融实验数据

#### 相比Meta-TGLink (标签稀缺对标)

- Meta-TGLink的"few-shot"使用跨域meta-learning(MAML)
- ProGATE_v2的"label scarcity"是在目标物种上减少训练标签
- **评估协议不可直接对齐,需新增跨物种迁移实验**

#### 相比MFGAC-PPI (结构特征缺口)

- MFGAC-PPI整合了AlphaFold结构特征,**本项目未使用结构模态**
- 需要在Discussion中论证序列+图特征在当前场景下足够

### 1.5 相关工作的后续引用 (2024-2026)

- **EPGAT**被Bingo等约25+篇论文引用
- **Bingo** (2024年初发表) 暂未发现将其方法直接移植到植物病原真菌的后续工作
- 暂未发现将ESM2+GNN+zero-shot应用于真菌植物病原体的同类工作 -> **本项目在该交叉点上具有先发优势**

---

## Phase 2: Figure叙事线核实

### 2.1 当前叙事线评估

| Figure | 当前Claim | 支撑数据状态 | 问题 |
|--------|----------|-------------|------|
| **Figure1** | 标签稀缺挑战 + 标签重建流程 + 来源构成 | 标签来源Sankey数据存在,标签重建文档详尽 | 1C的Sankey图可能需要主文精简版 |
| **Figure2** | GraphSAGE在多标签率下优于传统方法 | 7模型x9标签率x5种子已运行 | **缺少置信区间/显著性检验**;未与Meta-TGLink对齐 |
| **Figure3a** | ESM2 vs 无ESM2对比 | 结果存在 | 仅fgraminearum+scerevisiae |
| **Figure3b** | ESM2维度消融(160/320/640/1280) | 结果存在 | 维度仅覆盖pooling维度,未消融不同ESM2层 |
| **Figure3c系列** | 融合策略消融 | 结果存在且详尽 | **设计合理,这是论文最强消融证据之一** |
| **Figure4** | 图结构鲁棒性 | STRING阈值扫描+eFG多置信度 | 仅fgraminearum,未跨物种 |
| **Figure5a-d** | 可解释性(UMAP+特征归因+候选基因) | workflow定义完整 | 特征归因是group-level zero-out ablation,非SHAP |

### 2.2 更合理的Figure组织建议

```
Figure 1: 问题框架
  1A - 真菌病原体标签稀缺挑战示意图
  1B - 跨物种直系同源标签迁移流程图
  1C - 标签来源构成 (Sankey, 精简版)

Figure 2: 核心方法
  2A - ProGATE_v2整体架构 (GraphSAGE + gated residual fusion)
  2B - 融合策略对比 (gated vs residual_gated_concat vs concat)
  2C - 标签重建流程与oldlabel vs newlabel对比

Figure 3: 核心结果 - 标签稀缺
  3A - 标签稀缺benchmark主图 (当前Figure2)
  3B - ESM2贡献 + 维度消融 (当前Figure3a/3b)
  3C - 跨物种泛化 (需新增实验)

Figure 4: 鲁棒性与消融
  4A - 图结构鲁棒性 (当前Figure4)
  4B - GNN架构消融 (GAT vs GCN vs GIN vs GraphSAGE)
  4C - 特征消融 (当前Figure2b的GNN特征消融)

Figure 5: 可解释性与生物学验证
  5A - UMAP隐层表征 (当前Figure5a/5c)
  5B - 特征组归因 (当前Figure5d)
  5C - 新候选致病基因优先级排序

Supplemental:
  S1 - 多物种benchmark完整表
  S2 - 阈值扫描细节 (FigureS4C)
  S3 - 统计检验补充 (需新增)
```

### 2.3 逻辑跳跃/不一致

1. **Figure2与Figure3的连接**: Figure2用的是GraphSAGE_ORT_EXP_SUB_ESM2(早期拼接),但Figure3c才验证gated fusion更好 -> 存在"用次优配置做label scarcity benchmark"的可能
2. **物种覆盖不均**: Figure2仅fgraminearum, Figure3a含scerevisiae, Figure4仅fgraminearum -> 跨物种泛化证据不足
3. **Figure5与Figure3a的依赖**: Figure5基于Figure3a的run(早期拼接),如果gated residual更好,解释对象应该是更好的模型

---

## Phase 3: 缺口分析与补充实验清单

### 3.1 缺失对照组 (按优先级)

#### P0 (必须做,否则无法投稿)

1. **Bingo方法论对比/实验对比**
   - 当前状态: `src/train/run_bingo_440_replay.py` 是MLP+旧440标签回放,非Bingo实现
   - 方案A(推荐): 论证ProGATE_v2在目标领域(真菌)、直系同源图结构、gated residual融合、标签重建流程等方面与Bingo有本质区别
   - 方案B: 如Bingo代码可用,在fgraminearum数据上运行Bingo方法直接对比
   - **如Bingo代码未公开,方案A是唯一路径,但在Discussion中需明确标注"Bingo代码未公开"作为局限**

2. **GAT vs GraphSAGE vs GCN vs GIN架构消融**
   - 当前状态: Figure1 benchmark中含GAT和GraphSAGE但特征配置不同,非严格消融
   - 需补充: 在相同特征配置(ORT_EXP_SUB_ESM2)下,4个GNN架构的严格对比(5种子,含95%CI)
   - 成本: 低(仅需改config运行现有pipeline)
   - **此实验同时回应"为什么选GraphSAGE"和"与EPGAT的关系"两个问题**

#### P1 (应该做,显著增强说服力)

3. **跨物种泛化实验**
   - 需补充: 在至少1个额外物种(celegans或melanogaster)上验证,或模式生物train->真菌test

4. **统计显著性检验**
   - 当前状态: 仅计算mean和std(ddof=0),无CI,无显著性检验
   - 需补充: 95%置信区间 + 配对差异检验(Wilcoxon signed-rank test)

5. **纯序列baseline (无图结构)**
   - ESM2 embedding + MLP (仅序列,无图结构,无手工特征)
   - 目的: 量化图结构的边际贡献

#### P2 (锦上添花)

6. **与Meta-TGLink的直接对比** (代码可用但任务不同,可在Discussion中讨论)
7. **AlphaFold结构特征覆盖率分析** (如fgraminearum覆盖率<30%,可作为不纳入结构特征的客观理由)
8. **AraPathogen2.0对比** (如代码可用)

### 3.2 数据泄露风险评估

| 检查项 | 状态 | 风险 |
|--------|------|------|
| **ESM2特征预计算** | 预训练ESM2(非微调),提取蛋白质序列embedding | 无泄露 |
| **直系同源特征** | 基于InParanoid序列同源性 | 风险极低(同源检测不依赖必需性) |
| **标签迁移** | yeast essential -> fgraminearum ortholog | 有意信息迁移,论文需明确标注 |
| **Split划分** | 基因级别分层划分(`stratify=y`) | 标准做法 |
| **特征标准化** | Z-score using training split statistics | 标准做法(`maybe_normalize_features`在train_idx上计算) |
| **阈值选择** | Figure3c阈值调优 | **需确认阈值是否在验证集上调的(非测试集)** |

### 3.3 Figure5候选基因验证

当前new_candidate_prioritization依赖模型预测分数排序,验证渠道:
- 文献验证: PHI-base中对应基因的致病性记录(非lethal表型)
- 同源验证: 候选基因在yeast中的同源基因是否为必需基因
- 湿实验: 需独立实验验证(Discussion中作为future work)

---

## Phase 4: 代码审核

### 4.1 阻断性问题 (必须修复)

#### B1. `evaluate_graph_model.py`等核心评估文件源代码丢失

- **状态**: `src/eval/__pycache__/evaluate_graph_model.cpython-*.pyc`存在但`.py`源文件不存在
  - `evaluate_graph_model.py` - 缺失
  - `evaluate_baseline.py` - 缺失
  - `evaluate_support_graph_baseline.py` - 缺失
- **影响**: 无法审计评估逻辑正确性,无法复现结果,无法投稿
- **修复**: 从.pyc反编译或从git历史恢复

#### B2. 特征归因方法非标准SHAP

- **状态**: `build_figure5d_feature_group_attribution.py`使用group-level zero-out masking(设特征组为0,计算delta probability)
- **问题**: 这不是SHAP,不是GNNExplainer。论文中若声称使用了SHAP,则为错误陈述
- **修复**: 明确标注为"feature group ablation (zero-out masking)",不声称是SHAP

#### B3. 统计检验缺失

- **状态**: 所有evaluation仅报告mean和population std(ddof=0),无置信区间,无显著性检验
- **影响**: 生物信息学期刊(e.g., Briefings in Bioinformatics, Genome Biology)通常要求统计检验
- **修复**: 添加95% CI和配对差异检验

#### B4. `run_bingo_440_replay.py`命名误导

- **状态**: 文件名含"bingo"但实际是MLP+旧440标签回放,与Bingo(Ma et al., 2024)完全无关
- **修复**: 重命名为`run_old440_mlp_replay.py`

### 4.2 建议性改进

#### S1. Figure3c多配置间参数一致性

4个配置的差分已通过config记录,但缺少汇总表。建议在`docs/`下添加`figure3c_ablation_design.md`

#### S2. GraphSAGE实现质量

- `src/models/graph_models.py`中的自实现GraphSAGE(mean aggregator, index_add)实现正确
- PyG版本的SAGEConv使用正确
- Gated residual fusion模块(`multimodal_gated_fusion.py`)的两种模式设计合理:
  - `gated`: 双侧投影+元素级门控
  - `residual_gated_concat`: ESM2投影+门控+omics直通残差连接
- 无数据泄露风险(特征预计算与模型训练分离)

#### S3. Split种子控制

- `freeze_unified_protocol.py`的`assign_splits()`正确使用`random_state=seed`和`stratify=y`
- 5个种子(1029-1033)固定,可复现
- **但**: 跨种子结果报告为mean +/- sd,未做paired统计检验(见B3)

#### S4. workflow/*.smk路径一致性

- 各Figure的workflow共享`frozen_protocol.yaml`配置
- Figure3c系列有独立的`output_root`和`results_root`,确保输出隔离
- 未发现跨Figure硬编码路径冲突

---

## Phase 5: 投稿水平评估

### 5.1 期刊定位参考

| 同类工作 | 发表期刊 | IF | 评估标准特点 |
|---------|---------|-----|------------|
| EPGAT | IEEE/ACM TCBB | ~3.0 | 计算方法学+基准测试 |
| AraPathogen2.0 | J. Proteome Research | ~4.4 | 实验验证+方法对比 |
| Fungtion | J. Molecular Biology | ~5.6 | 方法+可视化+生物学insight |
| Bingo | Briefings in Bioinformatics | ~9.5 | 方法论创新+可解释性+zero-shot |
| Meta-TGLink | Genome Biology | ~12.3 | 方法论创新+小样本+跨域 |
| CLEF | Nature Communications | ~16.6 | 方法论创新+多模态+对比学习 |

### 5.2 推荐投稿目标

| 优先级 | 期刊 | 理由 |
|--------|------|------|
| **第一选择** | **Briefings in Bioinformatics** | Bingo同刊,接受方法学+benchmark+可解释性论文 |
| 第二选择 | **Bioinformatics** | 经典生信方法学期刊 |
| 第三选择 | **Genome Biology** | Meta-TGLink同刊,接受GNN+小样本方向 |
| 备选 | **PLOS Computational Biology** | 开放获取,方法学友好 |

### 5.3 模拟审稿意见(最可能的3-5个问题)

1. **"The authors use GraphSAGE instead of GAT (as in EPGAT). An ablation study comparing GAT, GCN, GIN, and GraphSAGE under the same feature configuration is needed."** -> P0-2

2. **"Bingo (Ma et al., 2024) also combines ESM-2 with GNNs and evaluates zero-shot transfer. A direct comparison or detailed methodological differentiation is required."** -> P0-1

3. **"The feature attribution analysis (Figure 5d) uses zero-out masking rather than an established method like SHAP or GNNExplainer. Please justify this choice."** -> B2

4. **"Statistical significance is not reported for any performance comparisons. Confidence intervals or hypothesis tests should be added throughout."** -> B3

5. **"The study focuses primarily on F. graminearum. Cross-species generalization experiments would strengthen the claims."** -> P1-3

### 5.4 整体录用概率判断

- **当前状态: 弱** — 方法论设计完整但缺少关键实验(Bingo对比/GAT消融/统计检验)会直接导致拒稿
- **补充P0后: 中** — 补齐后方法论叙事完整,真菌植物病原体场景有领域独特性
- **补充P0+P1后: 中强** — 有完整消融体系+跨物种泛化+统计显著性+生物学可解释性

---

## Phase 6: Agent时代定位讨论

### Discussion素材: 为什么专用GNN架构不可被Agentic AI替代

1. **数据稀缺场景下的架构设计需要领域知识**: LLM可生成GNN训练代码,但无法设计针对"真菌病原体标签仅1000+正例"场景的gated residual fusion策略。直系同源图的构建(Ontology-based canonical ID mapping, cross-species ortholog bridge)需要大量领域特定知识。

2. **严格消融验证不可被替代**: Figure3c的4种融合策略对比、Figure3b的4种ESM2维度消融、Figure4的9种图结构扰动——这种系统消融需要实验设计能力,而非LLM的文本生成能力。

3. **标签重建的protocolized provenance**: fgraminearum_newlabel的构建涉及PHI-base->lethal、yeast ortholog->transfer、XP_*->FGRAMPH1_* bridge等多层provenance,这种数据治理工作无法被Agent自动化。

4. **建议论文中的表述**:
   > "While recent advances in agentic AI and large language models have enabled automated literature review, hypothesis generation, and workflow orchestration, the core contributions of this work—designing a gated residual fusion architecture tailored to label-scarce fungal pathogens, systematically ablating each component, and establishing a protocolized label provenance pipeline—remain tasks that require domain-specific experimental design and rigorous empirical validation, which current agentic systems cannot yet perform autonomously."

---

## 补充实验优先级总结

```
P0 (必须做 -> 投稿门槛):
├── P0-1: Bingo方法论对比/实验对比 [约1-2周]
├── P0-2: GAT vs GraphSAGE vs GCN vs GIN架构消融 [约1-2天计算]
├── P0-3: 恢复丢失的评估源代码(evaluate_graph/evaluate_baseline) [约1天]
├── P0-4: 添加统计显著性检验(95%CI + 配对检验) [约2-3天]
└── P0-5: 确认Figure3c阈值在验证集上调(非测试集) [约1天审计]

P1 (应该做 -> 显著提升):
├── P1-1: 跨物种泛化实验(至少+1物种) [约1-2周计算]
├── P1-2: 纯序列baseline (ESM2+MLP, 无图结构) [约1天计算]
├── P1-3: 特征归因方法升级为GNNExplainer或明确标注 [约3-5天]
└── P1-4: 审计Figure2 label_scarcity是否用了最优融合策略 [约1天]

P2 (锦上添花):
├── P2-1: 与Meta-TGLink对比(如可适配)
├── P2-2: AlphaFold结构特征覆盖率分析
├── P2-3: Figure配色规范性检查
└── P2-4: 湿实验验证计划(Discussion future work)
```

---

## 附录

### A. 配置文件差分矩阵 (Figure3c系列)

| 配置 | 融合模式 | 损失函数 | 特征设置 |
|------|---------|---------|---------|
| Figure3c | `gated` | BCE (默认) | ORT_EXP_SUB_ESM2 vs ORT_EXP_SUB_ESM2_GATED |
| Figure3cA | `residual_gated_concat` | BCE (balanced) | + ORT_EXP_SUB_ESM2_GATED_RESIDUAL |
| Figure3cB | `residual_gated_concat` | **weighted BCE** | 同3cA |
| Figure3cC | `gated` | weighted BCE | 同3c, 旧版gate实现 |

### B. 已知代码问题清单

| 文件 | 问题 | 严重度 |
|------|------|--------|
| `src/eval/evaluate_graph_model.py` | 源代码丢失(仅.pyc) | 阻断 |
| `src/eval/evaluate_baseline.py` | 源代码丢失(仅.pyc) | 阻断 |
| `src/eval/evaluate_support_graph_baseline.py` | 源代码丢失(仅.pyc) | 阻断 |
| `src/train/run_bingo_440_replay.py` | 文件名误导(非Bingo模型) | 中 |
| `src/analysis/build_figure5d_feature_group_attribution.py` | 非标准SHAP方法 | 中 |
| `configs/frozen_protocol.yaml` | `env_name: EPGAT` 命名残留 | 低 |

### C. 数据规模快照

| 物种 | 节点数(约) | 边数(STRING300, 约) | 正标签 | 负标签 | ESM2维度 |
|------|-----------|-------------------|--------|--------|----------|
| Human | 18000 | 300K | OGEE标准 | OGEE标准 | 1280 |
| C. elegans | 14000 | 200K | OGEE标准 | OGEE标准 | 1280 |
| S. cerevisiae | 6000 | 200K | OGEE标准 | OGEE标准 | 1280 |
| D. melanogaster | 11000 | 100K | OGEE标准 | OGEE标准 | 1280 |
| **F. graminearum** (newlabel) | 12000 | 80K(STRING)/14K(eFG) | **1097** | **10868** | 1280 |
| F. graminearum (oldlabel) | 12000 | 同上 | 439 | 10750 | 1280 |
