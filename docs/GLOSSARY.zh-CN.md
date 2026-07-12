# EvoGATE 术语表

_科学、数据、模型和应用文档的受控术语。_

---

## 生物学术语

| 术语 | EvoGATE 中的定义 |
|---|---|
| Essential gene | 在明确生物学条件下维持生存所必需的基因。EvoGATE label 依据定义的 evidence 操作化该概念；prediction 不等于 essentiality 证明。 |
| Lethal phenotype | 记录 assay 中表示 non-viability 的 mutant 或 perturbation phenotype，是 PHI-supported positive component 使用的直接 phenotype category。 |
| Virulence | 病原体造成的损害或病害严重程度。Virulence 降低本身不意味着失去 viability。 |
| Pathogenicity | 导致疾病的能力。Loss of pathogenicity 不等于 lethality。 |
| Fitness | 在指定环境中的相对繁殖或生长表现。Reduced fitness 比 essentiality 更宽泛。 |
| Core essential | 在广泛条件下必需或跨 lineage 保守必需的基因。EvoGATE 未通过实验建立 universal core essentiality。 |
| Context-specific essentiality | 仅在特定 condition、host、tissue、developmental stage 或 environment 中成立的 essentiality。 |
| RNA target | 在 essentiality、accessibility、conservation 和 off-target assessment 后选择用于潜在 RNA-mediated suppression 的 transcript/gene。EvoGATE 尚未实现该阶段。 |
| dsRNA design | 选择 double-stranded RNA sequence，以产生有效 silencing fragment 并控制 off-target risk。状态：Planned。 |

## 进化与监督术语

| 术语 | EvoGATE 中的定义 |
|---|---|
| Evolution-aware supervision | 使用明确 evolutionary relationship 和 phenotype evidence 构建监督，而不是将稀缺 target-species annotation 视为唯一标签。 |
| Ortholog transfer | 在规定 confidence rule 下，将一个物种基因的证据迁移至 target-species ortholog。它是 derived evidence，不是 target-species direct validation。 |
| Single-copy ortholog | 在指定 scope 中每个相关物种有一个 copy 的 ortholog relationship/orthogroup。必须报告精确 scope。 |
| Orthogroup | 来自共同祖先基因的一组 genes，在 transfer artifact 中使用 `OG*` 等 ID。 |
| Positive label | 按规则接纳为 essential 的基因：PHI-supported lethal evidence 或 resolved high-confidence yeast-essentiality transfer。 |
| Negative label | 在移除 virulence/pathogenicity evidence 和 positive overlap 后保留的 resolved `none`-confidence gene。它是 operational negative，不是 universal non-essentiality proof。 |
| Exclusion set | 因 uncertain evidence、biological exclusion、unresolved mapping 或其他 protocol rule 而未进入 positive/negative 的基因。当前缺少完整 standalone manifest。 |
| Weak positive confidence | Yeast-transfer artifact 中已有 categorical field。下游使用已知；上游 producer 缺失，因此为 Unknown。 |

## Protocol 与评价术语

| 术语 | EvoGATE 中的定义 |
|---|---|
| Frozen protocol | 在正式比较中必须保持固定的 versioned label、split、graph/feature contract、seed、metric 和 setting。 |
| Split seed | `20260409`，用于物化 stratified 70/10/20 split。 |
| Training seed | `1029`-`1033` 之一，用于 frozen split 上的 model stochasticity。 |
| AUPRC | Precision-recall curve 下的面积，强调 positive-class ranking，是不平衡 essential-gene task 的主要指标。 |
| AUROC | Receiver operating characteristic curve 下的面积，测量 true-positive/false-positive rate 上的 ranking，是 AUPRC 的补充。 |
| MCC | Matthews correlation coefficient，使用 confusion matrix 四类结果的 threshold-dependent summary。 |
| Validation threshold | 只能从 validation prediction 选择的 decision threshold，不得来自 test outcome。 |
| Fixed threshold | Frozen runner 中标准 trainable-model 的 `0.5` decision threshold。 |

## 表征与模型术语

| 术语 | EvoGATE 中的定义 |
|---|---|
| PPI | 作为 graph structure 的 protein-protein interaction network。默认 processed source 为 STRING；eFG 用于 source-robustness analysis。 |
| ORT | 来自 processed cross-species orthology matrix 的 orthology feature block。 |
| EXP | Numeric expression feature block。 |
| SUB | Binary subcellular localization feature block。 |
| ESM2 embedding | 由 configured ESM2 model 生成的 mean-pooled protein representation；当前 full embedding dimension 为 1,280。 |
| GraphSAGE | 当前实验中的主要 message-passing implementation carrier，不是 EvoGATE 唯一创新。 |
| Gated fusion | Omics/ESM2 block 的 learned fusion variant；现有结果显示 metric-dependent trade-off。 |
| Group zero-out perturbation | 将 standardized feature column 在 inference 时置零以估计 group sensitivity。它不是 SHAP、GNNExplainer 或 causal attribution。 |

## 应用与状态术语

| 术语 | EvoGATE 中的定义 |
|---|---|
| Candidate prioritization | 使用 prediction score、cross-seed behavior、rank change 和 evidence profile 的计算排序，不表示 experimental validation。 |
| Predicted candidate | 由 computational model 优选的基因，应使用该词而不是 “validated target”。 |
| Implemented | 存在实现所述行为的 code 或 workflow。 |
| Validated | 已有证据在声明 scope 内验证所述 contract 或 result。 |
| Partially implemented | 部分所需 component 存在，但 end-to-end capability 不完整。 |
| Partially validated | 有证据，但不完整、冲突或不足以支持完整 claim。 |
| Planned | 未来计划，当前无 implementation evidence。 |
| Unknown | 现有证据无法确认。 |
| Blocked | 缺少 required dependency 或 artifact。 |
| Historical | 为 provenance、replay 或 migration 保留，而非 current mainline。 |

