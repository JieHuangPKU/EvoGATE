# EvoGATE 架构

_代码架构、科学架构、模块边界与执行关系。_

---

## 七层代码架构

```mermaid
flowchart TD
    accTitle: EvoGATE code layers
    accDescr: 七个仓库层次将证据转换为标签、特征、图、训练模型、评价和候选排序 artifact。

    knowledge_layer["Knowledge layer<br/>data/derived_labels, data/interim, data/manifests"]
    label_layer["Label layer<br/>data/processed/essential_gene, results/frozen_protocol"]
    feature_layer["Feature layer<br/>data/processed/OR, EXP, LC, ESM2"]
    graph_layer["Graph layer<br/>data/processed/PPI, src/data, src/graph"]
    model_layer["Model layer<br/>src/models, src/train"]
    evaluation_layer["Evaluation layer<br/>src/eval, src/analysis, src/plot"]
    application_layer["Application layer<br/>Figure5 candidate prioritization"]

    knowledge_layer --> label_layer
    knowledge_layer --> feature_layer
    label_layer --> graph_layer
    feature_layer --> graph_layer
    graph_layer --> model_layer --> evaluation_layer --> application_layer
```

| 层 | 职责 | 主要路径 | 规范入口 |
|---|---|---|---|
| Knowledge | Evidence、species scope、transfer artifact 和 ID provenance | `data/derived_labels/`、`data/interim/`、`data/manifests/` | `src.data.build_fgraminearum_newlabel_bridge` |
| Label | Positive/negative regime 和 frozen split | `data/processed/essential_gene/`、`results/frozen_protocol/` | `workflow/fgraminearum_label_materialization.smk`、`src.data.freeze_unified_protocol` |
| Feature | ORT、EXP、SUB、ESM2 feature block | `data/processed/OR/`、`data/processed/EXP/`、`data/processed/LC/`、`data/processed/ESM2/`、`src/features/` | 各模态 builder；无统一 feature workflow |
| Graph | PPI filtering、node universe、edge index、topology embedding | `data/processed/PPI/`、`src/data/frozen_protocol_loader.py`、`src/graph/` | `src.data.frozen_protocol_loader` |
| Model | Classical、topology、GNN 与 fusion model | `src/models/`、`src/train/` | `src.train.run_frozen_protocol_model` |
| Evaluation | Metric、aggregation、ablation、interpretation、plot | `src/eval/`、`src/analysis/`、`src/plot/`、`workflow/` | Figure workflow 与 evaluation module |
| Application | Candidate ranking 与未来 target discovery | `src/eval/build_figure5_candidate_prioritization.py`、`results/Figure5*` | Candidate module；RNA layer 为 Planned |

## 科学架构

```mermaid
flowchart TD
    accTitle: EvoGATE scientific architecture
    accDescr: 表型和分子证据经进化监督与图学习转化为候选排名，RNA 和 dsRNA 阶段明确标记为未来工作。

    evidence["Evidence<br/>PHI phenotypes, yeast essentiality, molecular data"]
    evolution["Evolution<br/>orthogroups, conservation, copy structure"]
    labels["Labels<br/>positive, negative, exclusion, frozen split"]
    representation["Representation<br/>ORT, EXP, SUB, ESM2"]
    graph_learning["Graph learning<br/>STRING/eFG and model families"]
    evaluation["Evaluation<br/>frozen metrics, ablation, robustness"]
    prioritization["Candidate prioritization<br/>scores, stability, rank shifts"]
    rna_target["RNA target discovery<br/>Planned"]
    dsrna_design["dsRNA design<br/>Planned"]

    evidence --> evolution --> labels
    evidence --> representation
    labels --> graph_learning
    representation --> graph_learning --> evaluation --> prioritization
    prioritization -. future .-> rna_target -. future .-> dsrna_design
```

| 科学阶段 | 状态 |
|---|---|
| Evidence assembly | Partially implemented |
| Evolutionary transfer artifact | Partially implemented |
| Label materialization | Validated |
| Multimodal representation | Validated |
| Graph learning | Validated |
| Evaluation | Partially validated |
| Candidate prioritization | Partially implemented |
| RNA target discovery | Planned |
| dsRNA design | Planned |

## 主要执行链

```mermaid
flowchart LR
    accTitle: Frozen benchmark call chain
    accDescr: Frozen protocol workflow 物化 manifest，加载对齐图数据，按配置与 seed 运行模型，再聚合 metric 并生成 plot。

    frozen_config["configs/frozen_protocol.yaml"]
    frozen_workflow["workflow/frozen_protocol_benchmark.smk"]
    freeze_protocol["src.data.freeze_unified_protocol"]
    frozen_loader["src.data.frozen_protocol_loader"]
    model_runner["src.train.run_frozen_protocol_model"]
    aggregate_runs["src.eval.aggregate_frozen_protocol_runs"]
    plot_rules["workflow/plots.smk"]

    frozen_config --> frozen_workflow
    frozen_workflow --> freeze_protocol
    frozen_workflow --> model_runner
    model_runner --> frozen_loader
    frozen_workflow --> aggregate_runs --> plot_rules
```

## 配置模型

`configs/frozen_protocol.yaml` 定义 repository-relative data root、protocol name、frozen runtime setting、feature root、ESM2 cache、label source、model family 和 hyperparameter。Figure-specific YAML 引用该 base config，并覆盖实验范围或模型 variant。

Per-run output directory 应保存 resolved runtime config。当前 `results/Figure3*/runtime/` 保留了部分 rendered config，但主要 `outputs/` tree 在本工作区缺失。

## 数据与 ID contract

Fusarium 主 canonical identifier 为 `fgraminearum::FGRAMPH1_*`；graph-facing file 可使用去掉前缀的 `FGRAMPH1_*`。`frozen_protocol_loader.py` 通过明确的 graph/canonical ID 连接 label、graph node、feature row 和 ESM2 key。

Loader 从 graph node 与 labeled node 的并集构建 node universe，将 frozen split 映射到 node index，使用 training node 归一化 numeric feature，并返回供所有 model family 消费的统一 bundle。

## 模型与输出 contract

标准运行写出 `predictions.tsv`、`metrics.tsv`、`feature_schema.tsv`、`edge_table.tsv`、`split_manifest.tsv`、`resolved_config.yaml`，以及 `best_model.pt`、`model.pkl`、`training_log.tsv` 或 ESM2 alignment audit 等 model-specific artifact。

当前工作区包含 aggregated result 与 Figure，但缺少主要 `outputs/` tree。因此 output contract 在代码中为 **Implemented**，但完整本地重建为 **Blocked**。

## 遗留边界

- `docs/epgat_migration/` 记录 EPGAT migration history
- `docs/protocol_refactor/` 记录 ProGATE_v2 protocol refactoring
- 若干 `scripts/run_*.sh` 硬编码历史 ProGATE_v2 path
- `src/` 中保留 legacy training/data adapter，用于 controlled replay
- 除非本文明确指定，historical artifact 不得作为当前 canonical entry point

非破坏性迁移策略见 [MIGRATION_GUIDE.zh-CN.md](MIGRATION_GUIDE.zh-CN.md)。

## 依赖与可移植性

代码 import Python、pandas、NumPy、PyYAML、scikit-learn、PyTorch、graph library、Snakemake、R 和 plotting package。当前没有 authoritative environment lock，部分 Python/cache path 绑定具体机器，因此 dependency 与 hardware 的发布级复现为 **Blocked**。
