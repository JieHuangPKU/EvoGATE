# Workflow Runtime Results Dependency Audit

Update: this document records a historical audit state. The current mainline no longer uses
`results/label_rebuild_experiments/` as a runtime input. That historical directory is now an
archive / provenance location, and mainline Fusarium label assets come from
`data/processed/essential_gene/fgraminearum/...` via the label materialization workflow.

## 1. 总结结论

本次审计同时检查了 `workflow/*.smk`、对应 `scripts/*.sh`、以及实际被调用的 Python 入口和其下游读取链。结论如下：

- `workflow/classical_baseline_benchmark.smk`：这里记录的是一次历史审计状态；当前主线标签来源已经迁移到 `data/processed/essential_gene/fgraminearum/...`，不再以 `results/label_rebuild_experiments/...` 为运行时输入。
- `workflow/fgraminearum_old_newlabel_benchmark.smk`：这里记录的是迁移前的历史依赖状态；当前应将 `results/label_rebuild_experiments/...` 理解为 archive / provenance，而不是主线输入。
- `workflow/frozen_protocol_benchmark.smk`：当前主线通过 `configs/frozen_protocol.yaml` 读取 materialized processed labels，再生成 `results/frozen_protocol/labels/*.tsv` 与 `results/frozen_protocol/splits/*.tsv`。

直接回答重点问题：

- `results/label_rebuild_experiments/`：**现状态是历史 archive / provenance 目录**。本文件保留的是旧审计结论，不应再把它当成当前 3 个 workflow 的主线输入。
- `results/frozen_protocol/`：**是当前 3 个 workflow 的协议化中间输出根**。其中 `labels/` 与 `splits/` 会被后续训练步骤读取，因此“再跑 workflow 时”它们会被当作当前生成链的上游输入；但其来源是本 workflow 的 `freeze_protocol` 规则，逻辑上属于当前 DAG 内部产物，不属于外部历史分析结果依赖。`summary/` 仅 summary 用。
- `results/phase2a_figures/`：**未发现被这 3 个 workflow 或其调用 Python 代码读取**。
- `results/phase2b_new_label/`：**未发现被这 3 个 workflow 直接读取**。当前主线 frozen protocol 读取的是 materialized processed labels，不是 `results/phase2b_new_label/...`。
- `results/fgraminearum_feature_combo_analysis/`：**未发现被这 3 个 workflow 读取**。
- 其他 `results/...` 历史 benchmark/figure/analysis 目录：**未发现参与这 3 个 workflow 当前运行**。

对“运行 snakemake 的依赖是不是不应该在 results 里面？”的直接回答：

- **原则上，不应把非协议化、历史分析型 `results/...` 目录当作 workflow 运行依赖。**
- 但在这个仓库当前实现里，`results/frozen_protocol/labels`、`results/frozen_protocol/splits` 这类“冻结标签 / 冻结切分 / 协议化 manifest”虽然放在 `results/` 下，实际上承担的是**协议化中间输入**角色；从工程语义上它们更像 `data/processed` 或 `artifacts/frozen_protocol`，放在 `results/` 里容易让人误判，但不是最严重的问题。
- 这份文档中描述的风险已经被后续 materialization 主线吸收；现在 `results/label_rebuild_experiments/...` 应被视为历史 archive，而不是 frozen protocol 的当前源标签输入。
- 像 `results/phase2a_figures/`、`results/phase2b_new_label/figures`、`results/fgraminearum_feature_combo_analysis/`、旧 benchmark summary 这类明显是图表/分析/历史输出的目录，**不应该再作为运行依赖**。本次审计未发现它们被当前 3 个 workflow 读取。

## 2. 逐 Workflow 审计

### 2.1 `workflow/classical_baseline_benchmark.smk`

#### 直接 Snakemake input 依赖

代码位置：

- `freeze_protocol` 产出 `results/frozen_protocol/labels/*.md` 与 `results/frozen_protocol/splits/*.md` 摘要作为下游 rule 输入：[workflow/classical_baseline_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/classical_baseline_benchmark.smk#L63)
- 训练 rule 只声明依赖这些 `freeze_protocol` 摘要文件，而未显式声明原始标签源：[workflow/classical_baseline_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/classical_baseline_benchmark.smk#L79)
- 聚合 rule 只依赖本次 `OUTPUT_ROOT` 下生成的 `metrics.tsv`：[workflow/classical_baseline_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/classical_baseline_benchmark.smk#L148)

审计判断：

- Snakemake 表面上**没有直接把** `results/label_rebuild_experiments/...` 列在 `input:` 里。
- 但 `freeze_protocol` shell 实际调用 `src.data.freeze_unified_protocol`，[workflow/classical_baseline_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/classical_baseline_benchmark.smk#L71) 到 [workflow/classical_baseline_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/classical_baseline_benchmark.smk#L76)，因此真实源依赖藏在 Python 里。

#### 间接 Python 读取依赖

- 这里保留的是一次历史审计记录；当前 `configs/frozen_protocol.yaml` 已把 Fusarium old/new label 源文件配置到 `data/processed/essential_gene/fgraminearum/...`。
- `src.data.freeze_unified_protocol` 用 `_read_tsv()` 读取这些源文件：[src/data/freeze_unified_protocol.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/freeze_unified_protocol.py#L36)
- old/new label 的正负样本读取发生在 `load_pair_labels()`：[src/data/freeze_unified_protocol.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/freeze_unified_protocol.py#L94)
- 训练入口 `src.train.run_frozen_protocol_model` 调用 `load_protocol_dataset()`：[src/train/run_frozen_protocol_model.py](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_frozen_protocol_model.py#L441)
- `load_protocol_dataset()` 运行时读取 `results/frozen_protocol/labels/*.tsv` 和 `results/frozen_protocol/splits/*.tsv`，再从 `data/processed/*` 读取 STRING 和特征表，见 [src/data/frozen_protocol_loader.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/frozen_protocol_loader.py#L49)、[src/data/frozen_protocol_loader.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/frozen_protocol_loader.py#L56)、[src/data/frozen_protocol_loader.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/frozen_protocol_loader.py#L67)、[src/data/frozen_protocol_loader.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/frozen_protocol_loader.py#L105)、[src/data/frozen_protocol_loader.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/frozen_protocol_loader.py#L171)
- 聚合脚本只 `glob("**/metrics.tsv")` 扫描当前 `OUTPUT_ROOT`，不读取历史 `results/...`：[src/eval/aggregate_frozen_protocol_runs.py](/home/jiehuang/software/fungi/ProGATE_v2/src/eval/aggregate_frozen_protocol_runs.py#L42)

#### 对重点目录的判断

- `results/label_rebuild_experiments/`：**历史上曾是隐藏在 Python 里的上游输入；现已降级为 archive / provenance**
- `results/phase2a_figures/`：**不依赖**
- `results/phase2b_new_label/`：**不依赖**
- `results/fgraminearum_feature_combo_analysis/`：**不依赖**
- `results/frozen_protocol/` 旧结果：
  - `labels/`、`splits/`：训练运行时会读取，但它们属于当前 workflow 的协议化中间产物
  - `summary/`：classical workflow 本身不读取 `results/frozen_protocol/summary`
  - 如果这些文件已存在，Snakemake 可能因时间戳判定跳过 `freeze_protocol`，从而复用现有 `results/frozen_protocol/labels` 与 `splits`；这属于“当前 output root 被复用为中间产物”，不是去读旧 benchmark 结果

### 2.2 `workflow/fgraminearum_old_newlabel_benchmark.smk`

#### 直接 Snakemake input 依赖

- `freeze_protocol` 与 classical/frozen benchmark 一样，先产出 frozen protocol 摘要：[workflow/fgraminearum_old_newlabel_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/fgraminearum_old_newlabel_benchmark.smk#L41)
- feature-combo 训练 rule 只声明依赖 frozen protocol 摘要文件：[workflow/fgraminearum_old_newlabel_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/fgraminearum_old_newlabel_benchmark.smk#L57)
- 聚合 rule 只依赖本次 `outputs/fgraminearum_old_newlabel_benchmark/.../metrics.tsv`：[workflow/fgraminearum_old_newlabel_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/fgraminearum_old_newlabel_benchmark.smk#L91)
- 分析 rule 只依赖本 workflow 当前生成的 summary TSV：[workflow/fgraminearum_old_newlabel_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/fgraminearum_old_newlabel_benchmark.smk#L114)

#### 间接 Python 读取依赖

- 这一步记录的是历史审计状态；当前 `freeze_protocol` 已改为经由 `configs/frozen_protocol.yaml` 读取 materialized processed labels，而不是 `results/label_rebuild_experiments/...`。
- 训练入口 `src.train.run_frozen_protocol_feature_combo_model` 同样调用 `load_protocol_dataset()`，因此读取 `results/frozen_protocol/labels/*.tsv`、`results/frozen_protocol/splits/*.tsv` 与 `data/processed/*`，不读取其他历史 `results/...`：[src/train/run_frozen_protocol_feature_combo_model.py](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_frozen_protocol_feature_combo_model.py#L36) 到 [src/train/run_frozen_protocol_feature_combo_model.py](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_frozen_protocol_feature_combo_model.py#L53)
- 分析脚本 `src.eval.analyze_fgraminearum_old_newlabel_feature_combo` 只读取本次 `summary_dir` 下的 `aggregated_metrics.tsv` 和 `per_run_metrics.tsv`：[src/eval/analyze_fgraminearum_old_newlabel_feature_combo.py](/home/jiehuang/software/fungi/ProGATE_v2/src/eval/analyze_fgraminearum_old_newlabel_feature_combo.py#L54)

#### 对重点目录的判断

- `results/label_rebuild_experiments/`：**当前应视为 historical archive / provenance，不是 frozen protocol 的现行标签源**
- `results/phase2a_figures/`：**不依赖**
- `results/phase2b_new_label/`：**不依赖**
- `results/fgraminearum_feature_combo_analysis/`：**不依赖**；当前 workflow 的 analysis 输出目录是 `results/fgraminearum_old_newlabel_feature_combo_analysis`
- `results/frozen_protocol/` 旧结果：
  - `labels/`、`splits/`：运行时读取，但属于当前 workflow 的协议化中间产物
  - `summary/`：当前 workflow 不读取 `results/frozen_protocol/summary`

### 2.3 `workflow/frozen_protocol_benchmark.smk`

#### 直接 Snakemake input 依赖

- `freeze_protocol` 是整个 benchmark 的起点：[workflow/frozen_protocol_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/frozen_protocol_benchmark.smk#L55)
- trainable 与 deterministic 两类训练 rule 均只声明依赖 frozen protocol 摘要文件：[workflow/frozen_protocol_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/frozen_protocol_benchmark.smk#L71) 与 [workflow/frozen_protocol_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/frozen_protocol_benchmark.smk#L104)
- 聚合 rule 只依赖本次 `outputs/frozen_protocol_benchmark_v2/.../metrics.tsv`：[workflow/frozen_protocol_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/frozen_protocol_benchmark.smk#L136)

#### 间接 Python 读取依赖

- 这里同样记录的是历史状态；当前 `freeze_protocol` 已从 `configs/frozen_protocol.yaml` 读取 materialized processed labels 作为 Fusarium labels 源。
- `src.train.run_frozen_protocol_model` 在运行时读取 `results/frozen_protocol/labels/*.tsv`、`results/frozen_protocol/splits/*.tsv`，并写出 `metrics.tsv` 中记录这些 manifest 路径：[src/train/run_frozen_protocol_model.py](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_frozen_protocol_model.py#L454) 到 [src/train/run_frozen_protocol_model.py](/home/jiehuang/software/fungi/ProGATE_v2/src/train/run_frozen_protocol_model.py#L558)
- 聚合脚本只扫描当前 benchmark `output_root`：[src/eval/aggregate_frozen_protocol_runs.py](/home/jiehuang/software/fungi/ProGATE_v2/src/eval/aggregate_frozen_protocol_runs.py#L42)

#### 对重点目录的判断

- `results/label_rebuild_experiments/`：**当前不是 frozen protocol 源标签；只保留为 historical archive / provenance**
- `results/phase2a_figures/`：**不依赖**
- `results/phase2b_new_label/`：**不依赖**
- `results/fgraminearum_feature_combo_analysis/`：**不依赖**
- `results/frozen_protocol/` 旧结果：
  - `labels/`、`splits/`：**会被当作输入读取**，但它们是本 workflow 自己生成的协议化中间件
  - `summary/`：聚合输出目录；不是训练输入
  - 已有历史文件存在时，Snakemake 可能直接复用现有 `labels/`、`splits/`，除非上游 config/代码/时间戳触发重建；因此它们是“当前 output root 中可被复用的协议化产物”，不是“完全无关的历史 benchmark 结果”

## 3. 依赖分类

### 3.1 必需上游输入

- `results/label_rebuild_experiments/old440/labels/positive_old440.tsv`
- `results/label_rebuild_experiments/old440/labels/negative_old440.tsv`
- `results/label_rebuild_experiments/labels/positive_set_P1.tsv`
- `results/label_rebuild_experiments/labels/negative_set.tsv`
- `results/frozen_protocol/labels/*.tsv`
- `results/frozen_protocol/splits/*.tsv`

分类说明：

- 前 4 个是**历史实验产物输入**。这份文档记录的是它们曾被 `freeze_unified_protocol` 当成源标签读取的历史状态；当前主线已迁移到 materialized processed labels。
- `results/frozen_protocol/labels/*.tsv` 与 `results/frozen_protocol/splits/*.tsv` 是**当前 workflow 生成链中的协议化 frozen labels / frozen splits**。它们作为运行输入是合理的，但目录位置放在 `results/` 下不理想，容易与“历史结果”混淆。

### 3.2 当前 workflow 运行生成的输出

- `results/classical_baseline_benchmark/`
- `results/fgraminearum_old_newlabel_benchmark/`
- `results/fgraminearum_old_newlabel_feature_combo_analysis/`
- `results/frozen_protocol/summary/`
- `results/frozen_protocol/labels/`
- `results/frozen_protocol/splits/`

说明：

- 这些目录在当前 3 个 workflow 中承担 output 或 workflow-internal artifact 角色。
- 其中 `results/frozen_protocol/labels/`、`splits/` 同时也是后续 rule 的输入，因此它们属于“当前 DAG 内部中间产物”，不是单纯终态输出。

### 3.3 历史分析产物但不参与当前运行

- `results/phase2a_figures/`
- `results/phase2b_new_label/figures/`
- `results/fgraminearum_feature_combo_analysis/`
- `results/phase2a_fixed300_multirun_analysis/`
- `results/tables/`

说明：

- 对这 3 个 workflow 及其调用 Python 代码，未发现读取这些目录的语句。
- 它们属于历史分析/画图/报告结果，不参与当前 benchmark 执行链。

### 3.4 可疑 / 需要迁移的残留依赖

- `results/label_rebuild_experiments/...`

说明：

- 这是本次审计最明确的“需要迁移/去依赖”问题点。
- 问题不在于路径名叫 `results`，而在于它承载的是**历史 label rebuild 实验输出**，却被 frozen protocol 当作当前主线 benchmark 的源输入。
- 更合理的归宿应是协议化数据目录，例如 `data/processed/frozen_labels_sources/`、`data/interim/protocol_refactor/` 或类似专用 artifact 根；至少不应继续依附于历史 benchmark/exploration 目录。

## 4. 风险判断

- 这份文档记录的是一次迁移前的风险审计：当时确实存在隐性依赖，且 frozen protocol 仍会在运行时读取这些历史结果文件。
- 当前状态下，这个风险已经通过 label materialization 主线被移除；`results/label_rebuild_experiments/...` 不再是 frozen protocol 的现行输入链。
- **`results/frozen_protocol/labels` 与 `splits` 属于可接受但需要明确语义的依赖。** 它们本质上是 frozen labels / frozen splits / rebuilt label manifests 这类协议化输入，作为运行依赖是合理的；问题是它们当前也放在 `results/` 下，和历史分析结果混在一起，容易掩盖真正的不当依赖。
- **未发现对旧 figures / analysis 目录的隐性读取。** 例如 `results/phase2a_figures/`、`results/phase2b_new_label/`、`results/fgraminearum_feature_combo_analysis/` 没有进入这 3 条 workflow 的当前运行时读取链。

## 审计依据

关键证据位置：

- `workflow/classical_baseline_benchmark.smk`：[workflow/classical_baseline_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/classical_baseline_benchmark.smk#L63)
- `workflow/fgraminearum_old_newlabel_benchmark.smk`：[workflow/fgraminearum_old_newlabel_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/fgraminearum_old_newlabel_benchmark.smk#L41)
- `workflow/frozen_protocol_benchmark.smk`：[workflow/frozen_protocol_benchmark.smk](/home/jiehuang/software/fungi/ProGATE_v2/workflow/frozen_protocol_benchmark.smk#L55)
- `configs/frozen_protocol.yaml` 的真实 label source 配置：[configs/frozen_protocol.yaml](/home/jiehuang/software/fungi/ProGATE_v2/configs/frozen_protocol.yaml#L42)
- `src.data.freeze_unified_protocol` 的运行时读取：[src/data/freeze_unified_protocol.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/freeze_unified_protocol.py#L36) 与 [src/data/freeze_unified_protocol.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/freeze_unified_protocol.py#L94)
- `src.data.frozen_protocol_loader` 的 manifest/feature 读取：[src/data/frozen_protocol_loader.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/frozen_protocol_loader.py#L49) 与 [src/data/frozen_protocol_loader.py](/home/jiehuang/software/fungi/ProGATE_v2/src/data/frozen_protocol_loader.py#L171)
- 聚合与分析脚本仅扫描当前 output/summary，而非历史 `results/...`：[src/eval/aggregate_frozen_protocol_runs.py](/home/jiehuang/software/fungi/ProGATE_v2/src/eval/aggregate_frozen_protocol_runs.py#L42)，[src/eval/analyze_fgraminearum_old_newlabel_feature_combo.py](/home/jiehuang/software/fungi/ProGATE_v2/src/eval/analyze_fgraminearum_old_newlabel_feature_combo.py#L54)
