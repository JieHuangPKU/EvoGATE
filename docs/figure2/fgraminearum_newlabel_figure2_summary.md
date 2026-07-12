# Figure 2: `fgraminearum_newlabel` 构建总结

## 1. 可直接用于稿件的描述

### 主文简述

当前仓库中的 `fgraminearum_newlabel` 不是直接沿用历史 `positive_set_P1.tsv`，而是在规范化 `fgraminearum::FGRAMPH1_*` 基因空间中重建得到：正类由 PHI 支持的 lethal 基因与高置信 yeast-essential ortholog transfer 基因并集组成，负类由 yeast-transfer 表中 `weak_positive_confidence == none` 的基因经 XP-to-canonical bridge 后，去除 virulence/pathogenicity 相关基因及全部正类后得到。

### Methods 描述

`fgraminearum_newlabel` 的当前 authoritative 构建流程由 `workflow/fgraminearum_label_materialization.smk` 调度，并由 `src/data/build_fgraminearum_newlabel_bridge.py`、`src/data/prepare_fgraminearum_label_materialization_sources.py` 和 `src/data/materialize_fgraminearum_label_regimes.py` 物化。首先，从 `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` 读取 PH-1 `XP_*` 蛋白层面的酵母 essential ortholog transfer 候选，并通过 protocolized bridge 显式映射到 canonical `fgraminearum::FGRAMPH1_*` 空间；该 bridge 优先使用 unified ID map 中唯一的 `XP_* -> FGRAMPH1_*` 对应关系，并结合 legacy proteome 的精确序列匹配与 NCBI header 中 `FGSG_*` 标记作为支持证据。随后，仅保留 `weak_positive_confidence == high` 且 bridge 成功解析的候选，并在 canonical gene 层面去重，得到 1045 个高置信 transfer-supported 基因，其中 25 个同时也出现在 lethal 集合中，故 transfer-only 贡献为 1020 个基因。其次，从 `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv` 读取 PHI 支持且满足 `evidence_term_raw == lethal`、`supports_gold_label == true` 的 77 个 lethal 基因。最终正类为两部分并集，共 1097 个基因。负类从同一 yeast-transfer 表中 `weak_positive_confidence == none` 且 bridge 成功解析的 11506 个基因开始，剔除 595 个带有 virulence/pathogenicity 证据的基因以及 43 个落入正类的基因，得到 10868 个最终负类。

## 2. 当前最准确的 stepwise 构建摘要

1. 证据源收集
   读取 PHI 镜像证据表、repo-local yeast-transfer 表和历史 old440 gene list。

2. ID 规范化与 canonical bridge
   将 PH-1 `XP_*` 蛋白 accession 显式映射到 `fgraminearum::FGRAMPH1_*`；这是当前 newlabel 重建中最关键的新 protocolized 层。

3. lethal 正类构建
   从 `master_evidence_table.preliminary.tsv` 中选取 `species == fgraminearum`、`evidence_source == phi-base_current.csv`、`evidence_term_raw == lethal`、`supports_gold_label == true` 且最终 canonical ID 可解析的基因，生成 `lethal_positive_gene_list.tsv`；当前计数为 77。

4. yeast-transfer 正类构建
   从 `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` 中取 `weak_positive_confidence == high` 的 PH-1 `XP_*` 行，经 bridge 成功解析后按 canonical gene 去重，得到 1045 个高置信 transfer-supported 基因；其中 8 行高置信候选仍 unresolved，被排除在外。

5. newlabel 正类合并
   最终正类定义为 `lethal_set ∪ high_confidence_transfer_set`。当前计数分解为：
   `52` lethal-only，`25` lethal+transfer overlap，`1020` transfer-only，总计 `1097`。

6. newlabel 负类构建
   从 yeast-transfer 表中取 `weak_positive_confidence == none` 且 bridge 成功解析的 canonical genes，形成 11506 基因的 none pool；再去掉 `595` 个 virulence/pathogenicity 相关基因和 `43` 个落入正类的基因，得到 `10868` 个最终负类。

7. oldlabel 回放
   `oldlabel` 不是当前 mainline，而是由 `gene_list.txt` 回放得到的 historical replay：正类来自 `old440_mapping_audit.tsv` 中 `Target == 1` 的已解析基因，共 `439`；负类复用 protocolized negative pool 后去除与 oldlabel 正类重叠的 `118` 个基因，得到 `10750`。

8. processed materialization 与 benchmark 关系
   mainline 构建结果写入 `data/processed/essential_gene/fgraminearum/newlabel/`；随后 `configs/frozen_protocol.yaml` 通过 `positive_genes.tsv` 与 `negative_genes.tsv` 冻结为 benchmark 实际消费的 `results/frozen_protocol/labels/fgraminearum_newlabel.tsv` 和 `results/frozen_protocol/splits/fgraminearum_newlabel_split.tsv`。

## 3. 当前 authoritative files 清单

### A. 历史中间产物，仅用于回顾，不应作为当前 newlabel 的 source of truth

- `results/label_rebuild_experiments/labels/positive_set_P1.tsv`
- `results/label_rebuild_experiments/labels/negative_set.tsv`
- `results/label_rebuild_experiments/old440/labels/positive_old440.tsv`
- `results/label_rebuild_experiments/old440/labels/negative_old440.tsv`
- `docs/protocol_refactor/fgraminearum_newlabel_provenance_reaudit.md`
- `docs/protocol_refactor/fgraminearum_newlabel_count_gap_audit.md`

这些文件对理解历史构建很重要，但当前 benchmark 已不再直接依赖它们。

### B. 当前 protocolized inputs / construction files

- `configs/fgraminearum_label_materialization.yaml`
- `workflow/fgraminearum_label_materialization.smk`
- `src/data/build_fgraminearum_newlabel_bridge.py`
- `src/data/prepare_fgraminearum_label_materialization_sources.py`
- `src/data/materialize_fgraminearum_label_regimes.py`
- `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
- `data/derived_labels/proteome_manifest.tsv`
- `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/old440_mapping_audit.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/high_confidence_yeast_transfer_candidates.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/unresolved_high_confidence_ids.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/bridge_summary.tsv`

### C. 当前 processed labels，且是 frozen benchmark 的上游 authoritative outputs

- `data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/negative_genes.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/labels.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/summary.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/source_manifest.tsv`
- `data/processed/essential_gene/fgraminearum/newlabel/label_construction_audit.tsv`
- `data/processed/essential_gene/fgraminearum/oldlabel/positive_genes.tsv`
- `data/processed/essential_gene/fgraminearum/oldlabel/negative_genes.tsv`
- `data/processed/essential_gene/fgraminearum/label_regime_comparison.tsv`

### D. benchmark 实际消费的冻结文件

- `results/frozen_protocol/labels/fgraminearum_newlabel.tsv`
- `results/frozen_protocol/splits/fgraminearum_newlabel_split.tsv`
- `results/frozen_protocol/labels/fgraminearum_oldlabel.tsv`
- `results/frozen_protocol/splits/fgraminearum_oldlabel_split.tsv`
- `results/Figure2a/Figure2a_old440_vs_newlabel_comparison.tsv`
- `results/Figure2a/Figure2a_best_combo_selection.tsv`

## 4. Figure 2-ready narrative

### Panel 2a: label sources

用 `newlabel` 正类的三段式构成直接说明来源：`52` 个 lethal-only、`25` 个 lethal+transfer overlap、`1020` 个 transfer-only。这个 panel 的核心信息是：newlabel 的主体不是旧 440 label 的简单扩容，而是“PHI lethal anchor + protocolized yeast-transfer expansion”。

### Panel 2b: construction / filtering logic

用 step-count 流程图展示从原始 transfer 候选到最终标签的过滤逻辑。高置信 yeast-transfer 候选共有 `1056` 行，其中 `1048` 行 bridge 成功解析，对应 `1045` 个唯一 canonical genes，另有 `8` 行 unresolved 被剔除。负类从 `11506` 个 bridged none genes 开始，经 virulence/pathogenicity 排除和 positive exclusion 后得到 `10868` 个最终负类。这个 panel 需要突出“bridge 显式化”和“virulence/positive exclusion”两个过滤节点。

### Panel 2c: label regime comparison and modeling consequence

先对比 regime 规模差异：`oldlabel` 为 `439` positive / `10750` negative，`newlabel` 为 `1097` positive / `10868` negative。然后用 Figure2a 已有 benchmark 结果连接下游后果：在 matched model-feature settings 上，`newlabel - oldlabel` 的平均 AUPRC 提升约 `0.245`，平均 MCC 提升约 `0.187`；若按每个模型各自最佳组合比较，`newlabel` 在 4 个 GNN backbone 上均优于 `oldlabel`。

## 5. 推荐的 2–3 个图

### 图 A: newlabel 正类来源构成

- 推荐图型：stacked bar 或 3-segment horizontal bar
- 直接驱动表：`docs/figure2/fgraminearum_newlabel_figure2_plot_data.tsv`
- 关键字段：`figure_panel == 2a_source_composition`
- 核心 message：newlabel 正类由 lethal anchor 与 yeast-transfer expansion 组成，且 transfer-only 是主要增量来源
- 适合位置：主文

### 图 B: label construction / filtering pipeline

- 推荐图型：stepwise count plot 或 provenance flow diagram
- 直接驱动表：`docs/figure2/fgraminearum_newlabel_figure2_plot_data.tsv`
- 关键字段：`figure_panel == 2b_pipeline`
- 核心 message：从 `XP_*` 候选到 canonical label 的 bridge、去重、virulence exclusion 和 positive exclusion 是 newlabel 可复现的关键步骤
- 适合位置：主文

### 图 C: oldlabel vs newlabel 及建模后果

- 推荐图型：双轴小面板，左侧 regime size bar，右侧 best-combo AUPRC/MCC bar 或 dot plot
- 直接驱动表：`docs/figure2/fgraminearum_newlabel_figure2_plot_data.tsv`
- 辅助原始表：`results/Figure2a/Figure2a_old440_vs_newlabel_comparison.tsv`、`results/Figure2a/Figure2a_best_combo_selection.tsv`
- 核心 message：标签制度变化不仅扩大了 biologically anchored positive space，也系统性改善了下游 GNN benchmark
- 适合位置：主文；若主文空间有限，可把 per-feature deltas 放补充材料

## 6. 最终推荐的 Figure 2 组合

如果标签构建部分只能给 `2` 个 panel，建议：

- `2a` 来源构成图
- `2b` stepwise pipeline 图

如果能给 `3` 个 panel，最佳主文组合是：

- `2a` newlabel 正类来源构成
- `2b` construction / filtering pipeline with counts
- `2c` oldlabel vs newlabel 规模与 benchmark consequence

这是最适合生信 / AI 交叉论文的组合，因为它先交代 biological provenance，再交代 protocol logic，最后把标签制度变化和性能收益闭环到模型评估。
