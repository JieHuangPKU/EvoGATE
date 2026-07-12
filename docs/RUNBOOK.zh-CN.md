# EvoGATE 运行手册

_只读检查、轻量验证、规范入口、阻断路径与执行安全。_

---

## 操作边界

本仓库不是可移植软件发布。运行前必须检查。未经用户明确批准，不运行 training、Snakemake workflow、ESM2 extraction、Figure regeneration 或大型数据处理。除非明确授权，不写入 `data/` 或现有 `results/`。

## 只读检查

```bash
find . -maxdepth 2 -type d
find configs docs scripts src workflow -maxdepth 2 -type f
du -sh data results data/* results/*
rg -n "ProGATE_v2|EPGAT|Bingo" configs docs scripts src workflow
head -n 5 data/processed/essential_gene/fgraminearum/newlabel/summary.tsv
head -n 5 results/Figure3a/data/Figure3a_final_summary.tsv
```

这些命令不会验证科学正确性，只检查仓库状态和小型文本 artifact。

## 轻量检查

以下检查仅解析 source/config，不启动 pipeline：

```bash
python -m src.data.freeze_unified_protocol --help
python -m src.train.run_frozen_protocol_model --help
python -m src.features.extract_esm2_pooled --help
python -m src.eval.build_figure5_candidate_prioritization --help
```

若 active environment 缺少 graph 或 machine-learning dependency，import check 可能失败。此类失败说明环境不完整；未经批准不得通过安装 package 自行修复。

## 规范入口

下列是 canonical code path，不代表可以不经批准执行。

| 任务 | 入口 | 预期输出 | 状态 |
|---|---|---|---|
| 物化 Fusarium label | `workflow/fgraminearum_label_materialization.smk` | `data/processed/essential_gene/fgraminearum/{oldlabel,newlabel}/` | Partially reproducible |
| Freeze label/split | `python -m src.data.freeze_unified_protocol --config configs/frozen_protocol.yaml` | `results/frozen_protocol/labels/`、`results/frozen_protocol/splits/` | Implemented；会写 frozen result |
| 运行单个 model task | `python -m src.train.run_frozen_protocol_model ...` | Per-run output directory | Implemented；需要完整环境 |
| 完整 frozen benchmark | `workflow/frozen_protocol_benchmark.smk` | `outputs/Figure1/`、`results/Figure1/` | 大型任务；需要批准 |
| 准备 ESM2 cache | `workflow/prepare_esm2_cache.smk` | `data/processed/ESM2/<species>/esm2_pooled.pt` | 大型任务；需要批准 |
| Candidate prioritization | `python -m src.eval.build_figure5_candidate_prioritization` | `results/Figure5_new_candidate_prioritization/` | 当前 Blocked |

## Dry-run 指南

Snakemake dry-run 是 workflow 的首选 preflight，但仍可能解析 configured path 与 import。当前项目规则排除了 workflow execution，因此 dry-run 也需要明确批准。

获批后，应从 repository root 直接使用 workflow file，而不是 historical wrapper。执行前确认所有 planned output 位于新目录。

## 历史且不可移植的入口

`scripts/run_Figure1_frozen_protocol_benchmark.sh`、`scripts/run_fgraminearum_label_materialization.sh`、`scripts/run_label_scarcity_benchmark.sh` 及多个 Figure2-Figure4 wrapper 硬编码 `/home/jiehuang/software/fungi/ProGATE_v2` 或特定机器环境。

状态：**Historical / non-portable**。不要将它们作为推荐 EvoGATE command。发现路径后先报告，不得批量替换。

部分 Figure5 wrapper 根据自身位置计算 repository root，结构上更可移植，但其上游 `outputs/` dependency 缺失。

这些旧执行约定整体属于 **Historical** 运行环境。

## 当前被阻断的操作

| 操作 | Blocker |
|---|---|
| 重建全部 published summary | 仓库内缺少 `outputs/` |
| 重建精确 software environment | 没有 environment lock 或 package manifest |
| 审计全部 evaluation implementation | 部分 evaluation source 缺失但 `.pyc` 保留 |
| 重建 yeast-transfer confidence | 上游 confidence generator 缺失 |
| 使用 Git history 追溯 provenance | `.git/` 为空 |
| 从当前 root 运行 historical wrapper | 硬编码 ProGATE_v2/macOS path |

## 需要确认的大型任务

- 完整 frozen benchmark
- 任何 Figure benchmark 或 ablation workflow
- 写入 `data/processed/` 的 label materialization
- ESM2 extraction
- Graph reconstruction 或 threshold sweep
- 覆盖 existing result 的 candidate regeneration
- 使用 GPU、多 CPU core 或大型 matrix 的任务

批准前应报告 exact command、input、output、overwrite behavior、environment、estimated scale 和 recovery plan。

## Preflight checklist

1. 确认当前工作目录为 EvoGATE root。
2. 阅读 `AGENTS.md`、`docs/INCONSISTENCIES.md` 和相关 config。
3. 验证所有 input path 存在。
4. 验证 output directory 是新的独立目录。
5. 记录 protocol version、split version、split seed 和 training seed。
6. 确认没有根据 test split 优化。
7. 确认 command 不删除或覆盖已有 artifact。
8. 为大型任务取得明确批准。

## 常见失败模式

| 症状 | 可能原因 | 响应 |
|---|---|---|
| Shell wrapper 切换到不存在目录 | Historical `PROJECT_ROOT` | 停止并报告；不得批量替换 |
| 找不到 ESM2 cache | 缺少 local model/cache path | 报告 configured path；不自动下载 |
| Graph import failure | 缺少 DGL/PyG/torch-scatter environment | 记录环境失败；不自动安装 |
| Candidate builder 找不到 prediction | 缺少 `outputs/Figure3a/` | 视为 Blocked |
| Summary 之间 metric 不同 | Conflicting historical artifact | 记录在 `INCONSISTENCIES.md`；不选择偏好值 |
| Snakemake 提议删除或覆盖 | Existing workflow cleanup logic | 执行前停止并请求确认 |

## 运行记录

对任何获批的未来运行，保留 command、timestamp、working directory、environment description、resolved config、input checksum、output directory、protocol version、split version、seed、exit status 和 log。缺少这些字段的结果不应成为 authoritative manuscript evidence。
