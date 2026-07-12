# EvoGATE 迁移指南

_未来从 ProGATE_v2 到 EvoGATE 的非破坏性迁移指导。_

---

## 范围

本文记录已知 migration 工作，不授权 path replacement、file movement、rename、symlink creation 或 result regeneration。当前文档阶段不执行迁移。

## 身份映射

| 项目 | 历史值 | 当前值 | 规范 |
|---|---|---|---|
| Project name | ProGATE / ProGATE_v2 | EvoGATE | 新文档和 claim 使用 EvoGATE |
| Method predecessor | EPGAT | EPGAT | 保留名称；不得改称 EvoGATE |
| External comparator/source | Bingo | Bingo | 保留名称和 external status |
| Historical Linux root | `/home/jiehuang/software/fungi/ProGATE_v2` | `/DATA/software/bioinfo/fungi/EvoGATE` | 先报告；不得批量替换 |
| Historical macOS root | `/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2` | 当前 repository root | 未来按 producer 逐个迁移 |

## 已知 hard-coded path 类型

### Shell wrapper

多个 `scripts/run_*.sh` wrapper 指定 `PROJECT_ROOT=/home/jiehuang/software/fungi/ProGATE_v2`，包括 Figure1、Figure2、Figure3、Figure4、label-scarcity 和 label-materialization wrapper。它们是 **Historical / non-portable**。

### Runtime config

`configs/frozen_protocol.yaml` 含历史 `mplconfigdir`、cache、Python environment 和 local ESM2 model path。这些路径记录曾使用的环境，但不是 portable default。

### Data builder

若干 `src/data/` builder 含用于 raw source、EPGAT asset 和 Fusarium evidence 的历史 macOS 或 external Linux path。现有 processed artifact 可使用，但其 exact raw rebuild 仍被阻断。

### Generated report

Result summary 和 migration document 含历史 absolute output path。若这些路径记录 provenance，应保持不变。新文档应使用 repository-relative path。

这些 generated path record 属于 **Historical** provenance。

## 历史入口

| Historical entry type | 当前解释 |
|---|---|
| 带 absolute root 的 ProGATE_v2 Shell wrapper | Non-portable；不推荐 |
| EPGAT replay runner | 仅用于 Historical comparison |
| Old ranking-only Fusarium workflow | Historical；不是 mainline frozen protocol |
| 使用 repository-relative input 的 Figure workflow | Current experiment entry，但受 missing outputs/environment 限制 |
| `workflow/frozen_protocol_benchmark.smk` | 当前 canonical benchmark workflow |
| `workflow/fgraminearum_label_materialization.smk` | 当前 canonical label workflow |

以上 legacy entry inventory 的用途是 **Historical** 追溯，而不是当前推荐命令。

## 未来迁移原则

1. 结构迁移前恢复或建立 version control。
2. 改变 path 前 inventory 所有 reference。
3. 分离 provenance text 与 executable config。
4. 每次只替换一类 executable path，并保留 validation record。
5. 保持 historical artifact immutable。
6. 不用 root-level symlink 重定向旧路径。
7. 不将 path migration 与 scientific protocol change 合并。
8. 将 oldlabel 和 legacy replay 保留为明确的 non-default branch。
9. 新运行写入 versioned output/result directory。
10. Deprecate old entry 前验证全部 consumer。

## 建议迁移顺序

| 阶段 | 变化 | 所需证据 |
|---|---|---|
| 1 | 建立 Git baseline 与 authoritative file inventory | Repository history 或 signed baseline manifest |
| 2 | 定义 tested environment | Locked Python/R/Snakemake dependencies |
| 3 | 创建 portable root-relative canonical command | Dry-run 与 lightweight check |
| 4 | 恢复 missing source 与 outputs | Checksum 与 provenance mapping |
| 5 | 将 data builder 从 absolute raw path 迁出 | Per-modality rebuild comparison |
| 6 | 合并 workflow entry point | Target/output equivalence audit |
| 7 | Archive deprecated wrapper | 确认没有 current consumer |

## 迁移停止条件

如果 path change 改变 sample count、label count、ID mapping、split、seed、graph edge set、feature dimension、metric 或 output destination，应停止并请求审核。这类变化属于科学或行为变化，不只是 path migration。

## 当前 blocker

Migration validation 因 empty Git metadata、missing environment lock、missing per-run `outputs/`、bytecode-only module 和 missing yeast-transfer confidence producer 而 **Blocked**。应在独立批准的工作中解决。
