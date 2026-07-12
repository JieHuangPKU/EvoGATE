# EvoGATE 可复现性

_当前复现 contract、已知缺口与发布级重建要求。_

---

## 当前可复现水平

EvoGATE 为 **Partially reproducible**。Frozen labels、splits、config、processed modality、aggregated metrics 和多个 Figure artifact 均存在。由于 software environment、部分上游 producer、部分 source module、Git history 和主要 per-run output tree 缺失，完整端到端复现为 **Blocked**。

## Frozen protocol

| 项目 | Authoritative value | 证据 |
|---|---|---|
| Protocol version | `frozen_protocol_v1` | `configs/frozen_protocol.yaml` |
| Split strategy | Stratified fixed split | `src/data/freeze_unified_protocol.py` |
| Train fraction | 70% | 由 test `0.20` 和 validation `0.10` 得到 |
| Validation fraction | 10% | Frozen config 和 split artifact |
| Test fraction | 20% | Frozen config 和 split artifact |
| Split seed | `20260409` | Frozen config 和 split version |
| Training seeds | `1029`-`1033` | Frozen config 和 result summary |
| Main Fusarium regime | `fgraminearum_newlabel` | Frozen protocol definition |
| Main newlabel counts | 1,097 positive；10,868 negative | Newlabel summary |

所有正式比较必须使用同一 frozen split。正式比较中的 model 不得自行生成私有 split。

## 随机性与模型选择

`src/train/run_frozen_protocol_model.py` 设置 Python、NumPy、PyTorch 和可用 CUDA device 的 seed。Neural checkpoint 在定义有效时按 validation AUPRC 选择，代码中 AUROC 或 accuracy 仅作为 fallback。

五个 training seed 量化同一 frozen split 上的优化随机性，不能替代独立 data split 或 biological replication。现有 summary 通常报告 mean 与 standard deviation，并非所有主要 comparison 都包含 paired confidence interval 或 hypothesis test。

## 阈值策略

标准 trainable-model output 使用固定 `0.5` threshold。Network heuristic 使用与 labeled positive count 相关的 top-k rule。Threshold-tuned Figure3c analysis 从 validation prediction 得到 F1/MCC optimal threshold，再应用于 test prediction。

Test split 不得用于选择 threshold、model、feature combination、epoch 或 hyperparameter。描述性 test-split threshold curve 必须标注为 descriptive，不能定义最终 operating point。

## 数据与特征确定性

Frozen label/split TSV 已物化。Numeric ORT/EXP/SUB/degree feature 与 ESM2 vector 使用 training-node statistics 归一化。ESM2 alignment 在缺少 node embedding 时失败，不会静默丢弃 node。

当前仓库包含 ESM2 `.pt` artifact 和 extraction log，但日常文档检查不应加载大型 binary。精确 ESM2 复现依赖 configured local model checkpoint 和 software environment。

## 环境状态

| 领域 | 当前证据 | 状态 |
|---|---|---|
| Python interpreter | Config/workflow 中的绝对 Python path | Historical / non-portable |
| Conda environment | Config/script 中的 `EPGAT` 和 `progate` 名称 | Historical |
| Python dependencies | Source import | Partially documented |
| R dependencies | 存在 R script，无 lock | Blocked |
| Snakemake version | 未锁定 | Blocked |
| PyTorch/graph backend | 多个 legacy environment/implementation | Blocked |
| ESM2 model | 有 local `esm2_t33_650M_UR50D` path 证据 | Partially documented |
| Hardware | ESM2 config 请求 CUDA；frozen benchmark config 使用 CPU | Partially documented |

没有 environment lock、package manifest、container 或 reproducible installation procedure。不得根据 import 推断 package installation。

## 缺失 artifact

- 仓库内包含 prediction、checkpoint 和 resolved config 的 `outputs/`
- 若干 bytecode-only evaluation module 的 source
- `weak_positive_confidence` 上游 producer
- 非空 Git metadata
- 完全可移植的 raw-data build input 与 command
- Environment lock file

这些缺口阻止 full reproducibility claim。

## 已有复现资产

- Frozen configuration 与 model hyperparameter
- Materialized label 与 split
- Label construction 的 source manifest 和部分 checksum
- Processed modality summary 与 mapping audit
- Per-Figure aggregated metric 与 report
- ESM2 extraction log 与 embedding metadata
- 若干 Figure3 experiment 的 runtime config snapshot

## 发布级要求

1. 恢复或建立 version control，并记录 baseline。
2. 锁定 Python、R、Snakemake、PyTorch、graph 和 plotting dependency。
3. 用经过审核的 portable config 替代 machine-specific runtime root。
4. 恢复 authoritative result 使用的全部 source module。
5. 恢复或重新生成 per-run output，且不覆盖 historical artifact。
6. 恢复 upstream evolutionary-confidence producer 和精确规则。
7. 创建显式 exclusion manifest。
8. 为命名 release 记录 input/output checksum。
9. 对 manuscript claim 进行 contract-matched paired validation。
10. 从命名 release 重建全部 manuscript table 与 Figure。

## 复现 claim 规范

在满足全部 release requirement 前，文档只能写 **Partially reproducible** 或 **Blocked for full reproduction**，不得写 fully reproducible、production-ready 或 independently reproduced。

