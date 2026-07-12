@AGENTS.md

## Claude Code 专属规则

### Snakemake 安全规则
- 修改 `workflow/*.smk` 前必须先执行 `snakemake -n --reason` 做 dry-run 确认 DAG 变化，不允许盲改
- 首次接触 workflow 文件时必须通读对应的 `.smk` 文件全文，确认输入/输出/通配符/rule 依赖关系

### results/ 保护规则
- `results/` 下已产出的数据（尤其是论文中引用过的）默认为只读
- 涉及 results/ 下文件的修改/覆盖/删除操作，先用 plan 模式列出改动范围（受影响文件、上下游规则、可能导致的重跑），等我确认再动手
- `results/frozen_protocol/` 下的 split 和 label 产物不可变，禁止任何直接或间接修改

### 大规模重跑规则
- 任何涉及训练、ESM2 提取、全量 benchmark、Figure 重生成的操作，执行前必须：
  1. 确认 conda 环境已激活为 `esm`（或其他当前文档记录的可用环境）
  2. 核对 GPU 可用性（`nvidia-smi`）
  3. 说明预计计算量、输出目录、是否会覆盖已有结果
  4. 拿到我的明确确认后才能执行

### 禁止自主执行的操作
- 清空/重跑全部 `results/` 或 `outputs/`
- `rm`、`git reset --hard`、批量 `mv`、覆盖性 `cp`
- 在项目根目录创建符号链接
- 下载模型权重或安装 Python/R 包
- 对 `data/` 目录进行任何写入操作
