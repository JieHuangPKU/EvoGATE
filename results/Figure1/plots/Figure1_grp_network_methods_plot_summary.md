# Figure1 Plot Summary: network_methods

- 组 ID：`grp_network_methods`
- 比较模型：GraphSAGE, N2V_MLP, CC, DC
- 跨目标按 AUROC / AUPRC / MCC / F1 / Recall 的平均表现看，`GraphSAGE` 综合均值最高。
- F. graminearum 新标签相对旧标签整体更强（平均差值 +0.166），变化最明显的是 `F1`。
- 存在明显权衡：`N2V_MLP` 的 Recall 最强，`DC` 的 Specificity 最高，`GraphSAGE` 的 Precision 最好。
- Barplot 仅调整展示层：y 轴下界固定为 0，负均值按 0 绘制，下误差条裁到 0；原始 TSV 数值不做修改。
