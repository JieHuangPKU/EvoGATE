import io
import os
import sys

import pandas as pd

from src.train.train_support_graph_baseline import (
    require_epgat_env,
    load_yaml,
    set_seed,
    train_one_model,
)


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def _write_metrics(path, result, species_set_label):
    pd.DataFrame(
        [
            {
                "model": "GraphSAGE",
                "species_set": species_set_label,
                "accuracy": result["metrics"]["accuracy"],
                "AUROC": result["metrics"]["auroc"],
                "AUPRC": result["metrics"]["auprc"],
                "F1": result["metrics"]["f1"],
                "node_count": result["nodes"],
                "edge_count": result["edges"],
                "run_status": "success",
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def main():
    enable_utf8_stdout()
    require_epgat_env()

    outdir = "outputs/support_graph_feature_ablation_v2"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    base_cfg = load_yaml("configs/support_graph_baseline.yaml")
    exp_cfg = load_yaml("configs/support_graph_experiments.yaml")
    config = dict(base_cfg)
    config.update(exp_cfg)
    config["feature_scope"] = {
        "embedding": True,
        "expression": False,
        "orthology": True,
        "localization": False,
        "prior_score": False,
    }
    config["species_sets"] = {
        "reduced": ["human", "scerevisiae"],
        "full": ["human", "scerevisiae", "celegans"],
    }

    print("【开始第二轮 species ablation（已接入 orthology 特征）】")
    print("配置1：")
    print("  • human + scerevisiae")
    print("配置2：")
    print("  • human + scerevisiae + celegans")

    results = []
    for set_name in ["reduced", "full"]:
        set_seed(int(config.get("seed", 20260403)))
        result = train_one_model("GraphSAGE", config, set_name)
        species_label = result["species_scope"]
        pred_path = os.path.join(
            outdir,
            "graphsage_{}_predictions.tsv".format(species_label.replace("+", "_")),
        )
        metrics_path = os.path.join(
            outdir,
            "graphsage_{}_metrics.tsv".format(species_label.replace("+", "_")),
        )
        result["predictions"].to_csv(pred_path, sep="\t", index=False)
        _write_metrics(metrics_path, result, species_label)
        results.append(
            {
                "model": "GraphSAGE",
                "species_set": species_label,
                "accuracy": result["metrics"]["accuracy"],
                "AUROC": result["metrics"]["auroc"],
                "AUPRC": result["metrics"]["auprc"],
                "F1": result["metrics"]["f1"],
                "node_count": result["nodes"],
                "edge_count": result["edges"],
                "run_status": "success",
            }
        )

    summary_df = pd.DataFrame(results)
    summary_df.to_csv(os.path.join(outdir, "feature_ablation_v2_summary.tsv"), sep="\t", index=False)

    reduced = summary_df[summary_df["species_set"] == "human+scerevisiae"].iloc[0]
    full = summary_df[summary_df["species_set"] == "human+scerevisiae+celegans"].iloc[0]
    delta_auroc = float(full["AUROC"]) - float(reduced["AUROC"])
    delta_auprc = float(full["AUPRC"]) - float(reduced["AUPRC"])
    delta_f1 = float(full["F1"]) - float(reduced["F1"])

    if delta_auroc > 0.01 and delta_auprc > 0.01:
        case = "CASE A: Positive contributor"
        cn_case = "正向"
        recommendation = "keep"
    elif delta_auroc < -0.01 and delta_auprc < -0.01:
        case = "CASE C: Negative contributor"
        cn_case = "负向"
        recommendation = "downweight"
    else:
        case = "CASE B: Neutral contributor"
        cn_case = "中性"
        recommendation = "auxiliary only"

    explanation = (
        "celegans 的 orthology coverage 只有 0.4587，missing rate 为 0.5413；"
        "这意味着即使加入了真实生物学特征，它仍然存在较高的特征稀疏性。"
        "最终表现取决于 sparse orthology signal 是否能抵消 partial-support graph 的结构噪音和跨物种不匹配。"
    )

    lines = [
        "# Feature Ablation V2 Summary",
        "",
        "## GraphSAGE 指标",
        "- human + scerevisiae: accuracy = {:.4f}, AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}".format(
            float(reduced["accuracy"]),
            float(reduced["AUROC"]),
            float(reduced["AUPRC"]),
            float(reduced["F1"]),
        ),
        "- human + scerevisiae + celegans: accuracy = {:.4f}, AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}".format(
            float(full["accuracy"]),
            float(full["AUROC"]),
            float(full["AUPRC"]),
            float(full["F1"]),
        ),
        "",
        "## Delta",
        "- delta_AUROC = {:.4f}".format(delta_auroc),
        "- delta_AUPRC = {:.4f}".format(delta_auprc),
        "- delta_F1 = {:.4f}".format(delta_f1),
        "",
        "## Classification",
        "- {}".format(case),
        "- celegans 贡献类型：{}".format(cn_case),
        "",
        "## Biological Explanation",
        "- {}".format(explanation),
        "- 需要同时考虑 feature sparsity impact、graph noise vs signal、以及 cross-species transfer mismatch。",
        "",
        "## Recommendation",
        "- {}".format(recommendation),
    ]
    with open(os.path.join(outdir, "feature_ablation_v2_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    with open("79_feature_upgraded_species_ablation.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    next_lines = [
        "# 80 Next Step After Feature-Upgraded Ablation",
        "",
        "建议下一步：",
        "- 如果 celegans 为正向贡献者，继续保留在默认 support graph 中，并做 GraphSAGE 调参与 species weighting 微调。",
        "- 如果 celegans 为中性或负向贡献者，优先测试 species weighting / graph normalization，而不是立即删除其图。",
        "- 继续保持 feature schema 不变，后续再考虑 expression / localization 的增量接入。",
    ]
    with open("80_next_step_after_feature_upgraded_ablation.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(next_lines))

    print("【训练完成】")
    print("【结果对比】")
    print("  • human + scerevisiae: AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}".format(
        float(reduced["AUROC"]), float(reduced["AUPRC"]), float(reduced["F1"])
    ))
    print("  • human + scerevisiae + celegans: AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}".format(
        float(full["AUROC"]), float(full["AUPRC"]), float(full["F1"])
    ))
    print("【结论】")
    print("  • celegans 是否有贡献：{}".format("是" if recommendation == "keep" else ("不明显" if recommendation == "auxiliary only" else "当前偏负面")))
    print("  • 类型：{}".format(cn_case))


if __name__ == "__main__":
    main()
