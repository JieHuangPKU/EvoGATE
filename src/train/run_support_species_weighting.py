import io
import os
import sys

import pandas as pd
import yaml

from src.train.train_support_graph_baseline import require_epgat_env, load_yaml, set_seed, train_one_model


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def main():
    enable_utf8_stdout()
    require_epgat_env()

    base_cfg = load_yaml("configs/support_graph_baseline.yaml")
    exp_cfg = load_yaml("configs/support_graph_experiments.yaml")
    weight_cfg = load_yaml("configs/support_species_weighting.yaml")

    outdir = "outputs/support_species_weighting"
    pred_dir = os.path.join(outdir, "per_config_predictions")
    os.makedirs(pred_dir, exist_ok=True)

    config_snapshot = {
        "base_cfg": base_cfg,
        "exp_cfg": exp_cfg,
        "weight_cfg": weight_cfg,
    }
    with open(os.path.join(outdir, "config_snapshot.yaml"), "w", encoding="utf-8") as handle:
        yaml.safe_dump(config_snapshot, handle, sort_keys=False)

    print("【开始 support species weighting 实验】")
    print("固定主物种：")
    print("  • human = 1.0")
    print("  • scerevisiae = 1.0")
    print("测试线虫权重：")
    for w in weight_cfg["celegans_weight_grid"]:
        print("  • {}".format(w))
    print("当前模型：")
    print("  • GraphSAGE")
    print("当前特征：")
    print("  • embedding + orthology/conservation + missing masks")

    rows = []
    for ce_weight in weight_cfg["celegans_weight_grid"]:
        run_cfg = dict(base_cfg)
        run_cfg.update(exp_cfg)
        run_cfg["feature_scope"] = dict(weight_cfg["feature_scope"])
        run_cfg["species_sets"] = dict(weight_cfg["species_sets"])
        run_cfg["species_loss_weights"] = {
            "human": float(weight_cfg["core_species_weights"]["human"]),
            "scerevisiae": float(weight_cfg["core_species_weights"]["scerevisiae"]),
            "celegans": float(ce_weight),
        }

        print("【开始训练】celegans_weight = {}".format(ce_weight))
        set_seed(int(run_cfg.get("seed", 20260403)))
        result = train_one_model("GraphSAGE", run_cfg, "default")
        pred_path = os.path.join(pred_dir, "graphsage_celegans_weight_{}.tsv".format(str(ce_weight).replace(".", "_")))
        result["predictions"].to_csv(pred_path, sep="\t", index=False)

        row = {
            "model": "GraphSAGE",
            "celegans_weight": float(ce_weight),
            "accuracy": result["metrics"]["accuracy"],
            "AUROC": result["metrics"]["auroc"],
            "AUPRC": result["metrics"]["auprc"],
            "F1": result["metrics"]["f1"],
            "node_count": result["nodes"],
            "edge_count": result["edges"],
            "run_status": "success",
            "species_scope": result["species_scope"],
            "weighting_mechanism": "training_loss_contribution_by_species",
            "prediction_path": pred_path,
        }
        rows.append(row)
        print("【训练完成】celegans_weight = {}".format(ce_weight))
        print("  • AUROC = {:.4f}".format(float(row["AUROC"])))
        print("  • AUPRC = {:.4f}".format(float(row["AUPRC"])))

    df = pd.DataFrame(rows)
    metrics_path = os.path.join(outdir, "species_weighting_metrics.tsv")
    df.to_csv(metrics_path, sep="\t", index=False)

    by_auroc = df.sort_values(["AUROC", "AUPRC", "F1"], ascending=[False, False, False], kind="stable").iloc[0]
    by_auprc = df.sort_values(["AUPRC", "AUROC", "F1"], ascending=[False, False, False], kind="stable").iloc[0]
    baseline = df[df["celegans_weight"] == 0.0].iloc[0]
    nonzero_best = df[df["celegans_weight"] > 0.0].sort_values(["AUROC", "AUPRC"], ascending=[False, False], kind="stable").iloc[0]
    any_nonzero_beats_baseline = (float(nonzero_best["AUROC"]) > float(baseline["AUROC"])) or (float(nonzero_best["AUPRC"]) > float(baseline["AUPRC"]))

    best_weight = float(by_auroc["celegans_weight"])
    if best_weight == 0.0:
        role = "排除"
        recommendation = "human + scerevisiae only"
    elif best_weight <= 0.25:
        role = "低权重辅助"
        recommendation = "human + scerevisiae + celegans with reduced weight"
    elif best_weight <= 0.75:
        role = "中权重辅助"
        recommendation = "human + scerevisiae + celegans with reduced weight"
    else:
        role = "全权重"
        recommendation = "human + scerevisiae + celegans with reduced weight"

    lines = [
        "# Species Weighting Summary",
        "",
        "## 设定",
        "- 模型：GraphSAGE",
        "- 特征：embedding + orthology/conservation + missing masks",
        "- weighting mechanism：training_loss_contribution_by_species",
        "",
        "## 关键回答",
        "- 是否存在非零 celegans 权重优于 no-celegans baseline：{}".format("是" if any_nonzero_beats_baseline else "否"),
        "- 最佳 AUROC 权重：{}".format(by_auroc["celegans_weight"]),
        "- 最佳 AUPRC 权重：{}".format(by_auprc["celegans_weight"]),
        "- celegans 最合适角色：{}".format(role),
        "- 推荐默认 support-species 设计：{}".format(recommendation),
        "",
        "## 全部结果",
    ]
    for _, row in df.sort_values("celegans_weight", kind="stable").iterrows():
        lines.append(
            "- celegans_weight={}: accuracy={:.4f}, AUROC={:.4f}, AUPRC={:.4f}, F1={:.4f}".format(
                row["celegans_weight"],
                float(row["accuracy"]),
                float(row["AUROC"]),
                float(row["AUPRC"]),
                float(row["F1"]),
            )
        )

    summary_path = os.path.join(outdir, "species_weighting_summary.md")
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    with open("81_species_weighting_results.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    next_lines = [
        "# 82 Next Step After Species Weighting",
        "",
        "建议下一步：",
        "- 固定当前最佳 celegans 权重后，继续做 GraphSAGE 调参与 graph normalization 比较。",
        "- 在确认权重策略后，再考虑 expression / localization 的增量接入。",
        "- 暂不切换到 GAT/GCN 作为主线。",
    ]
    with open("82_next_step_after_species_weighting.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(next_lines))

    print("【support species weighting 实验完成】")
    print("最佳线虫权重：")
    print("  • AUROC 最佳 = {}".format(by_auroc["celegans_weight"]))
    print("  • AUPRC 最佳 = {}".format(by_auprc["celegans_weight"]))
    print("当前结论：")
    print("  • 线虫最合适角色：{}".format(role))


if __name__ == "__main__":
    main()
