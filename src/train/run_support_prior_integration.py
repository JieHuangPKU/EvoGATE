import io
import os
import sys

import pandas as pd

from src.train.train_support_graph_baseline import require_epgat_env, load_yaml, set_seed, train_one_model


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def _metrics_row(result, feature_set):
    return {
        "model": "GraphSAGE",
        "feature_set": feature_set,
        "species_scope": result["species_scope"],
        "accuracy": result["metrics"]["accuracy"],
        "AUROC": result["metrics"]["auroc"],
        "AUPRC": result["metrics"]["auprc"],
        "F1": result["metrics"]["f1"],
        "node_count": result["nodes"],
        "edge_count": result["edges"],
        "run_status": "success",
    }


def main():
    enable_utf8_stdout()
    require_epgat_env()

    outdir = "outputs/support_graph_prior"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    base_cfg = load_yaml("configs/support_graph_baseline.yaml")
    exp_cfg = load_yaml("configs/support_graph_experiments.yaml")

    common = dict(base_cfg)
    common.update(exp_cfg)
    common["species_sets"] = {"default": ["human", "scerevisiae", "celegans"]}
    common["species_loss_weights"] = {"human": 1.0, "scerevisiae": 1.0, "celegans": 1.0}

    print("【开始训练】不含 prior_score")
    cfg_without = dict(common)
    cfg_without["feature_scope"] = {
        "embedding": True,
        "expression": False,
        "orthology": True,
        "localization": False,
        "prior_score": False,
    }
    set_seed(int(cfg_without.get("seed", 20260403)))
    without_result = train_one_model("GraphSAGE", cfg_without, "default")
    without_result["predictions"].to_csv(os.path.join(outdir, "graphsage_without_prior_predictions.tsv"), sep="\t", index=False)
    pd.DataFrame([_metrics_row(without_result, "embedding_plus_orthology")]).to_csv(
        os.path.join(outdir, "graphsage_without_prior_metrics.tsv"), sep="\t", index=False
    )

    print("【开始训练】含 prior_score")
    cfg_with = dict(common)
    cfg_with["feature_scope"] = {
        "embedding": True,
        "expression": False,
        "orthology": True,
        "localization": False,
        "prior_score": True,
    }
    set_seed(int(cfg_with.get("seed", 20260403)))
    with_result = train_one_model("GraphSAGE", cfg_with, "default")
    with_result["predictions"].to_csv(os.path.join(outdir, "graphsage_with_prior_predictions.tsv"), sep="\t", index=False)
    pd.DataFrame([_metrics_row(with_result, "embedding_plus_orthology_plus_prior")]).to_csv(
        os.path.join(outdir, "graphsage_with_prior_metrics.tsv"), sep="\t", index=False
    )

    summary = pd.DataFrame([
        _metrics_row(without_result, "embedding_plus_orthology"),
        _metrics_row(with_result, "embedding_plus_orthology_plus_prior"),
    ])
    summary.to_csv(os.path.join(outdir, "prior_feature_upgrade_summary.tsv"), sep="\t", index=False)

    auroc_delta = float(with_result["metrics"]["auroc"]) - float(without_result["metrics"]["auroc"])
    auprc_delta = float(with_result["metrics"]["auprc"]) - float(without_result["metrics"]["auprc"])
    f1_delta = float(with_result["metrics"]["f1"]) - float(without_result["metrics"]["f1"])
    keep_prior = (auroc_delta > 0) or (auprc_delta > 0)

    lines = [
        "# Prior Feature Upgrade Summary",
        "",
        "## 结果汇总",
        "- 不含 prior_score: AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}, accuracy = {:.4f}".format(
            float(without_result["metrics"]["auroc"]),
            float(without_result["metrics"]["auprc"]),
            float(without_result["metrics"]["f1"]),
            float(without_result["metrics"]["accuracy"]),
        ),
        "- 含 prior_score: AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}, accuracy = {:.4f}".format(
            float(with_result["metrics"]["auroc"]),
            float(with_result["metrics"]["auprc"]),
            float(with_result["metrics"]["f1"]),
            float(with_result["metrics"]["accuracy"]),
        ),
        "",
        "## 回答",
        "- 是否找到安全 numeric prior_score：是",
        "- prior_score 来源：outputs/support_graphs/*_edges_for_training.tsv 的 raw_edge_weight，经 per-node weighted degree sum -> log1p -> per-species min-max normalization 得到。",
        "- AUROC delta = {:.4f}".format(auroc_delta),
        "- AUPRC delta = {:.4f}".format(auprc_delta),
        "- F1 delta = {:.4f}".format(f1_delta),
        "- prior_score 是否值得保留：{}".format("是" if keep_prior else "否"),
    ]
    with open(os.path.join(outdir, "prior_feature_upgrade_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    with open("84_support_prior_integration_results.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    next_lines = [
        "# 85 Next Step After Support Prior Integration",
        "",
        "建议下一步：",
        "- 如果 prior_score 有增益，固定当前 GraphSAGE + embedding + orthology + prior 作为 support baseline。",
        "- 然后优先比较 graph normalization / dropout / hidden_dim 等轻量调参。",
        "- 暂不切换主模型，也不在本轮加入 expression / localization。",
    ]
    with open("85_next_step_after_support_prior_integration.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(next_lines))

    print("【support prior 特征集成完成】")
    print("结果汇总：")
    print("  • 不含 prior_score: AUROC = {:.4f}, AUPRC = {:.4f}".format(
        float(without_result["metrics"]["auroc"]), float(without_result["metrics"]["auprc"])
    ))
    print("  • 含 prior_score: AUROC = {:.4f}, AUPRC = {:.4f}".format(
        float(with_result["metrics"]["auroc"]), float(with_result["metrics"]["auprc"])
    ))
    print("当前结论：")
    print("  • prior_score 是否值得保留：{}".format("是" if keep_prior else "否"))


if __name__ == "__main__":
    main()
