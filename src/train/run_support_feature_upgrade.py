import io
import os
import sys

import pandas as pd
import yaml

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


def _write_metrics(path, result, feature_set):
    pd.DataFrame(
        [
            {
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
        ]
    ).to_csv(path, sep="\t", index=False)


def main():
    enable_utf8_stdout()
    require_epgat_env()

    outdir = "outputs/support_graph_features"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("【开始训练 embedding_only】")
    base_cfg = load_yaml("configs/support_graph_baseline.yaml")
    exp_cfg = load_yaml("configs/support_graph_experiments.yaml")
    cfg_embedding = dict(base_cfg)
    cfg_embedding.update(exp_cfg)
    cfg_embedding["feature_scope"] = {
        "embedding": True,
        "expression": False,
        "orthology": False,
        "localization": False,
        "prior_score": False,
    }
    set_seed(int(cfg_embedding.get("seed", 20260403)))
    embedding_result = train_one_model("GraphSAGE", cfg_embedding, "default")
    embedding_result["predictions"].to_csv(
        os.path.join(outdir, "graphsage_embedding_only_predictions.tsv"), sep="\t", index=False
    )
    _write_metrics(
        os.path.join(outdir, "graphsage_embedding_only_metrics.tsv"),
        embedding_result,
        "embedding_only",
    )

    print("【开始训练 embedding + orthology/conservation】")
    cfg_ortho = dict(base_cfg)
    cfg_ortho.update(exp_cfg)
    cfg_ortho["feature_scope"] = {
        "embedding": True,
        "expression": False,
        "orthology": True,
        "localization": False,
        "prior_score": False,
    }
    set_seed(int(cfg_ortho.get("seed", 20260403)))
    ortho_result = train_one_model("GraphSAGE", cfg_ortho, "default")
    ortho_result["predictions"].to_csv(
        os.path.join(outdir, "graphsage_embedding_plus_orthology_predictions.tsv"), sep="\t", index=False
    )
    _write_metrics(
        os.path.join(outdir, "graphsage_embedding_plus_orthology_metrics.tsv"),
        ortho_result,
        "embedding_plus_orthology",
    )

    celegans_matrix = pd.read_csv(
        os.path.join(outdir, "support_feature_matrix_celegans.tsv"), sep="\t", dtype=str
    ).fillna("")
    celegans_cov = pd.to_numeric(celegans_matrix["has_orthology_feature"], errors="coerce").fillna(0.0)
    celegans_coverage = float((celegans_cov > 0).mean()) if len(celegans_matrix) else 0.0
    celegans_missing = 1.0 - celegans_coverage

    embed_metrics = pd.read_csv(os.path.join(outdir, "graphsage_embedding_only_metrics.tsv"), sep="\t")
    ortho_metrics = pd.read_csv(os.path.join(outdir, "graphsage_embedding_plus_orthology_metrics.tsv"), sep="\t")
    merged = pd.concat([embed_metrics, ortho_metrics], ignore_index=True)
    merged.to_csv(os.path.join(outdir, "feature_upgrade_summary.tsv"), sep="\t", index=False)

    auroc_delta = float(ortho_result["metrics"]["auroc"]) - float(embedding_result["metrics"]["auroc"])
    auprc_delta = float(ortho_result["metrics"]["auprc"]) - float(embedding_result["metrics"]["auprc"])
    if auroc_delta > 0 and auprc_delta > 0:
        celegans_effect = "在加入 orthology 后，celegans 不再表现出明确拖累，值得继续保留并观察。"
    elif auroc_delta < 0 and auprc_delta < 0:
        celegans_effect = "即使加入 orthology 后，默认三物种集合仍未优于 embedding_only；celegans 是否拖累仍需进一步做 species weighting 或 graph normalization 检查。"
    else:
        celegans_effect = "加入 orthology 后结果混合，celegans 是否拖累仍不确定。"

    lines = [
        "# Feature Upgrade Summary",
        "",
        "## GraphSAGE 比较",
        "- embedding_only: accuracy = {:.4f}, AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}".format(
            float(embedding_result["metrics"]["accuracy"]),
            float(embedding_result["metrics"]["auroc"]),
            float(embedding_result["metrics"]["auprc"]),
            float(embedding_result["metrics"]["f1"]),
        ),
        "- embedding + orthology/conservation: accuracy = {:.4f}, AUROC = {:.4f}, AUPRC = {:.4f}, F1 = {:.4f}".format(
            float(ortho_result["metrics"]["accuracy"]),
            float(ortho_result["metrics"]["auroc"]),
            float(ortho_result["metrics"]["auprc"]),
            float(ortho_result["metrics"]["f1"]),
        ),
        "",
        "## celegans orthology 覆盖",
        "- coverage = {:.4f}".format(celegans_coverage),
        "- missing_rate = {:.4f}".format(celegans_missing),
        "",
        "## 回答",
        "- celegans orthology features coverage = {:.4f}".format(celegans_coverage),
        "- celegans missing rate = {:.4f}".format(celegans_missing),
        "- 加入 orthology 后 AUROC delta = {:.4f}".format(auroc_delta),
        "- 加入 orthology 后 AUPRC delta = {:.4f}".format(auprc_delta),
        "- 判断：{}".format(celegans_effect),
        "- 如果 orthology 升级后仍然波动，下一步优先检查 feature weighting / species weighting / graph normalization，而不是移除 celegans。",
    ]
    with open(os.path.join(outdir, "feature_upgrade_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    print("【本轮完成】")
    print("必须输出：")
    print("  • 是否成功接入 orthology/conservation：成功")
    print("  • celegans 是否仍拖累：{}".format("是" if auroc_delta < 0 and auprc_delta < 0 else "待定/不明显"))
    print("  • 下一步建议：feature weighting / species weighting / graph normalization")


if __name__ == "__main__":
    main()
