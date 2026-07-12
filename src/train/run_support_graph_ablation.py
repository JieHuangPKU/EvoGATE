import io
import json
import os
import sys

import pandas as pd
import yaml
from sklearn.metrics import accuracy_score, roc_auc_score, average_precision_score, f1_score

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
    config = {
        "models": ["GraphSAGE"],
        "species_sets": {
            "conservative": ["human", "scerevisiae"],
            "mixed": ["human", "scerevisiae", "celegans"],
        },
    }
    outdir = "outputs/support_graph_ablation"
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "ablation_config_snapshot.yaml"), "w", encoding="utf-8") as handle:
        yaml.safe_dump({"base": base_cfg, "ablation": config}, handle, sort_keys=False)

    print("【开始 support 物种消融实验】")
    print("对照组：")
    print("  • human + scerevisiae")
    print("实验组：")
    print("  • human + scerevisiae + celegans")
    print("当前模型：")
    print("  • GraphSAGE")

    rows = []
    for set_name in ["conservative", "mixed"]:
        set_seed(int(base_cfg.get("seed", 20260403)))
        result = train_one_model("GraphSAGE", dict(base_cfg, **config), set_name)
        species_label = result["species_scope"]
        pred_path = os.path.join(outdir, "graphsage_{}_predictions.tsv".format(species_label.replace("+", "_")))
        metrics_path = os.path.join(outdir, "graphsage_{}_metrics.tsv".format(species_label.replace("+", "_")))
        result["predictions"].to_csv(pred_path, sep="\t", index=False)
        metrics_df = pd.DataFrame(
            [
                {
                    "model": "GraphSAGE",
                    "species_set": species_label,
                    "accuracy": result["metrics"]["accuracy"],
                    "AUROC": result["metrics"]["auroc"],
                    "AUPRC": result["metrics"]["auprc"],
                    "F1": result["metrics"]["f1"],
                    "node_count": result["nodes"],
                    "edge_count": result["edges"],
                    "label_count": int((result["predictions"]["label"].astype(int) >= 0).sum()),
                    "run_status": "success",
                }
            ]
        )
        metrics_df.to_csv(metrics_path, sep="\t", index=False)
        rows.extend(metrics_df.to_dict(orient="records"))

    all_df = pd.DataFrame(rows)
    all_df.to_csv(os.path.join(outdir, "species_ablation_summary.tsv"), sep="\t", index=False)

    cons = all_df[all_df["species_set"] == "human+scerevisiae"].iloc[0]
    mixed = all_df[all_df["species_set"] == "human+scerevisiae+celegans"].iloc[0]
    if float(mixed["AUROC"]) > float(cons["AUROC"]) and float(mixed["AUPRC"]) > float(cons["AUPRC"]):
        effect = "celegans 带来提升"
    elif float(mixed["AUROC"]) < float(cons["AUROC"]) and float(mixed["AUPRC"]) < float(cons["AUPRC"]):
        effect = "celegans 可能带来噪音"
    else:
        effect = "celegans 无明显提升"

    lines = [
        "# Species Ablation Summary",
        "",
        "## GraphSAGE 消融结果",
        "- human + scerevisiae: AUROC = {:.4f}, AUPRC = {:.4f}, accuracy = {:.4f}, F1 = {:.4f}".format(
            float(cons["AUROC"]), float(cons["AUPRC"]), float(cons["accuracy"]), float(cons["F1"])
        ),
        "- human + scerevisiae + celegans: AUROC = {:.4f}, AUPRC = {:.4f}, accuracy = {:.4f}, F1 = {:.4f}".format(
            float(mixed["AUROC"]), float(mixed["AUPRC"]), float(mixed["accuracy"]), float(mixed["F1"])
        ),
        "",
        "## Interpretation",
        "- 当前 feature scope = embedding only",
        "- celegans 仍为 partial-support-graph，edge_weight = 0.8",
        "- 这一轮只回答 support-species graph baseline 中 celegans 是否值得保留，不涉及 Fusarium transfer。",
        "- 结论：{}".format(effect),
    ]
    with open(os.path.join(outdir, "species_ablation_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    with open("74_species_ablation_results.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    next_lines = [
        "# 75 Next Step After Species Ablation",
        "",
        "建议下一步：",
        "- 如果 celegans 带来提升或无明显负面影响，则保留 celegans 在默认 support graph 训练集合中。",
        "- 下一轮优先做更严谨的 support-graph 调参与 human+scerevisiae vs human+scerevisiae+celegans 的重复验证。",
    ]
    with open("75_next_step_after_species_ablation.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(next_lines))

    print("【GraphSAGE 消融结果】")
    print("  • human + scerevisiae: AUROC = {:.4f}, AUPRC = {:.4f}".format(float(cons["AUROC"]), float(cons["AUPRC"])))
    print("  • human + scerevisiae + celegans: AUROC = {:.4f}, AUPRC = {:.4f}".format(float(mixed["AUROC"]), float(mixed["AUPRC"])))
    print("【support 物种消融实验完成】")
    print("结论：")
    print("  • {}".format(effect))


if __name__ == "__main__":
    main()
