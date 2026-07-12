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
    tuning_cfg = load_yaml("configs/graphsage_tuning.yaml")
    outdir = "outputs/support_graph_tuning"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("【开始 GraphSAGE 调参与 Fusarium 推理输入构建】")
    print("当前固定 support 设计：")
    print("  • human")
    print("  • scerevisiae")
    print("  • celegans")
    print("当前固定特征：")
    print("  • embedding")
    print("  • orthology/conservation")
    print("  • true prior")
    print("  • missing masks")
    print("当前主模型：")
    print("  • GraphSAGE")

    rows = []
    for idx, grid in enumerate(tuning_cfg["tuning_grid"], start=1):
        cfg = dict(base_cfg)
        cfg["species_sets"] = dict(tuning_cfg["species_sets"])
        cfg["feature_scope"] = dict(tuning_cfg["feature_scope"])
        cfg["species_loss_weights"] = dict(tuning_cfg["species_loss_weights"])
        cfg["hidden_dim"] = int(grid["hidden_dim"])
        cfg["num_layers"] = int(grid["num_layers"])
        cfg["dropout"] = float(grid["dropout"])
        cfg["feature_normalization"] = bool(grid["feature_normalization"])
        print("【开始 GraphSAGE 调参】")
        print("  • hidden_dim = {}".format(cfg["hidden_dim"]))
        print("  • num_layers = {}".format(cfg["num_layers"]))
        print("  • dropout = {}".format(cfg["dropout"]))
        print("  • normalization = {}".format("on" if cfg["feature_normalization"] else "off"))
        set_seed(int(cfg.get("seed", 20260403)))
        try:
            result = train_one_model("GraphSAGE", cfg, "default")
            rows.append(
                {
                    "config_id": idx,
                    "hidden_dim": cfg["hidden_dim"],
                    "num_layers": cfg["num_layers"],
                    "dropout": cfg["dropout"],
                    "normalization": "on" if cfg["feature_normalization"] else "off",
                    "accuracy": result["metrics"]["accuracy"],
                    "AUROC": result["metrics"]["auroc"],
                    "AUPRC": result["metrics"]["auprc"],
                    "F1": result["metrics"]["f1"],
                    "run_status": "success",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "config_id": idx,
                    "hidden_dim": cfg["hidden_dim"],
                    "num_layers": cfg["num_layers"],
                    "dropout": cfg["dropout"],
                    "normalization": "on" if cfg["feature_normalization"] else "off",
                    "accuracy": "",
                    "AUROC": "",
                    "AUPRC": "",
                    "F1": "",
                    "run_status": "failed: {}".format(str(exc)),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, "graphsage_tuning_metrics.tsv"), sep="\t", index=False)
    success = df[df["run_status"] == "success"].copy()
    success["AUROC"] = pd.to_numeric(success["AUROC"])
    success["AUPRC"] = pd.to_numeric(success["AUPRC"])
    best_auroc = success.sort_values(["AUROC", "AUPRC"], ascending=[False, False], kind="stable").iloc[0]
    best_auprc = success.sort_values(["AUPRC", "AUROC"], ascending=[False, False], kind="stable").iloc[0]
    best_config = {
        "hidden_dim": int(best_auroc["hidden_dim"]),
        "num_layers": int(best_auroc["num_layers"]),
        "dropout": float(best_auroc["dropout"]),
        "feature_normalization": True if best_auroc["normalization"] == "on" else False,
    }
    with open(os.path.join(outdir, "best_graphsage_config.yaml"), "w", encoding="utf-8") as handle:
        yaml.safe_dump(best_config, handle, sort_keys=False)
    lines = [
        "# GraphSAGE Tuning Summary",
        "",
        "- best_by_AUROC: hidden_dim={}, num_layers={}, dropout={}, normalization={}, AUROC={:.4f}, AUPRC={:.4f}".format(
            int(best_auroc["hidden_dim"]),
            int(best_auroc["num_layers"]),
            float(best_auroc["dropout"]),
            best_auroc["normalization"],
            float(best_auroc["AUROC"]),
            float(best_auroc["AUPRC"]),
        ),
        "- best_by_AUPRC: hidden_dim={}, num_layers={}, dropout={}, normalization={}, AUROC={:.4f}, AUPRC={:.4f}".format(
            int(best_auprc["hidden_dim"]),
            int(best_auprc["num_layers"]),
            float(best_auprc["dropout"]),
            best_auprc["normalization"],
            float(best_auprc["AUROC"]),
            float(best_auprc["AUPRC"]),
        ),
        "- judgement: 当前最佳 config 已足够作为 support-side GraphSAGE 默认候选，但仍建议在进入 target ranking 前固定 seed 再做一次重复验证。",
    ]
    with open(os.path.join(outdir, "graphsage_tuning_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    with open("90_graphsage_tuning_results.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    print("【GraphSAGE 调参完成】")


if __name__ == "__main__":
    main()
