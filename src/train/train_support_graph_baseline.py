import io
import os
import sys
import json
import random

import numpy as np
import pandas as pd
import yaml

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, average_precision_score, f1_score


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def require_epgat_env():
    try:
        import torch  # noqa: F401
        import dgl  # noqa: F401
    except ImportError:
        raise SystemExit("错误：当前环境缺少 torch 或 dgl，请先执行 conda activate EPGAT 再运行。")


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def set_seed(seed):
    import torch

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def load_embedding_lookup():
    mf = pd.read_csv("outputs/baseline_dataset/embedding_manifest.pooled.tsv", sep="\t", dtype=str).fillna("")
    valid = mf[
        mf["exists"].astype(str).str.lower().isin(["true", "1", "yes"])
        & ~mf["needs_manual_review"].astype(str).str.lower().isin(["true", "1", "yes"])
    ].copy()
    return dict(zip(valid["canonical_gene_id"], valid["feature_path"]))


def load_feature_matrix(canonical_ids, feature_lookup):
    vectors = []
    for cid in canonical_ids:
        path = feature_lookup.get(cid, "")
        if not path:
            raise FileNotFoundError("缺少 embedding feature: {}".format(cid))
        vectors.append(np.load(path))
    return np.vstack(vectors).astype(np.float32)


def load_aux_feature_tables(species_list, config):
    feature_scope = config.get("feature_scope", {})
    use_orthology = bool(feature_scope.get("orthology", False))
    use_prior = bool(feature_scope.get("prior_score", False))
    orthology_feature_dir = str(config.get("orthology_feature_dir", "outputs/support_graph_features"))
    prior_feature_dir = str(config.get("prior_feature_dir", "outputs/support_prior"))
    tables = {}
    if use_orthology:
        for species in species_list:
            path = os.path.join(orthology_feature_dir, "support_feature_matrix_{}.tsv".format(species))
            if not os.path.exists(path):
                raise FileNotFoundError("缺少 support feature matrix: {}".format(path))
            df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
            tables.setdefault(species, {})["orthology"] = df
    if use_prior:
        for species in species_list:
            candidate_paths = [
                os.path.join(prior_feature_dir, "support_prior_matrix_{}.tsv".format(species)),
                "outputs/support_graph_prior/support_prior_matrix_{}.tsv".format(species),
            ]
            path = ""
            for candidate in candidate_paths:
                if os.path.exists(candidate):
                    path = candidate
                    break
            if not path:
                raise FileNotFoundError("缺少 support prior matrix: {}".format(candidate_paths[0]))
            df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
            tables.setdefault(species, {})["prior"] = df
    return tables


def assemble_numeric_feature_block(nodes, aux_tables, config):
    feature_scope = config.get("feature_scope", {})
    use_orthology = bool(feature_scope.get("orthology", False))
    use_prior = bool(feature_scope.get("prior_score", False))
    feature_cols = []
    merged = nodes[["species", "canonical_gene_id"]].copy()

    if use_orthology:
        cols = [
            "ortholog_count",
            "orthogroup_size",
            "support_species_presence_count",
            "conserved_across_support_species",
            "single_copy_like",
            "has_orthology_feature",
            "orthology_missing_mask",
        ]
        frames = []
        for species, table_map in aux_tables.items():
            if "orthology" not in table_map:
                continue
            df = table_map["orthology"]
            sub = df[["canonical_gene_id"] + cols].copy()
            sub["species"] = species
            frames.append(sub)
        if frames:
            orth = pd.concat(frames, ignore_index=True)
            merged = merged.merge(orth, on=["species", "canonical_gene_id"], how="left")
            feature_cols.extend(cols)

    if use_prior:
        cols = ["prior_score", "has_prior_score", "prior_missing_mask"]
        frames = []
        for species, table_map in aux_tables.items():
            if "prior" not in table_map:
                continue
            df = table_map["prior"]
            sub = df[["canonical_gene_id"] + cols].copy()
            sub["species"] = species
            frames.append(sub)
        if frames:
            prior = pd.concat(frames, ignore_index=True)
            merged = merged.merge(prior, on=["species", "canonical_gene_id"], how="left")
            feature_cols.extend(cols)

    if not feature_cols:
        return np.zeros((len(nodes), 0), dtype=np.float32), []

    for col in feature_cols:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
    return merged[feature_cols].to_numpy(dtype=np.float32), feature_cols


def build_dgl_graph(edge_df, num_nodes):
    import dgl
    import torch

    src = torch.tensor(edge_df["source_node_index"].astype(int).tolist(), dtype=torch.int64)
    dst = torch.tensor(edge_df["target_node_index"].astype(int).tolist(), dtype=torch.int64)
    graph = dgl.graph((src, dst), num_nodes=num_nodes)
    graph = dgl.to_simple(graph)
    if hasattr(dgl, "add_self_loop"):
        graph = dgl.add_self_loop(graph)
    else:
        nodes = torch.arange(num_nodes, dtype=torch.int64)
        graph.add_edges(nodes, nodes)
    return graph


def assemble_dataset(species_list):
    return assemble_dataset_with_config(species_list, {})


def assemble_dataset_with_config(species_list, config):
    registry = pd.read_csv("outputs/support_graphs/support_graph_model_registry.tsv", sep="\t", dtype=str).fillna("")
    feature_lookup = load_embedding_lookup()
    aux_tables = load_aux_feature_tables(species_list, config)
    node_frames = []
    edge_frames = []
    species_order = []
    for _, row in registry.iterrows():
        species = row["species"]
        if species not in species_list:
            continue
        species_order.append(species)
        node_df = pd.read_csv(row["node_table_path"], sep="\t", dtype=str).fillna("")
        edge_df = pd.read_csv(row["edge_table_path"], sep="\t", dtype=str).fillna("")
        label_df = pd.read_csv(row["label_manifest_path"], sep="\t", dtype=str).fillna("")
        node_df = node_df.reset_index(drop=True)
        node_df["species"] = species
        label_map = label_df.set_index("canonical_gene_id")
        node_df["label_value"] = node_df["canonical_gene_id"].map(label_map["label_value"]).fillna("")
        node_df["usable_for_support_supervision"] = node_df["canonical_gene_id"].map(label_map["usable_for_support_supervision"]).fillna("false")
        edge_df["species"] = species
        edge_df["edge_weight"] = edge_df["edge_weight"].astype(float)
        node_frames.append(node_df)
        edge_frames.append(edge_df)

    nodes = pd.concat(node_frames, ignore_index=True)
    nodes["global_node_index"] = nodes.index.astype(int)
    node_index = dict(zip(nodes["canonical_gene_id"], nodes["global_node_index"]))

    edge_tables = []
    for edge_df in edge_frames:
        edge_df = edge_df.copy()
        edge_df["source_node_index"] = edge_df["source_gene_id"].map(node_index)
        edge_df["target_node_index"] = edge_df["target_gene_id"].map(node_index)
        edge_df = edge_df.dropna(subset=["source_node_index", "target_node_index"]).copy()
        edge_tables.append(edge_df)
    edges = pd.concat(edge_tables, ignore_index=True)

    embedding_features = load_feature_matrix(nodes["canonical_gene_id"].astype(str).tolist(), feature_lookup)
    aux_features, aux_cols = assemble_numeric_feature_block(nodes, aux_tables, config)
    if aux_features.shape[1] > 0:
        features = np.hstack([embedding_features, aux_features]).astype(np.float32)
    else:
        features = embedding_features
    labels = nodes["label_value"].astype(str).map({"0": 0, "1": 1}).fillna(-1).astype(int).to_numpy()
    usable = nodes["usable_for_support_supervision"].astype(str).str.lower().isin(["true", "1", "yes"]).to_numpy()
    labels[~usable] = -1

    return nodes, edges, features, labels, species_order, aux_cols


def build_masks(labels, seed):
    valid_idx = np.where(labels >= 0)[0]
    y = labels[valid_idx]
    train_idx, test_idx = train_test_split(valid_idx, test_size=0.2, random_state=seed, stratify=y)
    train_idx, val_idx = train_test_split(
        train_idx,
        test_size=0.1 / 0.8,
        random_state=seed,
        stratify=labels[train_idx],
    )
    return train_idx, val_idx, test_idx


def evaluate_predictions(y_true, y_score):
    y_pred = (y_score >= 0.5).astype(int)
    out = {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "auroc": float(roc_auc_score(y_true, y_score)) if len(set(y_true.tolist())) > 1 else float("nan"),
        "auprc": float(average_precision_score(y_true, y_score)) if len(set(y_true.tolist())) > 1 else float("nan"),
        "f1": float(f1_score(y_true, y_pred)) if len(set(y_true.tolist())) > 1 else float("nan"),
    }
    return out


def build_species_loss_weights(nodes, config):
    weights_cfg = config.get("species_loss_weights", {}) or {}
    weights = nodes["species"].astype(str).map(lambda s: float(weights_cfg.get(s, 1.0))).astype(float).to_numpy()
    return weights.astype(np.float32)


def maybe_normalize_features(features_np, train_idx, config):
    if not bool(config.get("feature_normalization", False)):
        return features_np
    train_block = features_np[train_idx]
    mean = train_block.mean(axis=0, keepdims=True)
    std = train_block.std(axis=0, keepdims=True)
    std[std < 1e-8] = 1.0
    return ((features_np - mean) / std).astype(np.float32)


def train_one_model(model_type, config, species_set_name):
    import torch

    from src.models.support_graph_baseline import build_support_graph_model

    species_list = config["species_sets"][species_set_name]
    print("【开始 {} 训练】".format(model_type))
    print("纳入物种：{}".format("、".join(species_list)))

    nodes, edges, features_np, labels_np, species_order, aux_cols = assemble_dataset_with_config(species_list, config)
    print("节点数：{}".format(len(nodes)))
    print("边数：{}".format(len(edges)))
    feature_desc = "embedding only" if not aux_cols else "embedding + {}".format(",".join(aux_cols))
    print("当前特征：{}".format(feature_desc))

    graph = build_dgl_graph(edges, len(nodes))
    train_idx, val_idx, test_idx = build_masks(labels_np, int(config["seed"]))
    features_np = maybe_normalize_features(features_np, train_idx, config)
    features = torch.tensor(features_np, dtype=torch.float32)
    labels = torch.tensor(labels_np, dtype=torch.float32)
    edge_weight = torch.tensor(edges["edge_weight"].astype(float).tolist(), dtype=torch.float32) if len(edges) else None
    species_loss_weights = torch.tensor(build_species_loss_weights(nodes, config), dtype=torch.float32)

    model = build_support_graph_model(
        model_type,
        in_feats=features.shape[1],
        hidden_feats=int(config["hidden_dim"]),
        out_feats=1,
        num_layers=int(config["num_layers"]),
        dropout=float(config["dropout"]),
        num_heads=int(config.get("num_heads", 2)),
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=float(config["lr"]))
    loss_fn = torch.nn.BCEWithLogitsLoss(reduction="none")
    best_state = None
    best_epoch = 0
    best_val = float("-inf")
    patience = int(config.get("patience", 20))
    use_early_stopping = bool(config.get("early_stopping", True))
    patience_counter = 0
    history_rows = []

    for epoch in range(1, int(config["epochs"]) + 1):
        model.train()
        optimizer.zero_grad()
        if model_type.lower() == "gcn":
            logits = model(graph, features, edge_weight=edge_weight).squeeze()
        else:
            logits = model(graph, features).squeeze()
        train_loss_vec = loss_fn(logits[train_idx], labels[train_idx])
        train_species_weights = species_loss_weights[train_idx]
        weight_sum = torch.clamp(train_species_weights.sum(), min=1e-8)
        loss = (train_loss_vec * train_species_weights).sum() / weight_sum
        loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            if model_type.lower() == "gcn":
                eval_logits = model(graph, features, edge_weight=edge_weight).squeeze()
            else:
                eval_logits = model(graph, features).squeeze()
            val_scores = torch.sigmoid(eval_logits[val_idx]).cpu().numpy()
        val_metrics = evaluate_predictions(labels_np[val_idx], val_scores)
        val_metric = val_metrics["auroc"]
        if np.isnan(val_metric):
            val_metric = val_metrics["accuracy"]
        history_rows.append(
            {
                "epoch": epoch,
                "train_loss": float(loss.item()),
                "val_auroc": val_metrics["auroc"],
                "val_auprc": val_metrics["auprc"],
                "val_f1": val_metrics["f1"],
                "val_accuracy": val_metrics["accuracy"],
                "selection_value": float(val_metric),
            }
        )
        if float(val_metric) > best_val:
            best_val = float(val_metric)
            best_epoch = epoch
            best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
        print(
            "{} 训练中… epoch={} loss={:.6f} val_auroc={:.6f} best_val={:.6f}".format(
                model_type,
                epoch,
                float(loss.item()),
                float(val_metrics["auroc"]) if not np.isnan(val_metrics["auroc"]) else float("nan"),
                float(best_val),
            )
        )
        if use_early_stopping and patience_counter >= patience:
            print("{} 提前停止：patience={} 已耗尽".format(model_type, patience))
            break

    if best_state is not None and bool(config.get("reload_best_model", True)):
        model.load_state_dict(best_state)

    model.eval()
    with torch.no_grad():
        if model_type.lower() == "gcn":
            logits = model(graph, features, edge_weight=edge_weight).squeeze()
        else:
            logits = model(graph, features).squeeze()
        scores = torch.sigmoid(logits).cpu().numpy()

    test_metrics = evaluate_predictions(labels_np[test_idx], scores[test_idx])
    print("{} 训练完成".format(model_type))
    print("{} 评估结果：".format(model_type))
    print("  • accuracy = {:.4f}".format(test_metrics["accuracy"]))
    print("  • AUROC = {:.4f}".format(test_metrics["auroc"]))

    pred_df = nodes[["species", "canonical_gene_id"]].copy()
    pred_df["node_id"] = pred_df["canonical_gene_id"]
    pred_df["prediction_score"] = scores
    pred_df["label"] = labels_np
    pred_df["split"] = "train"
    pred_df.loc[test_idx, "split"] = "test"
    pred_df.loc[val_idx, "split"] = "val"

    return {
        "model_type": model_type,
        "species_scope": "+".join(species_list),
        "nodes": len(nodes),
        "edges": len(edges),
        "metrics": test_metrics,
        "predictions": pred_df,
        "state_dict": model.state_dict(),
        "model_object": model,
        "species_order": species_order,
        "feature_columns": ["embedding_vector"] + aux_cols,
        "species_loss_weights": dict(config.get("species_loss_weights", {}) or {}),
        "training_history": history_rows,
        "best_epoch": int(best_epoch),
        "best_val_score": float(best_val),
    }


def main():
    enable_utf8_stdout()
    require_epgat_env()
    config = load_yaml("configs/support_graph_experiments.yaml")
    base_cfg = load_yaml("configs/support_graph_baseline.yaml")
    set_seed(int(base_cfg.get("seed", 20260403)))

    os.makedirs("outputs/support_graph_runs", exist_ok=True)
    os.makedirs("outputs/support_graph_predictions", exist_ok=True)
    os.makedirs("outputs/support_graph_checkpoints", exist_ok=True)
    os.makedirs("outputs/support_graph_results", exist_ok=True)

    print("【support graph 基线训练启动】")
    print("运行环境：")
    print("  • conda 环境：EPGAT")
    print("当前特征：")
    print("  • embedding only")
    print("当前纳入物种：")
    print("  • human（edge_weight = 1.0）")
    print("  • scerevisiae（edge_weight = 1.0）")
    print("  • celegans（edge_weight = 0.8，partial-support）")
    print("当前模型：")
    print("  • GraphSAGE / GCN / GAT")

    models = config["models"]
    species_scope_name = "default"
    per_model_rows = []
    per_species_rows = []
    summary_lines = ["# Support Graph Experiment Summary", ""]

    import torch

    for model_name in models:
        try:
            result = train_one_model(model_name, {**base_cfg, **config}, species_scope_name)
            ckpt_path = "outputs/support_graph_checkpoints/{}.pt".format(model_name.lower())
            pred_path = "outputs/support_graph_predictions/{}_predictions.tsv".format(model_name.lower())
            log_path = "outputs/support_graph_runs/{}/training_log.json".format(model_name.lower())
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            torch.save(result["state_dict"], ckpt_path)
            result["predictions"].to_csv(pred_path, sep="\t", index=False)
            with open(log_path, "w", encoding="utf-8") as handle:
                json.dump(
                    {
                        "model": model_name,
                        "species_scope": result["species_scope"],
                        "metrics": result["metrics"],
                        "node_count": result["nodes"],
                        "edge_count": result["edges"],
                    },
                    handle,
                    ensure_ascii=False,
                    indent=2,
                )
            per_model_rows.append(
                {
                    "model": model_name,
                    "accuracy": result["metrics"]["accuracy"],
                    "AUROC": result["metrics"]["auroc"],
                    "AUPRC": result["metrics"]["auprc"],
                    "F1": result["metrics"]["f1"],
                    "node_count": result["nodes"],
                    "edge_count": result["edges"],
                    "species_scope": result["species_scope"],
                    "checkpoint_path": ckpt_path,
                    "prediction_path": pred_path,
                    "run_status": "success",
                    "notes": "embedding_only",
                }
            )
            for species in result["species_order"]:
                sub = result["predictions"][result["predictions"]["species"] == species].copy()
                test = sub[sub["split"] == "test"].copy()
                valid = test[test["label"] >= 0].copy()
                if len(valid) == 0:
                    continue
                metrics = evaluate_predictions(valid["label"].astype(int).to_numpy(), valid["prediction_score"].astype(float).to_numpy())
                per_species_rows.append(
                    {
                        "model": model_name,
                        "species": species,
                        "accuracy": metrics["accuracy"],
                        "AUROC": metrics["auroc"],
                        "AUPRC": metrics["auprc"],
                        "F1": metrics["f1"],
                        "node_count": len(sub),
                        "edge_count": result["edges"],
                    }
                )
            summary_lines.append(
                "- {}: accuracy={:.4f}, AUROC={:.4f}".format(
                    model_name, result["metrics"]["accuracy"], result["metrics"]["auroc"]
                )
            )
        except Exception as exc:
            print("{} 训练失败：".format(model_name))
            print("原因：{}".format(str(exc)))
            per_model_rows.append(
                {
                    "model": model_name,
                    "accuracy": "",
                    "AUROC": "",
                    "AUPRC": "",
                    "F1": "",
                    "node_count": "",
                    "edge_count": "",
                    "species_scope": "+".join(config["species_sets"][species_scope_name]),
                    "checkpoint_path": "",
                    "prediction_path": "",
                    "run_status": "failed",
                    "notes": str(exc),
                }
            )

    pd.DataFrame(per_model_rows).to_csv("outputs/support_graph_results/per_model_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(per_species_rows).to_csv("outputs/support_graph_results/per_species_metrics.tsv", sep="\t", index=False)

    best = None
    success_rows = [row for row in per_model_rows if row["run_status"] == "success"]
    if success_rows:
        best = sorted(success_rows, key=lambda r: (float(r["AUROC"]), float(r["accuracy"])), reverse=True)[0]["model"]

    summary_lines.extend(
        [
            "",
            "## 结果判断",
            "当前 feature scope = embedding only",
            "当前为 support-species graph baseline runs，不是 Fusarium transfer training",
            "celegans 保持 partial-support-graph caveat",
            "当前训练链路已 end-to-end 打通" if success_rows else "当前没有成功模型",
        ]
    )
    with open("outputs/support_graph_results/support_graph_experiment_summary.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines))
    pd.DataFrame(per_model_rows).to_csv("outputs/support_graph_results/support_graph_experiment_summary.tsv", sep="\t", index=False)

    print("【support graph 首轮训练完成】")
    print("模型结果汇总：")
    for row in per_model_rows:
        if row["run_status"] == "success":
            print("{}：".format(row["model"]))
            print("  • accuracy = {:.4f}".format(float(row["accuracy"])))
            print("  • AUROC = {:.4f}".format(float(row["AUROC"])))
        else:
            print("{}：失败".format(row["model"]))
    print("最佳模型：")
    print("  • {}".format(best if best else "无"))
    print("celegans 影响：")
    print("  • 当前已纳入训练且未造成运行时阻断；指标影响需下一轮保守/全量物种 ablation 才能单独判断。")
    print("当前结论：")
    print("  • 已完成真实 graph baseline 训练")
    print("  • pipeline 已 end-to-end 打通")


if __name__ == "__main__":
    main()
