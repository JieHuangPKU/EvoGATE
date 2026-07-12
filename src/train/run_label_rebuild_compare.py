"""Legacy label rebuild comparison script.

Historical experiment only. This script is not part of the current mainline
Fusarium label pipeline. Outputs written under
`results/label_rebuild_experiments/` are archival / provenance artifacts.

Mainline label assets are now produced by the label materialization pipeline:
- scripts/run_fgraminearum_label_materialization.sh
- workflow/fgraminearum_label_materialization.smk
"""

import io
import os
import sys

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, matthews_corrcoef, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import train_test_split

from src.train.train_support_graph_baseline import require_epgat_env, build_dgl_graph
from src.models.support_graph_baseline import build_support_graph_model


RESULT_ROOT = "results/label_rebuild_experiments"


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def _bool_series(series):
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def _load_feature_matrix(base_feature_path, orthology_path, use_ortholog):
    base = _read_tsv(base_feature_path)
    emb_cols = [c for c in base.columns if c.startswith("emb_")]
    out = base[["canonical_gene_id"] + emb_cols].copy()
    if use_ortholog:
        ortho = _read_tsv(orthology_path)
        patch = ortho[["canonical_gene_id", "presence_species_count", "occupancy", "single_copy_fraction", "pangenome_class", "orthology_join_status"]].copy()
        patch["presence_species_count"] = pd.to_numeric(patch["presence_species_count"], errors="coerce").fillna(0.0)
        patch["occupancy"] = pd.to_numeric(patch["occupancy"], errors="coerce").fillna(0.0)
        patch["single_copy_fraction"] = pd.to_numeric(patch["single_copy_fraction"], errors="coerce").fillna(0.0)
        patch["pangenome_core"] = (patch["pangenome_class"] == "core").astype(float)
        patch["pangenome_softcore"] = (patch["pangenome_class"] == "softcore").astype(float)
        patch["pangenome_shell"] = (patch["pangenome_class"] == "shell").astype(float)
        patch["pangenome_cloud"] = (patch["pangenome_class"] == "cloud").astype(float)
        patch["orthology_patch_has_feature"] = (patch["orthology_join_status"] == "exact").astype(float)
        patch["orthology_patch_missing_mask"] = (patch["orthology_join_status"] == "missing_in_orthofinder_or_unresolved").astype(float)
        patch["orthology_patch_ambiguous_mask"] = (patch["orthology_join_status"] == "ambiguous_multiple_orthogroups").astype(float)
        patch = patch[
            [
                "canonical_gene_id",
                "presence_species_count",
                "occupancy",
                "single_copy_fraction",
                "pangenome_core",
                "pangenome_softcore",
                "pangenome_shell",
                "pangenome_cloud",
                "orthology_patch_has_feature",
                "orthology_patch_missing_mask",
                "orthology_patch_ambiguous_mask",
            ]
        ]
        out = out.merge(patch, on="canonical_gene_id", how="left")
        for c in patch.columns:
            if c != "canonical_gene_id":
                out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
    return out


def _prepare_graph(edge_path, node_ids):
    edges = _read_tsv(edge_path)
    node_index = dict(zip(node_ids, range(len(node_ids))))
    edges["source_node_index"] = edges["source_canonical_gene_id"].map(node_index)
    edges["target_node_index"] = edges["target_canonical_gene_id"].map(node_index)
    edges = edges.dropna(subset=["source_node_index", "target_node_index"]).copy()
    graph = build_dgl_graph(edges, len(node_ids))
    return graph, edges


def _safe_metrics(y_true, y_score):
    y_pred = (y_score >= 0.5).astype(int)
    return {
        "AUROC": float(roc_auc_score(y_true, y_score)),
        "AUPRC": float(average_precision_score(y_true, y_score)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "Precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "Recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "F1": float(f1_score(y_true, y_pred, zero_division=0)),
        "MCC": float(matthews_corrcoef(y_true, y_pred)),
    }


def _make_label_table(cfg):
    strict_df = _read_tsv(cfg["paths"]["strict29"])
    broad_df = _read_tsv(cfg["paths"]["broad79"])
    conflict_df = _read_tsv(cfg["paths"]["conflict8"])

    merged_pos = pd.concat(
        [
            strict_df.assign(source_set="strict29"),
            broad_df.assign(source_set="broad79"),
            conflict_df.assign(source_set="conflict8"),
        ],
        ignore_index=True,
    )
    unique_pos = merged_pos.sort_values(["canonical_gene_id", "source_set"], kind="stable").drop_duplicates("canonical_gene_id", keep="first")
    lethal_df = unique_pos[unique_pos["primary_evidence_term"].astype(str) == "lethal"].copy()

    outdir = os.path.join(RESULT_ROOT, "labels")
    _mkdir(outdir)
    lethal_df.to_csv(os.path.join(outdir, "lethal_positive_gene_list.tsv"), sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "strict29_gene_count": int(strict_df["canonical_gene_id"].nunique()),
                "broad79_gene_count": int(broad_df["canonical_gene_id"].nunique()),
                "conflict8_gene_count": int(conflict_df["canonical_gene_id"].nunique()),
                "merged_unique_gene_count": int(unique_pos["canonical_gene_id"].nunique()),
                "lethal_hard_positive_count": int(lethal_df["canonical_gene_id"].nunique()),
            }
        ]
    ).to_csv(os.path.join(outdir, "lethal_positive_summary.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, "lethal_positive_summary.md"), "w", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "# Lethal Positive Summary",
                    "- strict29 gene count = {}".format(int(strict_df["canonical_gene_id"].nunique())),
                    "- broad79 gene count = {}".format(int(broad_df["canonical_gene_id"].nunique())),
                    "- conflict8 gene count = {}".format(int(conflict_df["canonical_gene_id"].nunique())),
                    "- merged unique gene count = {}".format(int(unique_pos["canonical_gene_id"].nunique())),
                    "- lethal hard positive count = {}".format(int(lethal_df["canonical_gene_id"].nunique())),
                ]
            )
        )

    weak = _read_tsv(cfg["paths"]["yeast_labels"])
    bridge = _read_tsv(cfg["paths"]["yeast_bridge"])
    weak = weak.merge(bridge[["old_gene_id", "canonical_gene_id", "mapping_status", "mapping_rule", "needs_manual_review"]], left_on="ph1_gene_id", right_on="old_gene_id", how="left")
    weak_exact = weak[weak["mapping_status"] == "exact"].copy()

    lethal_set = set(lethal_df["canonical_gene_id"].astype(str))
    high_set = set(weak_exact[weak_exact["weak_positive_confidence"] == "high"]["canonical_gene_id"].astype(str))
    high_both_set = set(
        weak_exact[
            (weak_exact["weak_positive_confidence"] == "high")
            & (weak_exact["yeast_essential_support_class"] == "both")
        ]["canonical_gene_id"].astype(str)
    )

    master_evidence = _read_tsv(cfg["paths"]["master_evidence"])
    fg_evidence = master_evidence[master_evidence["species"] == "fgraminearum"].copy()
    virulence_genes = set(
        fg_evidence[
            (fg_evidence["evidence_class"].astype(str) == "virulence_only")
            | fg_evidence["evidence_term_raw"].astype(str).str.contains("virulence|pathogenicity", case=False, regex=True)
        ]["canonical_gene_id"].astype(str)
    )

    none_set = set(weak_exact[weak_exact["weak_positive_confidence"] == "none"]["canonical_gene_id"].astype(str))
    negative_set = sorted(none_set - virulence_genes - lethal_set - high_set)

    p1_set = sorted(lethal_set | high_set)
    p2_set = sorted(lethal_set | high_both_set)

    def build_positive_df(name, genes, weak_source_set):
        rows = []
        for g in genes:
            tags = []
            if g in lethal_set:
                tags.append("lethal")
            if g in weak_source_set:
                tags.append("weak_positive")
            rows.append({"canonical_gene_id": g, "positive_sources": ";".join(tags), "positive_scheme": name})
        return pd.DataFrame(rows)

    p1_df = build_positive_df("P1", p1_set, high_set)
    p2_df = build_positive_df("P2", p2_set, high_both_set)
    neg_df = pd.DataFrame({"canonical_gene_id": negative_set, "label": 0})

    p1_df.to_csv(os.path.join(outdir, "positive_set_P1.tsv"), sep="\t", index=False)
    p2_df.to_csv(os.path.join(outdir, "positive_set_P2.tsv"), sep="\t", index=False)
    neg_df.to_csv(os.path.join(outdir, "negative_set.tsv"), sep="\t", index=False)

    pos_comp = pd.DataFrame(
        [
            {
                "lethal_only_count": len(lethal_set - high_set),
                "weak_high_only_count": len(high_set - lethal_set),
                "weak_high_both_only_count": len(high_both_set - lethal_set),
                "union_P1_count": len(p1_set),
                "union_P2_count": len(p2_set),
                "overlap_lethal_high": len(lethal_set & high_set),
                "overlap_lethal_high_both": len(lethal_set & high_both_set),
            }
        ]
    )
    pos_comp.to_csv(os.path.join(outdir, "positive_set_comparison.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, "positive_set_summary.md"), "w", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "# Positive Set Summary",
                    "- lethal-only = {}".format(len(lethal_set - high_set)),
                    "- weak-high-only = {}".format(len(high_set - lethal_set)),
                    "- weak-high-both-only = {}".format(len(high_both_set - lethal_set)),
                    "- P1 union = {}".format(len(p1_set)),
                    "- P2 union = {}".format(len(p2_set)),
                    "- overlap(lethal, high) = {}".format(len(lethal_set & high_set)),
                    "- overlap(lethal, high&both) = {}".format(len(lethal_set & high_both_set)),
                ]
            )
        )
    pd.DataFrame(
        [
            {
                "none_total": len(none_set),
                "virulence_excluded": len(none_set & virulence_genes),
                "negative_final_count": len(negative_set),
            }
        ]
    ).to_csv(os.path.join(outdir, "negative_set_summary.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, "negative_set_summary.md"), "w", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "# Negative Set Summary",
                    "- none total = {}".format(len(none_set)),
                    "- virulence excluded = {}".format(len(none_set & virulence_genes)),
                    "- final negative set = {}".format(len(negative_set)),
                ]
            )
        )
    return {
        "lethal_df": lethal_df,
        "weak_exact": weak_exact,
        "p1_df": p1_df,
        "p2_df": p2_df,
        "neg_df": neg_df,
        "virulence_genes": virulence_genes,
    }


def _run_one_experiment(experiment_name, positive_df, neg_df, feature_df, edge_df, cfg, detail_dir, lethal_df):
    import torch
    labels = pd.concat(
        [
            positive_df[["canonical_gene_id"]].assign(label=1),
            neg_df[["canonical_gene_id", "label"]],
        ],
        ignore_index=True,
    ).drop_duplicates(subset=["canonical_gene_id"], keep="first")

    work = feature_df.merge(labels, on="canonical_gene_id", how="inner")
    node_ids = feature_df["canonical_gene_id"].astype(str).tolist()
    node_index = dict(zip(node_ids, range(len(node_ids))))

    y = work["label"].astype(int).to_numpy()
    train_idx_raw, test_idx_raw = train_test_split(
        np.arange(len(work)),
        test_size=float(cfg["split"]["test_fraction"]),
        random_state=int(cfg["random_seed"]),
        stratify=y,
    )
    train_only_idx = train_idx_raw
    val_rel, test_rel = train_test_split(
        np.arange(len(train_only_idx)),
        test_size=float(cfg["split"]["val_fraction"]) / (1.0 - float(cfg["split"]["test_fraction"])),
        random_state=int(cfg["random_seed"]),
        stratify=y[train_only_idx],
    )
    val_idx_raw = train_only_idx[test_rel]
    train_idx_raw = train_only_idx[val_rel]

    all_x = feature_df.drop(columns=["canonical_gene_id"]).apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
    graph = build_dgl_graph(
        edge_df.assign(
            source_node_index=edge_df["source_canonical_gene_id"].map(node_index),
            target_node_index=edge_df["target_canonical_gene_id"].map(node_index),
        ).dropna(subset=["source_node_index", "target_node_index"]),
        len(node_ids),
    )

    label_lookup = dict(zip(work["canonical_gene_id"], work["label"]))
    labels_all = np.full(len(node_ids), -1, dtype=np.int64)
    labeled_node_positions = []
    for idx, cid in enumerate(node_ids):
        if cid in label_lookup:
            labels_all[idx] = int(label_lookup[cid])
            labeled_node_positions.append(idx)

    cid_to_rawidx = dict(zip(work["canonical_gene_id"], range(len(work))))
    train_nodes = np.array([node_index[cid] for cid in work.iloc[train_idx_raw]["canonical_gene_id"]], dtype=np.int64)
    val_nodes = np.array([node_index[cid] for cid in work.iloc[val_idx_raw]["canonical_gene_id"]], dtype=np.int64)
    test_nodes = np.array([node_index[cid] for cid in work.iloc[test_idx_raw]["canonical_gene_id"]], dtype=np.int64)

    x_tensor = torch.tensor(all_x, dtype=torch.float32)
    y_tensor = torch.tensor(labels_all.astype(np.float32), dtype=torch.float32)

    model = build_support_graph_model(
        "GraphSAGE",
        in_feats=all_x.shape[1],
        hidden_feats=int(cfg["graphsage"]["hidden_dim"]),
        out_feats=1,
        num_layers=int(cfg["graphsage"]["num_layers"]),
        dropout=float(cfg["graphsage"]["dropout"]),
        num_heads=2,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=float(cfg["graphsage"]["lr"]))
    pos = float((labels_all[train_nodes] == 1).sum())
    neg = float((labels_all[train_nodes] == 0).sum())
    pos_weight = torch.tensor([neg / pos if pos > 0 else 1.0], dtype=torch.float32)
    loss_fn = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight, reduction="mean")

    for epoch in range(int(cfg["graphsage"]["epochs"])):
        model.train()
        optimizer.zero_grad()
        logits = model(graph, x_tensor).squeeze()
        loss = loss_fn(logits[train_nodes], y_tensor[train_nodes])
        loss.backward()
        optimizer.step()

    model.eval()
    with torch.no_grad():
        logits = model(graph, x_tensor).squeeze()
        scores = torch.sigmoid(logits).cpu().numpy()

    test_true = labels_all[test_nodes]
    test_score = scores[test_nodes]
    metrics = _safe_metrics(test_true, test_score)

    all_scores = pd.DataFrame({"canonical_gene_id": node_ids, "prediction_score": scores}).sort_values(
        "prediction_score", ascending=False, kind="stable"
    ).reset_index(drop=True)
    all_scores["rank"] = np.arange(1, len(all_scores) + 1)

    pos_set = set(positive_df["canonical_gene_id"].astype(str))
    lethal_set = set(lethal_df["canonical_gene_id"].astype(str))
    weak_set = pos_set - lethal_set
    neg_set = set(neg_df["canonical_gene_id"].astype(str))
    rank_rows = []
    for label_name, gene_set in [("positive_set", pos_set), ("lethal_set", lethal_set)]:
        ranked = all_scores["canonical_gene_id"].astype(str).tolist()
        for k in [50, 100, 200, 500]:
            topk = ranked[:k]
            hit_count = len([g for g in topk if g in gene_set])
            rank_rows.append(
                {
                    "experiment": experiment_name,
                    "target_group": label_name,
                    "k": k,
                    "hit_count": hit_count,
                    "precision_at_k": hit_count / float(k),
                    "recall_at_k": hit_count / float(len(gene_set)) if gene_set else 0.0,
                }
            )
    rank_df = pd.DataFrame(rank_rows)

    dist_rows = []
    for group_name, gene_set in [("lethal", lethal_set), ("weak_positive", weak_set), ("negative", neg_set)]:
        sub = all_scores[all_scores["canonical_gene_id"].isin(gene_set)]
        dist_rows.append(
            {
                "experiment": experiment_name,
                "group_name": group_name,
                "count": len(sub),
                "score_mean": float(sub["prediction_score"].mean()) if len(sub) else float("nan"),
                "score_median": float(sub["prediction_score"].median()) if len(sub) else float("nan"),
                "score_q75": float(sub["prediction_score"].quantile(0.75)) if len(sub) else float("nan"),
            }
        )
    dist_df = pd.DataFrame(dist_rows)

    _mkdir(detail_dir)
    work.assign(
        split=np.where(
            np.isin(np.arange(len(work)), test_idx_raw),
            "test",
            np.where(np.isin(np.arange(len(work)), val_idx_raw), "val", "train"),
        )
    ).to_csv(os.path.join(detail_dir, "labeled_split_table.tsv"), sep="\t", index=False)
    all_scores.to_csv(os.path.join(detail_dir, "all_node_scores.tsv"), sep="\t", index=False)
    rank_df.to_csv(os.path.join(detail_dir, "ranking_metrics.tsv"), sep="\t", index=False)
    dist_df.to_csv(os.path.join(detail_dir, "score_distribution_summary.tsv"), sep="\t", index=False)

    metrics.update(
        {
            "train_count": int(len(train_idx_raw)),
            "val_count": int(len(val_idx_raw)),
            "test_count": int(len(test_idx_raw)),
            "node_count": int(len(node_ids)),
            "edge_count": int(len(edge_df)),
            "positive_count": int(len(pos_set)),
            "negative_count": int(len(neg_df)),
        }
    )
    return metrics, rank_df, dist_df


def _build_experiment_features(cfg):
    base = _load_feature_matrix(
        cfg["paths"]["fusarium_feature_matrix"],
        cfg["paths"]["fusarium_orthology_features"],
        use_ortholog=False,
    )
    emb_cols = ["canonical_gene_id"] + [c for c in base.columns if c.startswith("emb_")]
    f1 = base[emb_cols].copy()

    ortho = _read_tsv(cfg["paths"]["fusarium_orthology_features"])
    ortho["presence_species_count"] = pd.to_numeric(ortho["presence_species_count"], errors="coerce").fillna(0.0)
    ortho["occupancy"] = pd.to_numeric(ortho["occupancy"], errors="coerce").fillna(0.0)
    ortho["single_copy_fraction"] = pd.to_numeric(ortho["single_copy_fraction"], errors="coerce").fillna(0.0)
    ortho["pangenome_core"] = (ortho["pangenome_class"] == "core").astype(float)
    ortho["pangenome_softcore"] = (ortho["pangenome_class"] == "softcore").astype(float)
    ortho["pangenome_shell"] = (ortho["pangenome_class"] == "shell").astype(float)
    ortho["pangenome_cloud"] = (ortho["pangenome_class"] == "cloud").astype(float)
    ortho["has_fusarium_ortholog_feature"] = (ortho["orthology_join_status"] == "exact").astype(float)
    ortho["fusarium_ortholog_missing_mask"] = (ortho["orthology_join_status"] == "missing_in_orthofinder_or_unresolved").astype(float)
    ortho["fusarium_ortholog_ambiguous_mask"] = (ortho["orthology_join_status"] == "ambiguous_multiple_orthogroups").astype(float)
    ortho_patch = ortho[
        [
            "canonical_gene_id",
            "presence_species_count",
            "occupancy",
            "single_copy_fraction",
            "pangenome_core",
            "pangenome_softcore",
            "pangenome_shell",
            "pangenome_cloud",
            "has_fusarium_ortholog_feature",
            "fusarium_ortholog_missing_mask",
            "fusarium_ortholog_ambiguous_mask",
        ]
    ].copy()
    f2 = f1.merge(ortho_patch, on="canonical_gene_id", how="left")
    for c in ortho_patch.columns:
        if c != "canonical_gene_id":
            f2[c] = pd.to_numeric(f2[c], errors="coerce").fillna(0.0)
    return f1, f2


if __name__ == "__main__":
    enable_utf8_stdout()
    require_epgat_env()
    print(
        "WARNING: This script is legacy and not part of the current mainline label pipeline.",
        file=sys.stderr,
        flush=True,
    )
    cfg = load_yaml("configs/label_rebuild_compare.yaml")
    _mkdir(RESULT_ROOT)
    _mkdir(os.path.join(RESULT_ROOT, "metrics"))
    _mkdir(os.path.join(RESULT_ROOT, "metrics", "per_experiment_detailed"))
    _mkdir(os.path.join(RESULT_ROOT, "plots"))
    _mkdir(os.path.join(RESULT_ROOT, "labels"))

    labels_bundle = _make_label_table(cfg)
    print("【开始 label rebuild + GraphSAGE 对照实验】")
    print("已重建 lethal / weak positive / negative sets")

    f1_df, f2_df = _build_experiment_features(cfg)
    edge_df = _read_tsv(cfg["paths"]["fusarium_edges"])

    experiments = [
        ("L1_F1", labels_bundle["p1_df"], labels_bundle["neg_df"], f1_df, "lethal + high", "PPI + Fusarium embedding"),
        ("L1_F2", labels_bundle["p1_df"], labels_bundle["neg_df"], f2_df, "lethal + high", "PPI + Fusarium embedding + Fusarium ortholog"),
        ("L2_F1", labels_bundle["p2_df"], labels_bundle["neg_df"], f1_df, "lethal + high&both", "PPI + Fusarium embedding"),
        ("L2_F2", labels_bundle["p2_df"], labels_bundle["neg_df"], f2_df, "lethal + high&both", "PPI + Fusarium embedding + Fusarium ortholog"),
    ]

    metric_rows = []
    ranking_rows = []
    dist_rows = []

    for exp_name, pos_df, neg_df, feat_df, label_scheme, feature_cfg in experiments:
        print("【开始实验】{}".format(exp_name))
        metrics, rank_df, dist_df = _run_one_experiment(
            exp_name,
            pos_df,
            neg_df,
            feat_df,
            edge_df,
            cfg,
            os.path.join(RESULT_ROOT, "metrics", "per_experiment_detailed", exp_name),
            labels_bundle["lethal_df"],
        )
        metric_rows.append(
            {
                "experiment": exp_name,
                "label_scheme": label_scheme,
                "feature_config": feature_cfg,
                **metrics,
            }
        )
        ranking_rows.append(rank_df)
        dist_rows.append(dist_df)

    all_metrics = pd.DataFrame(metric_rows)
    all_metrics.to_csv(os.path.join(RESULT_ROOT, "metrics", "all_experiment_metrics.tsv"), sep="\t", index=False)
    pd.concat(ranking_rows, ignore_index=True).to_csv(
        os.path.join(RESULT_ROOT, "metrics", "per_experiment_detailed", "all_ranking_metrics.tsv"),
        sep="\t",
        index=False,
    )
    pd.concat(dist_rows, ignore_index=True).to_csv(
        os.path.join(RESULT_ROOT, "metrics", "per_experiment_detailed", "all_score_distributions.tsv"),
        sep="\t",
        index=False,
    )

    best = all_metrics.copy()
    best["AUROC"] = pd.to_numeric(best["AUROC"])
    best["AUPRC"] = pd.to_numeric(best["AUPRC"])
    best = best.sort_values(["AUPRC", "AUROC", "F1"], ascending=[False, False, False], kind="stable")
    best_row = best.iloc[0]
    best_label = best_row["label_scheme"]
    best_feature = best_row["feature_config"]

    md_lines = ["# All Experiment Metrics", ""]
    for _, row in all_metrics.iterrows():
        md_lines.append(
            "- {} | {} | {} | AUROC={:.4f} | AUPRC={:.4f} | F1={:.4f} | MCC={:.4f}".format(
                row["experiment"],
                row["label_scheme"],
                row["feature_config"],
                float(row["AUROC"]),
                float(row["AUPRC"]),
                float(row["F1"]),
                float(row["MCC"]),
            )
        )
    open(os.path.join(RESULT_ROOT, "metrics", "all_experiment_metrics.md"), "w", encoding="utf-8").write("\n".join(md_lines))

    final_lines = [
        "# Final Recommendation",
        "",
        "- recommended_label_scheme: {}".format(best_label),
        "- recommended_feature_config: {}".format(best_feature),
        "- best_experiment: {}".format(best_row["experiment"]),
        "- judgement: {}".format(
            "差异足够明确，可作为默认方案" if len(best) > 1 and float(best.iloc[0]["AUPRC"]) - float(best.iloc[1]["AUPRC"]) > 0.01 else "差异不显著/谨慎作为默认复杂化"
        ),
    ]
    open(os.path.join(RESULT_ROOT, "final_recommendation.md"), "w", encoding="utf-8").write("\n".join(final_lines))

    print("【实验完成】")
    print(all_metrics[["experiment", "label_scheme", "feature_config", "AUROC", "AUPRC", "F1", "MCC"]].to_string(index=False))
    print("最终推荐 label scheme: {}".format(best_label))
    print("最终推荐 feature config: {}".format(best_feature))
