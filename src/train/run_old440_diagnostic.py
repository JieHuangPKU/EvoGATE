"""Legacy old440 diagnostic script.

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

from src.train.run_label_rebuild_compare import (
    enable_utf8_stdout,
    require_epgat_env,
    load_yaml,
    _mkdir,
    _read_tsv,
    _load_feature_matrix,
    _run_one_experiment,
)


RESULT_ROOT = "results/label_rebuild_experiments"


def _build_old440_labels():
    outdir = os.path.join(RESULT_ROOT, "old440", "labels")
    _mkdir(outdir)

    old = pd.read_csv(
        "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/EssentialGenes/gene_list.txt",
        sep="\t",
        dtype=str,
    ).fillna("")
    raw_rows = len(old)
    dedup = old.drop_duplicates(subset=["Ensembl"], keep="first").copy()
    source_gene_count = dedup["Ensembl"].nunique()

    current_canonical = _read_tsv("data_registry/master_label_table.preliminary.tsv")
    fg = current_canonical[current_canonical["species"] == "fgraminearum"][
        ["canonical_gene_id", "raw_gene_id"]
    ].copy()
    raw_gene_to_canonical = dict(zip(fg["raw_gene_id"], fg["canonical_gene_id"]))
    bridge = _read_tsv("outputs/fusarium_orthology/fusarium_orthology_id_bridge.tsv")
    xp_to_canonical = dict(
        zip(
            bridge[bridge["mapping_status"] == "exact"]["old_gene_id"],
            bridge[bridge["mapping_status"] == "exact"]["canonical_gene_id"],
        )
    )

    rows = []
    for _, row in dedup.iterrows():
        src = row["Ensembl"]
        canonical = ""
        status = "unresolved"
        rule = "unresolved"
        source = ""
        if "::" in src and src in set(fg["canonical_gene_id"].astype(str)):
            canonical = src
            status = "matched"
            rule = "exact_canonical_match"
            source = "gene_list.txt -> canonical_gene_id"
        elif src in raw_gene_to_canonical:
            canonical = raw_gene_to_canonical[src]
            status = "matched"
            rule = "raw_gene_exact_match"
            source = "gene_list.txt -> raw_gene_id -> canonical_gene_id"
        elif src in xp_to_canonical:
            canonical = xp_to_canonical[src]
            status = "matched"
            rule = "bridge_match"
            source = "gene_list.txt -> orthology bridge -> canonical_gene_id"
        rows.append(
            {
                "source_gene_id": src,
                "canonical_gene_id": canonical,
                "mapping_status": status,
                "mapping_rule": rule,
                "mapping_source": source,
                "target_label": row["Target"],
            }
        )
    mapping = pd.DataFrame(rows)
    mapping.to_csv(os.path.join(outdir, "old440_mapping_audit.tsv"), sep="\t", index=False)

    matched = mapping[mapping["mapping_status"] == "matched"].copy()
    pos = matched[matched["target_label"] == "1"].copy()
    pos = pos.drop_duplicates(subset=["canonical_gene_id"], keep="first").copy()
    pos.to_csv(os.path.join(outdir, "positive_old440.tsv"), sep="\t", index=False)

    neg = _read_tsv(os.path.join(RESULT_ROOT, "labels", "negative_set.tsv"))
    overlap = set(pos["canonical_gene_id"].astype(str)).intersection(set(neg["canonical_gene_id"].astype(str)))
    neg = neg[~neg["canonical_gene_id"].astype(str).isin(overlap)].copy()
    neg.to_csv(os.path.join(outdir, "negative_old440.tsv"), sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "gene_list_raw_rows": raw_rows,
                "source_gene_count_deduplicated": source_gene_count,
                "matched_count": int((mapping["mapping_status"] == "matched").sum()),
                "unresolved_count": int((mapping["mapping_status"] == "unresolved").sum()),
                "ambiguous_count": int((mapping["mapping_status"] == "ambiguous").sum()),
                "final_positive_old440_count": int(len(pos)),
                "final_negative_old440_count": int(len(neg)),
                "negative_overlap_removed_count": int(len(overlap)),
            }
        ]
    )
    summary.to_csv(os.path.join(outdir, "old440_label_summary.tsv"), sep="\t", index=False)
    lines = [
        "# old440 Label Summary",
        "- gene_list.txt raw rows = {}".format(raw_rows),
        "- deduplicated source genes = {}".format(source_gene_count),
        "- matched = {}".format(int((mapping["mapping_status"] == "matched").sum())),
        "- unresolved = {}".format(int((mapping["mapping_status"] == "unresolved").sum())),
        "- ambiguous = {}".format(int((mapping["mapping_status"] == "ambiguous").sum())),
        "- final positive_old440 = {}".format(int(len(pos))),
        "- final negative_old440 = {}".format(int(len(neg))),
        "- negative overlap removed = {}".format(int(len(overlap))),
    ]
    open(os.path.join(outdir, "old440_label_summary.md"), "w", encoding="utf-8").write("\n".join(lines))
    return mapping, pos, neg, summary


if __name__ == "__main__":
    enable_utf8_stdout()
    require_epgat_env()
    print(
        "WARNING: This script is legacy and not part of the current mainline label pipeline.",
        file=sys.stderr,
        flush=True,
    )

    print("【开始 OLD440 诊断实验】")
    cfg = load_yaml("configs/label_rebuild_compare.yaml")
    mapping, pos_old440, neg_old440, summary = _build_old440_labels()

    old_metrics_dir = os.path.join(RESULT_ROOT, "old440", "metrics")
    detail_dir = os.path.join(old_metrics_dir, "per_experiment_detailed")
    _mkdir(old_metrics_dir)
    _mkdir(detail_dir)

    edge_df = _read_tsv(cfg["paths"]["fusarium_edges"])
    f1_df = _load_feature_matrix(cfg["paths"]["fusarium_feature_matrix"], cfg["paths"]["fusarium_orthology_features"], use_ortholog=False)
    f2_df = _load_feature_matrix(cfg["paths"]["fusarium_feature_matrix"], cfg["paths"]["fusarium_orthology_features"], use_ortholog=True)

    lethal_df = _read_tsv(os.path.join(RESULT_ROOT, "labels", "lethal_positive_gene_list.tsv"))

    experiments = [
        ("OLD440_F1", pos_old440[["canonical_gene_id"]].assign(positive_sources="old440", positive_scheme="OLD440"), neg_old440, f1_df, "old440", "PPI + Fusarium embedding"),
        ("OLD440_F2", pos_old440[["canonical_gene_id"]].assign(positive_sources="old440", positive_scheme="OLD440"), neg_old440, f2_df, "old440", "PPI + Fusarium embedding + Fusarium ortholog"),
    ]

    rows = []
    for exp_name, pos_df, neg_df, feat_df, label_scheme, feature_cfg in experiments:
        print("【开始实验】{}".format(exp_name))
        metrics, rank_df, dist_df = _run_one_experiment(
            exp_name,
            pos_df,
            neg_df,
            feat_df,
            edge_df,
            cfg,
            os.path.join(detail_dir, exp_name),
            lethal_df,
        )
        rows.append(
            {
                "experiment": exp_name,
                "label_scheme": label_scheme,
                "feature_config": feature_cfg,
                **metrics,
            }
        )

    old_metrics = pd.DataFrame(rows)
    old_metrics.to_csv(os.path.join(old_metrics_dir, "old440_experiment_metrics.tsv"), sep="\t", index=False)
    md_lines = ["# OLD440 Experiment Metrics", ""]
    for _, row in old_metrics.iterrows():
        md_lines.append(
            "- {} | {} | AUROC={:.4f} | AUPRC={:.4f} | F1={:.4f} | MCC={:.4f}".format(
                row["experiment"], row["feature_config"], float(row["AUROC"]), float(row["AUPRC"]), float(row["F1"]), float(row["MCC"])
            )
        )
    open(os.path.join(old_metrics_dir, "old440_experiment_metrics.md"), "w", encoding="utf-8").write("\n".join(md_lines))

    current_metrics = pd.read_csv(os.path.join(RESULT_ROOT, "metrics", "all_experiment_metrics.tsv"), sep="\t")
    combined = pd.concat([current_metrics, old_metrics], ignore_index=True, sort=False)
    combined.to_csv(os.path.join(RESULT_ROOT, "metrics", "combined_with_old440_metrics.tsv"), sep="\t", index=False)
    combined_md = ["# Combined Metrics with OLD440", ""]
    for _, row in combined.iterrows():
        combined_md.append(
            "- {} | {} | {} | AUROC={:.4f} | AUPRC={:.4f} | F1={:.4f}".format(
                row["experiment"], row["label_scheme"], row["feature_config"], float(row["AUROC"]), float(row["AUPRC"]), float(row["F1"])
            )
        )
    open(os.path.join(RESULT_ROOT, "metrics", "combined_with_old440_metrics.md"), "w", encoding="utf-8").write("\n".join(combined_md))

    l1f1 = current_metrics[current_metrics["experiment"] == "L1_F1"].iloc[0]
    old_f1 = old_metrics[old_metrics["experiment"] == "OLD440_F1"].iloc[0]
    old_f2 = old_metrics[old_metrics["experiment"] == "OLD440_F2"].iloc[0]

    if float(old_f1["AUROC"]) > float(l1f1["AUROC"]) and float(old_f1["AUPRC"]) > float(l1f1["AUPRC"]):
        diagnosis = "更像是 label 任务变了"
        detail = "OLD440_F1 在当前同一套 pipeline 上明显优于 L1_F1，支持当前性能下降主要来自 label 任务变窄/变硬，而不是当前图/embedding 完全失效。"
    else:
        diagnosis = "更像当前 pipeline / input 不一致"
        detail = "OLD440 在当前同一套 pipeline 上没有明显优于 L1_F1，因此更支持当前图 / embedding / mapping / 实现至少有一项与历史设置不一致。"

    old440_lines = [
        "# old440 Diagnostic Conclusion",
        "",
        "- OLD440_F1 AUROC = {:.4f}, AUPRC = {:.4f}".format(float(old_f1["AUROC"]), float(old_f1["AUPRC"])),
        "- OLD440_F2 AUROC = {:.4f}, AUPRC = {:.4f}".format(float(old_f2["AUROC"]), float(old_f2["AUPRC"])),
        "- L1_F1 AUROC = {:.4f}, AUPRC = {:.4f}".format(float(l1f1["AUROC"]), float(l1f1["AUPRC"])),
        "- 诊断结论：{}".format(diagnosis),
        "- 解释：{}".format(detail),
        "- 在 OLD440 下 ortholog feature 是否仍然无增益：{}".format(
            "是" if float(old_f2["AUROC"]) <= float(old_f1["AUROC"]) and float(old_f2["AUPRC"]) <= float(old_f1["AUPRC"]) else "否"
        ),
    ]
    open(os.path.join(RESULT_ROOT, "old440", "old440_diagnostic_conclusion.md"), "w", encoding="utf-8").write("\n".join(old440_lines))

    # optional graph coverage audit
    graph_df = edge_df.copy()
    all_nodes = set(graph_df["source_canonical_gene_id"].astype(str)) | set(graph_df["target_canonical_gene_id"].astype(str))
    degree = {}
    for col in ["source_canonical_gene_id", "target_canonical_gene_id"]:
        for gene, cnt in graph_df[col].astype(str).value_counts().items():
            degree[gene] = degree.get(gene, 0) + int(cnt)
    groups = {
        "lethal": set(lethal_df["canonical_gene_id"].astype(str)),
        "high": set(_read_tsv(os.path.join(RESULT_ROOT, "labels", "positive_set_P1.tsv"))["canonical_gene_id"].astype(str)),
        "old440": set(pos_old440["canonical_gene_id"].astype(str)),
        "negative": set(neg_old440["canonical_gene_id"].astype(str)),
    }
    gc_rows = []
    for name, genes in groups.items():
        mapped = [g for g in genes if g in all_nodes]
        degs = [degree.get(g, 0) for g in mapped]
        gc_rows.append(
            {
                "group_name": name,
                "total_genes": len(genes),
                "mapped_to_graph": len(mapped),
                "mapping_rate": len(mapped) / float(len(genes)) if genes else 0.0,
                "mean_degree": float(np.mean(degs)) if degs else 0.0,
                "median_degree": float(np.median(degs)) if degs else 0.0,
                "isolated_count": int(sum(1 for d in degs if d == 0)),
                "largest_component_fraction": "",
            }
        )
    pd.DataFrame(gc_rows).to_csv(os.path.join(RESULT_ROOT, "metrics", "graph_coverage_audit.tsv"), sep="\t", index=False)
    open(os.path.join(RESULT_ROOT, "metrics", "graph_coverage_audit.md"), "w", encoding="utf-8").write(
        "# Graph Coverage Audit\n\n" + "\n".join(
            [
                "- {}: mapped_to_graph = {}, mapping_rate = {:.4f}, mean_degree = {:.2f}".format(
                    row["group_name"], row["mapped_to_graph"], row["mapping_rate"], row["mean_degree"]
                )
                for row in gc_rows
            ]
        )
    )

    print("【OLD440 诊断实验完成】")
    print(old_metrics[["experiment", "feature_config", "AUROC", "AUPRC", "F1", "MCC"]].to_string(index=False))
    print("最终诊断结论：{}".format(diagnosis))
