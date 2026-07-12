import io
import os
import sys
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


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
        raise SystemExit("错误：当前环境缺少 torch 或 dgl，请先执行 conda activate EPGAT")


def load_embedding_lookup():
    mf = pd.read_csv("outputs/baseline_dataset/embedding_manifest.pooled.tsv", sep="\t", dtype=str).fillna("")
    valid = mf[
        mf["exists"].astype(str).str.lower().isin(["true", "1", "yes"])
        & ~mf["needs_manual_review"].astype(str).str.lower().isin(["true", "1", "yes"])
    ].copy()
    return dict(zip(valid["canonical_gene_id"], valid["feature_path"])), dict(zip(valid["canonical_gene_id"], valid["embedding_source"]))


def read_simple_xlsx(path):
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    rel_ns = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
    with zipfile.ZipFile(path, "r") as zf:
        shared_strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall(ns + "si"):
                texts = []
                for t in si.iter(ns + "t"):
                    texts.append(t.text or "")
                shared_strings.append("".join(texts))
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        sheets = workbook.find(ns + "sheets")
        first_sheet = sheets.findall(ns + "sheet")[0]
        rid = first_sheet.attrib[rel_ns + "id"]
        rel_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        target = None
        for rel in rel_root:
            if rel.attrib.get("Id") == rid:
                target = rel.attrib.get("Target")
                break
        if target is None:
            raise ValueError("无法定位 xlsx 首个 sheet")
        if not target.startswith("worksheets/"):
            target = target.replace("../", "")
        sheet_path = "xl/" + target
        root = ET.fromstring(zf.read(sheet_path))
        rows = []
        for row in root.find(ns + "sheetData").findall(ns + "row"):
            values = []
            cells = row.findall(ns + "c")
            for cell in cells:
                cell_type = cell.attrib.get("t", "")
                value_node = cell.find(ns + "v")
                if value_node is None:
                    values.append("")
                else:
                    raw = value_node.text or ""
                    if cell_type == "s":
                        values.append(shared_strings[int(raw)])
                    else:
                        values.append(raw)
            rows.append(values)
    header = rows[0]
    body = rows[1:]
    width = len(header)
    fixed = []
    for row in body:
        if len(row) < width:
            row = row + [""] * (width - len(row))
        fixed.append(row[:width])
    return pd.DataFrame(fixed, columns=header)


def load_vectors(canonical_ids, lookup):
    return np.vstack([np.load(lookup[cid]) for cid in canonical_ids]).astype(np.float32)


def build_benchmark_table():
    pos = read_simple_xlsx("/home/jiehuang/software/fungi/Bingo/data/fgraminearum/orig_sample_list/FG_Essential_genes.xlsx")
    neg = read_simple_xlsx("/home/jiehuang/software/fungi/Bingo/data/fgraminearum/orig_sample_list/FG_NonEssential_genes.xlsx")
    pos = pos[["Ensembl"]].copy()
    pos["label"] = 1
    neg = neg[["Ensembl"]].copy()
    neg["label"] = 0
    neg = neg.sample(n=len(pos), random_state=20260404).reset_index(drop=True)
    bench = pd.concat([pos, neg], ignore_index=True)
    bench = bench.rename(columns={"Ensembl": "raw_gene_id"})
    bench["canonical_gene_id"] = "fgraminearum::" + bench["raw_gene_id"].astype(str)

    emb_lookup, emb_source = load_embedding_lookup()
    bench["has_embedding"] = bench["canonical_gene_id"].isin(emb_lookup).astype(bool)
    bench["embedding_source"] = bench["canonical_gene_id"].map(emb_source).fillna("")
    prior = pd.read_csv("outputs/support_prior/fusarium_prior_scores.tsv", sep="\t", dtype=str).fillna("")
    prior = prior[prior["model_name"] == "mlp"][["canonical_gene_id", "model_name", "prior_score"]].copy()
    prior = prior.rename(columns={"model_name": "prior_model"})
    bench = bench.merge(prior, on="canonical_gene_id", how="left")
    bench["has_true_prior"] = bench["prior_score"].astype(str).ne("")
    bench["prior_score"] = pd.to_numeric(bench["prior_score"], errors="coerce").fillna(0.0)
    bench["prior_model"] = bench["prior_model"].fillna("mlp")
    return bench


def build_mlp():
    cfg = yaml.safe_load(open("configs/baseline.yaml", "r", encoding="utf-8"))
    model_cfg = cfg["model"]["mlp"]
    hidden = tuple(int(v) for v in model_cfg["hidden_layer_sizes"])
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "classifier",
                MLPClassifier(
                    hidden_layer_sizes=hidden,
                    activation=model_cfg.get("activation", "relu"),
                    solver=model_cfg.get("solver", "adam"),
                    batch_size=int(model_cfg.get("batch_size", 64)),
                    max_iter=min(int(model_cfg.get("max_iter", 300)), 200),
                    learning_rate_init=float(model_cfg.get("learning_rate_init", 1e-3)),
                    alpha=float(model_cfg.get("alpha", 1e-4)),
                    early_stopping=bool(model_cfg.get("early_stopping", True)),
                    random_state=20260404,
                ),
            ),
        ]
    )


def run_oof(bench, use_prior):
    emb_lookup, _ = load_embedding_lookup()
    ids = bench["canonical_gene_id"].astype(str).tolist()
    x = load_vectors(ids, emb_lookup)
    if use_prior:
        prior = bench["prior_score"].to_numpy(dtype=np.float32).reshape(-1, 1)
        has_prior = bench["has_true_prior"].astype(float).to_numpy(dtype=np.float32).reshape(-1, 1)
        prior_missing = (1.0 - has_prior).astype(np.float32)
        x = np.hstack([x, prior, has_prior, prior_missing]).astype(np.float32)
    y = bench["label"].astype(int).to_numpy()
    folds = StratifiedKFold(n_splits=5, shuffle=True, random_state=20260404)
    oof = np.zeros(len(y), dtype=np.float64)
    fold_id = np.zeros(len(y), dtype=np.int32)
    for fold, (tr, te) in enumerate(folds.split(x, y), start=1):
        model = build_mlp()
        model.fit(x[tr], y[tr])
        score = model.predict_proba(x[te])[:, 1]
        oof[te] = score
        fold_id[te] = fold
    pred = (oof >= 0.5).astype(int)
    metrics = {
        "AUROC": float(roc_auc_score(y, oof)),
        "AUPRC": float(average_precision_score(y, oof)),
        "accuracy": float(accuracy_score(y, pred)),
        "F1": float(f1_score(y, pred)),
    }
    out = bench.copy()
    out["prediction_score"] = oof
    out["prediction_label"] = pred
    out["fold_id"] = fold_id
    return metrics, out


def main():
    enable_utf8_stdout()
    require_epgat_env()
    outdir = "outputs/bingo_replay"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("【本轮任务：仅执行 Bingo 440 基准重演，不做 GraphSAGE / 不做 Fusarium 全基因组推理】")
    print("【开始审计 440-gene benchmark】")
    bench = build_benchmark_table()
    bench.to_csv(os.path.join(outdir, "bingo_440_benchmark_table.tsv"), sep="\t", index=False)

    audit_rows = [
        {
            "item": "benchmark_name",
            "value": "Bingo fgraminearum ESM replay",
            "notes": "User requested 440-gene setting.",
        },
        {
            "item": "old_positive_file",
            "value": "/home/jiehuang/software/fungi/Bingo/data/fgraminearum/orig_sample_list/FG_Essential_genes.xlsx",
            "notes": "Actual row count is 439, not 440.",
        },
        {
            "item": "old_negative_file",
            "value": "/home/jiehuang/software/fungi/Bingo/data/fgraminearum/orig_sample_list/FG_NonEssential_genes.xlsx",
            "notes": "Contains 13706 rows.",
        },
        {
            "item": "old_bingo_esm_setting",
            "value": "3-fold CV on 14145-gene whole benchmark",
            "notes": "Recovered from run_esm_model.log and results_esm.csv; this is not an exact 440-gene benchmark.",
        },
        {
            "item": "replay_match_level",
            "value": "approximate",
            "notes": "Approximate replay uses 439 essential + 439 sampled nonessential for a balanced small benchmark closest to the requested 440-gene setting.",
        },
        {
            "item": "embedding_source",
            "value": "bingo_pooled_embedding",
            "notes": "Current replay uses pooled ESM2-derived vectors from ProGATE_v2 embedding manifest; exact raw .pt replay is blocked in EPGAT torch 1.4 by old file-format incompatibility.",
        },
        {
            "item": "true_prior_source",
            "value": "outputs/support_prior/fusarium_prior_scores.tsv",
            "notes": "Uses best prior model = mlp.",
        },
    ]
    pd.DataFrame(audit_rows).to_csv(os.path.join(outdir, "bingo_440_benchmark_audit.tsv"), sep="\t", index=False)
    audit_lines = [
        "# Bingo 440 Benchmark Audit",
        "",
        "- The old Bingo fgraminearum ESM result recovered from logs is a 14145-gene 3-fold CV benchmark, not an exactly documented 440-gene benchmark.",
        "- The closest recoverable small benchmark source is `FG_Essential_genes.xlsx`, which has 439 positives.",
        "- This replay is therefore approximate, not exact.",
        "- Approximate replay design: 439 essential + 439 sampled nonessential, same label definition, same ESM-family embedding source lineage, same classifier family across both comparisons.",
    ]
    open(os.path.join(outdir, "bingo_440_benchmark_audit.md"), "w", encoding="utf-8").write("\n".join(audit_lines))
    open("93_bingo_440_replay_audit.md", "w", encoding="utf-8").write("\n".join(audit_lines))

    print("【开始运行】ESM2 embedding only")
    m0, p0 = run_oof(bench, use_prior=False)
    p0.to_csv(os.path.join(outdir, "embedding_only_predictions.tsv"), sep="\t", index=False)
    pd.DataFrame([m0]).to_csv(os.path.join(outdir, "embedding_only_metrics.tsv"), sep="\t", index=False)

    print("【开始运行】ESM2 embedding + true prior")
    m1, p1 = run_oof(bench, use_prior=True)
    p1.to_csv(os.path.join(outdir, "embedding_plus_true_prior_predictions.tsv"), sep="\t", index=False)
    pd.DataFrame([m1]).to_csv(os.path.join(outdir, "embedding_plus_true_prior_metrics.tsv"), sep="\t", index=False)

    comp = pd.DataFrame(
        [
            {"setting": "embedding_only", **m0},
            {"setting": "embedding_plus_true_prior", **m1},
        ]
    )
    comp["delta_AUROC_vs_embedding_only"] = comp["AUROC"] - float(m0["AUROC"])
    comp["delta_AUPRC_vs_embedding_only"] = comp["AUPRC"] - float(m0["AUPRC"])
    comp["delta_F1_vs_embedding_only"] = comp["F1"] - float(m0["F1"])
    comp.to_csv(os.path.join(outdir, "bingo_replay_comparison.tsv"), sep="\t", index=False)

    delta_auroc = float(m1["AUROC"]) - float(m0["AUROC"])
    delta_auprc = float(m1["AUPRC"]) - float(m0["AUPRC"])
    delta_f1 = float(m1["F1"]) - float(m0["F1"])
    lines = [
        "# Bingo Replay Comparison",
        "",
        "- match_level: approximate",
        "- embedding_only: AUROC = {:.4f}, AUPRC = {:.4f}, accuracy = {:.4f}, F1 = {:.4f}".format(
            m0["AUROC"], m0["AUPRC"], m0["accuracy"], m0["F1"]
        ),
        "- embedding_plus_true_prior: AUROC = {:.4f}, AUPRC = {:.4f}, accuracy = {:.4f}, F1 = {:.4f}".format(
            m1["AUROC"], m1["AUPRC"], m1["accuracy"], m1["F1"]
        ),
        "- delta_AUROC = {:.4f}".format(delta_auroc),
        "- delta_AUPRC = {:.4f}".format(delta_auprc),
        "- delta_F1 = {:.4f}".format(delta_f1),
        "- embedding-only replay close to old Bingo result? no; old recovered result is on a much larger 14145-gene 3-fold CV setting, so only qualitative comparison is honest.",
        "- true prior adds real gain on this approximate 440-style replay: {}".format("yes" if (delta_auroc > 0 and delta_auprc > 0) else "no"),
    ]
    open(os.path.join(outdir, "bingo_replay_comparison.md"), "w", encoding="utf-8").write("\n".join(lines))
    open("94_bingo_replay_results.md", "w", encoding="utf-8").write("\n".join(lines))
    next_lines = [
        "# 95 Next Step After Bingo Replay",
        "",
        "- Use the validated true-prior signal to move back into the target-side graph ranking stage.",
        "- Keep Bingo replay marked as approximate unless an exact historical 440-gene split file is found.",
    ]
    open("95_next_step_after_bingo_replay.md", "w", encoding="utf-8").write("\n".join(next_lines))

    print("【Bingo 440 基准重演完成】")
    print("结果汇总：")
    print("  • ESM2 only: AUROC = {:.4f}, AUPRC = {:.4f}".format(m0["AUROC"], m0["AUPRC"]))
    print("  • ESM2 + true prior: AUROC = {:.4f}, AUPRC = {:.4f}".format(m1["AUROC"], m1["AUPRC"]))
    print("  • delta_AUROC = {:.4f}".format(delta_auroc))
    print("  • delta_AUPRC = {:.4f}".format(delta_auprc))
    print("当前结论：")
    print("  • true prior 是否在原始 440-gene 任务上带来增益：{}".format("是" if (delta_auroc > 0 and delta_auprc > 0) else "否"))


if __name__ == "__main__":
    main()
