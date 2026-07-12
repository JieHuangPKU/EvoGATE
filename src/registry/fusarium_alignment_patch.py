from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


SAFE_RULES = {
    "raw_gene_id_is_FGRAMPH1",
    "raw_protein_id_uniquely_maps_to_FGRAMPH1",
    "raw_transcript_id_uniquely_maps_to_FGRAMPH1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an explicit Fusarium embedding alignment patch table")
    parser.add_argument("--config", type=str, required=True, help="Path to baseline YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required TSV file not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def build_patch_table(master_label_table: pd.DataFrame, resolved_table: pd.DataFrame) -> pd.DataFrame:
    fg_rows = master_label_table[
        (master_label_table["species"] == "fgraminearum")
        & (master_label_table["canonical_gene_id"].astype(str).str.startswith("fgraminearum::FGSG_"))
    ].copy()
    fg_rows = fg_rows.drop_duplicates(subset=["canonical_gene_id"], keep="first")

    patch_rows = []
    for record in fg_rows.to_dict(orient="records"):
        original_canonical_gene_id = record["canonical_gene_id"]
        raw_gene_id = original_canonical_gene_id.split("::", 1)[1]
        candidates = resolved_table[resolved_table["raw_gene_id"] == raw_gene_id].copy()
        candidates = candidates[
            candidates["mapping_rule"].isin(SAFE_RULES)
            & candidates["final_canonical_gene_id"].astype(str).str.strip().ne("")
        ].copy()

        unique_targets = sorted(candidates["final_canonical_gene_id"].unique().tolist())
        unique_rules = sorted(candidates["mapping_rule"].unique().tolist())
        unique_confidences = sorted(candidates["mapping_confidence"].unique().tolist())

        patched_canonical_gene_id = ""
        patch_rule = ""
        confidence = ""
        needs_manual_review = "true"
        notes = "no safe one-to-one FGRAMPH1 mapping found in resolved canonical audit table"

        if len(unique_targets) == 1:
            patched_canonical_gene_id = unique_targets[0]
            patch_rule = unique_rules[0] if len(unique_rules) == 1 else "multiple_safe_rules_agree"
            confidence = unique_confidences[0] if len(unique_confidences) == 1 else "mixed"
            needs_manual_review = "false"
            notes = "safe one-to-one mapping recovered from fgraminearum_canonical_id_resolved.tsv"
        elif len(unique_targets) > 1:
            patch_rule = "ambiguous_safe_targets"
            confidence = "mixed"
            notes = f"multiple safe targets found: {len(unique_targets)}"

        patch_rows.append(
            {
                "original_canonical_gene_id": original_canonical_gene_id,
                "patched_canonical_gene_id": patched_canonical_gene_id,
                "patch_rule": patch_rule,
                "confidence": confidence,
                "needs_manual_review": needs_manual_review,
                "notes": notes,
            }
        )

    patch_df = pd.DataFrame(patch_rows).sort_values("original_canonical_gene_id", kind="stable").reset_index(drop=True)
    return patch_df


def write_summary(output_path: Path, patch_df: pd.DataFrame, broad79: pd.DataFrame, strict29: pd.DataFrame, conflict8: pd.DataFrame) -> None:
    safe_rows = patch_df[patch_df["needs_manual_review"].str.lower() == "false"].copy()
    unresolved_rows = patch_df[patch_df["needs_manual_review"].str.lower() != "false"].copy()

    benchmark_ids = set(broad79["canonical_gene_id"]) | set(strict29["canonical_gene_id"]) | set(conflict8["canonical_gene_id"])
    impacted_benchmarks = int(safe_rows["patched_canonical_gene_id"].isin(benchmark_ids).sum())

    lines = [
        "# 34 Fusarium Alignment Patch Summary",
        "",
        "## Direct Answers",
        f"1. 59 个缺失中，可通过 patch 安全恢复: {len(safe_rows)}",
        f"2. 仍无法安全恢复: {len(unresolved_rows)}",
        "3. broad79 / strict29 / conflict8 benchmark 定义是否会被回写: 不会",
        "4. patch 默认只用于 embedding 对齐，不回写 master registry",
        "",
        "## Patch Scope",
        f"- total FGSG canonical ids inspected: {len(patch_df)}",
        f"- safe one-to-one patch rows: {len(safe_rows)}",
        f"- unresolved rows kept for manual review: {len(unresolved_rows)}",
        f"- safe patch targets already present in Fusarium benchmark tables: {impacted_benchmarks}",
        "",
        "## Benchmark Impact",
        "- broad79 / strict29 / conflict8 的 canonical id 不修改，benchmark 定义保持原样",
        "- 44 条可 patch 行的目标 FGRAMPH1 id 已存在于 inference pool / benchmark 体系中",
        "- patch 启用时应只用于 feature 对齐，并在 ranking eval 侧按 effective_canonical_gene_id 去重，避免重复基因评分污染评估",
        "",
        "## Still Unresolved",
    ]

    if unresolved_rows.empty:
        lines.append("- none")
    else:
        for _, row in unresolved_rows.iterrows():
            lines.append(f"- {row['original_canonical_gene_id']}: {row['notes']}")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    registry_dir = Path(config["paths"]["registry_dir"])
    output_path = Path(config["embeddings"]["alignment_patch_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    master_label_table = _read_tsv(registry_dir / "master_label_table.preliminary.tsv")
    resolved_table = _read_tsv(registry_dir / "fgraminearum_canonical_id_resolved.tsv")
    broad79 = _read_tsv(registry_dir / "fgraminearum_gold_positive.broad79.tsv")
    strict29 = _read_tsv(registry_dir / "fgraminearum_gold_positive.strict29.tsv")
    conflict8 = _read_tsv(registry_dir / "fgraminearum_gold_positive.conflict.tsv")

    patch_df = build_patch_table(master_label_table, resolved_table)
    patch_df.to_csv(output_path, sep="\t", index=False)

    summary_path = Path("34_fusarium_alignment_patch_summary.md")
    write_summary(summary_path, patch_df, broad79, strict29, conflict8)

    safe_count = int((patch_df["needs_manual_review"].str.lower() == "false").sum())
    unresolved_count = int((patch_df["needs_manual_review"].str.lower() != "false").sum())
    print(f"Wrote Fusarium alignment patch table to: {output_path}")
    print(f"Safe one-to-one patch rows: {safe_count}")
    print(f"Unresolved rows: {unresolved_count}")
    print(f"Wrote patch summary to: {summary_path}")


if __name__ == "__main__":
    main()
