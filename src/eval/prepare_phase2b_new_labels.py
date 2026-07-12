"""Prepare the canonical new-label gene sets for the Phase 2B companion benchmark."""

from pathlib import Path
import pandas as pd

SOURCE_POSITIVE = Path("data/processed/essential_gene/fgraminearum/newlabel/positive_genes.tsv")
SOURCE_NEGATIVE = Path("data/processed/essential_gene/fgraminearum/newlabel/negative_genes.tsv")
OUTPUT_DIR = Path("results/phase2b_new_label/labels")


def main():
    if not SOURCE_POSITIVE.exists():
        raise FileNotFoundError(f"Missing source positive set: {SOURCE_POSITIVE}")
    if not SOURCE_NEGATIVE.exists():
        raise FileNotFoundError(f"Missing source negative set: {SOURCE_NEGATIVE}")

    pos = pd.read_csv(SOURCE_POSITIVE, sep="\t", dtype=str).fillna("")
    neg = pd.read_csv(SOURCE_NEGATIVE, sep="\t", dtype=str).fillna("")
    if "canonical_gene_id" not in pos.columns:
        raise ValueError("newlabel positive genes file must contain canonical_gene_id")
    if "canonical_gene_id" not in neg.columns:
        raise ValueError("newlabel negative genes file must contain canonical_gene_id")

    pos = pos.drop_duplicates(subset=["canonical_gene_id"], keep="first").reset_index(drop=True)
    neg = neg.drop_duplicates(subset=["canonical_gene_id"], keep="first").reset_index(drop=True)
    overlap = sorted(set(pos["canonical_gene_id"]).intersection(set(neg["canonical_gene_id"])))
    if overlap:
        raise ValueError(f"Positive/negative overlap detected for {len(overlap)} genes; first examples: {overlap[:10]}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pos_out = OUTPUT_DIR / "new_positive.tsv"
    neg_out = OUTPUT_DIR / "new_negative.tsv"
    summary_tsv = OUTPUT_DIR / "new_label_summary.tsv"
    summary_md = OUTPUT_DIR / "new_label_summary.md"

    pos.to_csv(pos_out, sep="\t", index=False)
    neg.to_csv(neg_out, sep="\t", index=False)

    lethal_only_count = int(pos["positive_sources"].astype(str).eq("lethal").sum()) if "positive_sources" in pos.columns else ""
    yeast_transfer_count = int(pos["positive_sources"].astype(str).str.contains("weak_positive", regex=False).sum()) if "positive_sources" in pos.columns else ""
    both_count = int(pos["positive_sources"].astype(str).str.contains("lethal;weak_positive|weak_positive;lethal", regex=True).sum()) if "positive_sources" in pos.columns else ""
    summary = pd.DataFrame([
        {
            "label_regime": "new_label",
            "positive_source_file": str(SOURCE_POSITIVE),
            "negative_source_file": str(SOURCE_NEGATIVE),
            "positive_count": int(len(pos)),
            "negative_count": int(len(neg)),
            "lethal_only_count": lethal_only_count,
            "yeast_transfer_positive_count": yeast_transfer_count,
            "lethal_and_yeast_transfer_overlap": both_count,
            "description": "lethal-only plus genes with at least one yeast single-copy ortholog support, paired with the final 10270 negative set",
        }
    ])
    summary.to_csv(summary_tsv, sep="\t", index=False)

    lines = [
        "# New Label Summary",
        "",
        f"- positive source: `{SOURCE_POSITIVE}`",
        f"- negative source: `{SOURCE_NEGATIVE}`",
        f"- positive count = {len(pos)}",
        f"- negative count = {len(neg)}",
        f"- lethal-only = {lethal_only_count}",
        f"- entries carrying yeast-transfer support = {yeast_transfer_count}",
        f"- overlap(lethal, yeast-transfer) = {both_count}",
        "- regime description: lethal-only plus genes with at least one yeast single-copy ortholog support; negatives are the final filtered none-set.",
    ]
    summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
