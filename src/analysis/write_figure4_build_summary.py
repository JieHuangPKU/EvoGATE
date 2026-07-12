import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Write Figure4 build summary markdown")
    parser.add_argument("--output", required=True, type=str)
    parser.add_argument("--figure4-root", default="results/Figure4_representation", type=str)
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    root = Path(args.figure4_root).resolve()
    manifest_a = pd.read_csv(root / "tables" / "Figure4A_manifest.tsv", sep="\t")
    manifest_b = pd.read_csv(root / "tables" / "Figure4B_manifest.tsv", sep="\t")
    manifest_c = pd.read_csv(root / "tables" / "Figure4C_manifest.tsv", sep="\t")

    lines = [
        "# Figure4 Build Summary",
        "",
        "- species: `fgraminearum_newlabel`, `scerevisiae`",
        "- label regime: `newlabel` for Fusarium, `standard` for yeast",
        "- model: `GraphSAGE`",
        "- feature settings: `ORT_EXP_SUB`, `ORT_EXP_SUB_ESM2`",
        "- seed: `1029`",
        "- subset: `test`",
        "- runtime config: `{0}`".format(args.runtime_config),
        "- UMAP parameters: `n_components=2, n_neighbors=15, min_dist=0.1, metric=euclidean, random_state=1029`",
        "",
        "## Figure4A",
        "",
        "- object: `ORT_EXP_SUB_ESM2` GraphSAGE penultimate hidden embedding on the shared test-node set.",
        "- error transition definition: `TP_stable`, `TN_stable`, `FN_to_TP_rescued`, `FP_to_TN_corrected`, `FN_persistent`, `FP_persistent`, plus regression categories when present.",
        "",
        "## Figure4B",
        "",
        "- separation metrics definition: computed on high-dimensional hidden embeddings, not on UMAP coordinates.",
        "- metrics: `centroid distance`, `silhouette score`, `Davies-Bouldin index`.",
        "- rescue summary: `total TP`, `total TN`, `FN_to_TP_rescued`, `FP_to_TN_corrected`, `persistent FN`, `persistent FP`.",
        "",
        "## Figure4C",
        "",
        "- rationale: input-level and hidden-level manifolds are different objects and should not be conflated.",
        "",
        "## Why Figure4 does not use old-style ESM-only plots as the main figure",
        "",
        "- old-style ESM-only manifolds are input-level or shallow sequence-model representations, not GraphSAGE hidden space.",
        "- the current Figure4 focuses on representation-level diagnostics inside the final graph model, which is the object relevant to the paper's mainline claim.",
        "",
        "## Why dual-species layout is the default main output",
        "",
        "- dual-species output increases robustness for the current submission stage.",
        "- analysis and assembly are separated so later removal of yeast only requires re-assembly, not re-analysis.",
        "",
        "## Mainline input manifests",
        "",
    ]
    merged = manifest_a[["protocol", "species", "baseline_checkpoint_path", "esm2_checkpoint_path", "split_manifest_path", "label_manifest_path"]].copy()
    merged = merged.drop_duplicates(subset=["protocol"]).reset_index(drop=True)
    for _, row in merged.iterrows():
            lines.append(
                "- `{0}`: baseline checkpoint=`{1}`, +ESM2 checkpoint=`{2}`, split manifest=`{3}`, label manifest=`{4}`".format(
                    row["protocol"], row["baseline_checkpoint_path"], row["esm2_checkpoint_path"], row["split_manifest_path"], row["label_manifest_path"]
                )
            )
    lines.extend(
        [
            "",
            "## Supplementary rules reserved but not included in the main figure",
            "",
            "- yeast-only supplementary variants",
            "- t-SNE robustness",
            "- old Bingo-style ESM manifold",
            "- seed robustness",
            "- style comparison",
            "- UMAP sensitivity grid",
            "",
            "## Output manifests",
            "",
            "- panel A manifest: `{0}`".format(root / "tables" / "Figure4A_manifest.tsv"),
            "- panel B manifest: `{0}`".format(root / "tables" / "Figure4B_manifest.tsv"),
            "- panel C manifest: `{0}`".format(root / "tables" / "Figure4C_manifest.tsv"),
        ]
    )
    Path(args.output).write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
