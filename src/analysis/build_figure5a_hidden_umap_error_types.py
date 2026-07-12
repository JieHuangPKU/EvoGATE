import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from src.analysis.figure5_representation_common import (
    TRANSITION_ORDER,
    UMAP_PARAMS,
    compute_umap,
    fine_scatter,
    load_hidden_case,
    pair_cases,
    save_png,
    save_pdf,
    species_title,
    write_json,
    write_markdown,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build Figure5A single-species hidden UMAP with error transitions")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--model", default="GraphSAGE", type=str)
    parser.add_argument("--seed", default=1029, type=int)
    parser.add_argument("--subset", default="test", type=str)
    parser.add_argument("--plot-dir", default="results/Figure5/plots", type=str)
    parser.add_argument("--data-dir", default="results/Figure5/data", type=str)
    parser.add_argument("--table-dir", default="results/Figure5/tables", type=str)
    parser.add_argument("--summary-dir", default="results/Figure5/summary", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    plot_dir = Path(args.plot_dir).resolve()
    data_dir = Path(args.data_dir).resolve()
    table_dir = Path(args.table_dir).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    for path in [plot_dir, data_dir, table_dir, summary_dir]:
        path.mkdir(parents=True, exist_ok=True)

    baseline_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB", args.model, args.seed, args.subset)
    esm2_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB_ESM2", args.model, args.seed, args.subset)
    paired = pair_cases(baseline_case, esm2_case)
    coords = compute_umap(esm2_case["hidden_matrix"], UMAP_PARAMS)
    species = esm2_case["species"]
    stem = f"Figure5a_hidden_umap_error_transition_{species}_seed{args.seed}"

    coords_df = paired.copy().rename(columns={"graph_gene_id": "node_id"})
    coords_df["species"] = species
    coords_df["protocol"] = args.protocol
    coords_df["feature_setting"] = "ORT_EXP_SUB_ESM2"
    coords_df["model"] = args.model
    coords_df["seed"] = args.seed
    coords_df["subset"] = args.subset
    coords_df["umap1"] = coords[:, 0]
    coords_df["umap2"] = coords[:, 1]
    coords_path = data_dir / f"{stem}_coords.tsv"
    transition_path = data_dir / f"{stem}_transition_labels.tsv"
    counts_path = table_dir / f"{stem}_transition_counts.tsv"
    pdf_path = plot_dir / f"{stem}.pdf"
    png_path = plot_dir / f"{stem}.png"
    summary_path = summary_dir / f"{stem}.md"
    json_path = table_dir / f"{stem}.json"

    coords_df.to_csv(coords_path, sep="\t", index=False)
    coords_df[
        [
            "node_id",
            "species",
            "protocol",
            "label",
            "split",
            "baseline_pred_label",
            "esm2_pred_label",
            "baseline_pred_score",
            "esm2_pred_score",
            "transition",
        ]
    ].to_csv(transition_path, sep="\t", index=False)

    counts_df = (
        coords_df.groupby(["species", "protocol", "transition"], dropna=False)
        .size()
        .reset_index(name="value")
    )
    counts_df["n"] = 1
    counts_df["mean"] = counts_df["value"].astype(float)
    counts_df["std"] = 0.0
    counts_df["variance"] = 0.0
    counts_df["transition_order"] = counts_df["transition"].map({name: idx for idx, name in enumerate(TRANSITION_ORDER)}).fillna(999).astype(int)
    counts_df = counts_df.sort_values(["species", "transition_order", "transition"]).drop(columns=["transition_order"])
    counts_df.to_csv(counts_path, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(5.3, 4.8), facecolor="white")
    fine_scatter(
        ax,
        coords,
        transitions=coords_df["transition"].astype(str).to_numpy(),
        title="{0}\nGraphSAGE hidden UMAP by error transition".format(species_title(args.protocol)),
        by="transition",
        legend=True,
    )
    save_png(fig, png_path)
    save_pdf(fig, pdf_path)

    rescued = int((coords_df["transition"] == "FN_to_TP_rescued").sum())
    corrected = int((coords_df["transition"] == "FP_to_TN_corrected").sum())
    write_markdown(
        summary_path,
        [
            "# Figure5a Hidden UMAP Error Transition",
            "",
            f"- Species: `{species}`.",
            f"- Protocol: `{args.protocol}`.",
            f"- Plot PDF: `{pdf_path}`.",
            f"- Plot PNG: `{png_path}`.",
            f"- Coordinates table: `{coords_path}`.",
            f"- Transition label table: `{transition_path}`.",
            f"- Transition counts table: `{counts_path}`.",
            f"- Rescued false negatives: `{rescued}`.",
            f"- Corrected false positives: `{corrected}`.",
        ],
    )

    write_json(
        json_path,
        {
            "panel": "Figure5a",
            "protocol": args.protocol,
            "species": species,
            "model": args.model,
            "seed": args.seed,
            "subset": args.subset,
            "umap_params": UMAP_PARAMS,
            "baseline_checkpoint_path": baseline_case["checkpoint_path"],
            "esm2_checkpoint_path": esm2_case["checkpoint_path"],
            "split_manifest_path": esm2_case["split_manifest_path"],
            "label_manifest_path": esm2_case["label_manifest_path"],
            "coords_tsv": str(coords_path),
            "transition_tsv": str(transition_path),
            "transition_counts_tsv": str(counts_path),
            "plot_pdf": str(pdf_path),
            "plot_png": str(png_path),
            "summary_md": str(summary_path),
        },
    )


if __name__ == "__main__":
    main()
