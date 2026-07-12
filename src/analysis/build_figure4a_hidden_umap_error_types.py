import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from src.analysis.figure4_representation_common import (
    UMAP_PARAMS,
    compute_umap,
    fine_scatter,
    load_hidden_case,
    pair_cases,
    save_pdf,
    species_title,
    write_json,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build Figure4A single-species hidden UMAP with error transitions")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--model", default="GraphSAGE", type=str)
    parser.add_argument("--seed", default=1029, type=int)
    parser.add_argument("--subset", default="test", type=str)
    parser.add_argument("--output-dir", required=True, type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    baseline_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB", args.model, args.seed, args.subset)
    esm2_case = load_hidden_case(args.runtime_config, args.upstream_root, args.protocol, "ORT_EXP_SUB_ESM2", args.model, args.seed, args.subset)
    paired = pair_cases(baseline_case, esm2_case)
    coords = compute_umap(esm2_case["hidden_matrix"], UMAP_PARAMS)

    coords_df = paired.copy()
    coords_df = coords_df.rename(columns={"graph_gene_id": "node_id"})
    coords_df["species"] = esm2_case["species"]
    coords_df["feature_setting"] = "ORT_EXP_SUB_ESM2"
    coords_df["model"] = args.model
    coords_df["seed"] = args.seed
    coords_df["umap1"] = coords[:, 0]
    coords_df["umap2"] = coords[:, 1]
    coords_tsv = output_dir / "Figure4A_{0}_hidden_umap_error_transition_seed{1}_coords.tsv".format(esm2_case["species"], args.seed)
    transition_tsv = output_dir / "Figure4A_{0}_hidden_umap_error_transition_seed{1}_transition_labels.tsv".format(esm2_case["species"], args.seed)
    pdf_path = output_dir / "Figure4A_{0}_hidden_umap_error_transition_seed{1}.pdf".format(esm2_case["species"], args.seed)
    coords_df.to_csv(coords_tsv, sep="\t", index=False)
    coords_df[
        [
            "node_id",
            "label",
            "split",
            "baseline_pred_label",
            "esm2_pred_label",
            "baseline_pred_score",
            "esm2_pred_score",
            "transition",
        ]
    ].to_csv(transition_tsv, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(5.3, 4.8), facecolor="white")
    fine_scatter(
        ax,
        coords,
        transitions=coords_df["transition"].astype(str).to_numpy(),
        title="{0}\nGraphSAGE hidden UMAP by error transition".format(species_title(args.protocol)),
        by="transition",
        legend=True,
    )
    save_pdf(fig, pdf_path)

    write_json(
        output_dir / "Figure4A_{0}_hidden_umap_error_transition_seed{1}.json".format(esm2_case["species"], args.seed),
        {
            "protocol": args.protocol,
            "species": esm2_case["species"],
            "model": args.model,
            "seed": args.seed,
            "subset": args.subset,
            "umap_params": UMAP_PARAMS,
            "baseline_checkpoint_path": baseline_case["checkpoint_path"],
            "esm2_checkpoint_path": esm2_case["checkpoint_path"],
            "split_manifest_path": esm2_case["split_manifest_path"],
            "label_manifest_path": esm2_case["label_manifest_path"],
            "coords_tsv": str(coords_tsv),
            "transition_tsv": str(transition_tsv),
            "pdf_path": str(pdf_path),
        },
    )


if __name__ == "__main__":
    main()
