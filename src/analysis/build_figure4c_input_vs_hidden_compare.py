import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from src.analysis.figure4_representation_common import (
    UMAP_PARAMS,
    compute_umap,
    fine_scatter,
    load_hidden_case,
    save_pdf,
    species_title,
    write_json,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build Figure4C single-species input vs hidden UMAP comparison")
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

    if baseline_case["node_ids"] != esm2_case["node_ids"]:
        raise ValueError("Node order mismatch for Figure4C")

    labels = baseline_case["labels"]
    objects = [
        ("Input ORT_EXP_SUB", baseline_case["input_matrix"]),
        ("Input ORT_EXP_SUB_ESM2", esm2_case["input_matrix"]),
        ("Hidden ORT_EXP_SUB", baseline_case["hidden_matrix"]),
        ("Hidden ORT_EXP_SUB_ESM2", esm2_case["hidden_matrix"]),
    ]
    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.6), facecolor="white")
    for ax, (name, matrix) in zip(axes.ravel(), objects):
        coords = compute_umap(matrix, UMAP_PARAMS)
        fine_scatter(ax, coords, labels=labels, title=name, by="label", legend=True)
        for idx, node_id in enumerate(baseline_case["node_ids"]):
            rows.append(
                {
                    "species": baseline_case["species"],
                    "protocol": args.protocol,
                    "node_id": node_id,
                    "label": int(labels[idx]),
                    "feature_object": name,
                    "seed": args.seed,
                    "subset": args.subset,
                    "umap1": float(coords[idx, 0]),
                    "umap2": float(coords[idx, 1]),
                }
            )
    fig.suptitle("{0}\nInput-level vs hidden-level UMAP".format(species_title(args.protocol)), fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    pdf_path = output_dir / "Figure4C_{0}_input_vs_hidden_umap_seed{1}.pdf".format(baseline_case["species"], args.seed)
    tsv_path = output_dir / "Figure4C_{0}_input_vs_hidden_umap_seed{1}_coords.tsv".format(baseline_case["species"], args.seed)
    pd.DataFrame(rows).to_csv(tsv_path, sep="\t", index=False)
    save_pdf(fig, pdf_path)

    write_json(
        output_dir / "Figure4C_{0}_input_vs_hidden_umap_seed{1}.json".format(baseline_case["species"], args.seed),
        {
            "protocol": args.protocol,
            "species": baseline_case["species"],
            "model": args.model,
            "seed": args.seed,
            "subset": args.subset,
            "umap_params": UMAP_PARAMS,
            "coords_tsv": str(tsv_path),
            "pdf_path": str(pdf_path),
        },
    )


if __name__ == "__main__":
    main()
