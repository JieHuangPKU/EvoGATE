import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from src.analysis.figure5_representation_common import (
    UMAP_PARAMS,
    compute_umap,
    fine_scatter,
    load_hidden_case,
    save_png,
    save_pdf,
    separation_metrics,
    species_title,
    write_json,
    write_markdown,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build Figure5C single-species input vs hidden UMAP comparison")
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
    if baseline_case["node_ids"] != esm2_case["node_ids"]:
        raise ValueError("Node order mismatch for Figure5C")

    species = baseline_case["species"]
    stem = f"Figure5c_input_vs_hidden_compare_{species}_seed{args.seed}"
    labels = baseline_case["labels"]
    objects = [
        ("Input ORT_EXP_SUB", baseline_case["input_matrix"]),
        ("Input ORT_EXP_SUB_ESM2", esm2_case["input_matrix"]),
        ("Hidden ORT_EXP_SUB", baseline_case["hidden_matrix"]),
        ("Hidden ORT_EXP_SUB_ESM2", esm2_case["hidden_matrix"]),
    ]
    rows = []
    metric_rows = []
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.6), facecolor="white")
    for ax, (name, matrix) in zip(axes.ravel(), objects):
        coords = compute_umap(matrix, UMAP_PARAMS)
        fine_scatter(ax, coords, labels=labels, title=name, by="label", legend=True)
        metrics = separation_metrics(matrix, labels)
        for metric_name, metric_value in metrics.items():
            metric_rows.append(
                {
                    "species": species,
                    "protocol": args.protocol,
                    "feature_object": name,
                    "metric": metric_name,
                    "value": metric_value,
                    "n": 1,
                    "mean": metric_value,
                    "std": 0.0,
                    "variance": 0.0,
                }
            )
        for idx, node_id in enumerate(baseline_case["node_ids"]):
            rows.append(
                {
                    "species": species,
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

    pdf_path = plot_dir / f"{stem}.pdf"
    png_path = plot_dir / f"{stem}.png"
    coords_path = data_dir / f"{stem}_coords.tsv"
    metrics_path = table_dir / f"{stem}_metrics.tsv"
    summary_path = summary_dir / f"{stem}.md"
    json_path = table_dir / f"{stem}.json"
    pd.DataFrame(rows).to_csv(coords_path, sep="\t", index=False)
    pd.DataFrame(metric_rows).to_csv(metrics_path, sep="\t", index=False)
    save_png(fig, png_path)
    save_pdf(fig, pdf_path)

    write_markdown(
        summary_path,
        [
            "# Figure5c Input-Vs-Hidden Comparison",
            "",
            f"- Species: `{species}`.",
            f"- Plot PDF: `{pdf_path}`.",
            f"- Plot PNG: `{png_path}`.",
            f"- Coordinate table: `{coords_path}`.",
            f"- Metric table: `{metrics_path}`.",
            "- The metric table is computed on the original input/hidden objects; the coordinate table stores the plotted UMAP positions.",
        ],
    )

    write_json(
        json_path,
        {
            "panel": "Figure5c",
            "protocol": args.protocol,
            "species": species,
            "model": args.model,
            "seed": args.seed,
            "subset": args.subset,
            "umap_params": UMAP_PARAMS,
            "coords_tsv": str(coords_path),
            "metrics_tsv": str(metrics_path),
            "plot_pdf": str(pdf_path),
            "plot_png": str(png_path),
            "summary_md": str(summary_path),
        },
    )


if __name__ == "__main__":
    main()
