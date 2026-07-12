from pathlib import Path

import numpy as np
import pandas as pd

from src.analysis.figure4_representation_common import (
    LABEL_COLORS,
    TRANSITION_COLORS,
    TRANSITION_ORDER,
    UMAP_PARAMS,
    compute_umap,
    fine_scatter,
    load_hidden_case,
    pair_cases,
    save_pdf,
    separation_metrics,
    species_title,
    write_json,
)


def save_plot_bundle(fig, pdf_path, png_path=None):
    pdf_path = Path(pdf_path)
    png_path = Path(png_path) if png_path else pdf_path.with_suffix(".png")
    save_pdf(fig, pdf_path)
    fig = None
    # save_pdf closes the figure after writing the PDF, so PNG must be rendered separately upstream.
    return str(pdf_path), str(png_path)


def save_png(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path), format="png", dpi=300, transparent=False, bbox_inches="tight")


def write_stats_table(raw_df, group_cols, value_col, output_path):
    stats_df = raw_df.groupby(group_cols, dropna=False)[value_col].agg(["count", "mean", "std", "var"]).reset_index()
    stats_df = stats_df.rename(columns={"count": "n", "var": "variance"})
    stats_df["std"] = stats_df["std"].fillna(0.0)
    stats_df["variance"] = stats_df["variance"].fillna(0.0)
    stats_df.to_csv(output_path, sep="\t", index=False)
    return stats_df


def write_markdown(path, lines):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def concat_tables(left_path, right_path, output_path):
    combined = pd.concat(
        [
            pd.read_csv(left_path, sep="\t"),
            pd.read_csv(right_path, sep="\t"),
        ],
        ignore_index=True,
    )
    combined.to_csv(output_path, sep="\t", index=False)
    return combined
