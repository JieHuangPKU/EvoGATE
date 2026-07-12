import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.neighbors import NearestNeighbors

from src.analysis.figure4_representation_common import (
    LABEL_COLORS,
    TRANSITION_ORDER,
    UMAP_PARAMS,
    compute_umap,
    fine_scatter,
    load_hidden_case,
    pair_cases,
    save_pdf,
    separation_metrics,
    species_title,
)


FUSARIUM_PROTOCOL = "fgraminearum_newlabel"
SCEREVISIAE_PROTOCOL = "scerevisiae"
FIGURE5_PROTOCOLS = [FUSARIUM_PROTOCOL, SCEREVISIAE_PROTOCOL]
MODEL = "GraphSAGE"
BASE_FEATURE = "ORT_EXP_SUB"
ESM2_FEATURE = "ORT_EXP_SUB_ESM2"
SUBSET = "test"
STABLE_SEED_THRESHOLD = 2

OUTPUT_SPECIES_SLUGS = {
    FUSARIUM_PROTOCOL: "fgraminearum",
    SCEREVISIAE_PROTOCOL: "scerevisiae",
}

FOCUS_TRANSITIONS = [
    {
        "raw": "FN_to_TP_rescued",
        "display": "FN_to_TP_rescued",
        "slug": "fn_to_tp_rescued",
    },
    {
        "raw": "FP_to_TN_corrected",
        "display": "FP_to_TN_corrected",
        "slug": "fp_to_tn_corrected",
    },
    {
        "raw": "FN_persistent",
        "display": "persistent_FN",
        "slug": "persistent_fn",
    },
    {
        "raw": "FP_persistent",
        "display": "persistent_FP",
        "slug": "persistent_fp",
    },
]

CONFIDENCE_TIER_ORDER = [
    "high_confidence_rescued",
    "stable_rescued",
    "seed_specific_rescued",
    "not_rescued",
]

CONFIDENCE_TIER_LABELS = {
    "high_confidence_rescued": "High confidence (4-5 seeds)",
    "stable_rescued": "Stable (2-3 seeds)",
    "seed_specific_rescued": "Seed-specific (1 seed)",
    "not_rescued": "Not rescued (0 seeds)",
}


def load_runtime_seed_list(runtime_config_path):
    config = yaml.safe_load(Path(runtime_config_path).read_text(encoding="utf-8"))
    return [int(seed) for seed in config["runtime"]["seed_list"]]


def load_metrics_row(upstream_root, protocol, feature_setting, seed):
    metrics_path = Path(upstream_root) / protocol / MODEL / feature_setting / f"run_{seed}" / "metrics.tsv"
    df = pd.read_csv(metrics_path, sep="\t")
    row = df.iloc[0].to_dict()
    row["metrics_path"] = str(metrics_path.resolve())
    return row


def load_paired_cases(runtime_config_path, upstream_root, protocol, seed, subset=SUBSET):
    baseline_case = load_hidden_case(runtime_config_path, upstream_root, protocol, BASE_FEATURE, MODEL, seed, subset)
    esm2_case = load_hidden_case(runtime_config_path, upstream_root, protocol, ESM2_FEATURE, MODEL, seed, subset)
    paired = pair_cases(baseline_case, esm2_case)
    return baseline_case, esm2_case, paired


def write_markdown(path, lines):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def save_png(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path), format="png", dpi=300, transparent=False, bbox_inches="tight")


def save_plot_pair(fig, pdf_path, png_path):
    save_png(fig, png_path)
    save_pdf(fig, pdf_path)


def mad_scale(values):
    values = np.asarray(values, dtype=float)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    scale = mad if mad > 1e-12 else max(float(np.std(values, ddof=0)), 1.0)
    return median, scale


def protocol_output_slug(protocol, species=None):
    if protocol in OUTPUT_SPECIES_SLUGS:
        return OUTPUT_SPECIES_SLUGS[protocol]
    if species:
        return str(species)
    return str(protocol)


def focus_transition_columns(stable_seed_threshold=STABLE_SEED_THRESHOLD):
    columns = []
    for item in FOCUS_TRANSITIONS:
        columns.append(
            {
                **item,
                "count_col": f"{item['slug']}_count",
                "frequency_col": f"{item['slug']}_frequency",
                "stable_col": f"{item['slug']}_stable_ge{stable_seed_threshold}",
            }
        )
    return columns


def confidence_tier_from_count(rescue_count, total_seeds, stable_seed_threshold=STABLE_SEED_THRESHOLD):
    if rescue_count >= max(4, total_seeds - 1):
        return "high_confidence_rescued"
    if rescue_count >= stable_seed_threshold:
        return "stable_rescued"
    if rescue_count > 0:
        return "seed_specific_rescued"
    return "not_rescued"


def compute_neighbor_fraction(matrix, labels, k):
    matrix = np.asarray(matrix, dtype=np.float32)
    labels = np.asarray(labels, dtype=int)
    n_neighbors = min(k + 1, matrix.shape[0])
    model = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    model.fit(matrix)
    indices = model.kneighbors(matrix, return_distance=False)
    if indices.shape[1] > 1:
        neighbor_indices = indices[:, 1:]
    else:
        neighbor_indices = indices
    return labels[neighbor_indices].mean(axis=1)


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def representative_label(seed):
    return f"representative_seed_{seed}"


def write_manifest(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest_df = pd.DataFrame(rows)
    if not manifest_df.empty:
        manifest_df = manifest_df.sort_values(["category", "path"]).reset_index(drop=True)
    manifest_df.to_csv(path, sep="\t", index=False)
    return manifest_df
