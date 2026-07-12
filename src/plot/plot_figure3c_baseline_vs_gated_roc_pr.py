from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BOOL_TRUE = {"true", "t", "1", "yes", "y"}
CURVE_GRID = np.linspace(0.0, 1.0, 501)
COLORS = {
    "ORT_EXP_SUB_ESM2": "#332288",
    "ORT_EXP_SUB_ESM2_GATED": "#CC6677",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Figure3c ROC/PR comparison between baseline and gated predictions.")
    parser.add_argument("--figure3c-root", default="outputs/Figure3c/fgraminearum_newlabel/GraphSAGE", type=str)
    parser.add_argument("--output-dir", default="results/Figure3c_threshold_tuned/plots", type=str)
    parser.add_argument("--baseline-setting", default="ORT_EXP_SUB_ESM2", type=str)
    parser.add_argument("--gated-setting", default="ORT_EXP_SUB_ESM2_GATED", type=str)
    parser.add_argument("--mode", default="pooled", choices=["pooled", "mean_seed"], type=str)
    return parser.parse_args()


def parse_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(BOOL_TRUE)


def roc_curve_points(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    order = np.argsort(-scores, kind="mergesort")
    labels_sorted = labels[order]
    positives = int(np.sum(labels_sorted == 1))
    negatives = int(np.sum(labels_sorted == 0))
    if positives == 0 or negatives == 0:
        raise ValueError("ROC requires both positive and negative labels.")
    tp = np.cumsum(labels_sorted == 1)
    fp = np.cumsum(labels_sorted == 0)
    tpr = np.concatenate(([0.0], tp / positives, [1.0]))
    fpr = np.concatenate(([0.0], fp / negatives, [1.0]))
    auc = float(np.trapz(tpr, fpr))
    return fpr, tpr, auc


def pr_curve_points(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    order = np.argsort(-scores, kind="mergesort")
    labels_sorted = labels[order]
    positives = int(np.sum(labels_sorted == 1))
    if positives == 0:
        raise ValueError("PR requires at least one positive label.")
    tp = np.cumsum(labels_sorted == 1)
    fp = np.cumsum(labels_sorted == 0)
    recall = tp / positives
    precision = tp / np.maximum(tp + fp, 1)
    recall_curve = np.concatenate(([0.0], recall))
    precision_curve = np.concatenate(([1.0], precision))
    auprc = float(np.sum(np.diff(recall_curve) * precision_curve[1:]))
    return recall_curve, precision_curve, auprc


def interpolate_step_curve(x: np.ndarray, y: np.ndarray, grid: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(x, grid, side="right") - 1
    indices = np.clip(indices, 0, len(y) - 1)
    return y[indices]


def load_feature_setting_runs(figure_root: Path, feature_setting: str) -> list[pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    split_versions: set[str] = set()
    for predictions_path in sorted((figure_root / feature_setting).glob("run_*/predictions.tsv")):
        df = pd.read_csv(predictions_path, sep="\t")
        df = df.loc[
            (df["split"] == "test")
            & parse_bool_series(df["is_labeled"])
        ].copy()
        df["label"] = pd.to_numeric(df["label"], errors="coerce")
        df["pred_score"] = pd.to_numeric(df["pred_score"], errors="coerce")
        df = df.dropna(subset=["label", "pred_score"])
        if df.empty:
            raise ValueError(f"No labeled test rows found in {predictions_path}")
        split_versions.update(df["split_version"].dropna().astype(str).unique().tolist())
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No predictions.tsv files found for {feature_setting}")
    if len(split_versions) != 1:
        raise ValueError(f"Expected one split_version for {feature_setting}, found {sorted(split_versions)}")
    return frames


def pooled_curve(run_frames: list[pd.DataFrame], curve_kind: str) -> tuple[np.ndarray, np.ndarray, float]:
    pooled = pd.concat(run_frames, ignore_index=True)
    labels = pooled["label"].to_numpy(dtype=int)
    scores = pooled["pred_score"].to_numpy(dtype=float)
    if curve_kind == "roc":
        return roc_curve_points(labels, scores)
    return pr_curve_points(labels, scores)


def mean_seed_curve(run_frames: list[pd.DataFrame], curve_kind: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    curves = []
    metrics = []
    for frame in run_frames:
        labels = frame["label"].to_numpy(dtype=int)
        scores = frame["pred_score"].to_numpy(dtype=float)
        if curve_kind == "roc":
            x, y, metric = roc_curve_points(labels, scores)
        else:
            x, y, metric = pr_curve_points(labels, scores)
        curves.append(interpolate_step_curve(x, y, CURVE_GRID))
        metrics.append(metric)
    stacked = np.vstack(curves)
    return CURVE_GRID, stacked.mean(axis=0), stacked.std(axis=0), float(np.mean(metrics)), float(np.std(metrics, ddof=0))


def plot_curves(args: argparse.Namespace) -> None:
    figure_root = Path(args.figure3c_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_settings = [args.baseline_setting, args.gated_setting]
    loaded = {setting: load_feature_setting_runs(figure_root, setting) for setting in feature_settings}

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    roc_fig, roc_ax = plt.subplots(figsize=(6.2, 5.2))
    pr_fig, pr_ax = plt.subplots(figsize=(6.2, 5.2))

    for setting in feature_settings:
        color = COLORS.get(setting, "#333333")
        label_name = setting
        if args.mode == "pooled":
            roc_x, roc_y, roc_auc = pooled_curve(loaded[setting], "roc")
            pr_x, pr_y, pr_auc = pooled_curve(loaded[setting], "pr")
            roc_ax.plot(roc_x, roc_y, lw=2.2, color=color, label=f"{label_name} (AUROC={roc_auc:.3f})")
            pr_ax.plot(pr_x, pr_y, lw=2.2, color=color, label=f"{label_name} (AUPRC={pr_auc:.3f})")
        else:
            roc_x, roc_mean, roc_std, roc_auc_mean, roc_auc_std = mean_seed_curve(loaded[setting], "roc")
            pr_x, pr_mean, pr_std, pr_auc_mean, pr_auc_std = mean_seed_curve(loaded[setting], "pr")
            roc_ax.plot(roc_x, roc_mean, lw=2.2, color=color, label=f"{label_name} (AUROC={roc_auc_mean:.3f}±{roc_auc_std:.3f})")
            roc_ax.fill_between(roc_x, np.clip(roc_mean - roc_std, 0, 1), np.clip(roc_mean + roc_std, 0, 1), color=color, alpha=0.18)
            pr_ax.plot(pr_x, pr_mean, lw=2.2, color=color, label=f"{label_name} (AUPRC={pr_auc_mean:.3f}±{pr_auc_std:.3f})")
            pr_ax.fill_between(pr_x, np.clip(pr_mean - pr_std, 0, 1), np.clip(pr_mean + pr_std, 0, 1), color=color, alpha=0.18)

    roc_ax.plot([0, 1], [0, 1], linestyle="--", color="#777777", lw=1.2)
    roc_ax.set_xlim(0, 1)
    roc_ax.set_ylim(0, 1)
    roc_ax.set_xlabel("False Positive Rate")
    roc_ax.set_ylabel("True Positive Rate")
    roc_ax.set_title(f"Figure3c ROC ({args.mode})")
    roc_ax.legend(frameon=False, loc="lower right")
    roc_ax.grid(alpha=0.18, linewidth=0.6)

    pr_ax.set_xlim(0, 1)
    pr_ax.set_ylim(0, 1)
    pr_ax.set_xlabel("Recall")
    pr_ax.set_ylabel("Precision")
    pr_ax.set_title(f"Figure3c PR ({args.mode})")
    pr_ax.legend(frameon=False, loc="upper right")
    pr_ax.grid(alpha=0.18, linewidth=0.6)

    roc_stub = "Figure3c_baseline_vs_gated_roc" if args.mode == "pooled" else "Figure3c_baseline_vs_gated_roc_mean_seed"
    pr_stub = "Figure3c_baseline_vs_gated_pr" if args.mode == "pooled" else "Figure3c_baseline_vs_gated_pr_mean_seed"
    roc_fig.tight_layout()
    pr_fig.tight_layout()
    roc_fig.savefig(output_dir / f"{roc_stub}.pdf")
    roc_fig.savefig(output_dir / f"{roc_stub}.png", dpi=300)
    pr_fig.savefig(output_dir / f"{pr_stub}.pdf")
    pr_fig.savefig(output_dir / f"{pr_stub}.png", dpi=300)
    plt.close(roc_fig)
    plt.close(pr_fig)


def main() -> None:
    args = parse_args()
    plot_curves(args)


if __name__ == "__main__":
    main()
