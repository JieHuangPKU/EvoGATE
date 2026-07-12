from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from palettable.wesanderson import Darjeeling4_5

from src.eval.aggregate_frozen_protocol_runs import METRIC_COLUMNS, collect_metrics
from src.eval.publication_summary import build_publication_summary


FEATURE_ORDER = ["ESM2", "ORT_ESM2", "ORT_EXP_SUB", "ORT_EXP_SUB_ESM2"]
PROTOCOL_TO_SPECIES = {
    "fgraminearum_newlabel": "fgraminearum",
    "scerevisiae": "scerevisiae",
}
SPECIES_ORDER = ["fgraminearum", "scerevisiae"]
SPECIES_TITLES = {
    "fgraminearum": "F. graminearum",
    "scerevisiae": "S. cerevisiae",
}
METRIC_SPECS = [
    ("AUPRC", "AUPRC_mean", "AUPRC_std", "Figure3a publication comparison"),
    ("MCC", "MCC_mean", "MCC_std", "Figure3a publication comparison"),
]
FEATURE_COLORS = dict(zip(FEATURE_ORDER, Darjeeling4_5.hex_colors[:4]))
SUMMARY_GROUP_COLUMNS = [
    "protocol",
    "species",
    "regime",
    "model",
    "feature_setting",
    "label_regime",
    "split_version",
    "esm2_dim",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate audited Figure3a summaries and plots.")
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--results-root", required=True, type=str)
    parser.add_argument("--prefix", default="Figure3a", type=str)
    return parser.parse_args()


def configure_plotting() -> None:
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["font.size"] = 9
    plt.rcParams["axes.titlesize"] = 10
    plt.rcParams["axes.labelsize"] = 9
    plt.rcParams["xtick.labelsize"] = 8
    plt.rcParams["ytick.labelsize"] = 8
    plt.rcParams["legend.fontsize"] = 8
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1.0
    plt.rcParams["xtick.color"] = "black"
    plt.rcParams["ytick.color"] = "black"
    plt.rcParams["axes.grid"] = False
    sns.set_theme(style="white")


def apply_axis_style(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.grid(False)
    ax.xaxis.grid(False, which="both")
    ax.yaxis.grid(False, which="both")
    ax.minorticks_off()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["left"].set_linewidth(1.0)
    ax.tick_params(axis="both", colors="black", width=0.9, length=3.5)


def save_figure(fig: plt.Figure, output_stub: Path) -> None:
    output_stub.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stub.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    fig.savefig(output_stub.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def normalize_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out = out[out["protocol"].isin(PROTOCOL_TO_SPECIES)].copy()
    out["species"] = out["protocol"].map(PROTOCOL_TO_SPECIES)
    out["feature_setting"] = out["feature_setting"].astype(str).str.strip().str.upper()
    out["model"] = out["model"].astype(str).str.strip()
    out = out[out["model"].eq("GraphSAGE")].copy()
    out = out[out["feature_setting"].isin(FEATURE_ORDER)].copy()
    out["feature_setting"] = pd.Categorical(out["feature_setting"], categories=FEATURE_ORDER, ordered=True)
    out["species"] = pd.Categorical(out["species"], categories=SPECIES_ORDER, ordered=True)
    sort_columns = [column for column in ["species", "feature_setting", "seed"] if column in out.columns]
    return out.sort_values(sort_columns, kind="stable").reset_index(drop=True)


def validate_per_run_coverage(per_run: pd.DataFrame) -> None:
    problems: list[str] = []
    for species in SPECIES_ORDER:
        for feature in FEATURE_ORDER:
            subset = per_run[
                (per_run["species"].astype(str) == species)
                & (per_run["feature_setting"].astype(str) == feature)
            ].copy()
            if subset.empty:
                problems.append(f"missing per-run rows for {species}/{feature}")
                continue
            seeds = sorted(subset["seed"].astype(str).tolist())
            if len(seeds) != len(set(seeds)):
                problems.append(f"duplicated seeds for {species}/{feature}: {','.join(seeds)}")
    if problems:
        raise RuntimeError("Figure3a per-run coverage audit failed:\n- " + "\n- ".join(problems))


def build_run_manifest(per_run: pd.DataFrame) -> pd.DataFrame:
    keep_columns = [
        "protocol",
        "species",
        "regime",
        "model",
        "feature_setting",
        "label_regime",
        "split_version",
        "run_id",
        "seed",
        "test_auprc",
        "test_mcc",
        "metrics_path",
    ]
    manifest = per_run.copy()
    for column in keep_columns:
        if column not in manifest.columns:
            manifest[column] = pd.NA
    manifest = manifest[keep_columns].rename(
        columns={
            "protocol": "Target",
            "species": "Species",
            "regime": "Regime",
            "model": "Model",
            "feature_setting": "Feature_Setting",
            "label_regime": "Label_Regime",
            "split_version": "Split_Version",
            "run_id": "Run_ID",
            "seed": "Seed",
            "test_auprc": "AUPRC",
            "test_mcc": "MCC",
            "metrics_path": "Metrics_Path",
        }
    )
    return manifest.reset_index(drop=True)


def build_feature_audit(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species in SPECIES_ORDER:
        for feature in FEATURE_ORDER:
            subset = summary_df[
                (summary_df["Species"].astype(str) == species)
                & (summary_df["Feature_Setting"].astype(str) == feature)
            ].copy()
            present = not subset.empty
            duplicated = len(subset) > 1
            if present:
                row = subset.iloc[0]
                mean_available = pd.notna(row.get("AUPRC_mean")) and pd.notna(row.get("MCC_mean"))
                std_available = pd.notna(row.get("AUPRC_std")) and pd.notna(row.get("MCC_std"))
                runs_available = pd.notna(row.get("Runs")) and pd.notna(row.get("Seed_List")) and str(row.get("Seed_List")).strip() != ""
                na_metric = any(
                    pd.isna(row.get(column))
                    for column in ["AUPRC_mean", "AUPRC_std", "MCC_mean", "MCC_std"]
                )
            else:
                mean_available = False
                std_available = False
                runs_available = False
                na_metric = False

            if not present:
                status = "MISSING_ROW"
            elif duplicated:
                status = "DUPLICATED_ROW"
            elif not mean_available or not std_available or not runs_available:
                status = "MISSING_METRIC"
            elif na_metric:
                status = "NA_METRIC"
            else:
                status = "OK"

            rows.append(
                {
                    "Species": species,
                    "Feature_Setting": feature,
                    "present_in_summary": present,
                    "mean_available": mean_available,
                    "std_available": std_available,
                    "runs_available": runs_available,
                    "status": status,
                }
            )
    audit = pd.DataFrame(rows)
    bad_rows = audit[audit["status"] != "OK"]
    if not bad_rows.empty:
        raise RuntimeError(
            "Figure3a feature coverage audit failed:\n"
            + bad_rows.to_string(index=False)
        )
    return audit


def format_seed_list(value: object) -> str:
    parts = []
    for token in str(value).split(","):
        text = str(token).strip()
        if not text:
            continue
        try:
            numeric = float(text)
        except ValueError:
            parts.append(text)
            continue
        parts.append(str(int(numeric)) if numeric.is_integer() else text)
    return ",".join(parts)


def build_summary_tables(aggregated: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    publication = build_publication_summary(aggregated)
    if publication.empty:
        raise RuntimeError("Figure3a aggregation produced no publication rows")
    summary = publication.copy()
    summary["Species"] = summary["Species"].astype(str)
    summary["Feature_Setting"] = pd.Categorical(summary["Feature_Setting"], categories=FEATURE_ORDER, ordered=True)
    summary["Species"] = pd.Categorical(summary["Species"], categories=SPECIES_ORDER, ordered=True)
    summary["Seed_List"] = summary["Seed_List"].map(format_seed_list)
    summary = summary.sort_values(["Species", "Feature_Setting"], kind="stable").reset_index(drop=True)

    required_columns = [
        "Target",
        "Species",
        "Regime",
        "Model",
        "Feature_Setting",
        "Label_Regime",
        "Split_Version",
        "Runs",
        "Seed_List",
        "AUPRC_mean",
        "AUPRC_std",
        "MCC_mean",
        "MCC_std",
        "AUROC_mean",
        "AUROC_std",
        "F1_mean",
        "F1_std",
        "Precision_mean",
        "Precision_std",
        "Recall_mean",
        "Recall_std",
        "Specificity_mean",
        "Specificity_std",
        "ESM2_Dim",
    ]
    for column in required_columns:
        if column not in summary.columns:
            summary[column] = pd.NA
    plot_data = summary[required_columns].copy()
    return summary, plot_data


def aggregate_figure3a_runs(per_run: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group_key, group_df in per_run.groupby(SUMMARY_GROUP_COLUMNS, dropna=False, sort=True, observed=False):
        row = {column: value for column, value in zip(SUMMARY_GROUP_COLUMNS, group_key)}
        row["n_runs"] = int(len(group_df))
        row["seed_list"] = ",".join(str(value).strip() for value in group_df["seed"].tolist() if str(value).strip())
        for metric in METRIC_COLUMNS:
            values = pd.to_numeric(group_df.get(metric), errors="coerce")
            row[f"{metric}_mean"] = float(values.mean())
            row[f"{metric}_std"] = float(values.std(ddof=0))
        rows.append(row)
    aggregated = pd.DataFrame(rows)
    if aggregated.empty:
        return aggregated
    aggregated["feature_setting"] = pd.Categorical(aggregated["feature_setting"], categories=FEATURE_ORDER, ordered=True)
    aggregated["species"] = pd.Categorical(aggregated["species"], categories=SPECIES_ORDER, ordered=True)
    return aggregated.sort_values(["species", "feature_setting"], kind="stable").reset_index(drop=True)


def species_panel(plot_df: pd.DataFrame, species: str, metric_name: str, mean_col: str, std_col: str, output_stub: Path) -> None:
    subset = plot_df[plot_df["Species"].astype(str) == species].copy()
    if len(subset) != len(FEATURE_ORDER):
        raise RuntimeError(f"Expected 4 rows for {species} panel, found {len(subset)}")
    subset["Feature_Setting"] = pd.Categorical(subset["Feature_Setting"], categories=FEATURE_ORDER, ordered=True)
    subset = subset.sort_values("Feature_Setting", kind="stable").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    x = np.arange(len(FEATURE_ORDER))
    heights = subset[mean_col].to_numpy(dtype=float)
    errors = subset[std_col].to_numpy(dtype=float)
    colors = [FEATURE_COLORS[feature] for feature in FEATURE_ORDER]
    ax.bar(x, heights, yerr=errors, color=colors, edgecolor="black", linewidth=0.8, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(FEATURE_ORDER, rotation=25, ha="right")
    ax.set_ylabel(f"{metric_name} mean +/- seed SD")
    ax.set_title(SPECIES_TITLES[species])
    ax.set_ylim(bottom=0)
    apply_axis_style(ax)
    save_figure(fig, output_stub)


def main_panel(plot_df: pd.DataFrame, output_stub: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 6.8), sharex=False)
    for row_index, species in enumerate(SPECIES_ORDER):
        species_df = plot_df[plot_df["Species"].astype(str) == species].copy()
        species_df["Feature_Setting"] = pd.Categorical(species_df["Feature_Setting"], categories=FEATURE_ORDER, ordered=True)
        species_df = species_df.sort_values("Feature_Setting", kind="stable").reset_index(drop=True)
        if len(species_df) != len(FEATURE_ORDER):
            raise RuntimeError(f"Expected 4 rows for {species} in main panel, found {len(species_df)}")
        x = np.arange(len(FEATURE_ORDER))
        for col_index, (metric_name, mean_col, std_col, _) in enumerate(METRIC_SPECS):
            ax = axes[row_index, col_index]
            heights = species_df[mean_col].to_numpy(dtype=float)
            errors = species_df[std_col].to_numpy(dtype=float)
            colors = [FEATURE_COLORS[feature] for feature in FEATURE_ORDER]
            ax.bar(x, heights, yerr=errors, color=colors, edgecolor="black", linewidth=0.8, capsize=3)
            ax.set_xticks(x)
            ax.set_xticklabels(FEATURE_ORDER, rotation=25, ha="right")
            ax.set_ylabel(f"{metric_name} mean +/- seed SD")
            ax.set_title(f"{SPECIES_TITLES[species]} ({metric_name})")
            ax.set_ylim(bottom=0)
            apply_axis_style(ax)
    fig.tight_layout()
    save_figure(fig, output_stub)


def write_generation_report(
    report_path: Path,
    summary_df: pd.DataFrame,
    audit_df: pd.DataFrame,
    generated_files: list[Path],
) -> None:
    lines = [
        "# Figure3a Generation Report",
        "",
        "## Coverage",
        "",
        audit_df.to_markdown(index=False),
        "",
        "## Final Summary",
        "",
        summary_df[["Species", "Feature_Setting", "AUPRC_mean", "AUPRC_std", "MCC_mean", "MCC_std", "Runs", "Seed_List"]].to_markdown(index=False),
        "",
        "## Generated Files",
        "",
    ]
    lines.extend(f"- {path.as_posix()}" for path in generated_files)
    lines.append("")
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    results_root = Path(args.results_root).resolve()
    data_dir = results_root / "data"
    plots_dir = results_root / "plots"
    summary_dir = results_root / "summary"
    for directory in [data_dir, plots_dir, summary_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    configure_plotting()
    per_run = normalize_frame(collect_metrics(output_root))
    validate_per_run_coverage(per_run)
    aggregated = aggregate_figure3a_runs(per_run)

    summary_df, plot_df = build_summary_tables(aggregated)
    audit_df = build_feature_audit(summary_df)
    run_manifest = build_run_manifest(per_run)

    final_summary_path = data_dir / f"{args.prefix}_final_summary.tsv"
    plot_data_path = data_dir / f"{args.prefix}_plot_data.tsv"
    panel_a_path = data_dir / f"{args.prefix}_panelA_fgraminearum_plot_data.tsv"
    panel_b_path = data_dir / f"{args.prefix}_panelB_scerevisiae_plot_data.tsv"
    run_manifest_path = summary_dir / f"{args.prefix}_run_manifest.tsv"
    audit_path = summary_dir / f"{args.prefix}_feature_coverage_audit.tsv"
    report_path = summary_dir / f"{args.prefix}_generation_report.md"

    summary_df.to_csv(final_summary_path, sep="\t", index=False)
    plot_df.to_csv(plot_data_path, sep="\t", index=False)
    plot_df[plot_df["Species"].astype(str) == "fgraminearum"].to_csv(panel_a_path, sep="\t", index=False)
    plot_df[plot_df["Species"].astype(str) == "scerevisiae"].to_csv(panel_b_path, sep="\t", index=False)
    run_manifest.to_csv(run_manifest_path, sep="\t", index=False)
    audit_df.to_csv(audit_path, sep="\t", index=False)

    species_panel(plot_df, "fgraminearum", "AUPRC", "AUPRC_mean", "AUPRC_std", plots_dir / f"{args.prefix}_panelA_fgraminearum_auprc")
    species_panel(plot_df, "fgraminearum", "MCC", "MCC_mean", "MCC_std", plots_dir / f"{args.prefix}_panelA_fgraminearum_mcc")
    species_panel(plot_df, "scerevisiae", "AUPRC", "AUPRC_mean", "AUPRC_std", plots_dir / f"{args.prefix}_panelB_scerevisiae_auprc")
    species_panel(plot_df, "scerevisiae", "MCC", "MCC_mean", "MCC_std", plots_dir / f"{args.prefix}_panelB_scerevisiae_mcc")
    main_panel(plot_df, plots_dir / f"{args.prefix}_main_panels")

    generated_files = sorted(path for path in results_root.rglob("*") if path.is_file())
    write_generation_report(report_path, summary_df, audit_df, generated_files)


if __name__ == "__main__":
    main()
