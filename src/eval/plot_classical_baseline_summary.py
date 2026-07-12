from __future__ import annotations

import argparse
import json
from math import pi
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D


BAR_METRICS = [
    ("AUROC", "AUROC_mean", "AUROC_std"),
    ("AUPRC", "AUPRC_mean", "AUPRC_std"),
    ("MCC", "MCC_mean", "MCC_std"),
    ("Precision", "Precision_mean", "Precision_std"),
    ("Specificity", "Specificity_mean", "Specificity_std"),
]
RADAR_METRICS = [
    ("AUROC", "AUROC_mean"),
    ("AUPRC", "AUPRC_mean"),
    ("MCC", "MCC_mean"),
    ("F1", "F1_mean"),
    ("Recall", "Recall_mean"),
]
RADAR_SCORE_COLUMNS = [column for _, column in RADAR_METRICS]
REQUIRED_ID_COLUMNS = ["Target", "Species", "Regime", "Model", "Feature_Setting"]
REQUIRED_OPTIONAL_ID_COLUMNS = ["Label_Regime", "Split_Version", "Runs", "Seed_List"]

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"]
plt.rcParams["axes.titlecolor"] = "black"
plt.rcParams["axes.labelcolor"] = "black"
plt.rcParams["xtick.color"] = "black"
plt.rcParams["ytick.color"] = "black"

PANEL_TITLE_FONT = FontProperties(family="Arial", weight="bold", style="italic")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Figure1 publication-style plots and short summaries from the frozen protocol benchmark."
    )
    parser.add_argument("--summary-tsv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--group-config", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=str)
    return parser.parse_args()


def load_group_config(config_path: Path) -> dict[str, object]:
    with config_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_summary(summary_path: Path, target_order: list[str]) -> pd.DataFrame:
    df = pd.read_csv(summary_path, sep="\t")
    required_columns = REQUIRED_ID_COLUMNS + REQUIRED_OPTIONAL_ID_COLUMNS
    for column in required_columns:
        if column not in df.columns:
            raise ValueError(f"Missing required column '{column}' in {summary_path}")

    needed_columns = sorted(
        {mean_column for _, mean_column, _ in BAR_METRICS}
        | {std_column for _, _, std_column in BAR_METRICS}
        | {column for _, column in RADAR_METRICS}
    )
    missing_after_fill = [column for column in needed_columns if column not in df.columns]
    if missing_after_fill:
        raise ValueError(
            "Summary table is missing required plotting columns: " + ", ".join(missing_after_fill)
        )

    df = df[df["Target"].isin(target_order)].copy()
    df["Target"] = pd.Categorical(df["Target"], categories=target_order, ordered=True)
    return df.sort_values(["Target", "Species", "Regime", "Model", "Feature_Setting"], kind="stable").reset_index(drop=True)


def style_bar_axis(ax: plt.Axes, y_limits: tuple[float, float]) -> None:
    ax.set_facecolor("#E5E5E5")
    ax.grid(color="white", linestyle="-", linewidth=1.0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#B0B0B0")
    ax.spines["bottom"].set_color("#B0B0B0")
    ax.axhline(0.0, color="#6F6F6F", linewidth=0.8)
    ax.set_ylim(*y_limits)


def clip_bar_stat(mean_value: float, std_value: float) -> tuple[float, np.ndarray]:
    displayed_mean = min(max(mean_value, 0.0), 1.0)
    lower_end = min(max(mean_value - std_value, 0.0), 1.0)
    upper_end = min(max(mean_value + std_value, 0.0), 1.0)
    lower_err = max(displayed_mean - lower_end, 0.0)
    upper_err = max(upper_end - displayed_mean, 0.0)
    return displayed_mean, np.array([[lower_err], [upper_err]])


def bar_limits() -> tuple[float, float]:
    return 0.0, 1.0


def style_radar_axis(ax: plt.Axes, radial_min: float) -> np.ndarray:
    angles = np.linspace(0, 2 * pi, len(RADAR_METRICS), endpoint=False)
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels([label for label, _ in RADAR_METRICS], fontsize=11)
    ax.set_ylim(radial_min, 1.0)
    ax.set_facecolor("white")
    ax.grid(color="#D7D7D7", linewidth=0.8)
    ax.spines["polar"].set_color("#BFBFBF")
    yticks = [0.0, 0.25, 0.5, 0.75, 1.0] if radial_min < 0 else [0.2, 0.4, 0.6, 0.8, 1.0]
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{tick:.2f}" for tick in yticks], fontsize=8, color="#6A6A6A")
    return angles


def target_rows(group_df: pd.DataFrame, target: str) -> pd.DataFrame:
    return group_df[group_df["Target"] == target].copy()


def available_models_for_target(target_df: pd.DataFrame, ordered_models: list[str]) -> list[str]:
    available = set(target_df["Model"].astype(str).tolist())
    return [model for model in ordered_models if model in available]


def draw_missing_text(ax: plt.Axes, missing_models: list[str]) -> None:
    ax.text(
        0.5,
        0.5,
        "No data" if not missing_models else f"No data: {', '.join(missing_models)}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=10,
        color="#555555",
    )


def save_barplot(
    group_name: str,
    group_df: pd.DataFrame,
    output_path: Path,
    group_models: list[str],
    colors: dict[str, str],
    target_order: list[str],
    target_positions: dict[str, int],
    target_labels: dict[str, str],
) -> list[str]:
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    axes = axes.flatten()
    x_positions = np.arange(len(BAR_METRICS))
    missing_notes: list[str] = []

    for target in target_order:
        ax = axes[target_positions[target]]
        style_bar_axis(ax, bar_limits())
        target_df = target_rows(group_df, target)
        present_models = available_models_for_target(target_df, group_models)
        missing_models = [model for model in group_models if model not in present_models]
        if missing_models:
            missing_notes.append(f"{target}: {', '.join(missing_models)}")
        if not present_models:
            draw_missing_text(ax, missing_models)
            ax.set_xticks(x_positions)
            ax.set_xticklabels([label for label, _, _ in BAR_METRICS], rotation=18, ha="right")
            ax.set_title(target_labels[target], fontsize=12, fontproperties=PANEL_TITLE_FONT, pad=10)
            continue

        bar_width = 0.18 if len(present_models) >= 4 else 0.22
        ordered_df = target_df.set_index("Model").loc[present_models].reset_index()
        for index, model in enumerate(present_models):
            raw_means = [
                float(ordered_df.loc[ordered_df["Model"] == model, mean_column].iloc[0])
                for _, mean_column, _ in BAR_METRICS
            ]
            raw_stds = [
                float(ordered_df.loc[ordered_df["Model"] == model, std_column].iloc[0])
                for _, _, std_column in BAR_METRICS
            ]
            displayed_means = []
            lower_errors = []
            upper_errors = []
            for raw_mean, raw_std in zip(raw_means, raw_stds):
                displayed_mean, errors = clip_bar_stat(raw_mean, raw_std)
                displayed_means.append(displayed_mean)
                lower_errors.append(float(errors[0, 0]))
                upper_errors.append(float(errors[1, 0]))
            offsets = x_positions - (bar_width * (len(present_models) - 1) / 2.0) + index * bar_width
            ax.bar(
                offsets,
                displayed_means,
                width=bar_width,
                yerr=np.array([lower_errors, upper_errors]),
                capsize=3,
                color=colors[model],
                edgecolor="white",
                linewidth=0.8,
                label=model,
            )
        ax.set_xticks(x_positions)
        ax.set_xticklabels([label for label, _, _ in BAR_METRICS], rotation=18, ha="right")
        ax.set_title(target_labels[target], fontsize=12, fontproperties=PANEL_TITLE_FONT, pad=10)

    handles = [Line2D([0], [0], color=colors[model], lw=8, label=model) for model in group_models]
    legend = fig.legend(handles=handles, loc="upper center", ncol=len(group_models), bbox_to_anchor=(0.5, 0.985), frameon=True)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("#C9C9C9")
    fig.text(0.04, 0.5, "Performance", va="center", rotation="vertical", fontsize=12)
    fig.tight_layout(rect=(0.05, 0.04, 1.0, 0.92))
    fig.savefig(output_path.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return missing_notes


def save_radarplot(
    group_df: pd.DataFrame,
    output_path: Path,
    group_models: list[str],
    colors: dict[str, str],
    target_order: list[str],
    target_positions: dict[str, int],
    target_labels: dict[str, str],
) -> list[str]:
    metric_values = []
    for _, column in RADAR_METRICS:
        metric_values.extend(group_df[column].dropna().tolist())
    radial_min = min(-0.1, np.floor((min(metric_values) - 0.02) * 10.0) / 10.0) if metric_values else -0.1

    fig, axes = plt.subplots(2, 3, figsize=(16, 10), subplot_kw={"polar": True})
    axes = axes.flatten()
    missing_notes: list[str] = []

    for target in target_order:
        ax = axes[target_positions[target]]
        angles = style_radar_axis(ax, radial_min)
        target_df = target_rows(group_df, target)
        present_models = available_models_for_target(target_df, group_models)
        missing_models = [model for model in group_models if model not in present_models]
        if missing_models:
            missing_notes.append(f"{target}: {', '.join(missing_models)}")
        if not present_models:
            draw_missing_text(ax, missing_models)
            ax.set_title(target_labels[target], fontsize=12, fontproperties=PANEL_TITLE_FONT, pad=18)
            continue

        ordered_df = target_df.set_index("Model").loc[present_models].reset_index()
        closed_angles = np.append(angles, angles[0])
        for model in present_models:
            values = [float(ordered_df.loc[ordered_df["Model"] == model, column].iloc[0]) for _, column in RADAR_METRICS]
            closed_values = np.append(values, values[0])
            ax.plot(closed_angles, closed_values, linewidth=2.2, color=colors[model], label=model)
            ax.fill(closed_angles, closed_values, color=colors[model], alpha=0.08)
        if missing_models:
            ax.text(0.5, -0.15, f"Missing: {', '.join(missing_models)}", transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#555555")
        ax.set_title(target_labels[target], fontsize=12, fontproperties=PANEL_TITLE_FONT, pad=18)

    handles = [Line2D([0], [0], color=colors[model], lw=2.4, label=model) for model in group_models]
    legend = fig.legend(handles=handles, loc="upper center", ncol=len(group_models), bbox_to_anchor=(0.5, 0.985), frameon=True)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("#C9C9C9")
    fig.tight_layout(rect=(0.02, 0.02, 1.0, 0.92))
    fig.savefig(output_path.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return missing_notes


def mean_over_targets(group_df: pd.DataFrame, model: str, columns: list[str]) -> float:
    model_df = group_df[group_df["Model"] == model]
    if model_df.empty:
        return float("nan")
    return float(model_df[columns].mean(axis=1).mean())


def describe_label_shift(group_df: pd.DataFrame, columns: list[str]) -> str:
    old_df = group_df[group_df["Target"] == "fgraminearum_oldlabel"]
    new_df = group_df[group_df["Target"] == "fgraminearum_newlabel"]
    if old_df.empty or new_df.empty:
        return "F. graminearum 新旧标签对比数据不完整，未生成标签迁移描述。"
    old_score = float(old_df[columns].mean(axis=1).mean())
    new_score = float(new_df[columns].mean(axis=1).mean())
    delta = new_score - old_score
    direction = "整体更强" if delta > 0.02 else "整体更弱" if delta < -0.02 else "整体接近"
    changed_metric = max(columns, key=lambda column: abs(float(new_df[column].mean()) - float(old_df[column].mean())))
    changed_metric_label = changed_metric.replace("_mean", "")
    return f"F. graminearum 新标签相对旧标签{direction}（平均差值 {delta:+.3f}），变化最明显的是 `{changed_metric_label}`。"


def describe_tradeoff(group_df: pd.DataFrame, models: list[str]) -> str:
    precision_scores = {model: mean_over_targets(group_df, model, ["Precision_mean"]) for model in models if not group_df[group_df["Model"] == model].empty}
    recall_scores = {model: mean_over_targets(group_df, model, ["Recall_mean"]) for model in models if not group_df[group_df["Model"] == model].empty}
    specificity_scores = {model: mean_over_targets(group_df, model, ["Specificity_mean"]) for model in models if not group_df[group_df["Model"] == model].empty}
    if not precision_scores or not recall_scores or not specificity_scores:
        return "当前组缺少完整指标，未生成 Precision/Recall/Specificity 权衡解读。"
    best_precision = max(precision_scores, key=precision_scores.get)
    best_recall = max(recall_scores, key=recall_scores.get)
    best_specificity = max(specificity_scores, key=specificity_scores.get)
    if len({best_precision, best_recall, best_specificity}) == 1:
        return f"`{best_precision}` 在 Precision / Recall / Specificity 上都保持领先，整体最均衡。"
    return (
        f"存在明显权衡：`{best_recall}` 的 Recall 最强，`{best_specificity}` 的 Specificity 最高，"
        f"`{best_precision}` 的 Precision 最好。"
    )


def write_summary(
    group_id: str,
    group_label: str,
    group_df: pd.DataFrame,
    output_path: Path,
    group_models: list[str],
    missing_notes: list[str],
) -> None:
    observed_models = [model for model in group_models if not group_df[group_df["Model"] == model].empty]
    overall_scores = {model: mean_over_targets(group_df, model, RADAR_SCORE_COLUMNS) for model in observed_models}
    strongest_model = max(overall_scores, key=overall_scores.get) if overall_scores else "NA"
    lines = [
        f"# Figure1 Plot Summary: {group_label}",
        "",
        f"- 组 ID：`{group_id}`",
        f"- 比较模型：{', '.join(group_models)}",
        f"- 跨目标按 AUROC / AUPRC / MCC / F1 / Recall 的平均表现看，`{strongest_model}` 综合均值最高。",
        f"- {describe_label_shift(group_df, RADAR_SCORE_COLUMNS)}",
        f"- {describe_tradeoff(group_df, group_models)}",
        "- Barplot 仅调整展示层：y 轴下界固定为 0，负均值按 0 绘制，下误差条裁到 0；原始 TSV 数值不做修改。",
    ]
    if missing_notes:
        lines.append(f"- 缺失模型面板或局部缺失：{'; '.join(missing_notes)}")
    lines.append("")
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    group_config = load_group_config(args.group_config)
    target_order = list(group_config["target_order"])
    target_positions = {str(key): int(value) for key, value in group_config["target_positions"].items()}
    target_labels = {str(key): str(value) for key, value in group_config["target_labels"].items()}
    summary_df = load_summary(args.summary_tsv, target_order)

    for group in group_config["groups"]:
        group_id = str(group["id"])
        group_label = str(group["label"])
        group_models = [str(model) for model in group["models"]]
        colors = {str(model): str(color) for model, color in group["colors"].items()}
        group_df = summary_df[summary_df["Model"].isin(group_models)].copy()
        base_name = f"{args.output_prefix}_{group_id}"
        bar_missing = save_barplot(
            group_label,
            group_df,
            output_dir / f"{base_name}_barplot_2x3",
            group_models,
            colors,
            target_order,
            target_positions,
            target_labels,
        )
        radar_missing = save_radarplot(
            group_df,
            output_dir / f"{base_name}_radar_2x3",
            group_models,
            colors,
            target_order,
            target_positions,
            target_labels,
        )
        write_summary(
            group_id,
            group_label,
            group_df,
            output_dir / f"{base_name}_plot_summary.md",
            group_models,
            sorted(set(bar_missing + radar_missing)),
        )


if __name__ == "__main__":
    main()
