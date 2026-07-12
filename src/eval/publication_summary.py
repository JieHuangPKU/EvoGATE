from __future__ import annotations

from typing import Iterable

import pandas as pd


PUBLICATION_METRICS = [
    ("test_auroc", "AUROC"),
    ("test_auprc", "AUPRC"),
    ("test_mcc", "MCC"),
    ("test_f1", "F1"),
    ("test_precision", "Precision"),
    ("test_recall", "Recall"),
    ("test_specificity", "Specificity"),
]

DISPLAY_COLUMNS = [
    ("target", "Target"),
    ("protocol", "Target"),
    ("species", "Species"),
    ("regime", "Regime"),
    ("model", "Model"),
    ("feature_setting", "Feature_Setting"),
    ("esm2_dim", "ESM2_Dim"),
    ("label_regime", "Label_Regime"),
    ("split_version", "Split_Version"),
    ("n_runs", "Runs"),
    ("seed_list", "Seed_List"),
]

MARKDOWN_ID_COLUMNS = ["Target", "Model", "Feature_Setting", "Runs"]


def format_float(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def ensure_columns(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    df = frame.copy()
    for column in columns:
        if column not in df.columns:
            df[column] = pd.NA
    return df


def build_publication_summary(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    df = aggregated_df.copy()
    if df.empty:
        columns = [display for _, display in DISPLAY_COLUMNS]
        for _, public_name in PUBLICATION_METRICS:
            columns.extend([f"{public_name}_mean", f"{public_name}_std"])
        return pd.DataFrame(columns=columns)

    selected_columns: list[str] = []
    rename_map: dict[str, str] = {}
    seen_display_names: set[str] = set()
    for source_name, display_name in DISPLAY_COLUMNS:
        if source_name in df.columns and display_name not in seen_display_names:
            selected_columns.append(source_name)
            rename_map[source_name] = display_name
            seen_display_names.add(display_name)

    for technical_name, public_name in PUBLICATION_METRICS:
        mean_column = f"{technical_name}_mean"
        std_column = f"{technical_name}_std"
        if mean_column in df.columns:
            selected_columns.append(mean_column)
            rename_map[mean_column] = f"{public_name}_mean"
        if std_column in df.columns:
            selected_columns.append(std_column)
            rename_map[std_column] = f"{public_name}_std"

    publication_df = df[selected_columns].rename(columns=rename_map)
    sort_columns = [column for column in ["Target", "Model", "Feature_Setting"] if column in publication_df.columns]
    if sort_columns:
        publication_df = publication_df.sort_values(sort_columns, kind="stable").reset_index(drop=True)
    return publication_df


def publication_markdown(publication_df: pd.DataFrame, title: str, intro: str) -> str:
    lines = [f"# {title}", "", intro, ""]
    if publication_df.empty:
        lines.extend(["No benchmark rows available.", ""])
        return "\n".join(lines)

    display_df = publication_df.copy()
    metric_labels: list[str] = []
    for _, public_name in PUBLICATION_METRICS:
        mean_column = f"{public_name}_mean"
        std_column = f"{public_name}_std"
        if mean_column in display_df.columns:
            display_df[public_name] = display_df.apply(
                lambda row: f"{format_float(row.get(mean_column))} ± {format_float(row.get(std_column))}",
                axis=1,
            )
            metric_labels.append(public_name)

    drop_columns = [
        column
        for _, public_name in PUBLICATION_METRICS
        for column in [f"{public_name}_mean", f"{public_name}_std"]
        if column in display_df.columns
    ]
    display_df = display_df.drop(columns=drop_columns)

    ordered_columns = [column for column in MARKDOWN_ID_COLUMNS if column in display_df.columns]
    ordered_columns.extend(label for label in metric_labels if label in display_df.columns)
    trailing_columns = [column for column in display_df.columns if column not in ordered_columns]
    display_df = display_df[ordered_columns + trailing_columns]

    lines.append(display_df.to_markdown(index=False))
    lines.append("")
    return "\n".join(lines)
