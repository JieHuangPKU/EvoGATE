import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from utils_inparanoid_orthoxml import (
    build_binary_ortholog_matrix,
    download_orthoxml_files,
    extract_species_from_filename,
    find_labels_file,
    list_orthoxml_links,
    load_labels_table,
    load_string_aliases,
    map_focal_ids_to_gene_ids,
    parse_orthoxml_file,
)


PROGATE_ROOT = Path("/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2")
PROCESSED_ROOT = PROGATE_ROOT / "data" / "processed"
OR_ROOT = PROCESSED_ROOT / "OR"
DOWNLOAD_ROOT = OR_ROOT / "Download" / "InParanoid8"
PPI_ALIAS_ROOT = PROCESSED_ROOT / "PPI" / "stringDB"
SRC_ROOT = PROGATE_ROOT / "src" / "data"

FOCAL_SPECIES = {
    "fgraminearum": {"folder": "G.zeae", "taxid": "229533"},
    "celegans": {"folder": "C.elegans", "taxid": "6239"},
    "scerevisiae": {"folder": "S.cerevisiae", "taxid": "4932"},
    "melanogaster": {"folder": "D.melanogaster", "taxid": "7227"},
    "human": {"folder": "H.sapiens", "taxid": "9606"},
}
INPARANOID_BASE_URL = "https://inparanoid8.sbc.su.se/download/8.0_current/Orthologs_OrthoXML"


def load_species_taxonomy(or_root: Path) -> pd.DataFrame:
    candidate_paths = [
        or_root / "inparanoid_species_taxonomy.tsv",
        or_root / "Download" / "InParanoiDB9_species.txt",
    ]
    for path in candidate_paths:
        if not path.exists():
            continue
        if path.suffix == ".tsv":
            taxonomy_df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        else:
            taxonomy_df = pd.read_csv(path, sep=",", comment="#", dtype=str).fillna("")
            taxonomy_df.columns = ["taxid", "uniprot_id", "species_name", "domain"]
        taxonomy_df.columns = [column.strip() for column in taxonomy_df.columns]
        if "taxid" not in taxonomy_df.columns:
            continue
        taxonomy_df["taxid"] = taxonomy_df["taxid"].astype(str).str.strip()
        if "species_name" not in taxonomy_df.columns and "scientific_name" in taxonomy_df.columns:
            taxonomy_df["species_name"] = taxonomy_df["scientific_name"]
        for column in ("species_name", "scientific_name"):
            if column in taxonomy_df.columns:
                taxonomy_df[column] = taxonomy_df[column].astype(str).str.strip()
        return taxonomy_df
    return pd.DataFrame(columns=["taxid", "species_name", "scientific_name"])


def build_species_name_taxid_lookup(taxonomy_df: pd.DataFrame) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    for _, row in taxonomy_df.iterrows():
        taxid = str(row.get("taxid", "")).strip()
        if not taxid:
            continue
        for column in ("species_name", "scientific_name"):
            value = str(row.get(column, "")).strip()
            if not value:
                continue
            lookup[value] = taxid
            lookup[value.replace(" ", ".")] = taxid
            parts = value.split()
            if len(parts) >= 2:
                lookup[f"{parts[0][0]}.{parts[1]}"] = taxid
                lookup[f"{parts[0][0]}.{parts[1].split('_')[0]}"] = taxid
    for species, meta in FOCAL_SPECIES.items():
        lookup[meta["folder"]] = meta["taxid"]
        lookup[species] = meta["taxid"]
    return lookup


def resolve_target_column_id(target_species_name: str, target_species_taxid: str, taxonomy_lookup: Dict[str, str]) -> str:
    if target_species_taxid:
        return target_species_taxid
    name = target_species_name.strip()
    return taxonomy_lookup.get(name, name)


def species_outputs_complete(species: str) -> bool:
    outdir = OR_ROOT / species
    required_paths = [
        outdir / "orthologs.csv",
        outdir / "orthologs_summary.tsv",
    ]
    return all(path.exists() and path.stat().st_size > 0 for path in required_paths)


def build_species_ortholog_matrix(
    species: str,
    meta: Dict[str, str],
    taxonomy_lookup: Dict[str, str],
) -> Tuple[Dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    folder = meta["folder"]
    focal_taxid = meta["taxid"]
    folder_url = f"{INPARANOID_BASE_URL}/{folder}/"
    download_dir = DOWNLOAD_ROOT / folder
    print(f"[INFO] species={species} downloading orthoXML files from {folder}")
    listed_links = list_orthoxml_links(folder_url)
    xml_paths = download_orthoxml_files(folder_url, download_dir)
    print(f"[INFO] species={species} downloaded_or_cached={len(xml_paths)}")

    labels_path = find_labels_file(PROCESSED_ROOT, species)
    labels_df = load_labels_table(labels_path)
    valid_gene_ids = labels_df["gene_id"].tolist()

    alias_path = PPI_ALIAS_ROOT / f"{focal_taxid}.protein.aliases.v12.0.txt.gz"
    alias_df: Optional[pd.DataFrame] = None
    used_string_aliases = False
    if alias_path.exists():
        alias_df = load_string_aliases(alias_path)
        used_string_aliases = True

    per_file_rows: List[pd.DataFrame] = []
    target_species_seen = set()
    target_columns_seen = set()
    raw_id_values = set()

    for xml_path in xml_paths:
        focal_folder, target_folder = extract_species_from_filename(xml_path)
        if focal_folder != folder and target_folder != folder:
            print(f"[WARNING] species={species} xml={xml_path.name} skipped due to folder mismatch")
            continue
        parsed_df = parse_orthoxml_file(xml_path)
        if parsed_df.empty:
            print(f"[INFO] parsed xml={xml_path.name} rows=0")
            continue

        focal_mask = parsed_df["focal_species_taxid"].eq(focal_taxid)
        if not focal_mask.any():
            focal_mask = parsed_df["source_file"].eq(xml_path.name) & (
                parsed_df["focal_species_name"].str.contains(folder.split(".")[-1], regex=False)
            )
        focal_df = parsed_df[focal_mask].copy()
        if focal_df.empty:
            print(f"[WARNING] species={species} xml={xml_path.name} has no focal rows after filtering")
            continue

        focal_df["target_column_id"] = focal_df.apply(
            lambda row: resolve_target_column_id(row["target_species_name"], row["target_species_taxid"], taxonomy_lookup),
            axis=1,
        )
        raw_id_values.update(focal_df["focal_raw_id"].dropna().astype(str))
        target_species_seen.update(
            focal_df[["target_species_name", "target_species_taxid", "target_column_id"]].drop_duplicates().itertuples(index=False, name=None)
        )
        target_columns_seen.update(focal_df["target_column_id"].dropna().astype(str))
        print(f"[INFO] parsed xml={xml_path.name} rows={len(focal_df)}")
        per_file_rows.append(focal_df)

    all_pairs_df = pd.concat(per_file_rows, ignore_index=True) if per_file_rows else pd.DataFrame()
    mapping_df = map_focal_ids_to_gene_ids(
        all_pairs_df["focal_raw_id"] if not all_pairs_df.empty else pd.Series(dtype="string"),
        labels_df,
        alias_df,
        species,
    )
    mapping_summary = mapping_df["mapping_status"].value_counts().to_dict()
    print(
        "[INFO] species=%s mapped focal ids=%s unmapped=%s ambiguous=%s"
        % (
            species,
            mapping_summary.get("mapped", 0),
            mapping_summary.get("unmapped", 0),
            mapping_summary.get("ambiguous", 0),
        )
    )

    positive_long_df = pd.DataFrame(
        columns=[
            "species",
            "focal_gene_id",
            "target_species_name",
            "target_species_taxid",
            "target_column_id",
            "has_ortholog",
            "source_file",
        ]
    )
    if not all_pairs_df.empty:
        mapped_pairs_df = all_pairs_df.merge(mapping_df, left_on="focal_raw_id", right_on="raw_focal_id", how="left")
        mapped_pairs_df = mapped_pairs_df[mapped_pairs_df["mapping_status"] == "mapped"].copy()
        if not mapped_pairs_df.empty:
            positive_long_df = (
                mapped_pairs_df.assign(
                    species=species,
                    focal_gene_id=mapped_pairs_df["mapped_gene_id"],
                    has_ortholog=1,
                )[
                    [
                        "species",
                        "focal_gene_id",
                        "target_species_name",
                        "target_species_taxid",
                        "target_column_id",
                        "has_ortholog",
                        "source_file",
                    ]
                ]
                .drop_duplicates()
                .sort_values(["focal_gene_id", "target_column_id", "source_file"])
                .reset_index(drop=True)
            )

    column_ids = sorted(target_columns_seen, key=lambda value: (not str(value).isdigit(), int(value) if str(value).isdigit() else str(value)))
    matrix_df = build_binary_ortholog_matrix(positive_long_df, valid_gene_ids, column_ids)
    source_lookup_df = pd.DataFrame(columns=["focal_gene_id", "target_column_id", "source_file"])
    meta_lookup_df = pd.DataFrame(columns=["target_column_id", "target_species_name", "target_species_taxid"])
    if not positive_long_df.empty:
        source_lookup_df = (
            positive_long_df.groupby(["focal_gene_id", "target_column_id"], as_index=False)["source_file"]
            .agg(lambda values: "|".join(sorted(set(values))))
        )
        meta_lookup_df = (
            positive_long_df[["target_column_id", "target_species_name", "target_species_taxid"]]
            .drop_duplicates(subset=["target_column_id"])
            .sort_values("target_column_id")
        )
    elif column_ids:
        meta_lookup_df = pd.DataFrame(
            [{"target_column_id": column_id, "target_species_name": "", "target_species_taxid": column_id if str(column_id).isdigit() else ""} for column_id in column_ids]
        )
    long_df = (
        matrix_df.melt(id_vars=["Gene"], var_name="target_column_id", value_name="has_ortholog")
        .rename(columns={"Gene": "focal_gene_id"})
        .assign(species=species)
        .merge(meta_lookup_df, on="target_column_id", how="left")
        .merge(source_lookup_df, on=["focal_gene_id", "target_column_id"], how="left")
    )
    long_df = long_df[
        [
            "species",
            "focal_gene_id",
            "target_species_name",
            "target_species_taxid",
            "target_column_id",
            "has_ortholog",
            "source_file",
        ]
    ].sort_values(["focal_gene_id", "target_column_id"]).reset_index(drop=True)
    print(f"[INFO] species={species} target species count={len(target_species_seen)}")
    print(f"[INFO] species={species} final matrix shape={matrix_df.shape}")

    summary = {
        "species": species,
        "folder": folder,
        "labels_path": str(labels_path),
        "download_dir": str(download_dir),
        "downloaded_file_count": len(xml_paths),
        "target_species_count": len(target_species_seen),
        "raw_focal_id_count": len(raw_id_values),
        "mapped_focal_id_count": mapping_summary.get("mapped", 0),
        "unmapped_focal_id_count": mapping_summary.get("unmapped", 0),
        "ambiguous_focal_id_count": mapping_summary.get("ambiguous", 0),
        "label_gene_count": len(valid_gene_ids),
        "final_gene_count_in_matrix": len(matrix_df),
        "matrix_column_count": len(column_ids),
        "download_link_count": len(listed_links),
        "used_string_aliases": used_string_aliases,
        "notes": "",
    }
    if summary["unmapped_focal_id_count"] > summary["mapped_focal_id_count"]:
        summary["notes"] = "unmapped raw ids exceed mapped raw ids"
    return summary, matrix_df, mapping_df, long_df


def write_species_outputs(
    species: str,
    summary: Dict[str, object],
    matrix_df: pd.DataFrame,
    mapping_df: pd.DataFrame,
    long_df: pd.DataFrame,
) -> Dict[str, str]:
    outdir = OR_ROOT / species
    outdir.mkdir(parents=True, exist_ok=True)
    orthologs_csv_path = outdir / "orthologs.csv"
    mapping_path = outdir / "orthologs_id_mapping.tsv"
    long_path = outdir / "orthologs_long.tsv"
    summary_path = outdir / "orthologs_summary.tsv"

    matrix_df = matrix_df.rename(columns={"Gene": "Gene"})
    matrix_df.to_csv(orthologs_csv_path, index=False)
    pd.DataFrame([summary]).to_csv(summary_path, sep="\t", index=False)
    if mapping_path.exists():
        mapping_path.unlink()
    if long_path.exists():
        long_path.unlink()
    print(f"[DONE] wrote {orthologs_csv_path}")
    return {
        "orthologs_csv_path": str(orthologs_csv_path),
        "mapping_path": "",
        "long_path": "",
        "summary_path": str(summary_path),
    }


def build_download_summary(all_summaries: List[Dict[str, object]]) -> pd.DataFrame:
    df = pd.DataFrame(all_summaries)[
        ["species", "folder", "downloaded_file_count", "target_species_count", "download_dir"]
    ].copy()
    df.to_csv(OR_ROOT / "inparanoid_download_summary.tsv", sep="\t", index=False)
    return df


def build_global_summary(all_summaries: List[Dict[str, object]]) -> pd.DataFrame:
    summary_df = pd.DataFrame(all_summaries)[
        [
            "species",
            "downloaded_file_count",
            "target_species_count",
            "label_gene_count",
            "final_gene_count_in_matrix",
            "matrix_column_count",
            "output_csv_path",
        ]
    ].copy()
    summary_path = OR_ROOT / "orthologs_build_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    return summary_df


def write_readme(all_summaries: List[Dict[str, object]]) -> Path:
    readme_path = OR_ROOT / "README_inparanoid_orthologs.md"
    lines: List[str] = [
        "# InParanoid v8 Ortholog Matrix Build",
        "",
        "## Scope",
        "",
        "This workflow downloads pairwise InParanoid v8 `*.orthoXML` files for 5 focal species and builds EPGAT-style orthology binary matrices aligned to each species `labels.standard.tsv` final `gene_id` set.",
        "",
        "## Download Directories And File Counts",
        "",
        "| species | folder | downloaded_file_count | target_species_count | download_dir |",
        "| --- | --- | ---: | ---: | --- |",
    ]
    for item in all_summaries:
        lines.append(
            f"| {item['species']} | {item['folder']} | {item['downloaded_file_count']} | {item['target_species_count']} | `{item['download_dir']}` |"
        )
    lines.extend(
        [
            "",
            "## XML Parsing",
            "",
            "Each orthoXML file is parsed with `xml.etree.ElementTree`. The parser first reads every `<species>/<gene>` entry to map internal `geneRef id` values to `{species name, taxid, geneId, protId}`. It then walks each top-level `orthologGroup`, collects all nested `geneRef` members, and marks a focal raw identifier as ortholog-positive for a target species if they co-occur in the same ortholog group.",
            "",
            "## Focal Raw ID Mapping",
            "",
            "Mapping priority is: direct match to `labels.standard.tsv` `gene_id`, exact match against available labels auxiliary ID-like columns, exact match through STRING aliases, then a Fusarium-specific conservative transcript-to-gene conversion when aliases expose `FGRAMPH1_*T*` transcript IDs and the corresponding `FGRAMPH1_*G*` gene exists in labels.",
            "",
            "Ambiguous mappings are retained in `orthologs_id_mapping.tsv` as `mapping_status=ambiguous` and are excluded from final matrices. Unmapped raw IDs are also excluded from final matrices.",
            "",
            "## STRING Alias Usage",
            "",
            "STRING alias files under `data/processed/PPI/stringDB/` were used as auxiliary evidence for focal-side ID harmonization. The workflow does not use old EPGAT, Bingo orthology outputs, or `data_registry` files.",
            "",
            "## Matrix Definition",
            "",
            "For each focal species, `orthologs.csv` uses `Gene` as the first column, rows equal to all final label `gene_id` values, columns equal to target support species observed in downloaded orthoXML files, and 0/1 values indicating whether at least one ortholog exists in that target species. Column names prefer target taxid parsed from XML; species names are only used as fallback when taxid cannot be resolved.",
            "",
            "## High-Unmapped Species",
            "",
        ]
    )
    high_unmapped = [item for item in all_summaries if item["unmapped_focal_id_count"] > item["mapped_focal_id_count"]]
    if high_unmapped:
        for item in high_unmapped:
            lines.append(
                f"- `{item['species']}`: mapped={item['mapped_focal_id_count']}, unmapped={item['unmapped_focal_id_count']}, ambiguous={item['ambiguous_focal_id_count']}"
            )
    else:
        lines.append("- No focal species had unmapped raw IDs exceeding mapped raw IDs in the final audit.")
    lines.append("")
    readme_path.write_text("\n".join(lines) + "\n")
    return readme_path


def main() -> None:
    taxonomy_df = load_species_taxonomy(OR_ROOT)
    taxonomy_lookup = build_species_name_taxid_lookup(taxonomy_df)
    all_summaries: List[Dict[str, object]] = []

    for species, meta in FOCAL_SPECIES.items():
        if species_outputs_complete(species):
            summary_path = OR_ROOT / species / "orthologs_summary.tsv"
            summary_df = pd.read_csv(summary_path, sep="\t", dtype=str).fillna("")
            summary = summary_df.iloc[0].to_dict() if not summary_df.empty else {"species": species}
            summary["species"] = species
            summary["output_csv_path"] = str(OR_ROOT / species / "orthologs.csv")
            all_summaries.append(summary)
            print(f"[INFO] species={species} skip completed outputs at {OR_ROOT / species}")
            continue
        summary, matrix_df, mapping_df, long_df = build_species_ortholog_matrix(species, meta, taxonomy_lookup)
        outputs = write_species_outputs(species, summary, matrix_df, mapping_df, long_df)
        summary["output_csv_path"] = outputs["orthologs_csv_path"]
        all_summaries.append(summary)

    build_download_summary(all_summaries)
    build_global_summary(all_summaries)
    readme_path = write_readme(all_summaries)
    print(f"[DONE] wrote {readme_path}")


if __name__ == "__main__":
    pd.options.mode.copy_on_write = True
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
