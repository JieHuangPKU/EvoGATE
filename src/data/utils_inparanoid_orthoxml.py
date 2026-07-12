import gzip
import re
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


NS = {"ox": "http://orthoXML.org/2011/"}
ORTHOXML_SUFFIX = ".orthoXML"
MAPPING_CANDIDATE_COLUMNS = [
    "gene_id",
    "original_gene_id",
    "gene_symbol",
    "protein_id",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_protein_id",
]
TEXT_COLUMN_KEYWORDS = ("gene", "protein", "transcript", "symbol", "locus", "orf", "id")


class _HrefCollector(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: List[str] = []

    def handle_starttag(self, tag: str, attrs: List[Tuple[str, Optional[str]]]) -> None:
        if tag.lower() != "a":
            return
        attr_map = dict(attrs)
        href = attr_map.get("href")
        if href:
            self.hrefs.append(href)


def _fetch_url_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8", "ignore")


def _normalize_token(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _candidate_id_forms(raw_id: str) -> List[str]:
    value = _normalize_token(raw_id)
    if not value:
        return []
    candidates = {value}
    if "." in value:
        candidates.add(value.split(".", 1)[1])
    if "-" in value:
        candidates.add(value.split("-", 1)[0])
    if re.fullmatch(r"FGRAMPH1_\d{2}T\d{5}(?:-p\d+)?", value):
        tx = value.split("-")[0]
        candidates.add(tx)
        candidates.add(tx.replace("T", "G", 1))
    if re.fullmatch(r"FGRAMPH1_\d{2}T\d{5}", value):
        candidates.add(value.replace("T", "G", 1))
    return [item for item in candidates if item]


def _build_fgram_transcript_gene_map(valid_gene_ids: Iterable[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for gene_id in valid_gene_ids:
        gene_id = _normalize_token(gene_id)
        if not re.fullmatch(r"FGRAMPH1_\d{2}G\d{5}", gene_id):
            continue
        transcript_id = gene_id.replace("G", "T", 1)
        mapping[transcript_id] = gene_id
    return mapping


def list_orthoxml_links(folder_url: str) -> List[str]:
    html = _fetch_url_text(folder_url)
    parser = _HrefCollector()
    parser.feed(html)
    links = sorted(
        {
            href
            for href in parser.hrefs
            if href.endswith(ORTHOXML_SUFFIX) and "?" not in href and not href.startswith("/")
        }
    )
    return links


def _download_one_orthoxml(url: str, destination: Path) -> Path:
    with urllib.request.urlopen(url, timeout=120) as response:
        destination.write_bytes(response.read())
    return destination


def download_orthoxml_files(folder_url: str, outdir: Path, max_workers: int = 8) -> List[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    filenames = list_orthoxml_links(folder_url)
    downloaded_paths: List[Path] = []
    pending_jobs = []
    for filename in filenames:
        destination = outdir / filename
        if destination.exists() and destination.stat().st_size > 0:
            downloaded_paths.append(destination)
            continue
        pending_jobs.append((filename, folder_url.rstrip("/") + "/" + filename, destination))
    if pending_jobs:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_job = {}
            for filename, url, destination in pending_jobs:
                print(f"[INFO] queue download xml={filename}")
                future_to_job[executor.submit(_download_one_orthoxml, url, destination)] = filename
            for future in as_completed(future_to_job):
                filename = future_to_job[future]
                destination = future.result()
                print(f"[INFO] downloaded xml={filename}")
                downloaded_paths.append(destination)
    return sorted(downloaded_paths)


def find_labels_file(processed_root: Path, species: str) -> Path:
    candidates = sorted(processed_root.glob(f"**/{species}/labels.standard.tsv"))
    if not candidates:
        raise FileNotFoundError(f"No labels.standard.tsv found for species={species}")
    preferred = [path for path in candidates if "essential_gene" in str(path)]
    return preferred[0] if preferred else candidates[0]


def load_labels_table(labels_path: Path) -> pd.DataFrame:
    labels_df = pd.read_csv(labels_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in labels_df.columns:
        raise ValueError(f"Missing gene_id column in {labels_path}")
    labels_df["gene_id"] = labels_df["gene_id"].map(_normalize_token)
    labels_df = labels_df[labels_df["gene_id"].str.lower() != "gene_id"].copy()
    if "included_in_final" in labels_df.columns:
        included = labels_df["included_in_final"].astype(str).str.lower()
        labels_df = labels_df[included.isin({"true", "1", "yes"})].copy()
    labels_df = labels_df[labels_df["gene_id"] != ""].drop_duplicates(subset=["gene_id"]).reset_index(drop=True)
    return labels_df


def load_string_aliases(alias_path: Path) -> pd.DataFrame:
    records: List[Dict[str, str]] = []
    with gzip.open(alias_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            string_protein_id = _normalize_token(fields[0])
            alias = _normalize_token(fields[1])
            source = _normalize_token(fields[2]) if len(fields) > 2 else ""
            if not string_protein_id or not alias:
                continue
            records.append(
                {
                    "string_protein_id": string_protein_id,
                    "alias": alias,
                    "source": source,
                }
            )
    return pd.DataFrame(records)


def extract_species_from_filename(xml_path: Path) -> Tuple[str, str]:
    stem = xml_path.name[: -len(ORTHOXML_SUFFIX)] if xml_path.name.endswith(ORTHOXML_SUFFIX) else xml_path.stem
    if "-" not in stem:
        raise ValueError(f"Unexpected orthoXML filename: {xml_path.name}")
    left, right = stem.split("-", 1)
    return left, right


def _gene_attribute_candidates(gene_info: Dict[str, str]) -> List[str]:
    tokens: List[str] = []
    for key in ("geneId", "protId", "internal_id"):
        value = gene_info.get(key, "")
        tokens.extend(_candidate_id_forms(value))
    seen = set()
    ordered = []
    for token in tokens:
        if token not in seen:
            seen.add(token)
            ordered.append(token)
    return ordered


def _collect_group_gene_refs(element: ET.Element) -> List[str]:
    refs: List[str] = []
    tag = element.tag.split("}")[-1]
    if tag == "geneRef":
        ref_id = element.attrib.get("id")
        if ref_id:
            refs.append(ref_id)
        return refs
    for child in list(element):
        refs.extend(_collect_group_gene_refs(child))
    return refs


def parse_orthoxml_file(xml_path: Path) -> pd.DataFrame:
    tree = ET.parse(xml_path)
    root = tree.getroot()
    gene_index: Dict[str, Dict[str, str]] = {}
    species_rows: List[Dict[str, str]] = []

    for species_elem in root.findall("ox:species", NS):
        species_name = _normalize_token(species_elem.attrib.get("name"))
        species_taxid = _normalize_token(species_elem.attrib.get("NCBITaxId"))
        for gene_elem in species_elem.findall(".//ox:gene", NS):
            internal_id = _normalize_token(gene_elem.attrib.get("id"))
            if not internal_id:
                continue
            gene_index[internal_id] = {
                "species_name": species_name,
                "species_taxid": species_taxid,
                "internal_id": internal_id,
                "geneId": _normalize_token(gene_elem.attrib.get("geneId")),
                "protId": _normalize_token(gene_elem.attrib.get("protId")),
            }
        species_rows.append({"species_name": species_name, "species_taxid": species_taxid})

    records: List[Dict[str, str]] = []
    groups_elem = root.find("ox:groups", NS)
    if groups_elem is None:
        return pd.DataFrame(
            columns=[
                "source_file",
                "group_id",
                "focal_species_name",
                "focal_species_taxid",
                "focal_raw_id",
                "focal_prot_id",
                "target_species_name",
                "target_species_taxid",
                "target_raw_id",
                "target_prot_id",
            ]
        )

    for group_elem in groups_elem.findall("ox:orthologGroup", NS):
        ref_ids = _collect_group_gene_refs(group_elem)
        if not ref_ids:
            continue
        group_genes: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
        for ref_id in ref_ids:
            gene_info = gene_index.get(ref_id)
            if not gene_info:
                continue
            group_genes[(gene_info["species_taxid"], gene_info["species_name"])].append(gene_info)
        if len(group_genes) < 2:
            continue
        species_keys = sorted(group_genes.keys())
        for focal_key in species_keys:
            for target_key in species_keys:
                if focal_key == target_key:
                    continue
                for focal_gene in group_genes[focal_key]:
                    for target_gene in group_genes[target_key]:
                        records.append(
                            {
                                "source_file": xml_path.name,
                                "group_id": _normalize_token(group_elem.attrib.get("id")),
                                "focal_species_name": focal_gene["species_name"],
                                "focal_species_taxid": focal_gene["species_taxid"],
                                "focal_raw_id": focal_gene["geneId"] or focal_gene["protId"] or focal_gene["internal_id"],
                                "focal_prot_id": focal_gene["protId"],
                                "target_species_name": target_gene["species_name"],
                                "target_species_taxid": target_gene["species_taxid"],
                                "target_raw_id": target_gene["geneId"] or target_gene["protId"] or target_gene["internal_id"],
                                "target_prot_id": target_gene["protId"],
                            }
                        )
    return pd.DataFrame(records)


def _build_label_lookup(labels_df: pd.DataFrame) -> Dict[str, set]:
    lookup: Dict[str, set] = defaultdict(set)
    candidate_columns = [column for column in MAPPING_CANDIDATE_COLUMNS if column in labels_df.columns]
    for column in labels_df.columns:
        column_lower = column.lower()
        if column in candidate_columns:
            continue
        if labels_df[column].dtype == object and any(keyword in column_lower for keyword in TEXT_COLUMN_KEYWORDS):
            candidate_columns.append(column)
    for _, row in labels_df.iterrows():
        gene_id = _normalize_token(row["gene_id"])
        if not gene_id:
            continue
        for column in candidate_columns:
            value = _normalize_token(row.get(column, ""))
            if value:
                lookup[value].add(gene_id)
    return lookup


def _build_alias_gene_lookup(
    alias_df: Optional[pd.DataFrame],
    valid_gene_ids: Iterable[str],
    species: str,
) -> Dict[str, set]:
    if alias_df is None or alias_df.empty:
        return {}
    valid_gene_ids = {_normalize_token(gene_id) for gene_id in valid_gene_ids}
    transcript_gene_map = _build_fgram_transcript_gene_map(valid_gene_ids) if species == "fgraminearum" else {}
    protein_to_gene: Dict[str, set] = defaultdict(set)
    protein_to_aliases: Dict[str, List[Tuple[str, str]]] = defaultdict(list)

    for _, row in alias_df.iterrows():
        protein_id = _normalize_token(row["string_protein_id"])
        alias = _normalize_token(row["alias"])
        source = _normalize_token(row.get("source", ""))
        if not protein_id or not alias:
            continue
        protein_to_aliases[protein_id].append((alias, source))
        for candidate in _candidate_id_forms(alias):
            if candidate in valid_gene_ids:
                protein_to_gene[protein_id].add(candidate)
            if species == "fgraminearum" and candidate in transcript_gene_map:
                protein_to_gene[protein_id].add(transcript_gene_map[candidate])

    alias_lookup: Dict[str, set] = defaultdict(set)
    for protein_id, aliases in protein_to_aliases.items():
        gene_candidates = protein_to_gene.get(protein_id, set())
        protein_core = protein_id.split(".", 1)[1] if "." in protein_id else protein_id
        for token in _candidate_id_forms(protein_id) + _candidate_id_forms(protein_core):
            if token in valid_gene_ids:
                gene_candidates.add(token)
            if species == "fgraminearum" and token in transcript_gene_map:
                gene_candidates.add(transcript_gene_map[token])
        if not gene_candidates:
            continue
        observed_tokens = set()
        for alias, _source in aliases:
            observed_tokens.update(_candidate_id_forms(alias))
        observed_tokens.update(_candidate_id_forms(protein_id))
        observed_tokens.update(_candidate_id_forms(protein_core))
        for token in observed_tokens:
            alias_lookup[token].update(gene_candidates)
    return alias_lookup


def map_focal_ids_to_gene_ids(
    raw_ids: pd.Series,
    labels_df: pd.DataFrame,
    alias_df: Optional[pd.DataFrame],
    species: str,
) -> pd.DataFrame:
    unique_raw_ids = pd.Series(raw_ids, dtype="string").fillna("").map(_normalize_token)
    unique_raw_ids = unique_raw_ids[unique_raw_ids != ""].drop_duplicates().sort_values()
    label_lookup = _build_label_lookup(labels_df)
    valid_gene_ids = set(labels_df["gene_id"].map(_normalize_token))
    alias_lookup = _build_alias_gene_lookup(alias_df, valid_gene_ids, species)
    transcript_gene_map = _build_fgram_transcript_gene_map(valid_gene_ids) if species == "fgraminearum" else {}

    records: List[Dict[str, str]] = []
    for raw_id in unique_raw_ids:
        candidate_tokens = _candidate_id_forms(raw_id)
        direct_candidates = set()
        label_candidates = set()
        alias_candidates = set()
        transcript_candidates = set()

        for token in candidate_tokens:
            if token in valid_gene_ids:
                direct_candidates.add(token)
            label_candidates.update(label_lookup.get(token, set()))
            alias_candidates.update(alias_lookup.get(token, set()))
            if species == "fgraminearum" and token in transcript_gene_map:
                transcript_candidates.add(transcript_gene_map[token])

        candidates = set()
        method = ""
        evidence = ""
        notes = ""
        if direct_candidates:
            candidates = direct_candidates
            method = "direct_gene_id"
            evidence = "; ".join(sorted(candidate_tokens))
        elif label_candidates:
            candidates = label_candidates
            method = "labels_auxiliary_column"
            evidence = "; ".join(sorted(candidate_tokens))
        elif alias_candidates:
            candidates = alias_candidates
            method = "string_aliases"
            evidence = "; ".join(sorted(candidate_tokens))
        elif transcript_candidates:
            candidates = transcript_candidates
            method = "fgram_alias_transcript_to_gene"
            evidence = "; ".join(sorted(candidate_tokens))

        mapped_gene_id = ""
        status = "unmapped"
        if len(candidates) == 1:
            status = "mapped"
            mapped_gene_id = sorted(candidates)[0]
        elif len(candidates) > 1:
            status = "ambiguous"
            notes = "multiple candidate gene_ids: " + "; ".join(sorted(candidates))

        records.append(
            {
                "raw_focal_id": raw_id,
                "mapped_gene_id": mapped_gene_id,
                "mapping_status": status,
                "mapping_method": method,
                "evidence": evidence,
                "notes": notes,
            }
        )
    return pd.DataFrame(records)


def build_binary_ortholog_matrix(
    long_df: pd.DataFrame,
    valid_gene_ids: List[str],
    column_ids: List[str],
) -> pd.DataFrame:
    matrix_df = pd.DataFrame(0, index=pd.Index(valid_gene_ids, name="Gene"), columns=column_ids, dtype=int)
    if not long_df.empty:
        positive_df = long_df[long_df["has_ortholog"].astype(int) == 1].copy()
        positive_df = positive_df[
            positive_df["focal_gene_id"].isin(valid_gene_ids) & positive_df["target_column_id"].isin(column_ids)
        ]
        if not positive_df.empty:
            for gene_id, column_id in positive_df[["focal_gene_id", "target_column_id"]].drop_duplicates().itertuples(index=False):
                matrix_df.at[gene_id, column_id] = 1
    matrix_df = matrix_df.reset_index()
    return matrix_df
