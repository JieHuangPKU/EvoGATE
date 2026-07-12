import pandas as pd
from pathlib import Path
import gzip
import csv
import io
import re
from typing import Dict, Set, List, Optional, Tuple, Union

def find_labels_file(processed_root: Path, species: str) -> Path:
    target = processed_root / 'essential_gene' / species / 'labels.standard.tsv'
    if target.exists(): return target
    candidates = list(processed_root.glob(f"**/{species}/labels.standard.tsv"))
    if not candidates:
        all_labels = list(processed_root.glob("**/labels.standard.tsv"))
        candidates = [c for c in all_labels if f"/{species}/" in str(c)]
    if not candidates: raise FileNotFoundError(f"No labels.standard.tsv found for species: {species}")
    return candidates[0]

def load_labels_table(labels_path: Path) -> pd.DataFrame:
    return pd.read_csv(labels_path, sep='\t')

def build_valid_gene_id_index(labels_df: pd.DataFrame) -> Dict[str, str]:
    index = {}
    main_ids = labels_df['gene_id'].astype(str).tolist()
    for mid in main_ids: index[mid] = mid
    potential_cols = ['original_gene_id', 'gene_symbol', 'protein_id', 
                      'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_protein_id']
    for col in potential_cols:
        if col in labels_df.columns:
            for _, row in labels_df.iterrows():
                val = str(row[col]).strip()
                if val and val != 'nan' and val not in index:
                    index[val] = str(row['gene_id'])
    return index

def load_string_aliases(alias_path: Path) -> pd.DataFrame:
    data = []
    if not alias_path or not alias_path.exists():
        return pd.DataFrame(columns=['raw_string_protein_id', 'alias', 'source'])
    with gzip.open(alias_path, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        try: next(reader) 
        except StopIteration: return pd.DataFrame(columns=['raw_string_protein_id', 'alias', 'source'])
        for row in reader:
            if len(row) >= 2:
                data.append({'raw_string_protein_id': row[0], 'alias': row[1], 'source': row[2] if len(row) > 2 else "unknown"})
    return pd.DataFrame(data)

def extract_platform_annotation_from_family_soft(family_soft_path: Path) -> Dict[str, pd.DataFrame]:
    platforms = {}
    current_gpl = None
    in_table = False
    table_lines = []
    with gzip.open(family_soft_path, 'rt', errors='ignore') as f:
        for line in f:
            if line.startswith('^PLATFORM'):
                current_gpl = line.strip().split('=')[1].strip()
            elif line.startswith('!platform_table_begin'):
                in_table = True; table_lines = []
            elif line.startswith('!platform_table_end'):
                in_table = False
                if current_gpl and table_lines:
                    df = pd.read_csv(io.StringIO("".join(table_lines)), sep='\t', low_memory=False)
                    platforms[current_gpl] = df
            elif in_table:
                table_lines.append(line)
    return platforms

def load_geo_series_matrix(matrix_path: Path) -> Tuple[pd.DataFrame, List[str]]:
    header_line_idx = -1; platform_gpls = []
    with gzip.open(matrix_path, 'rt') as f:
        for i, line in enumerate(f):
            if line.startswith('!Series_platform_id'): platform_gpls.append(line.strip().split('\t')[1].strip('"'))
            if line.startswith('"ID_REF"'): header_line_idx = i; break
    if header_line_idx == -1: raise ValueError(f"Could not find ID_REF in {matrix_path}")
    df = pd.read_csv(matrix_path, sep='\t', compression='gzip', skiprows=header_line_idx)
    df = df[~df.iloc[:, 0].astype(str).str.startswith('!')]
    df.columns = [c.strip('"') for c in df.columns]
    if not df.empty: df.iloc[:, 0] = df.iloc[:, 0].astype(str).str.strip('"')
    return df, platform_gpls

def load_wormbase_id_mapping(path: Path) -> Dict[str, str]:
    mapping = {}
    if not path.exists(): return mapping
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                wbgene = parts[2]; public = parts[3]; mapping[public] = wbgene
                if len(parts) >= 5: mapping[parts[4]] = wbgene
    return mapping

def load_fusarium_tx2gene(path: Path) -> Dict[str, str]:
    mapping = {}
    if not path.exists(): return mapping
    df = pd.read_csv(path, sep='\t')
    for _, row in df.iterrows():
        tx = str(row['transcript']).strip(); gene = str(row['gene']).strip()
        if tx and gene: mapping[tx] = gene
    return mapping

def build_probe_to_gene_mapping(
    annotation_df: pd.DataFrame,
    valid_id_index: Dict[str, str],
    alias_df: Optional[pd.DataFrame],
    species: str
) -> pd.DataFrame:
    print(f"[INFO] Building probe-to-gene mapping for {species}...")
    
    id_candidates = ['ID', 'SPOT_ID', 'ID_REF']
    gene_candidates = ['ORF', 'Gene Symbol', 'GENE_SYMBOL', 'SYMBOL', 'ENTREZ_GENE_ID', 'Ensembl', 'WormBase', 'FlyBase', 'SGD accession number', 'description', 'representative public id', 'GB_ACC']
    
    col_map = {c.lower(): c for c in annotation_df.columns}
    actual_id_col = next((col_map[c.lower()] for c in id_candidates if c.lower() in col_map), annotation_df.columns[0])
    actual_gene_cols = [col_map[c.lower()] for c in gene_candidates if c.lower() in col_map]
    
    alias_to_genes = {}
    if alias_df is not None and not alias_df.empty:
        prot_to_genes = {}
        for prot_id, group in alias_df.groupby('raw_string_protein_id'):
            genes = set()
            core_id = prot_id.split('.', 1)[1] if '.' in prot_id else prot_id
            if core_id in valid_id_index: genes.add(valid_id_index[core_id])
            for a in group['alias'].astype(str):
                if a in valid_id_index: genes.add(valid_id_index[a])
            if genes: prot_to_genes[prot_id] = genes
        
        alias_to_prots = alias_df.groupby('alias')['raw_string_protein_id'].apply(set).to_dict()
        for alias, prots in alias_to_prots.items():
            all_genes = set()
            for p in prots:
                if p in prot_to_genes: all_genes.update(prot_to_genes[p])
            if all_genes: alias_to_genes[alias] = all_genes

    wb_mapping = {}
    if species == 'celegans':
        wb_mapping = load_wormbase_id_mapping(Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/PPI/wormbase.WS240.gene_ids.txt.gz'))

    fg_tx2gene = {}
    if species == 'fgraminearum':
        fg_tx2gene = load_fusarium_tx2gene(Path('/Users/jiehuang/work/2025禾谷镰刀菌/RNA农药联合攻关/小分子/FgraminearumPH-1.tsv'))

    results = []
    for _, row in annotation_df.iterrows():
        raw_probe_id = str(row[actual_id_col]).strip()
        ann_candidates = []
        for c in actual_gene_cols:
            val = str(row[c]).strip()
            if val and val != 'nan' and len(val) > 1:
                if '///' in val: ann_candidates.extend([v.strip() for v in val.split('///')])
                else: ann_candidates.append(val)
        ann_candidates.append(raw_probe_id)

        found_gene = None; method = "none"; status = "unmapped"
        for cand in ann_candidates:
            c = str(cand).strip().strip('"')
            c_base = c.split('_')[0] if '_' in c else c
            for test_c in [c, c_base]:
                if species == 'fgraminearum' and test_c.startswith('fg') and test_c[2:].isdigit():
                    test_c = f"FGRAMPH1_01G{test_c[2:].zfill(5)}"
                if test_c in valid_id_index:
                    found_gene = valid_id_index[test_c]; method = "mapped_direct"; status = "mapped"; break
                if species == 'fgraminearum' and test_c in fg_tx2gene:
                    g = fg_tx2gene[test_c]
                    if g in valid_id_index: found_gene = valid_id_index[g]; method = "fg_tx2gene"; status = "mapped"; break
                if species == 'celegans' and test_c in wb_mapping:
                    g = wb_mapping[test_c]
                    if g in valid_id_index: found_gene = valid_id_index[g]; method = "wb_mapping"; status = "mapped"; break
                if test_c in alias_to_genes:
                    genes = alias_to_genes[test_c]
                    if len(genes) == 1: found_gene = list(genes)[0]; method = "alias_bridge"; status = "mapped"; break
                    else: status = "ambiguous"
            if found_gene: break
        results.append({'raw_probe_id': raw_probe_id, 'mapped_gene_id': found_gene if found_gene else "", 'mapping_status': status, 'mapping_method': method})
    return pd.DataFrame(results)

def detect_expression_columns(df: pd.DataFrame, id_col: str) -> List[str]:
    expr_cols = []
    metadata_keywords = ['id', 'name', 'symbol', 'description', 'location', 'score', 'product', 'taxid', 'organism', 'source', 'title']
    for i, c in enumerate(df.columns):
        if str(c) == str(id_col): continue
        c_lower = str(c).lower()
        if any(k == c_lower or f"{k}_" in c_lower or f" {k}" in c_lower for k in metadata_keywords):
            if i < 5: continue
        if pd.api.types.is_numeric_dtype(df[c]): expr_cols.append(c)
        elif df[c].dtype == object:
            sample = df[c].dropna().head(10)
            if sample.empty: continue
            try: pd.to_numeric(sample); expr_cols.append(c)
            except: continue
    return expr_cols

def map_expression_ids_to_gene_ids(exp_df: pd.DataFrame, probe_map_df: pd.DataFrame, species: str) -> pd.DataFrame:
    id_col = exp_df.columns[0]
    return exp_df.merge(probe_map_df, left_on=id_col, right_on='raw_probe_id', how='left')

def aggregate_gene_level_expression(mapped_df: pd.DataFrame, expr_cols: List[str]) -> pd.DataFrame:
    valid_df = mapped_df[mapped_df['mapping_status'] == 'mapped'].copy()
    if valid_df.empty: return pd.DataFrame(columns=['Gene'] + expr_cols)
    for c in expr_cols: valid_df[c] = pd.to_numeric(valid_df[c], errors='coerce')
    grouped = valid_df.groupby('mapped_gene_id')[expr_cols].mean().reset_index()
    grouped = grouped.rename(columns={'mapped_gene_id': 'Gene'})
    return grouped

def calculate_tpm(counts_df, lengths):
    """
    counts_df: rows=genes, cols=samples
    lengths: series index=gene_id, values=bp
    """
    common_genes = counts_df.index.intersection(lengths.index)
    counts_df = counts_df.loc[common_genes]
    lengths = lengths.loc[common_genes]
    # RPK = count / (length / 1000)
    rpk = counts_df.divide(lengths / 1000, axis=0)
    # TPM = RPK / (sum(RPK) / 1,000,000)
    tpm = rpk.divide(rpk.sum(axis=0) / 1e6, axis=1)
    return tpm
