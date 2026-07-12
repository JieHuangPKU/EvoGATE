import pandas as pd
from pathlib import Path
import gzip
import csv
from typing import Dict, Set, List, Optional, Tuple

def find_labels_file(processed_root: Path, species: str) -> Path:
    """
    Recursively search for labels.standard.tsv for a given species.
    """
    candidates = list(processed_root.glob(f"**/{species}/labels.standard.tsv"))
    if not candidates:
        all_labels = list(processed_root.glob("**/labels.standard.tsv"))
        candidates = [c for c in all_labels if f"/{species}/" in str(c)]
        
    if not candidates:
        raise FileNotFoundError(f"No labels.standard.tsv found for species: {species}")
    
    if len(candidates) > 1:
        essential_candidates = [c for c in candidates if "essential_gene" in str(c)]
        if essential_candidates:
            return essential_candidates[0]
    
    return candidates[0]

def load_labels_table(labels_path: Path) -> pd.DataFrame:
    return pd.read_csv(labels_path, sep='\t')

def build_valid_gene_id_index(labels_df: pd.DataFrame) -> Dict[str, str]:
    """
    Builds an index of all valid IDs pointing to the main gene_id.
    """
    index = {}
    main_ids = labels_df['gene_id'].astype(str).tolist()
    for mid in main_ids:
        index[mid] = mid
        
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
    if not alias_path.exists():
        return pd.DataFrame(columns=['raw_string_protein_id', 'alias', 'source'])
        
    with gzip.open(alias_path, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        try:
            header = next(reader)
        except StopIteration:
            return pd.DataFrame(columns=['raw_string_protein_id', 'alias', 'source'])
            
        for row in reader:
            if len(row) >= 2:
                data.append({
                    'raw_string_protein_id': row[0],
                    'alias': row[1],
                    'source': row[2] if len(row) > 2 else "unknown"
                })
    return pd.DataFrame(data)

def load_compartments_table(input_path: Path, species: str) -> pd.DataFrame:
    """
    Loads raw compartments file and returns a long-format DataFrame with raw_id and location.
    """
    if species == 'fgraminearum':
        data = []
        with open(input_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split('\t')
                raw_id = parts[0]
                locs = [p.strip() for p in parts[1:] if p.strip()]
                for l in locs:
                    data.append({'raw_id': raw_id, 'location': l, 'score': 5.0, 'evidence': 'eFG'})
        return pd.DataFrame(data)
    
    elif species == 'scerevisiae':
        df = pd.read_csv(input_path, sep='\t', header=None, 
                         names=['raw_id', 'name', 'go_id', 'location', 'score'])
        df['evidence'] = 'integrated'
        return df
    
    else:
        df = pd.read_csv(input_path, sep='\t', header=None,
                         names=['raw_id', 'name', 'go_id', 'location', 'source', 'evidence', 'score'])
        return df

def build_fusarium_transcript_to_gene_map(ph1_tsv_path: Path) -> Dict[str, str]:
    if not ph1_tsv_path.exists():
        return {}
    df = pd.read_csv(ph1_tsv_path, sep='\t')
    
    tx_col, gene_col = None, None
    if 'transcript' in df.columns: tx_col = 'transcript'
    if 'gene' in df.columns: gene_col = 'gene'
    
    if not tx_col or not gene_col:
        for col in df.columns:
            col_lower = col.lower()
            if not tx_col and any(x == col_lower for x in ['transcript', 'mrna', 'tx']):
                tx_col = col
            if not gene_col and any(x == col_lower for x in ['gene', 'gene_id', 'canonical_gene_id']):
                gene_col = col
                
    if not tx_col or not gene_col:
        return {}
        
    mapping = {}
    for _, row in df.iterrows():
        tx = str(row[tx_col]).strip()
        gene = str(row[gene_col]).strip()
        if tx and gene and tx != 'nan' and gene != 'nan':
            mapping[tx] = gene
    return mapping

def normalize_location_term(raw_term: str) -> Optional[str]:
    """
    Controlled vocabulary mapping for subcellular locations.
    """
    term = str(raw_term).lower().strip()
    
    # Generic/Noise blacklist
    blacklist = [
        'cellular_component', 'intracellular', 'anatomical entity', 
        'anatomical structure', 'organelle', 'membrane-bounded',
        'membrane-enclosed', 'lumen', 'intrinsic component',
        'integral component', 'protein-containing complex', 'go:'
    ]
    if any(b in term for b in blacklist):
        return None
    
    # Mapping Rules (Priority matters)
    if 'plasma membrane' in term or 'cell membrane' in term or 'cell periphery' in term:
        return 'Cell membrane'
    if 'multi-pass membrane' in term or 'multipass membrane' in term:
        return 'Multi-pass membrane'
    if 'nucleus' in term or 'nuclear' in term or 'nucleoplasm' in term or 'nucleolus' in term or 'chromatin' in term or 'chromosome' in term:
        return 'Nucleus'
    if 'cytoplasm' in term or 'cytosol' in term:
        return 'Cytoplasm'
    if 'mitochondrion' in term or 'mitochondria' in term or 'mitochondrial' in term:
        return 'Mitochondrion'
    if 'endoplasmic reticulum' in term or ' er' in term or term == 'er' or 'sarcoplasmic reticulum' in term:
        return 'Endoplasmic reticulum'
    if 'golgi' in term:
        return 'Golgi'
    if 'endosome' in term or 'endosomal' in term:
        return 'Endosome'
    if 'vacuole' in term or 'vacuolar' in term:
        return 'Vacuole'
    if 'peroxisome' in term or 'peroxisomal' in term:
        return 'Peroxisome'
    if 'secreted' in term or 'extracellular' in term:
        return 'Secreted'
    if 'lysosome' in term or 'lysosomal' in term:
        return 'Lysosome'
    if 'cytoskeleton' in term or 'actin' in term or 'microtubule' in term or 'spindle' in term or 'kinetochore' in term:
        return 'Cytoskeleton'
    if 'ribosome' in term or 'ribosomal' in term:
        return 'Ribosome'
    if 'extracellular matrix' in term:
        return 'Extracellular matrix'
    
    return None

def map_subloc_ids_to_gene_ids(
    subloc_df: pd.DataFrame,
    valid_id_index: Dict[str, str],
    alias_df: Optional[pd.DataFrame],
    species: str,
    fusarium_tx2gene: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    print(f"[INFO] Mapping {species} IDs...")
    alias_to_genes = {}
    
    if alias_df is not None and not alias_df.empty:
        prot_to_genes = {}
        for prot_id, group in alias_df.groupby('raw_string_protein_id'):
            genes = set()
            core_id = prot_id.split('.', 1)[1] if '.' in prot_id else prot_id
            if core_id in valid_id_index:
                genes.add(valid_id_index[core_id])
            aliases = group['alias'].astype(str).tolist()
            for a in aliases:
                if a in valid_id_index:
                    genes.add(valid_id_index[a])
                if species == 'fgraminearum' and fusarium_tx2gene and a in fusarium_tx2gene:
                    g = fusarium_tx2gene[a]
                    if g in valid_id_index:
                        genes.add(valid_id_index[g])
            if genes:
                prot_to_genes[prot_id] = genes
        
        alias_to_prots = alias_df.groupby('alias')['raw_string_protein_id'].apply(set).to_dict()
        for alias, prots in alias_to_prots.items():
            all_genes = set()
            for p in prots:
                if p in prot_to_genes:
                    all_genes.update(prot_to_genes[p])
            if all_genes:
                alias_to_genes[alias] = all_genes

    unique_raw_ids = subloc_df['raw_id'].unique()
    id_to_final = {}
    for rid in unique_raw_ids:
        rid_str = str(rid).strip()
        if rid_str in valid_id_index:
            id_to_final[rid] = (valid_id_index[rid_str], 'mapped', 'direct_or_label_column', 'labels.standard.tsv')
            continue
        if rid_str in alias_to_genes:
            genes = alias_to_genes[rid_str]
            if len(genes) == 1:
                id_to_final[rid] = (list(genes)[0], 'mapped', 'alias_bridge', 'STRING aliases bridge')
            else:
                id_to_final[rid] = ("", 'ambiguous', 'alias_bridge', f'multiple matches: {";".join(list(genes))}')
            continue
        id_to_final[rid] = ("", 'unmapped', 'none', 'no evidence found')

    mapped_rows = []
    for _, row in subloc_df.iterrows():
        mapped_gene_id, status, method, evidence = id_to_final[row['raw_id']]
        mapped_rows.append({
            **row,
            'mapped_gene_id': mapped_gene_id,
            'mapping_status': status,
            'mapping_method': method,
            'mapping_evidence': evidence
        })
    return pd.DataFrame(mapped_rows)

def build_gene_location_long_table(mapped_df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes location terms and filters for mapped genes.
    """
    valid_df = mapped_df[mapped_df['mapping_status'] == 'mapped'].copy()
    if valid_df.empty:
        return pd.DataFrame(columns=['gene_id', 'raw_location', 'normalized_location', 'source_db', 'evidence'])
    
    valid_df['normalized_location'] = valid_df['location'].apply(normalize_location_term)
    
    long_df = valid_df.dropna(subset=['normalized_location']).copy()
    long_df = long_df.rename(columns={'mapped_gene_id': 'gene_id', 'location': 'raw_location'})
    
    # Deduplicate gene-normalized_location pairs
    cols = ['gene_id', 'raw_location', 'normalized_location', 'evidence']
    if 'source_db' in valid_df.columns:
        cols.append('source_db')
    
    return long_df[cols].drop_duplicates()

def summarize_location_frequency(gene_loc_long_df: pd.DataFrame) -> pd.DataFrame:
    if gene_loc_long_df.empty:
        return pd.DataFrame(columns=['location', 'gene_count', 'raw_term_count', 'example_raw_terms'])
    
    # Count genes per normalized location
    freq = gene_loc_long_df.groupby('normalized_location')['gene_id'].nunique().reset_index()
    freq.columns = ['location', 'gene_count']
    
    # Count raw terms per normalized location
    raw_counts = gene_loc_long_df.groupby('normalized_location')['raw_location'].nunique().reset_index()
    raw_counts.columns = ['location', 'raw_term_count']
    
    # Example raw terms
    examples = gene_loc_long_df.groupby('normalized_location')['raw_location'].apply(lambda x: ";".join(list(set(x))[:5])).reset_index()
    examples.columns = ['location', 'example_raw_terms']
    
    summary = freq.merge(raw_counts, on='location').merge(examples, on='location')
    return summary.sort_values('gene_count', ascending=False)

def select_final_location_columns(freq_df: pd.DataFrame, min_cols: int = 12, max_cols: int = 15) -> List[str]:
    if freq_df.empty:
        return []
    
    # Already sorted by gene_count in summarize_location_frequency
    all_locs = freq_df['location'].tolist()
    return all_locs[:max_cols]

def build_binary_location_matrix(
    gene_loc_long_df: pd.DataFrame,
    valid_gene_ids: List[str],
    final_location_cols: List[str],
) -> pd.DataFrame:
    if not final_location_cols:
        return pd.DataFrame({'Gene': valid_gene_ids})
    
    # Pivot
    matrix = gene_loc_long_df[gene_loc_long_df['normalized_location'].isin(final_location_cols)].copy()
    matrix['val'] = 1
    pivot = matrix.pivot_table(index='gene_id', columns='normalized_location', values='val', fill_value=0)
    
    # Ensure all valid_gene_ids are present
    final_df = pd.DataFrame({'Gene': valid_gene_ids})
    pivot = pivot.reset_index().rename(columns={'gene_id': 'Gene'})
    
    final_df = final_df.merge(pivot, on='Gene', how='left').fillna(0)
    
    # Reorder columns: Gene, then final_location_cols (which are frequency-sorted)
    cols = ['Gene'] + [c for c in final_location_cols if c in final_df.columns]
    return final_df[cols].astype({c: int for c in cols if c != 'Gene'})
