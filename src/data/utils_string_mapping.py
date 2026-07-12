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
        # Try without species name if it's already in the path
        candidates = list(processed_root.glob("**/labels.standard.tsv"))
        candidates = [c for c in candidates if f"/{species}/" in str(c)]
        
    if not candidates:
        raise FileNotFoundError(f"No labels.standard.tsv found for species: {species}")
    
    if len(candidates) > 1:
        # Pick the one under essential_gene if multiple exist
        best = [c for c in candidates if "essential_gene" in str(c)]
        if best:
            return best[0]
        return candidates[0]
    
    return candidates[0]

def load_valid_gene_ids(labels_path: Path) -> set:
    """
    Load gene_id column from labels.standard.tsv.
    """
    df = pd.read_csv(labels_path, sep='\t')
    if 'gene_id' not in df.columns:
        raise ValueError(f"Column 'gene_id' not found in {labels_path}")
    return set(df['gene_id'].astype(str))

def load_string_aliases(alias_path: Path) -> pd.DataFrame:
    """
    Load STRING aliases file.
    Columns: string_protein_id, alias, source
    """
    data = []
    with gzip.open(alias_path, 'rt') as f:
        # STRING aliases usually has a header
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        # Check if header is actually data or not
        # Usually it is: #string_protein_id	alias	source
        for row in reader:
            if len(row) >= 2:
                data.append({
                    'raw_string_protein_id': row[0],
                    'alias': row[1],
                    'source': row[2] if len(row) > 2 else "unknown"
                })
    return pd.DataFrame(data)

def build_fusarium_transcript_to_gene_map(ph1_tsv_path: Path) -> Dict[str, str]:
    """
    Build FGRAMPH1_01Txxxxx -> FGRAMPH1_01Gxxxxx map.
    """
    df = pd.read_csv(ph1_tsv_path, sep='\t')
    # ID	transcript	gene	...
    # FGRAMPH1_01T00001-p1	FGRAMPH1_01T00001	FGRAMPH1_01G00001
    
    tx_col = None
    gene_col = None
    
    # Priority 1: Exact matches
    if 'transcript' in df.columns: tx_col = 'transcript'
    if 'gene' in df.columns: gene_col = 'gene'
    
    # Priority 2: Fuzzy matches
    if not tx_col or not gene_col:
        for col in df.columns:
            col_lower = col.lower()
            if not tx_col and any(x == col_lower for x in ['transcript', 'mrna', 'tx']):
                tx_col = col
            if not gene_col and any(x == col_lower for x in ['gene', 'gene_id', 'canonical_gene_id']):
                gene_col = col
                
    if not tx_col or not gene_col:
        # Last resort fuzzy
        for col in df.columns:
            col_lower = col.lower()
            if not tx_col and any(x in col_lower for x in ['transcript', 'mrna', 'tx']):
                tx_col = col
            if not gene_col and any(x in col_lower for x in ['gene', 'gene_id']):
                gene_col = col
        
    if not tx_col or not gene_col:
        raise ValueError(f"Could not identify transcript and gene columns in {ph1_tsv_path}")
        
    mapping = {}
    for _, row in df.iterrows():
        tx = str(row[tx_col]).strip()
        gene = str(row[gene_col]).strip()
        if tx and gene and tx != 'nan' and gene != 'nan':
            mapping[tx] = gene
            
    return mapping

def map_string_proteins_to_gene_ids(
    alias_df: pd.DataFrame,
    valid_gene_ids: Set[str],
    species: str,
    fusarium_tx2gene: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Map each raw_string_protein_id to a single gene_id based on aliases.
    """
    print(f"[INFO] Mapping STRING proteins for {species} using aliases...")
    
    # Group aliases by protein
    prot_to_aliases = alias_df.groupby('raw_string_protein_id')
    
    results = []
    
    for prot_id, group in prot_to_aliases:
        candidates = set()
        evidence = {} # gene_id -> list of (alias, source)
        
        # 1. Collect all aliases
        aliases = group['alias'].astype(str).tolist()
        sources = group['source'].astype(str).tolist()
        
        for a, s in zip(aliases, sources):
            # Try direct match
            if a in valid_gene_ids:
                candidates.add(a)
                evidence.setdefault(a, []).append((a, s))
            
            # Try removing taxon prefix from alias if it has one (rare for aliases but just in case)
            if '.' in a:
                a_core = a.split('.', 1)[1]
                if a_core in valid_gene_ids:
                    candidates.add(a_core)
                    evidence.setdefault(a_core, []).append((a_core, s))
            
            # 2. Fusarium specific tx -> gene
            if species == 'fgraminearum' and fusarium_tx2gene:
                if a in fusarium_tx2gene:
                    g = fusarium_tx2gene[a]
                    if g in valid_gene_ids:
                        candidates.add(g)
                        evidence.setdefault(g, []).append((a, s + " (via ph1_tsv)"))
        
        # Also try protein ID core itself
        prot_core = prot_id.split('.', 1)[1] if '.' in prot_id else prot_id
        if prot_core in valid_gene_ids:
            candidates.add(prot_core)
            evidence.setdefault(prot_core, []).append((prot_core, "direct_protein_id_match"))

        # Decide status
        status = "unmapped"
        mapped_gene = ""
        method = ""
        matched_aliases = []
        
        if len(candidates) == 1:
            status = "mapped"
            mapped_gene = list(candidates)[0]
            matched_aliases = [f"{a}({s})" for a, s in evidence[mapped_gene]]
            method = "alias_match"
        elif len(candidates) > 1:
            status = "ambiguous"
            method = "multiple_aliases"
            # Try to see if they all point to same gene (already handled by set)
            # Actually if they are different, it's ambiguous.
            matched_aliases = [f"{g}:{[f'{a}({s})' for a, s in ev]}" for g, ev in evidence.items()]
            
        results.append({
            'raw_string_protein_id': prot_id,
            'alias_count': len(aliases),
            'matched_aliases': "; ".join(matched_aliases),
            'mapped_gene_id': mapped_gene,
            'mapping_status': status,
            'mapping_method': method,
            'evidence': "; ".join([s for a, s in zip(aliases, sources)][:5]) # just a sample
        })
        
    return pd.DataFrame(results)

def summarize_mapping(mapping_df: pd.DataFrame) -> Dict:
    counts = mapping_df['mapping_status'].value_counts().to_dict()
    return {
        'mapped': counts.get('mapped', 0),
        'unmapped': counts.get('unmapped', 0),
        'ambiguous': counts.get('ambiguous', 0),
        'total': len(mapping_df)
    }
