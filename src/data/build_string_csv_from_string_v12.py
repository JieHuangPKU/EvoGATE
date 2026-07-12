import os
import sys
import pandas as pd
from pathlib import Path
import gzip
import csv
from typing import Dict, Set, List, Optional
from utils_string_mapping import (
    find_labels_file,
    load_valid_gene_ids,
    load_string_aliases,
    build_fusarium_transcript_to_gene_map,
    map_string_proteins_to_gene_ids,
    summarize_mapping
)

# Constants
SPECIES_MAP = {
    "229533": "fgraminearum",
    "4932": "scerevisiae",
    "6239": "celegans",
    "7227": "melanogaster",
    "9606": "human",
}

PROGATE_ROOT = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2')
PROCESSED_ROOT = PROGATE_ROOT / 'data' / 'processed'
PPI_ROOT = PROCESSED_ROOT / 'PPI'
STRING_DB_ROOT = PPI_ROOT / 'stringDB'
FGRAM_PH1_TSV = Path('/Users/jiehuang/work/2025禾谷镰刀菌/RNA农药联合攻关/小分子/FgraminearumPH-1.tsv')

def load_string_links(links_path: Path) -> pd.DataFrame:
    """
    Load STRING links detailed file.
    """
    print(f"[INFO] Loading links from {links_path}...")
    # Use chunking if needed, but here we assume it fits in memory or we process it carefully.
    # STRING detailed links are space-separated.
    data = []
    with gzip.open(links_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter=' ')
        for row in reader:
            data.append({
                'protein1': row['protein1'],
                'protein2': row['protein2'],
                'combined_score': int(row['combined_score']),
                'coexpression': int(row['coexpression']),
                'experimental': int(row['experimental']),
                'database': int(row['database']),
            })
    return pd.DataFrame(data)

def process_species(taxon_id: str, species_name: str, fusarium_tx2gene: Optional[Dict[str, str]] = None):
    print(f"\n[INFO] >>> Processing species: {species_name} (Taxon: {taxon_id})")
    
    # 1. Find and load labels
    try:
        labels_path = find_labels_file(PROCESSED_ROOT, species_name)
        print(f"[INFO] Using labels file: {labels_path}")
        valid_gene_ids = load_valid_gene_ids(labels_path)
        print(f"[INFO] Loaded {len(valid_gene_ids)} valid gene IDs.")
    except Exception as e:
        print(f"[ERROR] Failed to load labels for {species_name}: {e}")
        return None

    # 2. Load and map aliases
    alias_path = STRING_DB_ROOT / f"{taxon_id}.protein.aliases.v12.0.txt.gz"
    if not alias_path.exists():
        print(f"[ERROR] Alias file not found: {alias_path}")
        return None
    
    alias_df = load_string_aliases(alias_path)
    print(f"[INFO] Loaded {len(alias_df)} alias entries.")
    
    mapping_df = map_string_proteins_to_gene_ids(
        alias_df, valid_gene_ids, species_name, fusarium_tx2gene
    )
    
    mapping_summary = summarize_mapping(mapping_df)
    print(f"[INFO] Mapping results: {mapping_summary}")
    
    # Create fast mapping dict
    prot_to_gene = mapping_df[mapping_df['mapping_status'] == 'mapped'].set_index('raw_string_protein_id')['mapped_gene_id'].to_dict()

    # 3. Load links and process edges
    links_path = STRING_DB_ROOT / f"{taxon_id}.protein.links.detailed.v12.0.txt.gz"
    if not links_path.exists():
        print(f"[ERROR] Links file not found: {links_path}")
        return None
    
    # We will process links line by line to save memory
    out_dir = PPI_ROOT / species_name
    out_dir.mkdir(parents=True, exist_ok=True)
    
    final_edges = {} # (A, B) -> scores
    edge_audit = []
    
    stats = {
        'species': species_name,
        'taxon_id': taxon_id,
        'labels_path': str(labels_path),
        'alias_path': str(alias_path),
        'links_path': str(links_path),
        'raw_edge_count': 0,
        'kept_edge_count': 0,
        'dropped_edge_count': 0,
        'self_loop_removed_count': 0,
        'duplicate_collapsed_count': 0,
        'mapped_label_gene_count': 0,
        'label_gene_count': len(valid_gene_ids)
    }

    print(f"[INFO] Processing links from {links_path}...")
    with gzip.open(links_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter=' ')
        for row in reader:
            stats['raw_edge_count'] += 1
            p1, p2 = row['protein1'], row['protein2']
            
            m1 = prot_to_gene.get(p1)
            m2 = prot_to_gene.get(p2)
            
            keep = True
            reason = ""
            
            if not m1:
                keep = False
                reason = "unmapped_A"
            elif not m2:
                keep = False
                reason = "unmapped_B"
            elif m1 == m2:
                keep = False
                reason = "self_loop"
                stats['self_loop_removed_count'] += 1
            
            if keep:
                # Canonicalize (undirected)
                a, b = sorted([m1, m2])
                edge_key = (a, b)
                
                scores = {
                    'combined_score': int(row['combined_score']),
                    'coexpression': int(row['coexpression']),
                    'experimental': int(row['experimental']),
                    'database': int(row['database']),
                }
                
                if edge_key in final_edges:
                    stats['duplicate_collapsed_count'] += 1
                    for k in scores:
                        final_edges[edge_key][k] = max(final_edges[edge_key][k], scores[k])
                else:
                    final_edges[edge_key] = scores
                    stats['kept_edge_count'] += 1
            else:
                stats['dropped_edge_count'] += 1
                
            # Audit first 500 edges
            if stats['raw_edge_count'] <= 500:
                edge_audit.append({
                    'raw_protein1': p1, 'raw_protein2': p2,
                    'mapped_A': m1 if m1 else "", 'mapped_B': m2 if m2 else "",
                    'keep_edge': keep, 'drop_reason': reason,
                    'combined_score': row['combined_score'],
                    'coexpression': row['coexpression'],
                    'experimental': row['experimental'],
                    'database': row['database']
                })

    # Save outputs
    # 1. string.csv
    final_edge_list = []
    final_nodes = set()
    for (a, b), scores in final_edges.items():
        final_edge_list.append({'A': a, 'B': b, **scores})
        final_nodes.add(a)
        final_nodes.add(b)
    
    pd.DataFrame(final_edge_list).to_csv(out_dir / 'string.csv', index=False)
    
    # 2. mapping and audit
    mapping_df.to_csv(out_dir / 'string_id_mapping.tsv', sep='\t', index=False)
    pd.DataFrame(edge_audit).to_csv(out_dir / 'string_edge_audit.tsv', sep='\t', index=False)
    
    # 3. Summary
    stats['final_unique_node_count'] = len(final_nodes)
    stats['final_unique_edge_count'] = len(final_edge_list)
    stats['mapped_label_gene_count'] = len(final_nodes)
    stats['label_gene_coverage_ratio'] = stats['final_unique_node_count'] / stats['label_gene_count'] if stats['label_gene_count'] > 0 else 0
    stats['aggregation_rule'] = "max"
    stats['mapped_unique_raw_protein_count'] = mapping_summary['mapped']
    stats['unmapped_unique_raw_protein_count'] = mapping_summary['unmapped']
    stats['ambiguous_unique_raw_protein_count'] = mapping_summary['ambiguous']
    stats['raw_unique_protein_count'] = mapping_summary['total']
    
    pd.DataFrame([stats]).to_csv(out_dir / 'string_summary.tsv', sep='\t', index=False)
    
    if stats['label_gene_coverage_ratio'] < 0.30:
        print(f"[WARNING] {species_name} coverage ratio low: {stats['label_gene_coverage_ratio']:.2%}")
        
    print(f"[DONE] {species_name} processed. Kept {stats['kept_edge_count']} edges.")
    return stats

def main():
    # 1. Prepare Fusarium tx2gene if available
    fusarium_tx2gene = None
    if FGRAM_PH1_TSV.exists():
        print(f"[INFO] Building Fusarium transcript-to-gene map from {FGRAM_PH1_TSV}")
        fusarium_tx2gene = build_fusarium_transcript_to_gene_map(FGRAM_PH1_TSV)
        print(f"[INFO] Loaded {len(fusarium_tx2gene)} transcript mappings.")
    else:
        print(f"[WARNING] Fusarium PH-1 TSV not found at {FGRAM_PH1_TSV}")

    all_stats = []
    for taxon_id, species_name in SPECIES_MAP.items():
        res = process_species(taxon_id, species_name, fusarium_tx2gene if species_name == 'fgraminearum' else None)
        if res:
            all_stats.append(res)
            
    # Global Summary
    if all_stats:
        df_global = pd.DataFrame(all_stats)
        df_global.to_csv(PPI_ROOT / 'string_build_summary.tsv', sep='\t', index=False)
        
        # README
        with open(PPI_ROOT / 'README_string_build.md', 'w') as f:
            f.write("# STRING PPI Dataset Build Report (v12.0)\n\n")
            f.write("## Global Statistics\n\n")
            f.write(df_global[['species', 'taxon_id', 'kept_edge_count', 'final_unique_node_count', 'label_gene_coverage_ratio']].to_markdown(index=False))
            f.write("\n\n## Methodology\n")
            f.write("1. **Mapping Source**: Strictly used STRING `.aliases.v12.0.txt.gz` files.\n")
            f.write("2. **Fusarium Rule**: For `fgraminearum`, used `FgraminearumPH-1.tsv` to map transcript aliases to gene IDs.\n")
            f.write("3. **Filtering**: Only edges where both nodes map to IDs in `labels.standard.tsv` are kept.\n")
            f.write("4. **Aggregation**: Duplicate edges are collapsed using the **MAX** score.\n")
            f.write("5. **Format**: Output follows the EPGAT `string.csv` format.\n")
            
    print("\n[ALL DONE]")

if __name__ == "__main__":
    main()
