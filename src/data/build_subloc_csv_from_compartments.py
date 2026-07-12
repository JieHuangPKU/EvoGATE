import os
import sys
import pandas as pd
from pathlib import Path
from typing import Dict, Set, List, Optional
from utils_subloc_mapping import (
    find_labels_file,
    load_labels_table,
    build_valid_gene_id_index,
    load_string_aliases,
    load_compartments_table,
    build_fusarium_transcript_to_gene_map,
    map_subloc_ids_to_gene_ids,
    build_gene_location_long_table,
    summarize_location_frequency,
    select_final_location_columns,
    build_binary_location_matrix
)

# Constants
SPECIES_TAXON = {
    "fgraminearum": "229533",
    "scerevisiae": "4932",
    "celegans": "6239",
    "melanogaster": "7227",
    "human": "9606",
}

SPECIES_INPUT_FILE = {
    "fgraminearum": "fgraminearum_eFG_sublocation.txt",
    "scerevisiae": "yeast_compartment_integrated_full.tsv",
    "celegans": "worm_compartment_knowledge_full.tsv",
    "melanogaster": "fly_compartment_knowledge_full.tsv",
    "human": "human_compartment_knowledge_full.tsv",
}

PROGATE_ROOT = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2')
PROCESSED_ROOT = PROGATE_ROOT / 'data' / 'processed'
LC_ROOT = PROCESSED_ROOT / 'LC'
COMPARTMENTS_DIR = LC_ROOT / 'COMPARTMENTS'
PPI_DIR = PROCESSED_ROOT / 'PPI'
STRING_DB_ROOT = PPI_DIR / 'stringDB'
FGRAM_PH1_TSV = Path('/Users/jiehuang/work/2025禾谷镰刀菌/RNA农药联合攻关/小分子/FgraminearumPH-1.tsv')

def build_species_subloc(species: str, fusarium_tx2gene: Optional[Dict[str, str]] = None):
    print(f"\n[INFO] >>> Processing Subcellular Localization for: {species}")
    
    # 1. Find and load labels
    try:
        labels_path = find_labels_file(PROCESSED_ROOT, species)
        print(f"[INFO] Using labels file: {labels_path}")
        labels_df = load_labels_table(labels_path)
        valid_id_index = build_valid_gene_id_index(labels_df)
        all_valid_gene_ids = labels_df['gene_id'].astype(str).tolist()
        label_gene_count = len(labels_df)
        print(f"[INFO] Indexed valid ID aliases for {label_gene_count} genes.")
    except Exception as e:
        print(f"[ERROR] Failed to load labels for {species}: {e}")
        return None

    # 2. Load raw localization table
    input_file = COMPARTMENTS_DIR / SPECIES_INPUT_FILE[species]
    if not input_file.exists():
        print(f"[ERROR] Input file not found: {input_file}")
        return None
    
    subloc_df = load_compartments_table(input_file, species)
    raw_record_count = len(subloc_df)
    source_db = "eFG" if species == 'fgraminearum' else "COMPARTMENTS"
    print(f"[INFO] Loaded {raw_record_count} raw records from {source_db}.")

    # 3. Load STRING aliases for mapping bridge
    taxon_id = SPECIES_TAXON[species]
    alias_path = STRING_DB_ROOT / f"{taxon_id}.protein.aliases.v12.0.txt.gz"
    alias_df = load_string_aliases(alias_path)
    print(f"[INFO] Loaded {len(alias_df)} alias entries for mapping bridge.")

    # 4. Map IDs
    mapped_df = map_subloc_ids_to_gene_ids(subloc_df, valid_id_index, alias_df, species, fusarium_tx2gene)
    
    # 5. Build Long Table & Normalize
    gene_loc_long_df = build_gene_location_long_table(mapped_df)
    gene_loc_long_df['source_db'] = source_db
    
    # 6. Frequency analysis
    freq_df = summarize_location_frequency(gene_loc_long_df)
    
    # 7. Select columns (top 12-15)
    final_cols = select_final_location_columns(freq_df)
    print(f"[INFO] Selected {len(final_cols)} location columns: {', '.join(final_cols)}")
    
    # 8. Build Binary Matrix
    binary_matrix_df = build_binary_location_matrix(gene_loc_long_df, all_valid_gene_ids, final_cols)
    
    final_gene_count_with_loc = gene_loc_long_df['gene_id'].nunique()
    coverage_ratio = final_gene_count_with_loc / label_gene_count if label_gene_count > 0 else 0
    print(f"[INFO] Final gene count with at least one location: {final_gene_count_with_loc} ({coverage_ratio:.2%})")

    # 9. Save outputs
    out_dir = LC_ROOT / species
    out_dir.mkdir(parents=True, exist_ok=True)
    
    binary_matrix_df.to_csv(out_dir / 'subloc.csv', index=False)
    gene_loc_long_df.to_csv(out_dir / 'gene_location_long.tsv', sep='\t', index=False)
    freq_df.to_csv(out_dir / 'location_frequency.tsv', sep='\t', index=False)
    
    # Mapping details
    id_mapping_df = mapped_df[['raw_id', 'mapped_gene_id', 'mapping_status', 'mapping_method', 'mapping_evidence']].drop_duplicates('raw_id')
    id_mapping_df.to_csv(out_dir / 'subloc_id_mapping.tsv', sep='\t', index=False)
    
    # Audit
    audit_df = mapped_df.copy()
    audit_df['keep_record'] = audit_df['mapping_status'] == 'mapped'
    audit_df['drop_reason'] = audit_df['mapping_status'].apply(lambda x: "" if x == 'mapped' else x)
    audit_df['source_file'] = str(input_file.name)
    audit_df['source_db'] = source_db
    audit_df = audit_df.rename(columns={'location': 'raw_location'})
    audit_cols = ['raw_id', 'raw_location', 'mapped_gene_id', 'keep_record', 'drop_reason', 'source_file', 'source_db']
    audit_df[audit_cols].to_csv(out_dir / 'subloc_audit.tsv', sep='\t', index=False)
    
    # Summary
    summary = {
        'species': species,
        'labels_path': str(labels_path),
        'input_file': str(input_file.name),
        'source_db': source_db,
        'raw_record_count': raw_record_count,
        'label_gene_count': label_gene_count,
        'mapped_label_gene_count': final_gene_count_with_loc,
        'label_gene_coverage_ratio': coverage_ratio,
        'final_location_column_count': len(final_cols),
        'final_location_columns': ";".join(final_cols),
        'top_location_by_gene_count': final_cols[0] if final_cols else "",
        'location_selection_rule': "top_15_by_frequency",
        'location_normalization_version': "v2_controlled_vocabulary",
        'notes': f"Binary matrix with {len(final_cols)} columns."
    }
    pd.DataFrame([summary]).to_csv(out_dir / 'subloc_summary.tsv', sep='\t', index=False)
    
    return summary

def main():
    # 1. Prepare Fusarium tx2gene if available
    fusarium_tx2gene = build_fusarium_transcript_to_gene_map(FGRAM_PH1_TSV)
    if fusarium_tx2gene:
        print(f"[INFO] Loaded {len(fusarium_tx2gene)} Fusarium transcript mappings.")

    all_summaries = []
    for species in SPECIES_TAXON.keys():
        res = build_species_subloc(species, fusarium_tx2gene if species == 'fgraminearum' else None)
        if res:
            all_summaries.append(res)
            
    if not all_summaries:
        print("[ERROR] No species processed successfully.")
        return

    # Global summary
    df_global = pd.DataFrame(all_summaries)
    df_global.to_csv(LC_ROOT / 'subloc_build_summary.tsv', sep='\t', index=False)
    
    # README
    with open(LC_ROOT / 'README_subloc_build.md', 'w') as f:
        f.write("# Subcellular Localization Data Build Report (Binary Matrix Version)\n\n")
        f.write("## Overview\n")
        f.write("This directory contains standardized subcellular localization data for 5 species in binary matrix format (0/1).\n\n")
        f.write("## Global Statistics\n\n")
        f.write(df_global[['species', 'source_db', 'final_location_column_count', 'mapped_label_gene_count', 'label_gene_coverage_ratio']].to_markdown(index=False))
        f.write("\n\n## Methodology\n")
        f.write("1. **Format**: Output `subloc.csv` is a binary matrix where rows are Genes and columns are standard Location categories.\n")
        f.write("2. **Normalization**: Raw terms are mapped to a controlled vocabulary (e.g., 'cytosol' -> 'Cytoplasm'). Synonymous terms are merged.\n")
        f.write("3. **Column Selection**: For each species, the top 12–15 standard categories (by gene count) are selected as columns.\n")
        f.write("4. **ID Mapping**: Aligned to `labels.standard.tsv` using direct matches and STRING aliases as a bridge.\n")
        f.write("5. **Output Columns**: Columns are sorted by frequency (most frequent first).\n\n")
        
        f.write("## Controlled Vocabulary (Selection Pool)\n")
        f.write("- Cytoplasm\n- Nucleus\n- Cell membrane\n- Multi-pass membrane\n- Mitochondrion\n- Endoplasmic reticulum\n- Golgi\n- Endosome\n- Vacuole\n- Peroxisome\n- Secreted\n- Lysosome\n- Cytoskeleton\n- Ribosome\n- Extracellular matrix\n\n")

    print(f"\n[ALL DONE] Global summary saved to {LC_ROOT / 'subloc_build_summary.tsv'}")

if __name__ == "__main__":
    main()
