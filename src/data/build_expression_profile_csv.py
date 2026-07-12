import os
import sys
import pandas as pd
from pathlib import Path
from typing import Dict, Set, List, Optional, Union
from utils_expression_mapping import (
    find_labels_file,
    load_labels_table,
    build_valid_gene_id_index,
    load_string_aliases,
    extract_platform_annotation_from_family_soft,
    load_geo_series_matrix,
    detect_expression_columns,
    build_probe_to_gene_mapping,
    map_expression_ids_to_gene_ids,
    aggregate_gene_level_expression,
    calculate_tpm
)

# Constants
SPECIES_TAXON = {
    "fgraminearum": "229533",
    "scerevisiae": "4932",
    "celegans": "6239",
    "melanogaster": "7227",
    "human": "9606",
}

SPECIES_FILE_MAP = {
    "celegans": {
        "matrix": "celegans_GSE31422-GPL14144_series_matrix.txt.gz",
        "soft": "GSE31422_family.soft.gz"
    },
    "fgraminearum": {
        "matrix": "fgraminearum_GSE292521_Fg_SW_counts.txt",
        "soft": None, # Use custom TPM logic
        "ref": "/Users/jiehuang/work/2025禾谷镰刀菌/程序/Fusarium_v2/data/raw/essential_gene/ncbi_gene_reference/Fgraminearum_Summary-v1.tsv"
    },
    "human": {
        "matrix": "human_GSE86354_GTEx_FPKM.csv",
        "soft": None
    },
    "melanogaster": {
        "matrix": "melanogaster_GSE67547_series_matrix.txt.gz",
        "soft": "GSE67547_family.soft.gz"
    },
    "scerevisiae": {
        "matrix": "scerevisiae_GSE3431_series_matrix.txt.gz",
        "soft": "GSE3431_family.soft.gz"
    },
}

PROGATE_ROOT = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2')
PROCESSED_ROOT = PROGATE_ROOT / 'data' / 'processed'
EXP_ROOT = PROCESSED_ROOT / 'EXP'
DOWNLOAD_DIR = EXP_ROOT / 'Download'
PPI_DIR = PROCESSED_ROOT / 'PPI'
STRING_DB_ROOT = PPI_DIR / 'stringDB'

def build_species_profile(species: str):
    print(f"\n[INFO] >>> Processing Expression Profile for: {species}")
    
    # 1. Find and load labels
    try:
        labels_path = find_labels_file(PROCESSED_ROOT, species)
        print(f"[INFO] Using labels file: {labels_path}")
        labels_df = load_labels_table(labels_path)
        valid_id_index = build_valid_gene_id_index(labels_df)
        all_valid_gene_ids = labels_df['gene_id'].astype(str).tolist()
        label_gene_count = len(labels_df)
        print(f"[INFO] Indexed valid IDs for {label_gene_count} genes.")
    except Exception as e:
        print(f"[ERROR] Failed to load labels for {species}: {e}"); return None

    # 2. Load Expression and Annotation
    file_info = SPECIES_FILE_MAP[species]
    matrix_file = DOWNLOAD_DIR / file_info['matrix']
    soft_file = DOWNLOAD_DIR / file_info['soft'] if file_info['soft'] else None
    
    if not matrix_file.exists():
        print(f"[ERROR] Matrix file not found: {matrix_file}"); return None

    # 3. Load STRING aliases for mapping bridge
    taxon_id = SPECIES_TAXON[species]
    alias_path = STRING_DB_ROOT / f"{taxon_id}.protein.aliases.v12.0.txt.gz"
    alias_df = load_string_aliases(alias_path)
    print(f"[INFO] Loaded {len(alias_df)} alias entries for mapping bridge.")

    platform_gpl = "NA"
    
    # 4. Handle Data Loading and Mapping
    if species == 'fgraminearum':
        print(f"[INFO] fgraminearum: switching source to {matrix_file.name}")
        counts_df = pd.read_csv(matrix_file, sep='\t', index_col=0)
        print(f"[INFO] Loaded count matrix rows={len(counts_df)} cols={len(counts_df.columns)}")
        
        print(f"[INFO] Loading transcript lengths from {file_info['ref']}")
        ref_df = pd.read_csv(file_info['ref'], sep='\t')
        ref_df.columns = [c.strip() for c in ref_df.columns]
        id_col_ref = next((c for c in ref_df.columns if 'gene' in c.lower()), ref_df.columns[0])
        len_col_ref = next((c for c in ref_df.columns if 'length' in c.lower() and 'transcript' in c.lower()), 'Transcript Length')
        print(f"[INFO] Using length column: {len_col_ref}")
        
        lengths = ref_df.set_index(id_col_ref)[len_col_ref]
        exp_df_numeric = calculate_tpm(counts_df, lengths)
        print(f"[INFO] Converted counts to TPM. Intersection with lengths: {len(exp_df_numeric)}")
        
        exp_df = exp_df_numeric.reset_index().rename(columns={'index': 'raw_id'})
        id_col = 'raw_id'
        expr_cols = list(exp_df_numeric.columns)
        
        id_mapping = []
        for rid in exp_df[id_col]:
            rid_str = str(rid).strip()
            found_gene = rid_str if rid_str in valid_id_index else ""
            status = "mapped" if found_gene else "unmapped"
            id_mapping.append({'raw_probe_id': rid, 'mapped_gene_id': found_gene, 'mapping_status': status, 'mapping_method': 'direct_gene_id'})
        probe_map_df = pd.DataFrame(id_mapping)
        mapped_df = exp_df.merge(probe_map_df, left_on=id_col, right_on='raw_probe_id', how='left')

    elif soft_file:
        exp_df, platform_gpls = load_geo_series_matrix(matrix_file)
        print(f"[INFO] Loaded matrix {matrix_file.name} with {len(exp_df)} rows. GPLs: {platform_gpls}")
        platforms = extract_platform_annotation_from_family_soft(soft_file)
        all_probe_maps = []
        for gpl, ann_df in platforms.items():
            pm = build_probe_to_gene_mapping(ann_df, valid_id_index, alias_df, species)
            all_probe_maps.append(pm)
        probe_map_df = pd.concat(all_probe_maps).drop_duplicates('raw_probe_id')
        mapped_df = map_expression_ids_to_gene_ids(exp_df, probe_map_df, species)
        id_col = exp_df.columns[0]
        expr_cols = detect_expression_columns(exp_df, id_col)
        platform_gpl = ";".join(platform_gpls)
    else:
        # Human CSV case
        print(f"[INFO] Loading matrix: {matrix_file.name}")
        exp_df = pd.read_csv(matrix_file)
        id_col = exp_df.columns[0]
        dummy_ann_df = pd.DataFrame({'ID': exp_df[id_col], 'ORF': exp_df[id_col]})
        probe_map_df = build_probe_to_gene_mapping(dummy_ann_df, valid_id_index, alias_df, species)
        mapped_df = map_expression_ids_to_gene_ids(exp_df, probe_map_df, species)
        expr_cols = detect_expression_columns(exp_df, id_col)
        platform_gpl = "NA"

    # 5. Aggregate to gene-level
    final_profile_df_mapped = aggregate_gene_level_expression(mapped_df, expr_cols)
    
    # Ensure all labels are present, fill missing with 0
    final_profile_df = pd.DataFrame({'Gene': all_valid_gene_ids})
    final_profile_df = final_profile_df.merge(final_profile_df_mapped, on='Gene', how='left').fillna(0)
    
    final_gene_count_mapped = final_profile_df_mapped['Gene'].nunique()
    coverage_ratio = final_gene_count_mapped / label_gene_count if label_gene_count > 0 else 0
    print(f"[INFO] Final gene count in profile.csv: {len(final_profile_df)} (Mapped: {final_gene_count_mapped}, {coverage_ratio:.2%})")

    # 6. Save outputs
    out_dir = EXP_ROOT / species; out_dir.mkdir(parents=True, exist_ok=True)
    final_profile_df.to_csv(out_dir / 'profile.csv', index=False)
    probe_map_df.to_csv(out_dir / 'exp_id_mapping.tsv', sep='\t', index=False)
    
    audit_cols = ['raw_probe_id', 'mapped_gene_id', 'mapping_status', 'mapping_method']
    available_cols = [c for c in audit_cols if c in mapped_df.columns]
    audit_df = mapped_df[available_cols].copy()
    audit_df['keep_record'] = mapped_df['mapping_status'] == 'mapped' if 'mapping_status' in mapped_df.columns else False
    audit_df['source_file'] = matrix_file.name
    audit_df['platform_gpl'] = platform_gpl
    audit_df.to_csv(out_dir / 'exp_audit.tsv', sep='\t', index=False)
    
    summary = {
        'species': species, 'labels_path': str(labels_path), 'input_source': matrix_file.name,
        'series_matrix_file': matrix_file.name if soft_file else "NA", 
        'family_soft_file': soft_file.name if soft_file else "NA",
        'platform_gpl': platform_gpl, 'raw_row_count': len(exp_df), 'raw_unique_probe_count': exp_df[id_col].nunique() if id_col in exp_df.columns else 0,
        'mapped_unique_id_count': final_gene_count_mapped, 'label_gene_count': label_gene_count,
        'mapped_label_gene_count': final_gene_count_mapped, 'label_gene_coverage_ratio': coverage_ratio,
        'experiment_count': 1, 'exp_group_count': len(expr_cols),
        'final_gene_count_in_profile': len(final_profile_df), 'duplicate_gene_aggregation': 'mean', 
        'notes': "Converted counts to TPM using reference lengths. Missing genes filled with 0." if species == 'fgraminearum' else "Missing genes filled with 0."
    }
    pd.DataFrame([summary]).to_csv(out_dir / 'exp_summary.tsv', sep='\t', index=False)
    return summary

def main():
    all_summaries = []
    for species in SPECIES_FILE_MAP.keys():
        res = build_species_profile(species)
        if res: all_summaries.append(res)
    if not all_summaries: return
    df_global = pd.DataFrame(all_summaries)
    df_global.to_csv(EXP_ROOT / 'profile_build_summary.tsv', sep='\t', index=False)
    with open(EXP_ROOT / 'README_profile_build.md', 'w') as f:
        f.write("# Expression Profile Data Build Report (Final Integrated Version)\n\n")
        f.write("## Global Statistics\n\n")
        f.write(df_global[['species', 'input_source', 'exp_group_count', 'final_gene_count_in_profile', 'label_gene_coverage_ratio']].to_markdown(index=False))
        f.write("\n\n## Methodology\n")
        f.write("1. **GPL Integration**: Parsed `GSE*_family.soft.gz` to extract platform annotations (GPL) for microarray datasets (celegans, melanogaster, scerevisiae).\n")
        f.write("2. **TPM Conversion**: For `fgraminearum`, replaced old microarray source with RNA-seq counts (`fgraminearum_GSE292521_Fg_SW_counts.txt`) and converted to TPM using reference transcript lengths.\n")
        f.write("3. **Heuristic Mapping**: Aligned IDs to standard Gene IDs using ORF/Symbol/Aliases/WormBase.\n")
        f.write("4. **Aggregation**: Multiple probes/records mapping to the same gene are averaged (**MEAN**).\n")
        f.write("5. **Completeness**: All genes in `labels.standard.tsv` are included in `profile.csv`. Missing data points are filled with **0**.\n")
        f.write("6. **Independence**: Process is self-contained and aligns strictly to `labels.standard.tsv`.\n\n")
    print(f"\n[ALL DONE] Global summary saved to {EXP_ROOT / 'profile_build_summary.tsv'}")

if __name__ == "__main__":
    main()
