import pandas as pd
from pathlib import Path
import numpy as np

# Paths
COUNTS_FILE = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/EXP/Download/fgraminearum_GSE292521_Fg_SW_counts.txt')
REF_FILE = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/Fusarium_v2/data/raw/essential_gene/ncbi_gene_reference/Fgraminearum_Summary-v1.tsv')
LABELS_FILE = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/essential_gene/fgraminearum/labels.standard.tsv')
OUTPUT_DIR = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/EXP/fgraminearum')
GLOBAL_SUMMARY = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/EXP/profile_build_summary.tsv')

def calculate_tpm(counts_df, lengths):
    """
    Calculates TPM from counts and lengths.
    counts_df: DataFrame where index is gene_id and columns are samples.
    lengths: Series where index is gene_id and values are lengths in base pairs.
    """
    # Align lengths with counts_df
    common_genes = counts_df.index.intersection(lengths.index)
    print(f"[INFO] Genes in counts: {len(counts_df)}, Genes in ref: {len(lengths)}, Intersection: {len(common_genes)}")
    
    counts_df = counts_df.loc[common_genes]
    lengths = lengths.loc[common_genes]
    
    # RPK = count / (length / 1000)
    rpk = counts_df.divide(lengths / 1000, axis=0)
    
    # TPM = RPK / (sum(RPK) / 1,000,000)
    tpm = rpk.divide(rpk.sum(axis=0) / 1e6, axis=1)
    
    return tpm

def main():
    print("--- Updating Fusarium Expression Data to TPM ---")
    
    # 1. Load Reference for lengths
    print(f"[INFO] Loading reference from {REF_FILE}...")
    ref_df = pd.read_csv(REF_FILE, sep='\t')
    # Column names might have spaces
    ref_df.columns = [c.strip() for c in ref_df.columns]
    lengths = ref_df.set_index('Gene ID')['Transcript Length']
    
    # 2. Load Counts
    print(f"[INFO] Loading counts from {COUNTS_FILE}...")
    # The first column is Gene ID but header might be tricky
    counts_df = pd.read_csv(COUNTS_FILE, sep='\t', index_col=0)
    
    # 3. Calculate TPM
    tpm_df = calculate_tpm(counts_df, lengths)
    
    # 4. Filter by valid labels
    print(f"[INFO] Loading valid labels from {LABELS_FILE}...")
    labels_df = pd.read_csv(LABELS_FILE, sep='\t')
    valid_genes = set(labels_df['gene_id'].astype(str))
    
    # Only keep genes present in labels
    final_genes = tpm_df.index.intersection(valid_genes)
    tpm_df = tpm_df.loc[final_genes]
    print(f"[INFO] Final gene count in profile: {len(tpm_df)} (out of {len(valid_genes)} labels)")
    
    # 5. Save Outputs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Profile.csv
    profile_df = tpm_df.reset_index().rename(columns={'index': 'Gene'})
    profile_df.to_csv(OUTPUT_DIR / 'profile.csv', index=False)
    print(f"[DONE] Wrote {OUTPUT_DIR / 'profile.csv'}")
    
    # Mapping table (trivial here as IDs are same)
    mapping_df = pd.DataFrame({
        'raw_id': final_genes,
        'mapped_gene_id': final_genes,
        'mapping_status': 'mapped',
        'mapping_method': 'direct_count_id'
    })
    mapping_df.to_csv(OUTPUT_DIR / 'exp_id_mapping.tsv', sep='\t', index=False)
    
    # Summary
    summary = {
        'species': 'fgraminearum',
        'labels_path': str(LABELS_FILE),
        'input_source': COUNTS_FILE.name,
        'series_matrix_file': 'NA',
        'family_soft_file': 'NA',
        'platform_gpl': 'NA',
        'raw_row_count': len(counts_df),
        'raw_unique_probe_count': len(counts_df),
        'mapped_unique_id_count': len(final_genes),
        'label_gene_count': len(valid_genes),
        'mapped_label_gene_count': len(final_genes),
        'label_gene_coverage_ratio': len(final_genes) / len(valid_genes),
        'experiment_count': 1,
        'exp_group_count': len(tpm_df.columns),
        'final_gene_count_in_profile': len(final_genes),
        'duplicate_gene_aggregation': 'mean',
        'notes': "Converted from counts to TPM using reference lengths."
    }
    pd.DataFrame([summary]).to_csv(OUTPUT_DIR / 'exp_summary.tsv', sep='\t', index=False)
    
    # Update Global Summary if it exists
    if GLOBAL_SUMMARY.exists():
        gs_df = pd.read_csv(GLOBAL_SUMMARY, sep='\t')
        # Remove old fgraminearum record
        gs_df = gs_df[gs_df['species'] != 'fgraminearum']
        # Add new one
        new_row = pd.DataFrame([summary])
        gs_df = pd.concat([gs_df, new_row], ignore_index=True)
        gs_df.to_csv(GLOBAL_SUMMARY, sep='\t', index=False)
        print(f"[INFO] Updated {GLOBAL_SUMMARY}")

if __name__ == "__main__":
    main()
