import os
import sys
import pandas as pd
from pathlib import Path
import json
import numpy as np

# Set directories
EPGAT_DIR = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/EPGAT/data/essential_genes')
BINGO_DIR = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/Bingo/essential_genes/data')
PROGATE_DIR = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data')

# Output directories
AUDIT_DIR = PROGATE_DIR / 'audits'
PROCESSED_DIR = PROGATE_DIR / 'processed' / 'essential_gene'
MANIFEST_DIR = PROGATE_DIR / 'manifests'

for d in [AUDIT_DIR, PROCESSED_DIR, MANIFEST_DIR]:
    d.mkdir(parents=True, exist_ok=True)

SPECIES_LIST = ['celegans', 'fgraminearum', 'human', 'melanogaster', 'scerevisiae']

def scan_directory(base_path):
    inventory = []
    if not base_path.exists():
        return inventory
    for root, dirs, files in os.walk(base_path):
        for f in files:
            if f.startswith('.'):
                continue
            file_path = Path(root) / f
            rel_path = file_path.relative_to(base_path)
            
            species = rel_path.parts[0] if len(rel_path.parts) > 0 else 'unknown'
            if species not in SPECIES_LIST:
                species = 'unknown'
                
            stat = file_path.stat()
            
            role = 'unknown'
            if 'split' in f.lower() or 'fold' in f.lower(): role = 'split'
            elif 'train' in f.lower(): role = 'split_train'
            elif 'test' in f.lower(): role = 'split_test'
            elif 'val' in f.lower(): role = 'split_val'
            elif f.lower() in ['essential_nr_list.txt', 'e.txt']: role = 'essential_list'
            elif f.lower() in ['nessential_nr_list.txt', 'ne.txt']: role = 'nonessential_list'
            elif f.lower() in ['gene_list.txt', 'ogee.csv', 'all.txt'] or f.endswith('_all.txt'): role = 'labels'
            
            inventory.append({
                'project': base_path.parent.parent.name if 'Bingo' in str(base_path) else base_path.parent.parent.name, # Adjust this depending on exact structure
                'species': species,
                'relative_path': str(rel_path),
                'file_name': f,
                'file_type': file_path.suffix,
                'size': stat.st_size,
                'modified_time': stat.st_mtime,
                'guessed_role': role,
                'full_path': str(file_path)
            })
    return inventory

print("Scanning directories...")
epgat_inv = scan_directory(EPGAT_DIR)
bingo_inv = scan_directory(BINGO_DIR)

for row in epgat_inv: row['project'] = 'epgat'
for row in bingo_inv: row['project'] = 'bingo'

pd.DataFrame(epgat_inv).drop(columns=['full_path']).to_csv(AUDIT_DIR / 'epgat_file_inventory.tsv', sep='\t', index=False)
pd.DataFrame(bingo_inv).drop(columns=['full_path']).to_csv(AUDIT_DIR / 'bingo_file_inventory.tsv', sep='\t', index=False)

def read_gene_list(path):
    try:
        with open(path, 'r') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    except:
        return []

def load_labels_bingo(species, bingo_inv):
    # Try gene_list.txt
    files = [r for r in bingo_inv if r['species'] == species and r['file_name'] == 'gene_list.txt']
    if files:
        df = pd.read_csv(files[0]['full_path'], sep='\t')
        if 'Target' in df.columns:
            id_col = 'Ensembl' if 'Ensembl' in df.columns else df.columns[0]
            res = df[[id_col, 'Target']].copy()
            res.columns = ['gene_id', 'label']
            res['source'] = files[0]['relative_path']
            return res
            
    # Try Essential_NR_list.txt / NEssential_NR_list.txt
    e_files = [r for r in bingo_inv if r['species'] == species and r['file_name'] == 'Essential_NR_list.txt']
    ne_files = [r for r in bingo_inv if r['species'] == species and r['file_name'] == 'NEssential_NR_list.txt']
    
    if e_files and ne_files:
        e_genes = read_gene_list(e_files[0]['full_path'])
        ne_genes = read_gene_list(ne_files[0]['full_path'])
        
        df_e = pd.DataFrame({'gene_id': e_genes, 'label': 1, 'source': e_files[0]['relative_path']})
        df_ne = pd.DataFrame({'gene_id': ne_genes, 'label': 0, 'source': ne_files[0]['relative_path']})
        return pd.concat([df_e, df_ne], ignore_index=True)
    return pd.DataFrame(columns=['gene_id', 'label', 'source'])

def load_labels_epgat(species, epgat_inv):
    # Try Essential_NR_list.txt / NEssential_NR_list.txt first as they are most explicit
    e_files = [r for r in epgat_inv if r['species'] == species and r['file_name'] == 'Essential_NR_list.txt']
    ne_files = [r for r in epgat_inv if r['species'] == species and r['file_name'] == 'NEssential_NR_list.txt']
    if e_files and ne_files:
        e_genes = read_gene_list(e_files[0]['full_path'])
        ne_genes = read_gene_list(ne_files[0]['full_path'])
        df_e = pd.DataFrame({'gene_id': e_genes, 'label': 1, 'source': e_files[0]['relative_path']})
        df_ne = pd.DataFrame({'gene_id': ne_genes, 'label': 0, 'source': ne_files[0]['relative_path']})
        return pd.concat([df_e, df_ne], ignore_index=True)
        
    # Try ogee.csv
    ogee_files = [r for r in epgat_inv if r['species'] == species and r['file_name'] == 'ogee.csv']
    if ogee_files:
        df = pd.read_csv(ogee_files[0]['full_path'])
        if 'Gene' in df.columns and 'Label' in df.columns:
            res = df[['Gene', 'Label']].copy()
            res.columns = ['gene_id', 'label']
            res['source'] = ogee_files[0]['relative_path']
            return res
            
    # Check for <taxID>_all.txt or other _all files
    all_files = [r for r in epgat_inv if r['species'] == species and ('_all' in r['file_name'] or 'genes.csv' in r['file_name'])]
    if all_files:
        try:
            df = pd.read_csv(all_files[0]['full_path'], sep=None, engine='python')
            if 'essential' in df.columns and 'locus' in df.columns:
                df['label'] = df['essential'].apply(lambda x: 1 if x == 'E' else 0 if x == 'NE' else None)
                df = df.dropna(subset=['label'])
                res = df[['locus', 'label']].copy()
                res.columns = ['gene_id', 'label']
                res['source'] = all_files[0]['relative_path']
                return res
            if 'essentiality' in df.columns and 'locus' in df.columns:
                df['label'] = df['essentiality'].apply(lambda x: 1 if x == 'E' else 0 if x == 'NE' else None)
                df = df.dropna(subset=['label'])
                res = df[['locus', 'label']].copy()
                res.columns = ['gene_id', 'label']
                res['source'] = all_files[0]['relative_path']
                return res
        except:
            pass

    # Check for E.txt and NE.txt or ALL.txt
    e_txt = [r for r in epgat_inv if r['species'] == species and r['file_name'] == 'E.txt']
    ne_txt = [r for r in epgat_inv if r['species'] == species and r['file_name'] == 'NE.txt']
    all_txt = [r for r in epgat_inv if r['species'] == species and r['file_name'] == 'ALL.txt']
    
    if e_txt:
        e_genes = read_gene_list(e_txt[0]['full_path'])
        if e_genes and e_genes[0] == 'gene_id': e_genes = e_genes[1:]
        df_e = pd.DataFrame({'gene_id': e_genes, 'label': 1, 'source': e_txt[0]['relative_path']})
        
        if ne_txt:
            ne_genes = read_gene_list(ne_txt[0]['full_path'])
            if ne_genes and ne_genes[0] == 'gene_id': ne_genes = ne_genes[1:]
            df_ne = pd.DataFrame({'gene_id': ne_genes, 'label': 0, 'source': ne_txt[0]['relative_path']})
            return pd.concat([df_e, df_ne], ignore_index=True)
        elif all_txt:
            a_genes = read_gene_list(all_txt[0]['full_path'])
            if a_genes and a_genes[0] == 'gene_id': a_genes = a_genes[1:]
            ne_genes = list(set(a_genes) - set(e_genes))
            df_ne = pd.DataFrame({'gene_id': ne_genes, 'label': 0, 'source': all_txt[0]['relative_path']})
            return pd.concat([df_e, df_ne], ignore_index=True)
            
    return pd.DataFrame(columns=['gene_id', 'label', 'source'])

print("Parsing files...")

audit_records = []
epgat_dfs = {}
bingo_dfs = {}

# Step 3 & 4
for sp in SPECIES_LIST:
    df_epgat = load_labels_epgat(sp, epgat_inv)
    df_bingo = load_labels_bingo(sp, bingo_inv)
    
    if not df_epgat.empty:
        df_epgat['gene_id'] = df_epgat['gene_id'].astype(str).str.strip()
        df_epgat = df_epgat.drop_duplicates(subset=['gene_id'])
    if not df_bingo.empty:
        df_bingo['gene_id'] = df_bingo['gene_id'].astype(str).str.strip()
        df_bingo = df_bingo.drop_duplicates(subset=['gene_id'])
        
    epgat_dfs[sp] = df_epgat
    bingo_dfs[sp] = df_bingo
    
    e_epgat = len(df_epgat[df_epgat['label'] == 1]) if not df_epgat.empty else 0
    ne_epgat = len(df_epgat[df_epgat['label'] == 0]) if not df_epgat.empty else 0
    
    e_bingo = len(df_bingo[df_bingo['label'] == 1]) if not df_bingo.empty else 0
    ne_bingo = len(df_bingo[df_bingo['label'] == 0]) if not df_bingo.empty else 0
    
    genes_epgat = set(df_epgat['gene_id']) if not df_epgat.empty else set()
    genes_bingo = set(df_bingo['gene_id']) if not df_bingo.empty else set()
    
    common_genes = genes_epgat.intersection(genes_bingo)
    only_epgat = genes_epgat - genes_bingo
    only_bingo = genes_bingo - genes_epgat
    
    conflicts = 0
    if not df_epgat.empty and not df_bingo.empty:
        merged = df_epgat.merge(df_bingo, on='gene_id', suffixes=('_epgat', '_bingo'))
        conflicts = len(merged[merged['label_epgat'] != merged['label_bingo']])
    
    audit_records.append({
        'species': sp,
        'epgat_essential': e_epgat,
        'epgat_nonessential': ne_epgat,
        'bingo_essential': e_bingo,
        'bingo_nonessential': ne_bingo,
        'common_genes': len(common_genes),
        'only_epgat': len(only_epgat),
        'only_bingo': len(only_bingo),
        'conflicts': conflicts
    })

df_audit = pd.DataFrame(audit_records)
df_audit.to_csv(AUDIT_DIR / 'epgat_vs_bingo_comparison.tsv', sep='\t', index=False)
with open(AUDIT_DIR / 'epgat_vs_bingo_comparison.md', 'w') as f:
    f.write("# EPGAT vs Bingo Essential Gene Data Comparison\n\n")
    f.write(df_audit.to_markdown(index=False))

# Step 5 & 6 & 7: Standardize, Rebuild, Split
print("Rebuilding standard datasets...")
from sklearn.model_selection import train_test_split

manifest_records = []
recommendations = []

for sp in SPECIES_LIST:
    df_epgat = epgat_dfs.get(sp, pd.DataFrame())
    df_bingo = bingo_dfs.get(sp, pd.DataFrame())
    
    if df_bingo.empty and df_epgat.empty:
        recommendations.append({'species': sp, 'recommended_version': 'None', 'reason': 'No data found'})
        continue
    
    if not df_bingo.empty and df_epgat.empty:
        df_final = df_bingo.copy()
        df_final['label_binary'] = df_final['label']
        df_final['label_source_project'] = 'bingo'
        df_final['label_source_file'] = df_final['source']
        df_final['evidence_note'] = ''
        rec = 'Bingo'
        reason = 'Only Bingo data available'
    elif df_bingo.empty and not df_epgat.empty:
        df_final = df_epgat.copy()
        df_final['label_binary'] = df_final['label']
        df_final['label_source_project'] = 'epgat'
        df_final['label_source_file'] = df_final['source']
        df_final['evidence_note'] = ''
        rec = 'EPGAT'
        reason = 'Only EPGAT data available'
    else:
        genes_epgat = set(df_epgat['gene_id'])
        genes_bingo = set(df_bingo['gene_id'])
        common_genes = genes_epgat.intersection(genes_bingo)
        
        if len(common_genes) == 0:
            df_final = df_bingo.copy()
            df_final['label_binary'] = df_final['label']
            df_final['label_source_project'] = 'bingo'
            df_final['label_source_file'] = df_final['source']
            df_final['evidence_note'] = 'EPGAT used a different ID namespace. Using Bingo as primary.'
            rec = 'Bingo'
            reason = 'ID namespace mismatch between EPGAT and Bingo (0 overlap). Defaulted to Bingo.'
        else:
            df_epgat = df_epgat.rename(columns={'label': 'label_epgat', 'source': 'source_epgat'})
            df_bingo = df_bingo.rename(columns={'label': 'label_bingo', 'source': 'source_bingo'})
            
            merged = pd.merge(df_epgat, df_bingo, on='gene_id', how='outer')
            
            def resolve_label(row):
                if pd.notna(row['label_bingo']) and pd.notna(row['label_epgat']):
                    if row['label_bingo'] == row['label_epgat']:
                        return pd.Series([row['label_bingo'], 'merged', f"bingo:{row['source_bingo']}|epgat:{row['source_epgat']}", ""])
                    else:
                        return pd.Series([row['label_bingo'], 'merged_conflict_bingo_wins', f"bingo:{row['source_bingo']}", f"Conflict: EPGAT={row['label_epgat']}, Bingo={row['label_bingo']}"])
                elif pd.notna(row['label_bingo']):
                    return pd.Series([row['label_bingo'], 'bingo', row['source_bingo'], ""])
                else:
                    return pd.Series([row['label_epgat'], 'epgat', row['source_epgat'], ""])
                    
            merged[['label_binary', 'label_source_project', 'label_source_file', 'evidence_note']] = merged.apply(resolve_label, axis=1)
            df_final = merged[['gene_id', 'label_binary', 'label_source_project', 'label_source_file', 'evidence_note']].copy()
            
            rec = 'Bingo + EPGAT (Bingo priority)'
            reason = 'Merged both sources. Prioritized Bingo in case of conflict due to newer project state.'
    
    recommendations.append({'species': sp, 'recommended_version': rec, 'reason': reason})
    
    # Save standardized labels
    sp_dir = PROCESSED_DIR / sp
    sp_dir.mkdir(exist_ok=True)
    
    df_final['label_text'] = df_final['label_binary'].map({1.0: 'essential', 0.0: 'non-essential', 1: 'essential', 0: 'non-essential'})
    df_final['gene_id_source'] = 'Ensembl/Original'
    df_final['gene_symbol'] = df_final['gene_id'] # Fallback
    df_final['species'] = sp
    df_final['included_in_final'] = True
    df_final['exclusion_reason'] = ""
    
    df_final['label_binary'] = df_final['label_binary'].astype(int)
    
    cols = ['gene_id', 'gene_id_source', 'gene_symbol', 'species', 'label_binary', 'label_text', 'label_source_project', 'label_source_file', 'evidence_note', 'included_in_final', 'exclusion_reason']
    for c in cols:
        if c not in df_final.columns: df_final[c] = ""
    df_final = df_final[cols]
    
    out_labels = sp_dir / 'labels.standard.tsv'
    df_final.to_csv(out_labels, sep='\t', index=False)
    
    # Check for original splits - simplified rule: Rebuild splits to be uniform
    # 70/15/15 stratify
    
    train_ids, test_val_ids, y_train, y_test_val = train_test_split(
        df_final['gene_id'], df_final['label_binary'], test_size=0.3, random_state=42, stratify=df_final['label_binary']
    )
    val_ids, test_ids = train_test_split(
        test_val_ids, test_size=0.5, random_state=42, stratify=y_test_val
    )
    
    df_splits = pd.DataFrame({'gene_id': df_final['gene_id']})
    df_splits['split'] = 'infer'
    df_splits.loc[df_splits['gene_id'].isin(train_ids), 'split'] = 'train'
    df_splits.loc[df_splits['gene_id'].isin(val_ids), 'split'] = 'val'
    df_splits.loc[df_splits['gene_id'].isin(test_ids), 'split'] = 'test'
    df_splits['split_source'] = 'ProGATE_v2_stratified_rebuild'
    df_splits['split_strategy'] = '70/15/15'
    df_splits['random_seed'] = 42
    
    out_splits = sp_dir / 'splits.standard.tsv'
    df_splits.to_csv(out_splits, sep='\t', index=False)
    
    # Summary
    e_count = len(df_final[df_final['label_binary'] == 1])
    ne_count = len(df_final[df_final['label_binary'] == 0])
    tr_c = len(train_ids)
    va_c = len(val_ids)
    te_c = len(test_ids)
    
    df_summary = pd.DataFrame([{
        'species': sp,
        'final_essential_count': e_count,
        'final_nonessential_count': ne_count,
        'final_total_labeled': e_count + ne_count,
        'split_train_count': tr_c,
        'split_val_count': va_c,
        'split_test_count': te_c,
        'notes': f"Rebuilt from {rec}"
    }])
    out_summary = sp_dir / 'summary.tsv'
    df_summary.to_csv(out_summary, sep='\t', index=False)
    
    manifest_records.append({
        'species': sp,
        'labels_path': str(out_labels.relative_to(PROGATE_DIR)),
        'splits_path': str(out_splits.relative_to(PROGATE_DIR)),
        'summary_path': str(out_summary.relative_to(PROGATE_DIR)),
        'final_essential_count': e_count,
        'final_nonessential_count': ne_count,
        'source_strategy': rec,
        'notes': reason
    })

df_rec = pd.DataFrame(recommendations)
df_rec.to_csv(AUDIT_DIR / 'final_recommendation.tsv', sep='\t', index=False)

with open(AUDIT_DIR / 'final_recommendation.md', 'w') as f:
    f.write("# Final Dataset Recommendations\n\n")
    f.write(df_rec.to_markdown(index=False))

df_manifest = pd.DataFrame(manifest_records)
df_manifest.to_csv(MANIFEST_DIR / 'essential_gene_dataset_manifest.tsv', sep='\t', index=False)

print("Finished generating standardized datasets.")
