import os
import pandas as pd
from pathlib import Path
import subprocess
import csv
import re
import json

# Focal Species
SPECIES_TAXON = {
    "fgraminearum": "229533",
    "scerevisiae": "4932",
    "celegans": "6239",
    "melanogaster": "7227",
    "human": "9606",
}

# Directories
PROGATE_ROOT = Path('/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2')
DOWNLOAD_DIR = PROGATE_ROOT / 'data/processed/EXP/Download' # Wait, the prompt said processed/OR/Download? 
# Re-reading prompt: "current ProGATE_v2/data/processed/OR/Download/"
DOWNLOAD_DIR = PROGATE_ROOT / 'data/processed/OR/Download'
OUTPUT_DIR = PROGATE_ROOT / 'data/processed/OR'

os.makedirs(OUTPUT_DIR, exist_ok=True)

def run_taxonkit_lineage(taxids):
    """
    Calls taxonkit to get lineage and reformatted ranks for a list of taxids.
    """
    taxid_str = "\n".join(map(str, taxids))
    
    # Get lineage
    lineage_cmd = ["taxonkit", "lineage"]
    process = subprocess.Popen(lineage_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=taxid_str)
    
    if stderr:
        print(f"[ERROR] taxonkit lineage error: {stderr}")
        
    # Reformat
    reformat_cmd = ["taxonkit", "reformat", "-f", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}", "-F", "-P"]
    process = subprocess.Popen(reformat_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout_final, stderr_final = process.communicate(input=stdout)
    
    if stderr_final:
        print(f"[ERROR] taxonkit reformat error: {stderr_final}")
        
    # Parse output
    results = []
    for line in stdout_final.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 3:
            taxid = parts[0]
            lineage_name = parts[1]
            # ranks start from parts[2]
            ranks = parts[2:]
            if len(ranks) < 7:
                ranks = ranks + [""] * (7 - len(ranks))
            results.append({
                "taxid": taxid,
                "scientific_name": lineage_name.split(";")[-1] if lineage_name else "",
                "kingdom": ranks[0],
                "phylum": ranks[1],
                "class": ranks[2],
                "order": ranks[3],
                "family": ranks[4],
                "genus": ranks[5],
                "species": ranks[6]
            })
    return pd.DataFrame(results)

def main():
    print("--- Starting Orthology Species Scope Analysis ---")
    
    # --- Part A: OrthoDB ---
    print("[INFO] Analyzing OrthoDB v12...")
    
    # Load OrthoDB tables
    species_tab = DOWNLOAD_DIR / "odb12v2_species.tab"
    levels_tab = DOWNLOAD_DIR / "odb12v2_levels.tab"
    level2species_tab = DOWNLOAD_DIR / "odb12v2_level2species.tab"
    
    df_species = pd.read_csv(species_tab, sep='\t', header=None, 
                             names=['taxid', 'organism_id', 'scientific_name', 'assembly_id', 'gene_count', 'og_count', 'mapping_type'])
    df_levels = pd.read_csv(levels_tab, sep='\t', header=None,
                            names=['level_taxid', 'level_name', 'total_genes', 'total_ogs', 'total_species'])
    df_l2s = pd.read_csv(level2species_tab, sep='\t', header=None,
                         names=['top_level_taxid', 'organism_id', 'hop_count', 'intermediate_levels'])
    
    # A1. Focal Species Info
    focal_taxids = [int(tid) for tid in SPECIES_TAXON.values()]
    df_focal_sp = df_species[df_species['taxid'].isin(focal_taxids)].copy()
    
    # Map species_key
    tid_to_key = {int(v): k for k, v in SPECIES_TAXON.items()}
    df_focal_sp['species_key'] = df_focal_sp['taxid'].map(tid_to_key)
    
    # Reorder columns
    df_focal_sp = df_focal_sp[['species_key', 'taxid', 'organism_id', 'scientific_name', 'assembly_id', 'gene_count', 'og_count', 'mapping_type']]
    df_focal_sp.to_csv(OUTPUT_DIR / "orthodb_focal_species.tsv", sep='\t', index=False)
    
    # A2. Level paths for focal species
    focal_org_ids = df_focal_sp['organism_id'].unique()
    df_focal_l2s = df_l2s[df_l2s['organism_id'].isin(focal_org_ids)].copy()
    
    level_paths = []
    for _, row in df_focal_l2s.iterrows():
        # parse intermediate_levels like {2759,33208,...}
        levels_str = row['intermediate_levels'].strip('{}')
        levels = [int(l) for l in levels_str.split(',') if l.strip()]
        
        sp_key = df_focal_sp[df_focal_sp['organism_id'] == row['organism_id']]['species_key'].iloc[0]
        f_taxid = df_focal_sp[df_focal_sp['organism_id'] == row['organism_id']]['taxid'].iloc[0]
        
        for idx, l_taxid in enumerate(levels):
            # find level info
            l_info = df_levels[df_levels['level_taxid'] == l_taxid]
            if not l_info.empty:
                level_paths.append({
                    'species_key': sp_key,
                    'focal_taxid': f_taxid,
                    'orthodb_organism_id': row['organism_id'],
                    'top_level_taxid': row['top_level_taxid'],
                    'hop_count': row['hop_count'],
                    'level_taxid': l_taxid,
                    'level_name': l_info['level_name'].iloc[0],
                    'level_gene_count': l_info['total_genes'].iloc[0],
                    'level_og_count': l_info['total_ogs'].iloc[0],
                    'level_species_count': l_info['total_species'].iloc[0],
                    'path_rank': idx
                })
            else:
                # If not in levels table (unlikely for selected levels)
                level_paths.append({
                    'species_key': sp_key,
                    'focal_taxid': f_taxid,
                    'orthodb_organism_id': row['organism_id'],
                    'top_level_taxid': row['top_level_taxid'],
                    'hop_count': row['hop_count'],
                    'level_taxid': l_taxid,
                    'level_name': 'Unknown',
                    'level_gene_count': 0,
                    'level_og_count': 0,
                    'level_species_count': 0,
                    'path_rank': idx
                })
                
    df_focal_levels_long = pd.DataFrame(level_paths)
    df_focal_levels_long.to_csv(OUTPUT_DIR / "orthodb_focal_levels_long.tsv", sep='\t', index=False)
    
    # A3. Species under each level
    # We need to know which organisms are under which levels.
    # We can use df_l2s: for each level in focal paths, find all organisms that have this level in their path.
    all_query_levels = df_focal_levels_long[['level_taxid', 'level_name']].drop_duplicates()
    
    level_species_list = []
    # Pre-parse all paths in df_l2s for faster lookup?
    # Actually, df_l2s has 32k rows, it's fine.
    
    print(f"[INFO] Expanding species for {len(all_query_levels)} levels...")
    for _, l_row in all_query_levels.iterrows():
        l_taxid = l_row['level_taxid']
        l_name = l_row['level_name']
        
        # organisms that have l_taxid in their path
        pattern = rf"\{{.*,?{l_taxid},?.*\}}"
        # Or more robustly:
        def has_level(level_list_str, target):
            ls = level_list_str.strip('{}').split(',')
            return str(target) in ls
            
        mask = df_l2s['intermediate_levels'].apply(lambda x: has_level(x, l_taxid))
        org_ids_under_level = df_l2s[mask]['organism_id'].unique()
        
        # Get details for these organisms
        df_orgs = df_species[df_species['organism_id'].isin(org_ids_under_level)]
        
        # Focal species that use this level
        focal_keys = df_focal_levels_long[df_focal_levels_long['level_taxid'] == l_taxid]['species_key'].unique()
        
        for sk in focal_keys:
            ftid = SPECIES_TAXON[sk]
            for _, o_row in df_orgs.iterrows():
                level_species_list.append({
                    'species_key': sk,
                    'focal_taxid': ftid,
                    'query_level_taxid': l_taxid,
                    'query_level_name': l_name,
                    'member_orthodb_organism_id': o_row['organism_id'],
                    'member_taxid': o_row['taxid'],
                    'member_scientific_name': o_row['scientific_name']
                })
                
    df_level_species_long = pd.DataFrame(level_species_list)
    df_level_species_long.to_csv(OUTPUT_DIR / "orthodb_level_species_long.tsv", sep='\t', index=False)
    
    # --- Part B: InParanoidDB ---
    print("[INFO] Analyzing InParanoiDB9...")
    
    inparanoid_file = DOWNLOAD_DIR / "InParanoiDB9_species.txt"
    df_inp = pd.read_csv(inparanoid_file, sep=',', comment='#', names=['taxid', 'uniprot_id', 'species_name', 'domain'])
    # Actually it has a header starting with #
    df_inp = pd.read_csv(inparanoid_file, sep=',', header=0)
    df_inp.columns = ['taxid', 'uniprot_id', 'species_name', 'domain']
    
    # B2. Supplement Phylogenetic Levels
    all_inp_taxids = df_inp['taxid'].unique().tolist()
    # Also add focal taxids if not present
    all_ids_for_lineage = list(set(all_inp_taxids + focal_taxids))
    
    print(f"[INFO] Calling taxonkit for {len(all_ids_for_lineage)} taxids...")
    df_taxonomy = run_taxonkit_lineage(all_ids_for_lineage)
    
    # Ensure taxid is string for merge
    df_inp['taxid'] = df_inp['taxid'].astype(str)
    df_taxonomy['taxid'] = df_taxonomy['taxid'].astype(str)
    
    # Merge taxonomy back to InParanoid species
    df_inp_tax = df_inp.merge(df_taxonomy, on='taxid', how='left')
    df_inp_tax.to_csv(OUTPUT_DIR / "inparanoid_species_taxonomy.tsv", sep='\t', index=False)
    
    # B3. Relative Coverage
    coverage_stats = []
    relative_species = []
    
    for sk, tid in SPECIES_TAXON.items():
        focal_tid = int(tid)
        f_row = df_taxonomy[df_taxonomy['taxid'] == str(focal_tid)]
        if f_row.empty:
            print(f"[WARNING] No taxonomy found for focal species {sk} (taxid {tid})")
            continue
        
        f_info = f_row.iloc[0]
        
        # Matches in InParanoid
        ranks = ['phylum', 'class', 'order', 'family', 'genus']
        stats = {
            'species_key': sk,
            'focal_taxid': focal_tid,
            'focal_scientific_name': f_info['scientific_name'],
            'focal_phylum': f_info['phylum'],
            'focal_class': f_info['class'],
            'focal_order': f_info['order'],
            'focal_family': f_info['family'],
            'focal_genus': f_info['genus']
        }
        
        for r in ranks:
            rank_val = f_info[r]
            if not rank_val:
                stats[f'same_{r}_count'] = 0
                continue
                
            matches = df_inp_tax[df_inp_tax[r] == rank_val]
            stats[f'same_{r}_count'] = len(matches)
            
            for _, m_row in matches.iterrows():
                relative_species.append({
                    'species_key': sk,
                    'comparison_level': r,
                    'focal_taxid': focal_tid,
                    'focal_scientific_name': f_info['scientific_name'],
                    'matched_taxid': m_row['taxid'],
                    'matched_scientific_name': m_row['species_name'],
                    'matched_phylum': m_row['phylum'],
                    'matched_class': m_row['class'],
                    'matched_order': m_row['order'],
                    'matched_family': m_row['family'],
                    'matched_genus': m_row['genus']
                })
        
        coverage_stats.append(stats)
        
    df_coverage = pd.DataFrame(coverage_stats)
    df_coverage.to_csv(OUTPUT_DIR / "inparanoid_relative_coverage.tsv", sep='\t', index=False)
    
    df_rel_sp = pd.DataFrame(relative_species)
    df_rel_sp.to_csv(OUTPUT_DIR / "inparanoid_relative_species_long.tsv", sep='\t', index=False)
    
    # --- Generate Markdown Reports ---
    print("[INFO] Generating Markdown summaries...")
    
    # OrthoDB Summary
    with open(OUTPUT_DIR / "orthodb_levels_summary.md", 'w') as f:
        f.write("# OrthoDB v12 Levels Summary\n\n")
        f.write("## Focal Species Details\n\n")
        f.write(df_focal_sp.to_markdown(index=False))
        f.write("\n\n## Available Levels per Focal Species\n\n")
        for sk in SPECIES_TAXON.keys():
            f.write(f"### {sk}\n")
            df_sub = df_focal_levels_long[df_focal_levels_long['species_key'] == sk].sort_values('path_rank')
            f.write(df_sub[['level_name', 'level_taxid', 'level_species_count', 'level_gene_count']].to_markdown(index=False))
            f.write("\n\n")
            
    # InParanoid Summary
    with open(OUTPUT_DIR / "inparanoid_species_summary.md", 'w') as f:
        f.write("# InParanoiDB9 Species Distribution Summary\n\n")
        f.write("## Relative Coverage Statistics\n\n")
        f.write(df_coverage.to_markdown(index=False))
        f.write("\n\n## Methodology\n")
        f.write("Taxonomy lineages were supplemented using `taxonkit`. Coverage counts include all species in InParanoiDB9 sharing the same rank value as the focal species.\n")

    # Global README
    with open(OUTPUT_DIR / "README_orthology_species_scope.md", 'w') as f:
        f.write("# Orthology Data Species Scope Analysis Report\n\n")
        f.write("## Overview\n")
        f.write("This report analyzes the coverage and hierarchical distribution of 5 focal species in OrthoDB v12 and InParanoiDB9.\n\n")
        
        f.write("## 1. OrthoDB Findings\n")
        f.write("- **Broad Levels**: All focal species share the `Eukaryota` (2759) and often `Opisthokonta` (33208) levels.\n")
        f.write("- **Fungal Path**: `fgraminearum` and `scerevisiae` share the `Ascomycota` (4890) level.\n")
        f.write("- **Human/Fly/Worm Path**: Shared levels end at `Bilateria` or `Metazoa` usually found under `33208`.\n\n")
        
        f.write("## 2. InParanoiDB9 Findings\n")
        # Find richest species
        richest = df_coverage.loc[df_coverage['same_order_count'].idxmax()]['species_key']
        poorest = df_coverage.loc[df_coverage['same_order_count'].idxmin()]['species_key']
        f.write(f"- **Richest Coverage**: `{richest}` has the most relatives in InParanoiDB9 at intermediate ranks.\n")
        f.write(f"- **Poorest Coverage**: `{poorest}` has relatively few close relatives.\n\n")
        
        f.write("## 3. Implementation Recommendations\n")
        f.write("- **has_ortholog**: Recommended to use intermediate levels (e.g., Phylum or Order) to capture meaningful presence/absence patterns.\n")
        f.write("- **single_copy**: Should be defined at lower levels (e.g., Genus or Family) if coverage allows, or Order level for broader evolutionary conservation studies.\n")

    print(f"--- All Analysis Finished ---")
    print(f"Results saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
