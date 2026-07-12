# Subcellular Localization Data Build Report (Binary Matrix Version)

## Overview
This directory contains standardized subcellular localization data for 5 species in binary matrix format (0/1).

## Global Statistics

| species      | source_db    |   final_location_column_count |   mapped_label_gene_count |   label_gene_coverage_ratio |
|:-------------|:-------------|------------------------------:|--------------------------:|----------------------------:|
| fgraminearum | eFG          |                            12 |                     11508 |                    0.813516 |
| scerevisiae  | COMPARTMENTS |                            13 |                      5504 |                    0.977967 |
| celegans     | COMPARTMENTS |                            13 |                      5616 |                    0.410466 |
| melanogaster | COMPARTMENTS |                            13 |                      4304 |                    0.658104 |
| human        | COMPARTMENTS |                            13 |                     16602 |                    0.899447 |

## Methodology
1. **Format**: Output `subloc.csv` is a binary matrix where rows are Genes and columns are standard Location categories.
2. **Normalization**: Raw terms are mapped to a controlled vocabulary (e.g., 'cytosol' -> 'Cytoplasm'). Synonymous terms are merged.
3. **Column Selection**: For each species, the top 12–15 standard categories (by gene count) are selected as columns.
4. **ID Mapping**: Aligned to `labels.standard.tsv` using direct matches and STRING aliases as a bridge.
5. **Output Columns**: Columns are sorted by frequency (most frequent first).

## Controlled Vocabulary (Selection Pool)
- Cytoplasm
- Nucleus
- Cell membrane
- Multi-pass membrane
- Mitochondrion
- Endoplasmic reticulum
- Golgi
- Endosome
- Vacuole
- Peroxisome
- Secreted
- Lysosome
- Cytoskeleton
- Ribosome
- Extracellular matrix

