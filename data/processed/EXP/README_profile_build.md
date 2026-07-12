# Expression Profile Data Build Report (Final Integrated Version)

## Global Statistics

| species      | input_source                                    |   exp_group_count |   final_gene_count_in_profile |   label_gene_coverage_ratio |
|:-------------|:------------------------------------------------|------------------:|------------------------------:|----------------------------:|
| celegans     | celegans_GSE31422-GPL14144_series_matrix.txt.gz |                30 |                         13682 |                    0.917775 |
| fgraminearum | fgraminearum_GSE292521_Fg_SW_counts.txt         |                60 |                         14146 |                    0.933338 |
| human        | human_GSE86354_GTEx_FPKM.csv                    |                64 |                         18458 |                    0.967873 |
| melanogaster | melanogaster_GSE67547_series_matrix.txt.gz      |                66 |                          6540 |                    0.902752 |
| scerevisiae  | scerevisiae_GSE3431_series_matrix.txt.gz        |                36 |                          5628 |                    0.969083 |

## Methodology
1. **GPL Integration**: Parsed `GSE*_family.soft.gz` to extract platform annotations (GPL) for microarray datasets (celegans, melanogaster, scerevisiae).
2. **TPM Conversion**: For `fgraminearum`, replaced old microarray source with RNA-seq counts (`fgraminearum_GSE292521_Fg_SW_counts.txt`) and converted to TPM using reference transcript lengths.
3. **Heuristic Mapping**: Aligned IDs to standard Gene IDs using ORF/Symbol/Aliases/WormBase.
4. **Aggregation**: Multiple probes/records mapping to the same gene are averaged (**MEAN**).
5. **Completeness**: All genes in `labels.standard.tsv` are included in `profile.csv`. Missing data points are filled with **0**.
6. **Independence**: Process is self-contained and aligns strictly to `labels.standard.tsv`.

