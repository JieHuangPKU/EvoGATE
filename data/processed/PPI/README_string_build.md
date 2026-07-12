# STRING PPI Dataset Build Report (v12.0)

## Global Statistics

| species      |   taxon_id |   kept_edge_count |   final_unique_node_count |   label_gene_coverage_ratio |
|:-------------|-----------:|------------------:|--------------------------:|----------------------------:|
| fgraminearum |     229533 |           1439790 |                     11122 |                    0.786229 |
| scerevisiae  |       4932 |           1281599 |                      5517 |                    0.980277 |
| celegans     |       6239 |           1937760 |                     12703 |                    0.928446 |
| melanogaster |       7227 |            517640 |                      6497 |                    0.993425 |
| human        |       9606 |           6321712 |                     18276 |                    0.99014  |

## Methodology
1. **Mapping Source**: Strictly used STRING `.aliases.v12.0.txt.gz` files.
2. **Fusarium Rule**: For `fgraminearum`, used `FgraminearumPH-1.tsv` to map transcript aliases to gene IDs.
3. **Filtering**: Only edges where both nodes map to IDs in `labels.standard.tsv` are kept.
4. **Aggregation**: Duplicate edges are collapsed using the **MAX** score.
5. **Format**: Output follows the EPGAT `string.csv` format.
