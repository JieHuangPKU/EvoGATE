# InParanoid v8 Ortholog Matrix Build

## Scope

This workflow downloads pairwise InParanoid v8 `*.orthoXML` files for 5 focal species and builds EPGAT-style orthology binary matrices aligned to each species `labels.standard.tsv` final `gene_id` set.

## Download Directories And File Counts

| species | folder | downloaded_file_count | target_species_count | download_dir |
| --- | --- | ---: | ---: | --- |
| fgraminearum | G.zeae | 167 | 167 | `/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/OR/Download/InParanoid8/G.zeae` |
| celegans | C.elegans | 233 | 233 | `/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/OR/Download/InParanoid8/C.elegans` |
| scerevisiae | S.cerevisiae | 62 | 62 | `/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/OR/Download/InParanoid8/S.cerevisiae` |
| melanogaster | D.melanogaster | 203 | 203 | `/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/OR/Download/InParanoid8/D.melanogaster` |
| human | H.sapiens | 162 | 162 | `/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2/data/processed/OR/Download/InParanoid8/H.sapiens` |

## XML Parsing

Each orthoXML file is parsed with `xml.etree.ElementTree`. The parser first reads every `<species>/<gene>` entry to map internal `geneRef id` values to `{species name, taxid, geneId, protId}`. It then walks each top-level `orthologGroup`, collects all nested `geneRef` members, and marks a focal raw identifier as ortholog-positive for a target species if they co-occur in the same ortholog group.

## Focal Raw ID Mapping

Mapping priority is: direct match to `labels.standard.tsv` `gene_id`, exact match against available labels auxiliary ID-like columns, exact match through STRING aliases, then a Fusarium-specific conservative transcript-to-gene conversion when aliases expose `FGRAMPH1_*T*` transcript IDs and the corresponding `FGRAMPH1_*G*` gene exists in labels.

Ambiguous mappings are retained in `orthologs_id_mapping.tsv` as `mapping_status=ambiguous` and are excluded from final matrices. Unmapped raw IDs are also excluded from final matrices.

## STRING Alias Usage

STRING alias files under `data/processed/PPI/stringDB/` were used as auxiliary evidence for focal-side ID harmonization. The workflow does not use old EPGAT, Bingo orthology outputs, or `data_registry` files.

## Matrix Definition

For each focal species, `orthologs.csv` uses `Gene` as the first column, rows equal to all final label `gene_id` values, columns equal to target support species observed in downloaded orthoXML files, and 0/1 values indicating whether at least one ortholog exists in that target species. Column names prefer target taxid parsed from XML; species names are only used as fallback when taxid cannot be resolved.

## High-Unmapped Species

- No focal species had unmapped raw IDs exceeding mapped raw IDs in the final audit.

