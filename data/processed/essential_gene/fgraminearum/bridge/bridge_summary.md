# Fusarium Newlabel Protein-to-Canonical Bridge Summary

## Purpose
This bridge protocolizes the previously implicit XP-to-canonical mapping step required to convert the PH-1 yeast-transfer table into the final canonical Fusarium gene ID space.

## Inputs
- repo-local proteome manifest: `data/derived_labels/proteome_manifest.tsv`
- PH-1 NCBI protein FASTA: `/data276/jiehuang/fungi/Fusarium/fusarium-protein/ncbi_dataset/data/GCF_000240135.3/protein.faa`
- PH-1 legacy 229533 proteome: `/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/Orthologs/229533.fasta`
- legacy accession-to-canonical table: `/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/PPI/mapping.tab`
- unified XP/FGSG-to-PH1 map: `/data276/jiehuang/fungi/Fusarium/Evidence/00_idmap/FG_gene_id_unified_map.tsv`
- repo-local supplemental FGSG-to-canonical evidence: `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`, `data/processed/PPI/fgraminearum/string_id_mapping.tsv`, `data/processed/EXP/fgraminearum/exp_id_mapping.tsv`, `data/processed/LC/fgraminearum/subloc_id_mapping.tsv`

## Bridge Logic
The primary path is now the direct `XP_*` to `FGRAMPH1_*` mapping preserved in the shared unified ID map.
That route is preferred when it provides a unique canonical PH-1 gene because it directly records the NCBI protein accession lineage used by the historical Fusarium orthology workflow.
Sequence-based matching against the legacy `229533` proteome remains a second protocolized evidence path.
When the NCBI header itself contains an `FGSG_*` locus tag, a repo-local FGSG-to-canonical bridge is also evaluated and recorded as supporting evidence.

## Output Space
The final canonical ID space is `fgraminearum::FGRAMPH1_*`.

## Current Coverage
- total PH-1 NCBI proteins parsed: 13312
- resolved protein-to-canonical mappings: 13092
- unresolved mappings: 220
- ambiguous mappings: 0
- high-confidence yeast-transfer XP rows: 1056
- protocolized high-confidence transfer genes after canonical deduplication: 1045
- unresolved high-confidence transfer rows: 8
- high-confidence transfer rows resolved via the unified XP/FGSG mapping path: 111

## Caveats
This bridge removes the direct dependence on `positive_set_P1.tsv` by making the XP-to-canonical mapping explicit.
However, not every high-confidence XP protein currently resolves through the repo-local bridge evidence. Any unresolved XP rows are preserved in the audit outputs and are excluded from the protocolized rebuilt newlabel positive component until additional annotation evidence is protocolized.
