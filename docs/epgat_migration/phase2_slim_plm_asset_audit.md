# Phase 2 Slim PLM Asset Audit

## ESM2 assets

- Source format: HDF5 (`.h5`).
- Key space: legacy node keys stored as top-level HDF5 dataset names.
- Human key type: UniProt-like legacy gene/protein IDs used directly by the legacy replay node manifests.
- Celegans key type: legacy gene IDs used directly by the legacy replay node manifests.
- Fusarium key type: raw `FGRAMPH1_*` legacy gene IDs; mapped into canonical space by stripping the `fgraminearum::` prefix from canonical gene IDs.
- Embedding dimension: 1280.

## ProtT5 assets

- Source format: HDF5 (`.h5`).
- Key space: same as ESM2 for each species.
- Embedding dimension: 1024.

## Coverage by species

- human: ESM2 `17874/18822 = 0.9496`, ProtT5 `17874/18822 = 0.9496`.
- celegans: ESM2 `5423/5766 = 0.9405`, ProtT5 `5423/5766 = 0.9405`.
- fgraminearum_canonical: ESM2 `14143/14145 = 0.9999`, ProtT5 `14143/14145 = 0.9999`.

## Recommended loader paths and mapping keys

- human ESM2: `/home/jiehuang/software/fungi/EPGAT/scripts/embeddings/human_esm2_embeddings.h5`, key=`legacy_gene_id`.
- human ProtT5: `/home/jiehuang/software/fungi/EPGAT/scripts/embeddings/human_prott5_embeddings.h5`, key=`legacy_gene_id`.
- celegans ESM2: `/home/jiehuang/software/fungi/EPGAT/scripts/embeddings/celegans_esm2_embeddings.h5`, key=`legacy_gene_id`.
- celegans ProtT5: `/home/jiehuang/software/fungi/EPGAT/scripts/embeddings/celegans_prott5_embeddings.h5`, key=`legacy_gene_id`.
- fgraminearum ESM2: `/home/jiehuang/software/fungi/EPGAT/scripts/embeddings/fgraminearum_esm2_embeddings.h5`, key=`legacy_fgraminearum_gene_id`.
- fgraminearum ProtT5: `/home/jiehuang/software/fungi/EPGAT/scripts/embeddings/fgraminearum_prott5_embeddings.h5`, key=`legacy_fgraminearum_gene_id`.

## Availability conclusion

- All six slim Phase 2 assets required for this round already exist; no new embedding generation was needed.
- The main mapping difference is only on `fgraminearum_canonical`, where canonical node IDs must be projected back to raw `FGRAMPH1_*` keys for HDF5 lookup.