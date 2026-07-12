# Orthology Data Species Scope Analysis Report

## Overview
This report analyzes the coverage and hierarchical distribution of 5 focal species in OrthoDB v12 and InParanoiDB9.

## 1. OrthoDB Findings
- **Broad Levels**: All focal species share the `Eukaryota` (2759) and often `Opisthokonta` (33208) levels.
- **Fungal Path**: `fgraminearum` and `scerevisiae` share the `Ascomycota` (4890) level.
- **Human/Fly/Worm Path**: Shared levels end at `Bilateria` or `Metazoa` usually found under `33208`.

## 2. InParanoiDB9 Findings
- **Richest Coverage**: `scerevisiae` has the most relatives in InParanoiDB9 at intermediate ranks.
- **Poorest Coverage**: `fgraminearum` has relatively few close relatives.

## 3. Implementation Recommendations
- **has_ortholog**: Recommended to use intermediate levels (e.g., Phylum or Order) to capture meaningful presence/absence patterns.
- **single_copy**: Should be defined at lower levels (e.g., Genus or Family) if coverage allows, or Order level for broader evolutionary conservation studies.
