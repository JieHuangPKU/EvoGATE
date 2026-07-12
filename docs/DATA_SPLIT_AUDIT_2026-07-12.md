# EvoGATE Data Split Audit — 2026-07-12

_Focus: Homology-driven data leakage in train/val/test splits for essential gene prediction._

---

## 1. Split Method Identification

### 1.1 Primary split code

**File**: `src/data/freeze_unified_protocol.py`, function `assign_splits()` (line 184–212)

```python
train_val_idx, test_idx = train_test_split(
    working.index, test_size=test_fraction,
    random_state=seed, stratify=y,
)
```

**Method**: `sklearn.model_selection.train_test_split` with `stratify=y` (binary label only).

**Parameters**:
| Parameter | Value | Source |
|-----------|-------|--------|
| Split seed | `20260409` | `configs/frozen_protocol.yaml` |
| Test fraction | `0.20` | `configs/frozen_protocol.yaml` |
| Val fraction | `0.10` (after train split: `0.10 / (1 - 0.20) ≈ 0.125`) | `configs/frozen_protocol.yaml` |
| Train fraction | `0.70` (implicit) | Derived |
| Stratify | Binary label (0/1) | Hard-coded |

### 1.2 Label-scarcity splits

**File**: `src/data/build_label_scarcity_split_manifests.py`

- Uses the base frozen split as input, preserves the test set intact
- Subsamples train+val pool to target fractions (10%–90%)
- Seeds: `1029,1030,1031,1032,1033` (same as training seeds)
- Stratification: `stratify=pool_labels` — binary label only

### 1.3 Species scope

Each species (human, celegans, scerevisiae, dmelanogaster, fgraminearum) receives its **own independent split**. No cross-species split coordination exists.

---

## 2. Orthogroup (Homology) Leakage Analysis

### 2.1 Fusarium graminearum

**Data sources**:
- `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`: 13,152 PH-1 proteins → 12,573 unique orthogroups (OG*)
- `data/processed/.../bridge/high_confidence_yeast_transfer_candidates.tsv`: 1,045 high-confidence positive genes with orthogroup annotations (col 32: `canonical_gene_id`, col 2: `orthogroup_id`)
- `results/frozen_protocol/splits/fgraminearum_newlabel_split.tsv`: 11,965 genes with split assignments (col 5: `canonical_gene_id`, col 18: `split`)

**Key finding**: After the protein-to-canonical-gene bridge resolves multiple proteins mapping to the same gene:

| Metric | Count |
|--------|-------|
| Orthogroups with ≥1 canonical gene | 1,044 |
| Orthogroups with ≥2 canonical genes | **1** (OG0001039) |
| Orthogroups spanning train+test | **0** |
| Orthogroups spanning train+val | **1** (OG0001039: 2 genes, train+val) |

**Verdict**: Fusarium orthogroup leakage through the label transfer pipeline is **negligible**. The protein-to-gene bridge effectively deduplicates within-orthogroup redundancy.

### 2.2 Model organisms (human, C. elegans, S. cerevisiae, D. melanogaster)

**Data source**: `data/processed/essential_gene/<species>/labels.standard.tsv` — standard labels from OGEE/Bingo/EPGAT (exact provenance disputed, see INC-005).

| Species | Labeled genes | Positives | Positive rate | Split counts (train/val/test) |
|---------|-------------|-----------|---------------|-------------------------------|
| Human | 18,458 | 1,827 | 9.9% | 12,920 / 1,846 / 3,692 |
| S. cerevisiae | 5,628 | 1,050 | 18.7% | 3,939 / 563 / 1,126 |
| C. elegans | 13,682 | 578 | 4.2% | 9,576 / 1,369 / 2,737 |
| D. melanogaster | 6,540 | 402 | 6.1% | 4,578 / 654 / 1,308 |

**Key concern**: These species have extensive gene families (paralogs from duplication events). The split is **random stratified by label only** — no gene family, orthogroup, or sequence identity constraint.

**Untested risk**: Paralogous genes within model organisms may be split across train and test. This could inflate performance metrics, especially for:
- Human: large gene families (e.g., olfactory receptors, zinc finger proteins)
- S. cerevisiae: whole-genome duplication paralogs

**No orthogroup-level audit was performed for model organism data** because orthogroup annotations (OrthoFinder OG*) are available only for Fusarium (via the yeast transfer pipeline). Model organisms lack the equivalent orthogroup-to-gene mapping in processed data.

---

## 3. Sequence-Level Homology Leakage (ESM2)

### 3.1 Mechanism

ESM2 (`esm2_t33_650M_UR50D`, 650M parameters, embedding dim 1280) produces protein-level embeddings. Paralogous genes with similar amino acid sequences will produce similar ESM2 vectors.

If paralog A (train) and paralog B (test) have ≥60% sequence identity, their ESM2 embeddings will be highly correlated. The model can "memorize" the embedding pattern rather than learn generalizable features.

### 3.2 Absence of sequence deduplication

**No evidence found** of:
- CD-HIT (cluster at sequence identity threshold)
- MMseqs2 clustering
- BLAST-based homology reduction
- Pfam/InterPro domain-based grouping for split constraints

### 3.3 Literature context

Published work (e.g., essential gene prediction benchmarks) documents that splitting without sequence homology control can inflate AUROC by **7+ percentage points**. This is a known failure mode in the field.

### 3.4 Practical severity for EvoGATE

- **Fusarium**: 12,000 nodes, 1,097 positives. After protein-to-gene bridge, within-orthogroup paralogy is negligible. However, more distant paralogs (not in the same OrthoFinder OG) could still have sequence similarity through domain conservation.
- **Model organisms**: Higher risk due to larger gene families.
- **Feature-level**: The ORT features directly encode orthology information. For Fusarium, the ORT features indicate cross-species orthology (e.g., "this Fusarium gene has a yeast ortholog"), NOT within-Fusarium paralogy. So ORT features are NOT a direct leakage vector for within-species splits.

---

## 4. Random Seed and Replication Analysis

### 4.1 Split seed

- Single fixed split seed: `20260409`
- One train/val/test partition per species
- NOT repeated holdout / NOT multiple random splits

### 4.2 Training seeds

- 5 seeds: `1029`, `1030`, `1031`, `1032`, `1033`
- Applied per model configuration
- Results reported as mean ± standard deviation

### 4.3 What this measures

The 5 training seeds quantify **optimization variability** (random init, SGD noise) on a single fixed split. This is NOT:
- Independent biological replication
- Cross-split generalization assessment
- Estimate of split-induced variance

### 4.4 Missing statistical machinery

- No 95% confidence intervals
- No paired hypothesis tests (Wilcoxon signed-rank, paired t-test)
- No multiple comparison correction across models
- Std uses `ddof=0` (population std), inflating apparent precision

---

## 5. Leakage Risk Assessment

| Risk vector | Fusarium | Model organisms | Overall severity |
|-------------|----------|-----------------|------------------|
| Orthogroup label leakage (train↔test) | **LOW** — 0 OG crossing | **UNKNOWN** — not audited | Fusarium: safe; model orgs: investigate |
| ESM2 sequence similarity leakage | **MODERATE** — paralogs untested | **MODERATE** — gene families untested | Theoretical risk, literature-supported |
| ORT feature leakage | **LOW** — cross-species only | **LOW** — cross-species only | Not a within-species leakage vector |
| Label-stratified split without OG constraint | **MODERATE** — methodology limitation | **MODERATE** — methodology limitation | Standard practice but suboptimal for this domain |
| Single split (no repeated holdout) | **MODERATE** — no split-variance estimate | **MODERATE** — no split-variance estimate | Standard practice, but limits robustness claims |

---

## 6. Recommendations (for review before implementation)

### 6.1 P0 — Audit model organism orthogroup leakage

Use the existing InParanoid ortholog matrices (`data/processed/OR/<species>/orthologs.csv`) to identify within-species paralog pairs (genes with identical ORT target-species patterns) and check their split assignments.

```python
# Pseudocode for model organism audit
for species in [human, celegans, scerevisiae, dmelanogaster]:
    ort = pd.read_csv(f"data/processed/OR/{species}/orthologs.csv")
    split = pd.read_csv(f"results/frozen_protocol/splits/{species}_split.tsv", sep="\t")
    # Find gene pairs with identical ORT vectors
    # Check if any cross train/test boundary
```

### 6.2 P1 — Orthogroup-stratified splitting (future protocol version)

For `frozen_protocol_v2`:
1. Assign each gene to its primary orthogroup (from OrthoFinder or InParanoid)
2. Use `GroupShuffleSplit` or iterative constrained allocation to ensure all genes from the same orthogroup go to the SAME split
3. This prevents within-orthogroup leakage for both Fusarium and model organisms
4. Document this as an explicit protocol upgrade

### 6.3 P1 — Sequence identity deduplication

1. Run CD-HIT or MMseqs2 on each species' proteome at 40% sequence identity
2. Map cluster representatives to gene IDs
3. Verify no cluster spans train+test boundary
4. If clusters span splits, use cluster-stratified splitting

### 6.4 P1 — Statistical rigor

1. Add 95% confidence intervals (via bootstrap over test predictions, stratified by label)
2. Add Wilcoxon signed-rank test for paired model comparisons across the 5 seeds
3. Switch from `ddof=0` to `ddof=1` for sample standard deviation
4. Report per-model per-seed metrics, not just aggregated mean±std

### 6.5 P2 — Repeated holdout

Add 3–5 random split seeds (not just training seeds) to estimate split-induced variance. This quantifies how much the fixed split choice affects reported metrics.

### 6.6 P2 — Cross-species leakage check

If cross-species training is ever implemented (e.g., train on yeast, test on Fusarium), orthogroup-stratified splitting becomes mandatory — Fusarium orthologs of yeast training genes would be in the test set, creating direct label leakage through the ESM2+ORT feature combination.

---

## 7. What NOT to do

- Do NOT rerun splits or regenerate frozen protocol without explicit versioning as `frozen_protocol_v2`
- Do NOT modify `assign_splits()` without preserving the current split for comparison
- Do NOT delete or overwrite `results/frozen_protocol/splits/`
- Do NOT run CD-HIT/MMseqs2 without confirming input proteome completeness
- Do NOT claim "no data leakage" without completing the model organism orthogroup audit (6.1)
