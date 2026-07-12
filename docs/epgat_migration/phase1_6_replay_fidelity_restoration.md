# Phase 1.6 Replay Fidelity Restoration

## 1. Why Phase 1.6

- `why_replay_diff_report.md` showed that celegans and human replay gaps were dominated by fidelity mismatch, not by label inversion, eval-column inversion, or obvious merge-key breakage.
- For that reason, Phase 2 should not start until the legacy-compatible replay is pushed closer to the old conditions.

## 2. Conditions Restored in This Round

- Restored old-style seed semantics: command seed base retained, but per-run seeds follow the old `run_idx` loop behavior.
- Restored old validation split semantics: 20 percent test, then 5 percent of the remaining train used as validation, implemented as 4 percent of the total labeled set.
- Restored old final-state testing behavior: replay now evaluates the final model state rather than reloading the best validation checkpoint.
- Restored multi-run aggregation partially: celegans completed 3 runs; human completed 1 restored run in this practical restoration step.
- Restored epoch budget partially: celegans now uses 100 epochs; human was raised from 2 to 10 epochs, which is still below the old 100-epoch condition.

## 3. Celegans Restored Replay

- Previous AUROC/AUPRC/MCC: 0.2418 / 0.0589 / -0.0358
- Restored AUROC/AUPRC/MCC: 0.9169 / 0.5554 / 0.4643
- The restored celegans replay is now very close to the old raw full-contract GAT AUROC reference (gap -0.0067).
- This strongly supports the diagnosis that the original celegans failure was driven by replay fidelity mismatch rather than a label/probability/merge bug.
- Current remaining gap is mostly about which old reference should be treated as authoritative: the raw full-contract row or the higher threshold-sweep/processed figure proxy.

## 4. Human Restored Replay

- Previous AUROC/AUPRC/MCC: 0.7965 / 0.3610 / 0.3069
- Restored AUROC/AUPRC/MCC: 0.8516 / 0.4296 / 0.3100
- The restored human replay is clearly improved versus the previous replay, but still remains below the old raw full-contract GAT AUROC reference by 0.0389.
- The main unresolved difference is that human still does not fully match the old condition: it is still only a partial restoration on epochs and run aggregation, and it still is not reproducing the processed-summary provenance used by the old figure.

## 5. Conclusion

- Celegans replay fidelity is now largely explained and mostly restored.
- Human replay fidelity is only partially restored; the gap is reduced, but not yet small enough to claim the old condition has been faithfully reproduced.
- The overall replay fidelity question is therefore not fully closed.
- Recommended next step: fix replay fidelity further before Phase 2.