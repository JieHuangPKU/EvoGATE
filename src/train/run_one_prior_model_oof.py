import io
import os
import sys

import pandas as pd

from src.train.run_cross_species_prior_oof import (
    enable_utf8_stdout,
    require_epgat_env,
    load_yaml,
    load_datasets,
    alignment_audit_rows,
    run_oof,
)


def main():
    enable_utf8_stdout()
    require_epgat_env()
    model_name = os.environ.get("PRIOR_MODEL_NAME", "").strip()
    if not model_name:
        raise SystemExit("缺少 PRIOR_MODEL_NAME")
    os.environ["PRIOR_MODELS"] = model_name

    outdir = "outputs/support_prior/intermediate"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    baseline_config = load_yaml("configs/baseline.yaml")
    raw_support, support_ready, x_support, y_support, raw_fus, fus_ready, x_fus = load_datasets()
    cv_df, oof_df, fus_df, validation_rows, best_model = run_oof(
        baseline_config,
        support_ready,
        x_support,
        y_support,
        fus_ready,
        x_fus,
    )
    cv_df.to_csv(os.path.join(outdir, "{}_cv.tsv".format(model_name)), sep="\t", index=False)
    oof_df.to_csv(os.path.join(outdir, "{}_oof.tsv".format(model_name)), sep="\t", index=False)
    fus_df.to_csv(os.path.join(outdir, "{}_fusarium.tsv".format(model_name)), sep="\t", index=False)
    pd.DataFrame(validation_rows).to_csv(os.path.join(outdir, "{}_validation.tsv".format(model_name)), sep="\t", index=False)
    print("单模型 OOF 完成：{} -> {}".format(model_name, outdir))


if __name__ == "__main__":
    main()
