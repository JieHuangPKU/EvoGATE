import argparse
from pathlib import Path

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render a Figure3-specific frozen protocol config with ESM2 overrides")
    parser.add_argument("--base-config", required=True, type=str)
    parser.add_argument("--output-config", required=True, type=str)
    parser.add_argument("--esm2-cache-root", default=None, type=str)
    parser.add_argument("--mock-embedding-dim", default=None, type=int)
    parser.add_argument("--esm2-projection-dim", default=None, type=int)
    parser.add_argument("--fusion-mode", default=None, type=str)
    parser.add_argument("--fusion-hidden-dim", default=None, type=int)
    parser.add_argument("--fusion-dropout", default=None, type=float)
    parser.add_argument("--loss-type", default=None, type=str)
    parser.add_argument("--pos-weight-mode", default=None, type=str)
    parser.add_argument("--pos-weight-scale", default=None, type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(Path(args.base_config).read_text(encoding="utf-8"))
    config.setdefault("esm2", {})
    if args.esm2_cache_root is not None:
        config["esm2"]["cache_root"] = args.esm2_cache_root
    if args.mock_embedding_dim is not None:
        config["esm2"]["mock_embedding_dim"] = int(args.mock_embedding_dim)
    if args.esm2_projection_dim is not None:
        config.setdefault("runtime", {})
        config["runtime"]["esm2_projection_dim"] = int(args.esm2_projection_dim)

    models = config.setdefault("models", {})
    if "GraphSAGE" not in models:
        template_key = next(
            (key for key in ["GraphSAGE_ORT_EXP_SUB", "GraphSAGE_ESM2", "GraphSAGE_ORT_EXP_SUB_ESM2"] if key in models),
            None,
        )
        if template_key is None:
            raise ValueError("Could not infer a GraphSAGE model template from the base config")
        models["GraphSAGE"] = dict(models[template_key])
        models["GraphSAGE"]["feature_setting"] = "ORT_EXP_SUB"
    if args.fusion_mode is not None:
        models["GraphSAGE"]["fusion_mode"] = str(args.fusion_mode)
    if args.fusion_hidden_dim is not None:
        models["GraphSAGE"]["fusion_hidden_dim"] = int(args.fusion_hidden_dim)
    if args.fusion_dropout is not None:
        models["GraphSAGE"]["fusion_dropout"] = float(args.fusion_dropout)
    if args.loss_type is not None:
        models["GraphSAGE"]["loss_type"] = str(args.loss_type)
    if args.pos_weight_mode is not None:
        models["GraphSAGE"]["pos_weight_mode"] = str(args.pos_weight_mode)
    if args.pos_weight_scale is not None:
        models["GraphSAGE"]["pos_weight_scale"] = float(args.pos_weight_scale)

    output_path = Path(args.output_config)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


if __name__ == "__main__":
    main()
