import argparse
import json

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate JSON records into a TSV manifest")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    rows = [json.load(open(path)) for path in args.inputs]
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
