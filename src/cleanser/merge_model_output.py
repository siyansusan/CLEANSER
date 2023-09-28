import argparse
import glob
import os
import sys


def get_args():
    parser = argparse.ArgumentParser(
        "merge_model_output",
        description="merge per-guide guide mixture model output files into a single file ",
    )
    parser.add_argument(
        "posterior_probs_dir",
        help="The directory containing the model output files",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=argparse.FileType("w", encoding="utf-8"),
        default=sys.stdout,
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    output = args.output_file

    for filepath in glob.iglob(f"{args.posterior_probs_dir}/*.txt"):
        num = os.path.basename(filepath).strip().split(".")[0]
        with open(filepath, "r") as txt_file:
            for line in txt_file:
                cell, umi = line.strip().split("\t")
                print(f"{num}\t{cell}\t{umi}", file=output)
