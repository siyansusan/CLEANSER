import argparse
import gzip
import sys


def get_guide_barcodes(feature_filename):
    barcodes = set()
    with gzip.open(feature_filename, "r") as feature:
        for i, line in enumerate(feature, start=1):
            _, _, flag = line.decode("utf8").strip().split("\t")
            if flag == "CRISPR Guide Capture":
                barcodes.add(i)

    return barcodes


def process_mm(mm_filename, output_file, barcodes):
    with gzip.open(mm_filename, "r") as all_m:
        output_file.write(all_m.readline().decode("utf8"))
        output_file.write(all_m.readline().decode("utf8"))
        output_file.write(all_m.readline().decode("utf8"))

        for line in all_m:
            line = line.decode("utf8")
            guide, _, _ = line.split()
            if int(guide) in barcodes:
                output_file.write(line)


def get_args():
    parser = argparse.ArgumentParser(description="Generate a MM file suitable for CLEANSER from Cell Ranger outputs")
    parser.add_argument("-m", "--matrix-market", help="Cell Ranger matrix market output file", required=True)
    parser.add_argument("-f", "--features", help="Cell Ranger features output", required=True)
    parser.add_argument(
        "-o", "--output", help="output file for use by CLEANSER", type=argparse.FileType("w"), default=sys.stdout
    )

    return parser.parse_args()


def run_cli():
    args = get_args()
    barcodes = get_guide_barcodes(args.features)
    process_mm(args.matrix_market, args.output, barcodes)
