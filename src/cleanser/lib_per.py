# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author:Susan Liu
# =========================================================================

import argparse
import gzip
from collections import defaultdict


def mm_counts(mtx_file, col):
    counts = defaultdict(lambda: 0)

    for line in mtx_file:
        # Skip market matrix header/comments
        if line.startswith("%"):
            continue

        # skip the first non-comment line. It's just dimension info we
        # are ignoring
        break

    for line in mtx_file:
        mm_items = line.strip().split()
        counts[mm_items[col]] += int(mm_items[2])

    return counts


def get_args():
    parser = argparse.ArgumentParser("mm_counts", description="")
    parser.add_argument("matrix")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--cell", dest="col", action="store_const", const=1)
    group.add_argument("-g", "--guide", dest="col", action="store_const", const=0)

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    with gzip.open(args.matrix, "rt", encoding="utf8") as matr:
        counted_columns = mm_counts(matr, args.col)

    for key, value in counted_columns.items():
        print(f"{key}\t{value}")
