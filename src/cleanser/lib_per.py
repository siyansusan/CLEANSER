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

    # Skip market matrix headers
    mtx_file.readline()
    mtx_file.readline()
    mtx_file.readline()

    for line in mtx_file:
        mm_items = line.strip().split()
        counts[mm_items[col]] += int(mm_items[2])

    # lib_size=list(d.values())
    # med_lib_size=statistics.median(lib_size)

    # for key in d:
    #    d[key]=med_lib_size/d[key]

    # TO FIND TOP
    # sorted_d=sorted(d,key=d.get)
    # print(sorted_d)
    # print(d["44118"])

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
