import os
import sys


def normalize(lib_file):
    cell_counts = []
    norm_cell_counts = []
    count = 0
    total_size = 0

    for line in lib_file:
        cell_num, size = line.strip().split("\t")
        size = int(size)
        count += 1
        total_size += size
        cell_counts.append((cell_num, size))

    print(total_size, count, total_size // count)
    avg_size = total_size // count

    for cell_num, size in cell_counts:
        norm_size = size / avg_size
        norm_cell_counts.append((cell_num, norm_size))

    return norm_cell_counts


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"{os.path.basename(sys.argv[0])} <mRNA_lib.txt>\n")
        sys.exit(1)

    with os.open(sys.argv[1], "r") as lib_file:
        normalized_sizes = normalize(lib_file)

    for cell_num, norm_size in normalized_sizes:
        print(f"{cell_num}\t{norm_size}")
