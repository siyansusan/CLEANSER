# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2023 Siyan Liu (siyan.liu432@duke.edu)
# =========================================================================

import argparse
import asyncio
import concurrent.futures
import multiprocessing as mp
import sys
from collections import defaultdict
from itertools import islice

from cmdstanpy import CmdStanModel
import numpy as np

MODEL_FILE = "./guide-mixture.stan"

MMLines = list[tuple[int, int, int]]
CountData = dict[int, float]


def read_mm_file(mtx_file) -> MMLines:
    mm_lines = []
    for line in mtx_file:
        # Skip market matrix header/comments
        if line.startswith("%"):
            continue

        # skip the first non-comment line. It's just dimension info we
        # are ignoring
        break

    for line in mtx_file:
        guide, cell, count = line.strip().split()
        mm_lines.append((int(guide), int(cell), int(count)))

    return mm_lines


def mm_counts(mtx_lines: MMLines) -> tuple[dict[int, int], dict[int, list[tuple[int, int]]]]:
    cumulative_counts = defaultdict(lambda: 0)
    per_guide_counts = defaultdict(lambda: [])

    for guide, cell_id, guide_count in mtx_lines:
        cumulative_counts[cell_id] += guide_count
        per_guide_counts[guide].append((cell_id, guide_count))

    return cumulative_counts, per_guide_counts


def normalize(count_data: dict[int, int]) -> CountData:
    norm_cell_counts = {}
    count = 0
    total_size = 0

    for cell_id, size in count_data.items():
        count += 1
        total_size += size

    avg_size = total_size // count

    for cell_id, size in count_data.items():
        norm_size = size / avg_size
        norm_cell_counts[cell_id] = norm_size

    return norm_cell_counts


async def go():
    await asyncio.sleep(1)
    return 1, 1


def run_stan(stan_args):
    model, guide_id, X, L, num_warmup, num_samples = stan_args

    fit = model.sample(data={"N": len(X), "X": X, "L": L}, iter_warmup=num_warmup, iter_sampling=num_samples)
    return guide_id, fit


def output_posteriors(stan_output, output_file):
    if output_file is None:
        output_file = sys.stdout

    for x, y, z in stan_output:
        print(f"{x} {y} {z}", file=output_file)


def build_model():
    return CmdStanModel(stan_file=MODEL_FILE)


async def run(input_file, output_file, num_warmup, num_samples, num_cores):
    mm_lines = read_mm_file(input_file)
    sorted_mm_lines = sorted(mm_lines, key=lambda x: x[0])
    cumulative_counts, per_guide_counts = mm_counts(sorted_mm_lines)
    normalized_counts = normalize(cumulative_counts)

    stan_model = build_model()

    stan_params = [
        # X, L, num_warmup, num_samples
        (
            stan_model,
            guide_id,
            [guide_count for _, guide_count in guide_counts],
            [normalized_counts[cell_id] for cell_id, _ in guide_counts],
            num_warmup,
            num_samples,
        )
        for guide_id, guide_counts in per_guide_counts.items()
        # for guide_id, guide_counts in islice(per_guide_counts.items(), 3)
    ]

    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        for guide_id, samples in executor.map(run_stan, stan_params):
            results[guide_id] = samples

    for samples in results.values():
        print(samples.stan_variables().keys())
        print(
            f"r={np.mean(samples.stan_variable('r'))}\tmu={np.mean(samples.stan_variable('nbMean'))}\tlambda={np.mean(samples.stan_variable('lambda'))}"
        )

    # output_posteriors(stan_output, output_file)
    return results


def get_args():
    parser = argparse.ArgumentParser(
        "guide-mixture",
        description="Crispr Library Evaluation and Ambient Noise Suppression for Enhanced scRNA-seq",
    )

    parser.add_argument(
        "-i",
        "--input",
        help="Matrix Market file of guide library information",
        type=argparse.FileType("r", encoding="utf-8"),
        required=True,
    )
    parser.add_argument("-o", "--output", help="output file name of per-guide posterior probabilities")
    parser.add_argument("-n", "--num-samples", type=int, default=1000)
    parser.add_argument("-w", "--num-warmup", type=int, default=300)
    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        default=mp.cpu_count(),
        help="Number of CPUs to use for running the model",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    try:
        asyncio.run(run(args.input, args.output, args.num_warmup, args.num_samples, args.cores))
        print("done!")
    except KeyboardInterrupt:
        sys.exit(1)
