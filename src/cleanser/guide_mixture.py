# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2023 Siyan Liu (siyan.liu432@duke.edu)
# =========================================================================

import concurrent.futures
from collections import defaultdict

from cmdstanpy import CmdStanModel

MODEL_FILE = "./guide-mixture.stan"
MAX_SEED_INT = 4_294_967_295  # 2^32 - 1, the largest seed allowed by STAN

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


def run_stan(stan_args):
    model, guide_id, X, L, num_warmup, num_samples, chains, seed = stan_args

    fit = model.sample(
        data={"N": len(X), "X": X, "L": L},
        iter_warmup=num_warmup,
        iter_sampling=num_samples,
        chains=chains,
        seed=seed,
        show_progress=False,
    )
    return guide_id, fit


async def run(input_file, num_warmup, num_samples, num_parallel_runs, chains, seed):
    mm_lines = read_mm_file(input_file)
    sorted_mm_lines = sorted(mm_lines, key=lambda x: x[0])
    cumulative_counts, per_guide_counts = mm_counts(sorted_mm_lines)
    normalized_counts = normalize(cumulative_counts)

    stan_model = CmdStanModel(stan_file=MODEL_FILE)

    stan_params = [
        (
            stan_model,
            guide_id,
            [guide_count for _, guide_count in guide_counts],  # X
            [normalized_counts[cell_id] for cell_id, _ in guide_counts],  # L
            num_warmup,
            num_samples,
            chains,
            (seed + guide_id) % MAX_SEED_INT,
        )
        for guide_id, guide_counts in per_guide_counts.items()
    ]

    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel_runs) as executor:
        for guide_id, samples in executor.map(run_stan, stan_params):
            results[guide_id] = (samples, per_guide_counts[guide_id])

    return results
