# =========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2023 Siyan Liu (siyan.liu432@duke.edu)
# =========================================================================

import concurrent.futures
from collections import defaultdict
from importlib.resources import files
from operator import itemgetter

from cmdstanpy import CmdStanModel

from .constants import (
    DEFAULT_CHAINS,
    DEFAULT_NORM_LPF,
    DEFAULT_RUNS,
    DEFAULT_SAMPLE,
    DEFAULT_SEED,
    DEFAULT_WARMUP,
    MAX_SEED_INT,
)

MMLines = list[tuple[int, int, int]]
CountData = dict[int, float]


def mm_counts(mtx_lines: MMLines, norm_lpf: int) -> tuple[dict[int, int], dict[int, list[tuple[int, int]]]]:
    cumulative_counts = {}
    per_guide_counts = defaultdict(lambda: [])

    for guide, cell_id, guide_count in mtx_lines:
        if norm_lpf:
            if cell_id not in cumulative_counts:
                cumulative_counts[cell_id] = 0

            if guide_count <= norm_lpf:
                cumulative_counts[cell_id] += guide_count
        else:
            if cell_id not in cumulative_counts:
                cumulative_counts[cell_id] = guide_count
            else:
                cumulative_counts[cell_id] += guide_count

        per_guide_counts[guide].append((cell_id, guide_count))

    for key, value in cumulative_counts.items():
        if value == 0:
            cumulative_counts[key] = 1

    return cumulative_counts, per_guide_counts


def normalize(count_data: dict[int, int]) -> CountData:
    count = len(count_data)
    total_size = sum(count_data.values())
    avg_size = total_size / count

    norm_cell_counts = {cell_id: lib_size / avg_size for cell_id, lib_size in count_data.items()}

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


async def run(
    mm_lines,
    model,
    chains=DEFAULT_CHAINS,
    normalization_lpf=DEFAULT_NORM_LPF,
    num_parallel_runs=DEFAULT_RUNS,
    num_samples=DEFAULT_SAMPLE,
    num_warmup=DEFAULT_WARMUP,
    seed=DEFAULT_SEED,
):
    sorted_mm_lines = sorted(mm_lines, key=itemgetter(0, 1))
    cumulative_counts, per_guide_counts = mm_counts(sorted_mm_lines, normalization_lpf)
    normalized_counts = normalize(cumulative_counts)

    stan_model = CmdStanModel(stan_file=files("cleanser").joinpath(model))

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
