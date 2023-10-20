import argparse
import asyncio
import multiprocessing as mp
import sys
from random import randint

import numpy as np

from .guide_mixture import MAX_SEED_INT, run


def output_posteriors(stan_results, output_file):
    for guide_id, (samples, cell_info) in stan_results.items():
        for i, (cell_id, _) in enumerate(cell_info):
            with open(f"post/{guide_id}_{cell_id}.txt", "w") as post_out:
                post_out.write(f"{sorted(samples.stan_variable('PZi')[i])}")
            output_file.write(f"{guide_id}\t{cell_id}\t{np.median(samples.stan_variable('PZi')[i])}\n")


def output_samples(stan_results, output_file):
    output_file.write("guide id\tr\tmu\tDisp\tlambda\n")
    for guide_id, (samples, _) in stan_results.items():
        r = samples.stan_variable("r")
        mu = samples.stan_variable("nbMean")
        disp = samples.stan_variable("nbDisp")
        lamb = samples.stan_variable("lambda")
        for i, r_samp in enumerate(r):
            output_file.write(f"{guide_id}\t{r_samp}\t{mu[i]}\t{disp[i]}\t{lamb[i]}\n")


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
    parser.add_argument(
        "-o",
        "--posteriors-output",
        help="output file name of per-guide/cell posterior probabilities",
        type=argparse.FileType("w", encoding="utf-8"),
        default=sys.stdout,
    )
    parser.add_argument(
        "--so",
        "--samples-output",
        help="output file name of per-guide/cell posterior probabilities",
        type=argparse.FileType("w", encoding="utf-8"),
        default=sys.stdout,
    )
    parser.add_argument("-n", "--num-samples", type=int, default=1000)
    parser.add_argument("-w", "--num-warmup", type=int, default=300)
    parser.add_argument("-s", "--seed", type=int, default=randint(0, MAX_SEED_INT))
    parser.add_argument("-c", "--chains", type=int, default=4)
    parser.add_argument(
        "-p",
        "--parallel-runs",
        type=int,
        default=mp.cpu_count(),
        help="Number of guide models to run in parallel",
    )

    return parser.parse_args()


def run_cli():
    args = get_args()
    try:
        results = asyncio.run(
            run(args.input, args.num_warmup, args.num_samples, args.parallel_runs, args.chains, args.seed)
        )
        output_posteriors(results, args.posteriors_output)
        output_samples(results, args.so)
        for _, (samples, _) in results.items():
            print(
                f"r={np.median(samples.stan_variable('r'))}\tmu={np.median(samples.stan_variable('nbMean'))}\tlambda={np.median(samples.stan_variable('lambda'))}"
            )

        print(f"Random seed: {args.seed}")
    except KeyboardInterrupt:
        sys.exit(1)
