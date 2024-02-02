import argparse
import asyncio
import os.path
import sys

import numpy as np

from .constants import (
    CS_MODEL_FILE,
    DC_MODEL_FILE,
    DEFAULT_CHAINS,
    DEFAULT_NORM_LPF,
    DEFAULT_RUNS,
    DEFAULT_SAMPLE,
    DEFAULT_SEED,
    DEFAULT_WARMUP,
)
from .guide_mixture import MMLines, run


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


def output_posteriors(stan_results, output_file):
    output_all_posteriors = True
    if not os.path.isdir("posteriors"):
        print("Please create a 'posteriors' directory if you want all posterior values saved.")
        output_all_posteriors = False

    for guide_id, (samples, cell_info) in stan_results.items():
        pzi = np.transpose(samples.stan_variable("PZi"))
        for i, (cell_id, _) in enumerate(cell_info):
            if output_all_posteriors:
                with open(f"posteriors/{guide_id}_{cell_id}.txt", "w", encoding="ascii") as post_out:
                    post_out.write(f"{', '.join(str(n) for n in sorted(pzi[i]))}")

            output_file.write(f"{guide_id}\t{cell_id}\t{np.median(pzi[i])}\n")


def output_cs_samples(stan_results, output_file):
    output_file.write("guide id\tr\tmu\tDisp\tlambda\n")
    for guide_id, (samples, _) in stan_results.items():
        r = samples.stan_variable("r")
        mu = samples.stan_variable("nbMean")
        disp = samples.stan_variable("nbDisp")
        lamb = samples.stan_variable("lambda")
        for i, r_samp in enumerate(r):
            output_file.write(f"{guide_id}\t{r_samp}\t{mu[i]}\t{disp[i]}\t{lamb[i]}\n")


def output_dc_samples(stan_results, output_file):
    output_file.write("guide id\tr\tmu\tDisp\tn_nbMean\tn_nbDisp\n")
    for guide_id, (samples, _) in stan_results.items():
        r = samples.stan_variable("r")
        mu = samples.stan_variable("nbMean")
        disp = samples.stan_variable("nbDisp")
        n_mean = samples.stan_variable("n_nbMean")
        n_disp = samples.stan_variable("n_nbDisp")
        for i, r_samp in enumerate(r):
            output_file.write(f"{guide_id}\t{r_samp}\t{mu[i]}\t{disp[i]}\t{n_mean[i]}\t{n_disp[i]}\n")


def output_cs_stats(results):
    for _, (samples, _) in results.items():
        print(
            f"r={np.median(samples.stan_variable('r'))}\tmu={np.median(samples.stan_variable('nbMean'))}\tlambda={np.median(samples.stan_variable('lambda'))}"
        )


def output_dc_stats(results):
    for _, (samples, _) in results.items():
        print(
            f"r={np.median(samples.stan_variable('r'))}\tmu={np.median(samples.stan_variable('nbMean'))}\tn_nbMean={np.median(samples.stan_variable('n_nbMean'))}\tn_nbDisp={np.median(samples.stan_variable('n_nbDisp'))}"
        )


def get_args():
    parser = argparse.ArgumentParser(
        "cleanser",
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
        help="output file name of sample data",
        type=argparse.FileType("w", encoding="utf-8"),
        default=sys.stdout,
    )
    parser.add_argument(
        "-n", "--num-samples", type=int, default=DEFAULT_SAMPLE, help="The number of samples to take of the model"
    )
    parser.add_argument(
        "-w", "--num-warmup", type=int, default=DEFAULT_WARMUP, help="The number of warmup iterations per chain"
    )
    parser.add_argument("-s", "--seed", type=int, default=DEFAULT_SEED, help="The seed for the random number generator")
    parser.add_argument("-c", "--chains", type=int, default=DEFAULT_CHAINS, help="The number of Markov chains")
    parser.add_argument(
        "-p",
        "--parallel-runs",
        type=int,
        default=DEFAULT_RUNS,
        help="Number of guide models to run in parallel",
    )
    parser.add_argument(
        "--lpf",
        "--normalization-lpf",
        type=int,
        default=DEFAULT_NORM_LPF,
        help="The upper limit for including the guide counts in guide count normalization. Set to 0 for no limit.",
        dest="normalization_lpf",
    )
    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument(
        "--dc",
        "--direct-capture",
        action="store_true",
        help="Use direct capture mixture model",
    )
    model_group.add_argument(
        "--cs",
        "--crop-seq",
        action="store_true",
        help="Use crop-seq mixture model",
    )

    return parser.parse_args()


def run_cli():
    args = get_args()

    if args.dc:
        model_file = DC_MODEL_FILE
    elif args.cs:
        model_file = CS_MODEL_FILE

    try:
        mm_lines = read_mm_file(args.input)
        results = asyncio.run(
            run(
                mm_lines,
                model_file,
                chains=args.chains,
                normalization_lpf=args.normalization_lpf,
                num_parallel_runs=args.parallel_runs,
                num_samples=args.num_samples,
                num_warmup=args.num_warmup,
                seed=args.seed,
            )
        )
        output_posteriors(results, args.posteriors_output)
        if args.dc:
            output_dc_samples(results, args.so)
            output_dc_stats(results)
        elif args.cs:
            output_cs_samples(results, args.so)
            output_cs_stats(results)

        print(f"Random seed: {args.seed}")
    except KeyboardInterrupt:
        sys.exit(1)
