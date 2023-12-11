import argparse
import asyncio
import multiprocessing as mp
import sys
import os.path
from random import randint

import numpy as np

from .guide_mixture import CS_MODEL_FILE, DC_MODEL_FILE, MAX_SEED_INT, run


def output_posteriors(stan_results, output_file):
    output_all_posteriors = True
    if not os.path.isdir("posteriors"):
        print("Please create a 'posteriors' directory if you want all posterior values saved.")
        output_all_posteriors = False

    for guide_id, (samples, cell_info) in stan_results.items():
        for i, (cell_id, _) in enumerate(cell_info):
            if output_all_posteriors:
                with open(f"posteriors/{guide_id}_{cell_id}.txt", "w", encoding="ascii") as post_out:
                    post_out.write(f"{', '.join(str(n) for n in sorted(samples.stan_variable('PZi')[i]))}")

            output_file.write(f"{guide_id}\t{cell_id}\t{np.median(samples.stan_variable('PZi')[i])}\n")


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
    parser.add_argument("-w", "--num-warmup", type=int, default=300, help="The number of warmup iterations per chain")
    parser.add_argument(
        "-s", "--seed", type=int, default=randint(0, MAX_SEED_INT), help="The seed for the random number generator"
    )
    parser.add_argument("-c", "--chains", type=int, default=4, help="The number of Markov chains")
    parser.add_argument(
        "-p",
        "--parallel-runs",
        type=int,
        default=mp.cpu_count(),
        help="Number of guide models to run in parallel",
    )
    parser.add_argument(
        "--nt",
        "--normalization-threshold",
        type=int,
        default=2,
        help="The upper threshold for including the guide counts in guide count normalization",
        dest="normalization_threshold",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--dc",
        "--direct-capture",
        action="store_true",
        help="Use direct capture mixture model",
    )
    group.add_argument(
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
        results = asyncio.run(
            run(
                args.input,
                model_file,
                args.num_warmup,
                args.num_samples,
                args.parallel_runs,
                args.chains,
                args.normalization_threshold,
                args.seed,
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
