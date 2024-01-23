# CLEANSER1.0

**C**rispr **L**ibrary **E**valuation and **A**mbient **N**oise **S**uppression for **E**nhanced sc**R**NA-seq
CLEANSER is a gRNA-cell assignment method that uses a mixture of two distinct distributions to model ambient and native gRNA presence in perturb-seq CRISPR libraries. CLEANSER takes into account gRNA-specific and cell-specific biases and generates a probability value of whether or not a gRNA is expressed natively in a cell.

## Installation

**Note**: These steps have been tested on macOS and Linux, but not Windows.

tl;dr:

    git clone https://github.com/siyansusan/CLEANSER1.0.git
    cd CLEANSER1.0
    pip install .
    install_cmdstan

First, clone the repository using `git` using either the https or ssh url:

Https: `git clone https://github.com/siyansusan/CLEANSER1.0.git`
SSH: `git clone git@github.com:siyansusan/CLEANSER1.0.git`

Next change into the CLEANSER1.0 directory and run `pip install .`.

CLEANSER depends on something called [CmdStan](https://mc-stan.org/docs/cmdstan-guide/index.html) which, if you don't have, you'll need to install. Fortunately a script was installed as part of CLEANSER to make this easier. To install CmdStan run `install_cmdstan` which will download and install the latest version of CmdStan for CLEANSER to use.

## Usage

    usage: cleanser [-h] -i INPUT [-o POSTERIORS_OUTPUT] [--so SO] [-n NUM_SAMPLES] [-w NUM_WARMUP] [-s SEED] [-c CHAINS]
                        [-p PARALLEL_RUNS] [--lpf NORMALIZATION_LPF] (--dc | --cs)

    Crispr Library Evaluation and Ambient Noise Suppression for Enhanced scRNA-seq

    options:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
                            Matrix Market file of guide library information
    -o POSTERIORS_OUTPUT, --posteriors-output POSTERIORS_OUTPUT
                            output file name of per-guide/cell posterior probabilities
    --so SO, --samples-output SO
                            output file name of per-guide/cell posterior probabilities
    -n NUM_SAMPLES, --num-samples NUM_SAMPLES
    -w NUM_WARMUP, --num-warmup NUM_WARMUP
                            The number of warmup iterations per chain
    -s SEED, --seed SEED  The seed for the random number generator
    -c CHAINS, --chains CHAINS
                            The number of Markov chains
    -p PARALLEL_RUNS, --parallel-runs PARALLEL_RUNS
                            Number of guide models to run in parallel
    --lpf NORMALIZATION_LPF, --normalization-lpf NORMALIZATION_LPF
                            The upper limit for including the guide counts in guide count normalization. Set to 0 for no
                            limit.
    --dc, --direct-capture
                            Use direct capture mixture model
    --cs, --crop-seq      Use crop-seq mixture model
