import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from io import StringIO
from math import log2
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from cleanser.run import MMLines, read_mm_file


@dataclass
class DCSampleMetadata:
    "Direct capture sample information"

    guide_id: str
    r: float
    mu: float
    dispersion: float
    neg_binom_mean: float
    neg_binom_dispersion: float


@dataclass
class CSSampleMetadata:
    "CROP-seq sample information"

    guide_id: str
    r: float
    mu: float
    dispersion: float
    lamb: float  # lambda


@dataclass
class Prediction:
    "CLEANSER Prediction for a guide/cell pair"
    guide_id: str
    cell_id: str
    prediction: float


def per_guide_sample_averages(
    samples: list[CSSampleMetadata | DCSampleMetadata],
) -> list[CSSampleMetadata | DCSampleMetadata]:

    sample_count = len(samples)

    if sample_count == 0:
        raise ValueError("No samples")

    if isinstance(samples[0], CSSampleMetadata):
        sample_tallies = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0])
        for sample in samples:
            running_tally = sample_tallies[sample.guide_id]
            running_tally[0] += 1
            running_tally[1] += sample.r
            running_tally[2] += sample.mu
            running_tally[3] += sample.dispersion
            running_tally[4] += sample.lamb

        return [
            CSSampleMetadata(
                guide_id, tally[1] / tally[0], tally[2] / tally[0], tally[3] / tally[0], tally[4] / tally[0]
            )
            for guide_id, tally in sample_tallies.items()
        ]

    if isinstance(samples[0], DCSampleMetadata):
        sample_tallies = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])
        for sample in samples:
            running_tally = sample_tallies[sample.guide_id]
            running_tally[0] += 1
            running_tally[1] += sample.r
            running_tally[2] += sample.mu
            running_tally[3] += sample.dispersion
            running_tally[4] += sample.neg_binom_mean
            running_tally[5] += sample.neg_binom_dispersion

        return [
            DCSampleMetadata(
                guide_id,
                tally[1] / tally[0],
                tally[2] / tally[0],
                tally[3] / tally[0],
                tally[4] / tally[0],
                tally[5] / tally[0],
            )
            for guide_id, tally in sample_tallies.items()
        ]

    raise ValueError("Invalid sample type")


def plot_hist(values, bin_count, title="", x_label="", y_label="", density=False) -> Figure:
    fig, ax = plt.subplots()

    # the histogram of the data
    ax.hist(values, bins=bin_count, density=density)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()

    return fig


def plot_scatter(x, y, title="", x_label="", y_label="") -> Figure:
    fig, ax = plt.subplots()

    # the histogram of the data
    ax.scatter(x, y, s=0.5)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()

    return fig


def sample_average_output(averages: list[CSSampleMetadata | DCSampleMetadata]) -> str:
    average_text = StringIO()

    if isinstance(averages[0], CSSampleMetadata):
        average_text.write("guide\tr\tmu\tdisp\tlambda\n")
        for ave in averages:
            average_text.write(f"{ave.guide_id}\t{ave.r}\t{ave.mu}\t{ave.dispersion}\t{ave.lamb}\n")

    if isinstance(averages[0], DCSampleMetadata):
        average_text.write("guide\tr\tmu\tdisp\tn_mu\tn_disp\n")
        for ave in averages:
            average_text.write(
                f"{ave.guide_id}\t{ave.r}\t{ave.mu}\t{ave.dispersion}\t{ave.neg_binom_mean}\t{ave.neg_binom_dispersion}\n"
            )

    return average_text.getvalue()


def sample_average_histogram(averages: list[CSSampleMetadata | DCSampleMetadata]) -> Figure:
    average_count = len(averages)

    if average_count == 0:
        raise ValueError("No averages")

    avg_r = [a.r for a in averages]
    fig = plot_hist(avg_r, 30, title="cleanser", x_label="sample average", y_label="count")
    return fig


def assigned_counts_histogram(predictions: list[Prediction], mm_lines: MMLines, threshold: float) -> Figure:
    thresholded_preds = {(p.guide_id, p.cell_id) for p in predictions if p.prediction >= threshold}
    umis = [umi for guide_id, cell_id, umi in mm_lines if (guide_id, cell_id) in thresholded_preds]
    fig = plot_hist(umis, 100, title="cleanser", x_label="UMI", y_label="count")
    return fig


def posterior_umi_scatterplot(predictions: list[Prediction], mm_lines: MMLines) -> Figure:
    umi_keys = [(p.guide_id, p.cell_id) for p in predictions]
    posteriors = [p.prediction for p in predictions]

    # order the umis to match the posteriors
    umi_dict = {(guide_id, cell_id): umi for guide_id, cell_id, umi in mm_lines}
    umis = [umi_dict[k] for k in umi_keys]

    fig = plot_scatter(posteriors, umis, title="cleanser", x_label="Posterior", y_label="UMI")

    return fig


def posterior_umi_scatterplot_log2(predictions: list[Prediction], mm_lines: MMLines) -> Figure:
    umi_keys = [(p.guide_id, p.cell_id) for p in predictions]
    posteriors = [p.prediction for p in predictions]

    # order the umis to match the posteriors
    umi_dict = {(guide_id, cell_id): log2(umi) for guide_id, cell_id, umi in mm_lines}
    umis = [umi_dict[k] for k in umi_keys]

    fig = plot_scatter(posteriors, umis, title="cleanser", x_label="Posterior", y_label="log2(UMI)")

    return fig


def calc_moi(predictions: list[Prediction], threshold: float) -> float:
    cells = set()
    predict_count = 0
    for prediction in predictions:
        cells.add(prediction.cell_id)
        if prediction.prediction >= threshold:
            predict_count += 1

    return predict_count / len(cells)


def calc_coverage(predictions: list[Prediction], threshold: float) -> float:
    guides = set()
    predict_count = 0
    for prediction in predictions:
        guides.add(prediction.guide_id)
        if prediction.prediction >= threshold:
            predict_count += 1

    return predict_count / len(guides)


def read_predictions(input_file) -> list[Prediction]:
    reader = csv.reader(input_file, delimiter="\t")

    return [Prediction(guide_id=l[0], cell_id=l[1], prediction=float(l[2])) for l in reader]


def read_sample_data(sample_data_file) -> list[CSSampleMetadata | DCSampleMetadata]:
    reader = csv.DictReader(sample_data_file, delimiter="\t")
    if len(reader.fieldnames) == 5:
        return [
            CSSampleMetadata(l["guide id"], float(l["r"]), float(l["mu"]), float(l["Disp"]), float(l["lambda"]))
            for l in reader
        ]

    if len(reader.fieldnames) == 6:
        return [
            DCSampleMetadata(
                l["guide id"],
                float(l["r"]),
                float(l["mu"]),
                float(l["Disp"]),
                float(l["n_nbMean"]),
                float(l["n_nbDisp"]),
            )
            for l in reader
        ]

    raise ValueError("Invalid sample data file")


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate CLEANSER QC information")
    parser.add_argument("-i", "--input", type=argparse.FileType(), help="Cleanser posterior output file", required=True)
    parser.add_argument("-o", "--output-directory", help="Cleanser QC output directory", required=True)
    parser.add_argument("-g", "--guide-counts", type=argparse.FileType(), help="Guide count file")
    parser.add_argument("-s", "--samples", type=argparse.FileType(), help="Cleanser sampling data")
    parser.add_argument(
        "-t", "--threshold", type=float, default=0.0, help="Disregard assignment probabilities below this value"
    )

    return parser.parse_args()


def write(output_dir, filename, contents):
    output_file = output_dir / Path(filename)
    with open(output_file, "w", encoding="utf-8") as file:
        file.write(contents)


def run_cli():
    args = get_args()

    output_dir = Path(args.output_directory)
    if not output_dir.exists():
        output_dir.mkdir()
    elif not output_dir.is_dir():
        raise ValueError(f"Output directory {output_dir} must be a directory")

    predictions = read_predictions(args.input)
    moi = calc_moi(predictions, args.threshold)
    moi_text = f"MOI: {moi}"
    print(moi_text)
    write(output_dir, "moi.txt", moi_text + "\n")

    coverage = calc_coverage(predictions, args.threshold)
    coverage_text = f"Coverage: {coverage}"
    print(coverage_text)
    write(output_dir, "coverage.txt", coverage_text + "\n")

    if args.guide_counts is not None:
        mm_lines = read_mm_file(args.guide_counts)
        assigned_counts_histogram(predictions, mm_lines, args.threshold)
        plt.savefig(output_dir / Path("umi_hist.png"))

        posterior_umi_scatterplot(predictions, mm_lines)
        plt.savefig(output_dir / Path("umi_count_scatter.png"))

        posterior_umi_scatterplot_log2(predictions, mm_lines)
        plt.savefig(output_dir / Path("umi_count_scatter_log2.png"))

    if args.samples is not None:
        samples = read_sample_data(args.samples)
        assert len(samples) > 0
        averages = per_guide_sample_averages(samples)
        write(output_dir, "sample_avg.txt", sample_average_output(averages))
        sample_average_histogram(averages)
        plt.savefig(output_dir / Path("sample_average_hist.png"))
