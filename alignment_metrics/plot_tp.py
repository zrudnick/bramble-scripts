#!/usr/bin/env python3
"""
plot_tp.py

Usage:
    python plot_tp.py --config config/short_read.yaml [--metrics precision recall f1]

For each tool label listed in config.yaml's `tools:` (see ../tools.yaml for
what each label means), loads every per-sample, per-transcript TSV produced
by bam_to_tsv.py/run.py (one row per transcript, one file per sample) from
tsv/<config's name>/<label>/*.tsv, pools them by exact isoform count,
computes every derived metric from those pooled sums, then bins by
isoform-count range for plotting.

Expected per-transcript TSV columns (see bam_to_tsv.py):
    n_isoforms        - number of transcripts overlapping the locus
    n_reads_total      - true read count for this transcript, from the FASTQ
                         (blank/NaN if no truth cache was available)
    n_reads_mapped     - reads that got >=1 alignment in this method's BAM
    n_reads_hit        - of those, how many included the true transcript
    precision_sum      - sum of per-read fractional credit (1/n_alignments
                         when correct, else 0)
    n_alignments_sum   - sum of n_alignments over every mapped read (lets us
                         reconstruct average multimapping degree)

Derived metrics (computed after pooling, never cached):
    precision          precision_sum / n_reads_mapped  (per-read: 1/n_alignments
                        credit when correct, so ambiguous-but-correct reads
                        count for less than a unique correct read)
    recall             n_reads_hit / n_reads_total     (per-read: binary hit
                        or not, since a read has exactly one true transcript —
                        multimapping doesn't discount it the way it does precision)
    hit_rate           n_reads_hit / n_reads_mapped    (recall restricted to
                        mapped reads only, i.e. ignoring reads that never
                        mapped anywhere)
    avg_n_alignments   n_alignments_sum / n_reads_mapped
    f1                 harmonic mean of precision and recall

Output:
    plots/<name>/<metric>_isoforms.png, one per requested metric
"""

import glob
from pathlib import Path
import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from itertools import cycle

HERE = Path(__file__).resolve().parent
TOOLS_YAML = HERE / "tools.yaml"

PALETTE = [
    "#1E90FF", "#64B5CD", "#D8D8D8", "#AEAEAE",
]

METRIC_LABELS = {
    "precision": "Average per-read precision",
    "recall": "Average per-read recall",
    "f1": "F1",
    "hit_rate": "Hit rate (mapped reads)",
    "avg_n_alignments": "Avg. alignments / mapped read",
}
# Metrics that aren't already a 0-1 fraction shouldn't get a percent axis
NOT_A_FRACTION = {"avg_n_alignments"}

RAW_COLUMNS = ["n_isoforms", "n_reads_total", "n_reads_mapped", "n_reads_hit",
               "precision_sum", "n_alignments_sum"]

OVERFLOW_BIN = "80+"
OVERFLOW_THRESHOLD = 80


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def load_tool_registry() -> dict:
    with open(TOOLS_YAML) as f:
        return yaml.safe_load(f)["tools"]


def load_method_data(name: str, label: str) -> pd.DataFrame:
    """Load every per-sample per-transcript TSV run.py wrote for this tool
    (tsv/<name>/<label>/*.tsv) and pool (sum) by exact n_isoforms, across
    every transcript and every sample. Cheap — it's working off the cached
    BAM-parse results, not the BAM itself."""
    pattern = str(HERE / "tsv" / name / label / "*.tsv")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matched pattern {pattern!r}")

    frames = []
    for path in files:
        df = pd.read_csv(path, sep="\t")
        missing = set(RAW_COLUMNS) - set(df.columns)
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        frames.append(df[RAW_COLUMNS])
    raw = pd.concat(frames, ignore_index=True)

    pooled = (
        raw.groupby("n_isoforms", as_index=False)[RAW_COLUMNS[1:]]
        .sum(min_count=1)
    )
    return pooled, len(files)


def _safe_div(numer: pd.Series, denom: pd.Series) -> pd.Series:
    """numer/denom, with 0/0 and x/0 both left as NaN (renders as a gap in
    the plot rather than a misleading 0)."""
    return (numer / denom).replace([np.inf, -np.inf], np.nan)


def compute_metrics(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["precision"] = _safe_div(df["precision_sum"], df["n_reads_mapped"])
    df["recall"] = _safe_div(df["n_reads_hit"], df["n_reads_total"])
    df["hit_rate"] = _safe_div(df["n_reads_hit"], df["n_reads_mapped"])
    df["avg_n_alignments"] = _safe_div(df["n_alignments_sum"], df["n_reads_mapped"])
    denom = df["precision"] + df["recall"]
    df["f1"] = np.where(denom > 0, 2 * df["precision"] * df["recall"] / denom, np.nan)
    return df


def compute_all_bins(bin_size: int) -> list:
    """Return bin labels from 0 up to (but not including) OVERFLOW_THRESHOLD,
    followed by a single OVERFLOW_BIN label that captures everything above."""
    bin_edges = range(0, OVERFLOW_THRESHOLD, bin_size)
    regular = [f"{lo}-{lo + bin_size - 1}" for lo in bin_edges]
    return regular + [OVERFLOW_BIN]


def bin_dataframe(df: pd.DataFrame, bin_size: int, all_bins: list) -> pd.DataFrame:
    """Aggregate pooled-by-exact-isoform-count rows into n_isoforms bins of
    `bin_size` width, with all loci at or above OVERFLOW_THRESHOLD collapsed
    into a single overflow bin. Raw sums are combined first, then every
    metric is recomputed from the binned totals (a weighted combination, not
    an average of averages). Bins with no data are left as NaN so they
    render as gaps rather than zeros. Reindexes to `all_bins` so every
    method shares the same x-axis.
    """
    df = df.copy()

    bin_edges = list(range(0, OVERFLOW_THRESHOLD + bin_size, bin_size))
    labels = [f"{lo}-{lo + bin_size - 1}" for lo in bin_edges[:-1]]

    df["bin"] = pd.cut(
        df["n_isoforms"].clip(upper=OVERFLOW_THRESHOLD - 1),
        bins=bin_edges,
        labels=labels,
        right=False,
        include_lowest=True,
    ).astype(str)
    df.loc[df["n_isoforms"] >= OVERFLOW_THRESHOLD, "bin"] = OVERFLOW_BIN

    agg = (
        df.groupby("bin", observed=True)[RAW_COLUMNS[1:]]
        .sum(min_count=1)
        .reset_index()
    )
    agg = compute_metrics(agg)

    agg = agg.set_index("bin").reindex(all_bins).reset_index()
    agg.columns = ["bin"] + list(agg.columns[1:])
    return agg


def plot_metric(binned: list, labels: list, display_names: dict, metric: str, name: str, bin_size: int) -> None:
    """Plot one metric by isoform-count bin for each method and save to disk."""
    all_bins = list(binned[0]["bin"])
    n = len(binned)
    colors = [c for c, _ in zip(cycle(PALETTE), range(n))]

    total_group_width = 0.8
    bar_width = total_group_width / n
    x = np.arange(len(all_bins))
    offsets = np.linspace(-(n - 1) / 2, (n - 1) / 2, n) * bar_width

    fig_width = max(8, 0.8 * len(all_bins))
    fig, ax = plt.subplots(figsize=(fig_width, 8))

    for b, label, color, offset in zip(binned, labels, colors, offsets):
        ax.bar(
            x + offset,
            b[metric],
            width=bar_width,
            label=display_names.get(label, label),
            color=color,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.4,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(all_bins, rotation=45, ha="right", fontsize=21)
    ax.set_xlabel("Number of gene isoforms", fontsize=27)
    ax.set_ylabel(METRIC_LABELS.get(metric, metric), fontsize=27)
    max_val = max(b[metric].max(skipna=True) for b in binned)
    if pd.notna(max_val) and max_val > 0:
        ax.set_ylim(0, max_val * 1.05)
    if metric not in NOT_A_FRACTION:
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
    ax.tick_params(axis="y", labelsize=21, length=6, width=1.5)
    ax.legend(fontsize=24, loc="lower left")
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out_dir = HERE / "plots" / name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{metric}_isoforms.png"
    fig.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument(
        "--metrics", nargs="+", default=["precision", "recall", "f1"],
        help="Which metrics to plot: precision, recall, f1, hit_rate, "
             "avg_n_alignments (default: precision recall f1)",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    registry = load_tool_registry()
    name = cfg["name"]
    bin_size = cfg.get("bin_size", 10)
    display_names = {label: t["display_name"] for label, t in registry.items()}

    pooled, labels = [], []
    for label in cfg["tools"]:
        df, n_files = load_method_data(name, label)
        pooled.append(df)
        labels.append(label)
        total_reads = df["n_reads_mapped"].sum()
        has_truth = df["n_reads_total"].notna().any()
        print(
            f"Loaded {label!r} from {n_files} sample file(s): "
            f"{total_reads:,.0f} mapped reads"
            f"{'' if has_truth else ' (no truth counts available — recall will be blank)'}"
        )

    all_bins = compute_all_bins(bin_size)
    binned = [bin_dataframe(df, bin_size, all_bins) for df in pooled]

    print("Plotting...")
    for metric in args.metrics:
        plot_metric(binned, labels, display_names, metric, name, bin_size)


if __name__ == "__main__":
    main()
