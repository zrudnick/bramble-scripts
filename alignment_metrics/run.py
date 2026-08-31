#!/usr/bin/env python3
"""
run.py: end-to-end alignment-precision analysis for one platform's
quantification_benchmarking output.

Usage:
    python3 run.py --config config/short_read.yaml
    python3 run.py --config config/ont.yaml --metrics precision recall f1 hit_rate

Driven entirely by:
  - the analysis config (config/*.yaml): platform, where that platform's
    quantification_benchmarking pipeline wrote its BAMs (results_dir), and
    which tools (see tools.yaml) to compute precision/recall for
  - tools.yaml: the registry of what each tool's BAM is called
  - the repo-root config.sh: gtf_file, data_short_dir/data_ont_dir/
    data_pacbio_dir, fastq_extensions — sourced directly so these stay
    defined in exactly one place rather than duplicated into Python

An analysis config no longer hardcodes bam_pattern/read_format per tool —
this is quantification_benchmarking's own downstream analysis step, so it
follows quantification_benchmarking's actual output layout and BAM naming
(tools.yaml's bam_name templates mirror config.sh's *_bam_name() functions).

For each sample directory found under results_dir (any t<N>_s<N> directory
— no fixed tissue/sample grid assumed) and each requested tool, this:
  1. builds that tool's BAM path from tools.yaml's bam_name template
  2. parses it into a per-transcript TSV cache via bam_to_tsv.py's functions
     (skipping BAMs already cached under the current TSV schema)
  3. caches true per-transcript read totals per sample (shared across tools)
Then hands off to plot_tp.py to pool everything and render the plots.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import yaml

import bam_to_tsv

HERE = Path(__file__).resolve().parent
ROOT_CONFIG_SH = HERE.parent / "config.sh"
TOOLS_YAML = HERE / "tools.yaml"

PLATFORM_FASTQ_VAR = {
    "short": "data_short_dir",
    "ont": "data_ont_dir",
    "pacbio": "data_pacbio_dir",
}

SAMPLE_DIR_RE = re.compile(r"t(\d+)_s(\d+)")


def load_yaml(path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def source_root_config_scalars(*names: str) -> dict:
    """Source the repo-root config.sh and read back the requested scalar
    shell variables, so gtf_file/data_*_dir stay defined in exactly one
    place instead of being duplicated into this script."""
    script = f"source '{ROOT_CONFIG_SH}' >/dev/null 2>&1\n"
    script += "\n".join(f'echo "${{{n}}}"' for n in names)
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    values = out.stdout.splitlines()
    if len(values) != len(names):
        raise SystemExit(
            f"Expected {len(names)} value(s) from config.sh ({names}), got {values!r} "
            f"— check that every name is actually set there."
        )
    return dict(zip(names, values))


def source_root_config_array(name: str) -> list:
    script = f"source '{ROOT_CONFIG_SH}' >/dev/null 2>&1\necho \"${{{name}[@]}}\""
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    return out.stdout.split()


def discover_samples(results_dir: Path) -> list:
    """Every t<N>_s<N> subdirectory under results_dir — no fixed
    tissue/sample grid assumed, so this adapts to however many samples
    quantification_benchmarking actually ran."""
    return sorted(
        (d.name for d in results_dir.iterdir() if d.is_dir() and re.fullmatch(r"t\d+_s\d+", d.name)),
        key=lambda s: tuple(int(x) for x in SAMPLE_DIR_RE.match(s).groups()),
    )


def build_fastq_index(fastq_dir: Path, extensions: list) -> dict:
    """Map (tissue, sample) -> path for every FASTQ/FASTA anywhere under
    fastq_dir whose filename contains t<N>_s<N> (as a substring, not
    necessarily the whole filename — matches quantification_benchmarking's
    own pipeline.sh discovery convention). Built once and reused across
    every tool/sample lookup, since fastq_dir may be on network storage."""
    index = {}
    for p in sorted(fastq_dir.rglob("*")):
        if not p.is_file() or not any(p.name.endswith("." + ext) for ext in extensions):
            continue
        m = SAMPLE_DIR_RE.search(p.name)
        if m:
            key = (int(m.group(1)), int(m.group(2)))
            index.setdefault(key, p)  # first match wins
    return index


def is_cached(out_tsv: Path) -> bool:
    """True if out_tsv exists and matches bam_to_tsv.py's current schema —
    stale (old-format) caches get silently regenerated instead of crashing
    plot_tp.py downstream."""
    if not out_tsv.exists():
        return False
    with open(out_tsv) as f:
        header = f.readline()
    return header == bam_to_tsv.TSV_HEADER


def process_config(config_path: str) -> dict:
    """Run the BAM -> TSV step for every requested tool/sample. Returns the
    loaded config (for run() to hand to plot_tp.py's --config)."""
    cfg = load_yaml(config_path)
    registry = load_yaml(TOOLS_YAML)["tools"]

    name = cfg["name"]
    platform = cfg["platform"]
    if platform not in PLATFORM_FASTQ_VAR:
        raise SystemExit(f"config.platform must be one of {sorted(PLATFORM_FASTQ_VAR)}, got {platform!r}")
    results_dir = (Path(config_path).resolve().parent / cfg["results_dir"]).resolve()
    tool_labels = cfg["tools"]

    unknown = [t for t in tool_labels if t not in registry]
    if unknown:
        raise SystemExit(f"Unknown tool label(s) {unknown} in {config_path} — see tools.yaml")
    wrong_platform = [t for t in tool_labels if platform not in registry[t]["platforms"]]
    if wrong_platform:
        raise SystemExit(
            f"Tool(s) {wrong_platform} aren't valid for platform {platform!r} "
            f"(see each tool's 'platforms' list in tools.yaml)"
        )

    if not results_dir.is_dir():
        raise SystemExit(
            f"results_dir {results_dir} does not exist — has quantification_benchmarking's "
            f"'{platform}' pipeline been run yet?"
        )
    samples = discover_samples(results_dir)
    if not samples:
        raise SystemExit(f"No t<N>_s<N> sample directories found under {results_dir}")
    print(f"=== {name} ({platform}) ===")
    print(f"Found {len(samples)} sample(s) under {results_dir}: {', '.join(samples)}")

    shell_vars = source_root_config_scalars("gtf_file", PLATFORM_FASTQ_VAR[platform])
    gtf_path = shell_vars["gtf_file"]
    fastq_dir = Path(shell_vars[PLATFORM_FASTQ_VAR[platform]])
    fastq_extensions = source_root_config_array("fastq_extensions")

    print(f"Building isoform counts from {gtf_path!r}...")
    isoform_counts = bam_to_tsv.build_isoform_counts(gtf_path)
    print(f"  {len(isoform_counts):,} transcripts loaded.")

    read_fmt = bam_to_tsv.PLATFORM_READ_FORMAT[platform]
    fastq_index = None  # built lazily, on first sample that actually needs it

    for label in tool_labels:
        tool = registry[label]
        print(f"--- {name} / {label} ---")
        out_dir = HERE / "tsv" / name / label
        truth_dir = HERE / "tsv" / name / "truth"
        out_dir.mkdir(parents=True, exist_ok=True)
        truth_dir.mkdir(parents=True, exist_ok=True)

        for sample_id in samples:
            bam_path = results_dir / sample_id / tool["bam_name"].format(sample_id=sample_id)
            out_tsv = out_dir / f"{sample_id}.tsv"
            truth_tsv = truth_dir / f"{sample_id}.tsv"

            if is_cached(out_tsv):
                print(f"  [skip] {out_tsv} already cached")
                continue
            if not bam_path.exists():
                print(f"  WARNING: {bam_path} not found — skipping {sample_id}")
                continue

            true_totals = bam_to_tsv.load_truth_tsv(truth_tsv) if truth_tsv.exists() else None
            if true_totals is None:
                if fastq_index is None:
                    fastq_index = build_fastq_index(fastq_dir, fastq_extensions)
                t, s = (int(x) for x in SAMPLE_DIR_RE.match(sample_id).groups())
                fastq = fastq_index.get((t, s))
                if fastq is None:
                    print(f"  WARNING: no FASTQ/FASTA found for {sample_id} under {fastq_dir} "
                          f"— n_reads_total will be blank")
                else:
                    print(f"  Counting true read totals from {fastq}...")
                    true_totals = bam_to_tsv.count_true_totals(fastq, read_fmt)
                    bam_to_tsv.write_truth_tsv(true_totals, isoform_counts, truth_tsv)

            print(f"  Processing {bam_path} -> {out_tsv}")
            results = bam_to_tsv.process_bam(str(bam_path), read_fmt, isoform_counts)
            bam_to_tsv.write_tsv(results, isoform_counts, true_totals, out_tsv)

    return cfg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to an analysis config (see config/*.yaml)")
    parser.add_argument(
        "--metrics", nargs="+", default=["precision", "recall", "f1"],
        help="Which metrics to plot: precision, recall, f1, hit_rate, "
             "avg_n_alignments (default: precision recall f1)",
    )
    args = parser.parse_args()

    process_config(args.config)

    print("\nRunning plot_tp.py...")
    subprocess.run(
        [sys.executable, str(HERE / "plot_tp.py"), "--config", args.config, "--metrics", *args.metrics],
        check=True,
    )
    print("Done.")


if __name__ == "__main__":
    main()
