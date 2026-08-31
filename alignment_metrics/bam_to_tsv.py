#!/usr/bin/env python3
"""
bam_to_tsv.py: parse one BAM (expensive) into a per-transcript TSV cache
(cheap to re-read). Every column downstream metrics might ever need is kept
at full per-transcript granularity — binning/averaging happens later, in
plot_tp.py, off the cache, so a new question never requires re-parsing a BAM.
"""

import re
import gzip
import argparse
import pysam
from pathlib import Path
from collections import defaultdict

# Which read_format (see parse_true_transcript) each simulated platform's
# reads use — the single source of truth run.py and this file's own CLI
# both draw from, so "which simulator produced this platform's reads" is
# never duplicated per-analysis.
PLATFORM_READ_FORMAT = {"short": 1, "pacbio": 2, "ont": 3}


def parse_true_transcript(read_name, read_format):
    """Extract the true transcript ID from a read name."""
    if read_format == 1:
        # e.g. read349/CHS.2.3;mate1:2892-2992
        return read_name.split("/")[1].split(";")[0]
    elif read_format == 2:
        # e.g. CHS.2.3_PBSIM_simulated_read_6
        return read_name.rsplit("_PBSIM_", 1)[0]
    elif read_format == 3:
        # e.g. CHS.2.3_ONT_simulated_read_6
        return read_name.rsplit("_ONT_", 1)[0]
    else:
        raise ValueError(f"Unknown read_format: {read_format}")


def build_isoform_counts(gtf_path):
    """Return {transcript_id: n_isoforms} by counting transcripts per gene."""
    gene_to_transcripts = defaultdict(set)
    transcript_to_gene = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]

            gene_match = re.search(r'gene_id "([^"]+)"', attrs)
            tx_match = re.search(r'transcript_id "([^"]+)"', attrs)
            if not gene_match or not tx_match:
                continue

            gene_id = gene_match.group(1)
            tx_id = tx_match.group(1)
            gene_to_transcripts[gene_id].add(tx_id)
            transcript_to_gene[tx_id] = gene_id

    return {
        tx: len(gene_to_transcripts[gene])
        for tx, gene in transcript_to_gene.items()
    }


def _open_maybe_gzip(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _is_fasta(path):
    return str(path).rstrip(".gz").endswith((".fa", ".fasta"))


def count_true_totals(reads_path, read_format):
    """Return {transcript_id: n_reads_total} by counting read names directly
    in the FASTQ/FASTA that was fed to the aligner — independent of any
    aligner's behavior, since every simulated read has a record here whether
    or not it ever mapped anywhere."""
    counts = defaultdict(int)
    fasta = _is_fasta(reads_path)

    with _open_maybe_gzip(reads_path) as f:
        if fasta:
            for line in f:
                if not line.startswith(">"):
                    continue
                read_name = line[1:].strip().split()[0]
                if "ALL" in read_name:
                    continue
                counts[parse_true_transcript(read_name, read_format)] += 1
        else:
            for i, line in enumerate(f):
                if i % 4 != 0:
                    continue
                read_name = line[1:].strip().split()[0]
                if "ALL" in read_name:
                    continue
                counts[parse_true_transcript(read_name, read_format)] += 1

    return counts


def write_truth_tsv(true_totals, isoform_counts, out_path):
    with open(out_path, "w") as f:
        f.write("transcript_id\tn_isoforms\tn_reads_total\n")
        for tx in sorted(true_totals):
            if tx not in isoform_counts:
                continue
            f.write(f"{tx}\t{isoform_counts[tx]}\t{true_totals[tx]}\n")


TRUTH_HEADER = "transcript_id\tn_isoforms\tn_reads_total\n"


def load_truth_tsv(path):
    totals = {}
    with open(path) as f:
        header = f.readline()
        if header != TRUTH_HEADER:
            return None  # stale/foreign schema — caller should regenerate
        for line in f:
            tx, _n_iso, n_total = line.rstrip("\n").split("\t")
            totals[tx] = int(n_total)
    return totals


def process_bam(bam_path, read_format, isoform_counts):
    """For each read, record every distinct transcript it aligned to."""
    read_hits = defaultdict(set)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam:
            if aln.is_unmapped or "ALL" in aln.query_name:
                continue
            read_hits[aln.query_name].add(aln.reference_name)

    # true_tx -> {n_reads_mapped, n_reads_hit, precision_sum, n_alignments_sum}
    results = defaultdict(lambda: {
        "n_reads_mapped": 0,
        "n_reads_hit": 0,
        "precision_sum": 0.0,
        "n_alignments_sum": 0,
    })
    for read_name, aligned_refs in read_hits.items():
        true_tx = parse_true_transcript(read_name, read_format)
        if true_tx not in isoform_counts:
            continue

        n_alignments = len(aligned_refs)
        hit = true_tx in aligned_refs
        score = (1 / n_alignments) if hit else 0.0

        r = results[true_tx]
        r["n_reads_mapped"] += 1
        r["n_reads_hit"] += 1 if hit else 0
        r["precision_sum"] += score
        r["n_alignments_sum"] += n_alignments

    return results


TSV_HEADER = (
    "transcript_id\tn_isoforms\tn_reads_total\tn_reads_mapped\t"
    "n_reads_hit\tprecision_sum\tn_alignments_sum\n"
)


def write_tsv(results, isoform_counts, true_totals, out_path):
    all_tx = set(results) | (set(true_totals) if true_totals else set())
    empty = {"n_reads_mapped": 0, "n_reads_hit": 0, "precision_sum": 0.0, "n_alignments_sum": 0}

    with open(out_path, "w") as f:
        f.write(TSV_HEADER)
        for tx in sorted(all_tx):
            if tx not in isoform_counts:
                continue
            r = results.get(tx, empty)
            n_total = true_totals.get(tx, "") if true_totals is not None else ""
            f.write(
                f"{tx}\t{isoform_counts[tx]}\t{n_total}\t"
                f"{r['n_reads_mapped']}\t{r['n_reads_hit']}\t"
                f"{r['precision_sum']:.6f}\t{r['n_alignments_sum']}\n"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--out", required=True, help="Output TSV file")
    parser.add_argument("--gtf", required=True, help="Full annotation GTF")
    fmt_group = parser.add_mutually_exclusive_group(required=True)
    fmt_group.add_argument("--platform", choices=sorted(PLATFORM_READ_FORMAT),
                            help="Simulated platform this BAM's reads came from "
                                 "(determines read_format automatically)")
    fmt_group.add_argument("--read-format", type=int, choices=[1, 2, 3],
                            help="Read-name format directly, if not using --platform")
    parser.add_argument("--fastq", help="This sample's FASTQ/FASTA (for true read totals)")
    parser.add_argument("--truth-tsv", help="Path to read/write the per-sample truth-count cache")
    args = parser.parse_args()

    gtf_path = args.gtf
    read_fmt = args.read_format or PLATFORM_READ_FORMAT[args.platform]

    print(f"Building isoform counts from {gtf_path!r}...")
    isoform_counts = build_isoform_counts(gtf_path)
    print(f"  {len(isoform_counts):,} transcripts loaded.")

    true_totals = None
    if args.truth_tsv:
        truth_path = Path(args.truth_tsv)
        true_totals = load_truth_tsv(truth_path) if truth_path.exists() else None
        if true_totals is not None:
            print(f"Loaded cached truth counts from {truth_path!r} ({len(true_totals):,} transcripts).")
        elif args.fastq:
            print(f"Counting true read totals from {args.fastq!r}...")
            true_totals = count_true_totals(args.fastq, read_fmt)
            truth_path.parent.mkdir(parents=True, exist_ok=True)
            write_truth_tsv(true_totals, isoform_counts, truth_path)
            print(f"  {len(true_totals):,} transcripts. Cached -> {truth_path!r}")
        else:
            print("WARNING: no --fastq given and no cached --truth-tsv found; n_reads_total will be blank.")

    print(f"Processing {args.bam!r} (read_format={read_fmt})...")
    results = process_bam(args.bam, read_fmt, isoform_counts)
    print(f"  {len(results):,} transcripts observed.")

    write_tsv(results, isoform_counts, true_totals, args.out)
    print(f"Written -> {args.out!r}")


if __name__ == "__main__":
    main()
