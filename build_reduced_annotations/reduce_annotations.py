#!/usr/bin/env python3
"""Generate reduced-annotation GTFs by randomly removing a percentage of genes
or transcripts from a full GTF file, with bootstrapping and a TSV manifest of
every keep/remove decision. Run with --help for CLI usage."""

import argparse
import random
import sys
import os
import csv
from collections import defaultdict
from pathlib import Path
 
 
# ---------------------------------------------------------------------------
# GTF I/O
# ---------------------------------------------------------------------------
 
def parse_gtf_attributes(attr_string: str) -> dict:
    """Parse a GTF attribute string into a {key: value} dict."""
    attrs = {}
    for item in attr_string.strip().rstrip(';').split(';'):
        item = item.strip()
        if item:
            key, _, value = item.partition(' ')
            attrs[key] = value.strip('"')
    return attrs
 
 
def read_gtf(filename: str):
    """
    Read a GTF file.
 
    Returns
    -------
    genes               : dict  gene_id  -> [lines]
    transcripts         : dict  transcript_id -> [lines]
    gene_to_transcripts : dict  gene_id  -> set of transcript_ids
    header_lines        : list  comment / header lines (preserved verbatim)
    """
    genes = defaultdict(list)
    transcripts = defaultdict(list)
    gene_to_transcripts = defaultdict(set)
    header_lines = []
 
    with open(filename, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                header_lines.append(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get('gene_id')
            transcript_id = attrs.get('transcript_id')
 
            if gene_id:
                genes[gene_id].append(line)
                if transcript_id:
                    transcripts[transcript_id].append(line)
                    gene_to_transcripts[gene_id].add(transcript_id)
 
    return genes, transcripts, gene_to_transcripts, header_lines
 
 
# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------
 
def sample_genes(genes: dict, keep_pct: float, seed: int):
    """Return (kept_gene_ids, removed_gene_ids) sets."""
    random.seed(seed)
    gene_ids = sorted(genes.keys())          # sorted for determinism
    n_keep = round(len(gene_ids) * keep_pct / 100)
    kept = set(random.sample(gene_ids, n_keep))
    removed = set(gene_ids) - kept
    return kept, removed
 
 
def sample_transcripts(transcripts: dict, gene_to_transcripts: dict,
                       keep_pct: float, seed: int):
    """
    Return (kept_transcript_ids, removed_transcript_ids,
            kept_gene_ids, removed_gene_ids) sets.
 
    A gene is 'kept' if ≥1 of its transcripts is kept.
    A gene is 'removed' if ALL of its transcripts are removed.
    """
    random.seed(seed)
    tx_ids = sorted(transcripts.keys())
    n_keep = round(len(tx_ids) * keep_pct / 100)
    kept_tx = set(random.sample(tx_ids, n_keep))
    removed_tx = set(tx_ids) - kept_tx
 
    kept_genes = set()
    removed_genes = set()
    for gene_id, tx_set in gene_to_transcripts.items():
        if tx_set & kept_tx:
            kept_genes.add(gene_id)
        else:
            removed_genes.add(gene_id)
 
    return kept_tx, removed_tx, kept_genes, removed_genes
 
 
# ---------------------------------------------------------------------------
# GTF writing
# ---------------------------------------------------------------------------
 
def write_filtered_gtf(input_file: str, output_file: str,
                       kept_genes: set, kept_transcripts: set,
                       header_lines: list):
    """Write a filtered GTF keeping only entries in kept_genes / kept_transcripts."""
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(input_file, 'r') as inf, open(output_file, 'w') as outf:
        # Write preserved headers
        for hl in header_lines:
            outf.write(hl)
        for line in inf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get('gene_id')
            transcript_id = attrs.get('transcript_id')
            if gene_id in kept_genes:
                if not transcript_id or transcript_id in kept_transcripts:
                    outf.write(line)
 
 
# ---------------------------------------------------------------------------
# Manifest helpers
# ---------------------------------------------------------------------------
 
MANIFEST_FIELDS = [
    'mode',          # gene | transcript
    'keep_pct',      # target keep percentage (e.g. 80)
    'replicate',     # bootstrap replicate index (1-based)
    'seed',          # actual random seed used
    'entity_type',   # gene | transcript
    'entity_id',     # gene_id or transcript_id
    'status',        # kept | removed
    'gene_id',       # gene_id (= entity_id for gene mode; parent gene for tx mode)
    'n_transcripts_in_gene',   # total transcripts belonging to this gene
    'n_kept_transcripts_in_gene',  # how many of those were kept
    'output_gtf',    # path to the GTF this record belongs to
    'n_total_genes',
    'n_total_transcripts',
    'n_kept_genes',
    'n_kept_transcripts',
    'actual_kept_gene_pct',
    'actual_kept_tx_pct',
]
 
 
def build_manifest_rows(mode, keep_pct, replicate, seed,
                        genes, transcripts, gene_to_transcripts,
                        kept_genes, removed_genes,
                        kept_transcripts, removed_transcripts,
                        output_gtf):
    """Yield one dict per entity (gene or transcript) for the TSV manifest."""
    n_total_genes = len(genes)
    n_total_tx = len(transcripts)
    n_kept_genes = len(kept_genes)
    n_kept_tx = len(kept_transcripts)
    actual_gene_pct = (n_kept_genes / n_total_genes * 100) if n_total_genes else 0
    actual_tx_pct = (n_kept_tx / n_total_tx * 100) if n_total_tx else 0
 
    # Pre-compute per-gene transcript keep counts
    gene_total_tx = {g: len(txs) for g, txs in gene_to_transcripts.items()}
    gene_kept_tx = {g: len(txs & kept_transcripts)
                    for g, txs in gene_to_transcripts.items()}
 
    common = dict(
        mode=mode,
        keep_pct=keep_pct,
        replicate=replicate,
        seed=seed,
        output_gtf=output_gtf,
        n_total_genes=n_total_genes,
        n_total_transcripts=n_total_tx,
        n_kept_genes=n_kept_genes,
        n_kept_transcripts=n_kept_tx,
        actual_kept_gene_pct=f'{actual_gene_pct:.4f}',
        actual_kept_tx_pct=f'{actual_tx_pct:.4f}',
    )
 
    if mode == 'gene':
        for gene_id in sorted(kept_genes | removed_genes):
            yield {**common,
                   'entity_type': 'gene',
                   'entity_id': gene_id,
                   'status': 'kept' if gene_id in kept_genes else 'removed',
                   'gene_id': gene_id,
                   'n_transcripts_in_gene': gene_total_tx.get(gene_id, 0),
                   'n_kept_transcripts_in_gene': gene_kept_tx.get(gene_id, 0)}
    else:  # transcript
        all_tx = kept_transcripts | removed_transcripts
        # Build transcript->gene map
        tx_to_gene = {}
        for g, txs in gene_to_transcripts.items():
            for t in txs:
                tx_to_gene[t] = g
        for tx_id in sorted(all_tx):
            g = tx_to_gene.get(tx_id, '')
            yield {**common,
                   'entity_type': 'transcript',
                   'entity_id': tx_id,
                   'status': 'kept' if tx_id in kept_transcripts else 'removed',
                   'gene_id': g,
                   'n_transcripts_in_gene': gene_total_tx.get(g, 0),
                   'n_kept_transcripts_in_gene': gene_kept_tx.get(g, 0)}
 
 
# ---------------------------------------------------------------------------
# Core orchestration
# ---------------------------------------------------------------------------
 
def run_bootstrap(input_file, output_dir, mode, percentages,
                  n_bootstrap, seed_base, prefix, manifest_path, verbose):
    """
    Run all (percentage × replicate) combinations.
 
    Seeds are derived deterministically:
        seed = seed_base + replicate_index * 997 + int(percentage * 10)
    This ensures seeds are unique per (pct, rep) tuple but fully reproducible
    given the same seed_base.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
 
    print(f"Reading GTF: {input_file}", file=sys.stderr)
    genes, transcripts, gene_to_transcripts, header_lines = read_gtf(input_file)
    print(f"  {len(genes):,} genes | {len(transcripts):,} transcripts", file=sys.stderr)
 
    manifest_rows = []
    mode_char = 'g' if mode == 'gene' else 't'
 
    total_runs = len(percentages) * n_bootstrap
    run_idx = 0
 
    for pct in percentages:
        for rep in range(1, n_bootstrap + 1):
            run_idx += 1
            seed = seed_base + rep * 997 + int(pct * 10)
 
            # Build output filename  e.g.  chess/g_95_1.gtf  or  g_95_1.gtf
            fname = f"{mode_char}_{int(pct)}_{rep}.gtf"
            if prefix:
                out_subdir = os.path.join(output_dir, prefix)
                Path(out_subdir).mkdir(parents=True, exist_ok=True)
                out_path = os.path.join(out_subdir, fname)
            else:
                out_path = os.path.join(output_dir, fname)
 
            if verbose:
                print(f"[{run_idx}/{total_runs}] mode={mode} pct={pct} "
                      f"rep={rep} seed={seed} -> {out_path}", file=sys.stderr)
 
            # Sample
            if mode == 'gene':
                kept_genes, removed_genes = sample_genes(genes, pct, seed)
                # All transcripts of kept genes are kept
                kept_tx = set()
                for g in kept_genes:
                    kept_tx.update(gene_to_transcripts[g])
                removed_tx = set(transcripts.keys()) - kept_tx
            else:
                kept_tx, removed_tx, kept_genes, removed_genes = \
                    sample_transcripts(transcripts, gene_to_transcripts, pct, seed)
 
            # Write GTF
            write_filtered_gtf(input_file, out_path,
                                kept_genes, kept_tx, header_lines)
 
            # Accumulate manifest rows
            if manifest_path:
                manifest_rows.extend(
                    build_manifest_rows(
                        mode=mode, keep_pct=pct, replicate=rep, seed=seed,
                        genes=genes, transcripts=transcripts,
                        gene_to_transcripts=gene_to_transcripts,
                        kept_genes=kept_genes, removed_genes=removed_genes,
                        kept_transcripts=kept_tx, removed_transcripts=removed_tx,
                        output_gtf=out_path,
                    )
                )
 
    # Write manifest
    if manifest_path and manifest_rows:
        with open(manifest_path, 'w', newline='') as mf:
            writer = csv.DictWriter(mf, fieldnames=MANIFEST_FIELDS,
                                    delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(manifest_rows)
        print(f"\nManifest written: {manifest_path} "
              f"({len(manifest_rows):,} rows)", file=sys.stderr)
 
    print(f"\nDone. {total_runs} GTF files written to: {output_dir}", file=sys.stderr)
 
 
# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
 
def main():
    parser = argparse.ArgumentParser(
        description='Subsample genes or transcripts from a GTF with bootstrapping.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Single run (classic usage)
  %(prog)s -i input.gtf -o output.gtf -m gene -p 80 -s 42

  # Bootstrap 10 replicates, gene mode, output to reduced_annotations/ with chess/ prefix
  %(prog)s -i chess2.2_assembly.gtf -o reduced_annotations/ -m gene \\
      -p 95 90 85 80 75 70 65 60 55 50 \\
      --bootstrap 10 --prefix chess \\
      --manifest chess_gene_manifest.tsv

  # Same for transcripts
  %(prog)s -i chess2.2_assembly.gtf -o reduced_annotations/ -m transcript \\
      -p 95 90 85 80 75 70 65 60 55 50 \\
      --bootstrap 10 --prefix chess \\
      --manifest chess_tx_manifest.tsv
        """,
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input GTF file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output GTF file (single run) or output directory (bootstrap)')
    parser.add_argument('-m', '--mode', choices=['gene', 'transcript'], required=True,
                        help='Subsample at the gene or transcript level')
    parser.add_argument('-p', '--percentage', type=float, nargs='+', required=True,
                        metavar='PCT',
                        help='Percentage(s) to keep (0–100). Multiple values trigger batch mode.')
    parser.add_argument('-s', '--seed', type=int, default=42,
                        help='Base random seed (default: 42). '
                             'In bootstrap mode, per-run seeds are derived from this.')
    parser.add_argument('--bootstrap', type=int, default=1, metavar='N',
                        help='Number of bootstrap replicates per percentage (default: 1)')
    parser.add_argument('--prefix', type=str, default='',
                        help='Optional subdirectory prefix inside output dir '
                             '(e.g. "chess" -> reduced_annotations/chess/g_95_1.gtf)')
    parser.add_argument('--manifest', type=str, default='',
                        metavar='TSV',
                        help='Path for output TSV manifest of all kept/removed entities. '
                             'Only written in batch/bootstrap mode.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress for each run')
    parser.add_argument('--batch', action='store_true',
                        help="Force batch/bootstrap mode (-o is always a directory of "
                             "t_{pct}_{rep}.gtf/g_{pct}_{rep}.gtf files) even when only one "
                             "percentage and --bootstrap 1 are given. Without this, that "
                             "combination is ambiguous with classic single-file mode (-o a "
                             "literal output file) and single-file mode wins.")

    args = parser.parse_args()

    for pct in args.percentage:
        if not 0 <= pct <= 100:
            parser.error(f"Percentage must be between 0 and 100, got {pct}")

    single_run = (len(args.percentage) == 1 and args.bootstrap == 1 and not args.batch)
 
    if single_run:
        # Classic single-file mode
        genes, transcripts, gene_to_transcripts, header_lines = read_gtf(args.input)
        pct = args.percentage[0]
        seed = args.seed
 
        print(f"Input: {len(genes):,} genes | {len(transcripts):,} transcripts",
              file=sys.stderr)
 
        if args.mode == 'gene':
            kept_genes, removed_genes = sample_genes(genes, pct, seed)
            kept_tx = set()
            for g in kept_genes:
                kept_tx.update(gene_to_transcripts[g])
            removed_tx = set(transcripts.keys()) - kept_tx
        else:
            kept_tx, removed_tx, kept_genes, removed_genes = \
                sample_transcripts(transcripts, gene_to_transcripts, pct, seed)
 
        print(f"Kept {len(kept_genes):,} / {len(genes):,} genes "
              f"({len(kept_genes)/len(genes)*100:.1f}%)", file=sys.stderr)
        print(f"Kept {len(kept_tx):,} / {len(transcripts):,} transcripts "
              f"({len(kept_tx)/len(transcripts)*100:.1f}%)", file=sys.stderr)
 
        write_filtered_gtf(args.input, args.output,
                           kept_genes, kept_tx, header_lines)
        print(f"Output: {args.output}", file=sys.stderr)
 
        if args.manifest:
            rows = list(build_manifest_rows(
                mode=args.mode, keep_pct=pct, replicate=1, seed=seed,
                genes=genes, transcripts=transcripts,
                gene_to_transcripts=gene_to_transcripts,
                kept_genes=kept_genes, removed_genes=removed_genes,
                kept_transcripts=kept_tx, removed_transcripts=removed_tx,
                output_gtf=args.output,
            ))
            with open(args.manifest, 'w', newline='') as mf:
                writer = csv.DictWriter(mf, fieldnames=MANIFEST_FIELDS,
                                        delimiter='\t', extrasaction='ignore')
                writer.writeheader()
                writer.writerows(rows)
            print(f"Manifest: {args.manifest} ({len(rows):,} rows)", file=sys.stderr)
    else:
        # Bootstrap / batch mode
        run_bootstrap(
            input_file=args.input,
            output_dir=args.output,
            mode=args.mode,
            percentages=args.percentage,
            n_bootstrap=args.bootstrap,
            seed_base=args.seed,
            prefix=args.prefix,
            manifest_path=args.manifest,
            verbose=args.verbose,
        )
 
 
if __name__ == '__main__':
    main()