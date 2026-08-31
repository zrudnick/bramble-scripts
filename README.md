# bramble-scripts

Benchmarking pipelines that compare [Bramble](https://github.com/zrudnick/bramble)'s
genomic-to-transcriptomic BAM projection against native aligners/quantifiers, across short-read
and long-read (ONT, PacBio) data, at both full and reduced (subsampled) annotation completeness.

## Overview

Reproducing an experiment is four steps, run in order:

1. **[Set up](#1-set-up)** — point `config.sh` at your tools, add your reference files, and drop
   in your reads.
2. **[Build indices](#2-build-indices)** — one-time index build for every aligner/quantifier.
3. **[Run a benchmarking pipeline](#3-run-a-benchmarking-pipeline)** — align, project through
   Bramble, and quantify one platform's reads (short-read, ONT, or PacBio), optionally repeated
   across reduced annotations.
4. **[Analyze the results](#4-analyze-the-results)** — score quantification accuracy against
   ground truth, and/or alignment precision/recall.

## Repository layout

```
bramble_scripts/
├── config.sh                        shared configuration — edit this first
├── data/
│   ├── ref/                         reference genome/annotation (you provide this)
│   ├── short-read/
│   │   ├── FASTQ/                   short-read FASTQ files (you provide this)
│   │   └── truth/                   ground-truth per-transcript counts (for step 4)
│   └── long-read/
│       ├── FASTQ-ont/               ONT FASTQ files (you provide this)
│       ├── FASTQ-pacbio/            PacBio FASTA files (you provide this)
│       └── truth/                   ground-truth per-transcript counts (for step 4)
├── quantification_benchmarking/     the six benchmark pipelines
├── build_reduced_annotations/       generates the reduced annotations the *_reduced_annotation pipelines use
├── alignment_metrics/               alignment precision/recall analysis
└── README.md
```

```
quantification_benchmarking/
├── build_indices.sh                 step 2: base indices
├── build_reduced_annotation_indices.sh   step 2: reduced-annotation indices
├── scripts/                         one script per tool, shared across pipelines
├── short-read/                      short-read benchmark (HISAT2, STAR, Bowtie2, Bramble, Salmon)
├── ont/                             ONT long-read benchmark (minimap2, Bramble, Oarfish, TranSigner)
├── pacbio/                          PacBio long-read benchmark (same tools as ont/, different presets)
├── short-read_reduced_annotation/   short-read benchmark repeated across reduced annotations
├── ont_reduced_annotation/          ont/ repeated across reduced annotations
└── pacbio_reduced_annotation/       pacbio/ repeated across reduced annotations
```

---

## 1. Set up

### `config.sh`

Open `config.sh` at the repo root and set:

- **`env_*`** — the conda/mamba environment each tool runs in (e.g. `env_STAR="star_env"`).
  Leave a value as `""` to run that tool in whatever environment is already active.
- **`*_bin`** — the command or absolute path for each tool binary (`bramble_bin`, `HISAT2_bin`,
  `STAR_bin`, `bowtie2_bin`, `salmon_bin`, `minimap2_bin`, `oarfish_bin`, `transigner_bin`,
  `gffread_bin`, `samtools_bin`). These default to paths on the machine this repo was developed
  on — you will need to change most of them.
- **`threads`** — CPU threads for every tool invocation.

### Reference files

Place your reference genome and annotation in `data/ref/`, matching the filenames `config.sh`
expects (or edit `config.sh` to match your own filenames):

| `config.sh` variable | What it is |
|---|---|
| `genome_fa` | Reference genome FASTA |
| `gtf_file` | Full annotation GTF |
| `transcript_fa` | Transcript sequences extracted from `genome_fa`/`gtf_file` |
| `gentrome_fa` | `transcript_fa` + `genome_fa` concatenated, for Salmon's decoy-aware index |
| `decoys_txt` | Genome sequence names, for Salmon's decoy-aware index |

### Reads

Drop your FASTQ/FASTA files into the matching folder:

- Short-read: `data/short-read/FASTQ/`
- ONT: `data/long-read/FASTQ-ont/`
- PacBio: `data/long-read/FASTQ-pacbio/`

Every file just needs `t{tissue}_s{sample}` somewhere in its filename (e.g. `t0_s1.fastq`) — no
other naming convention or subdirectory structure is required. Accepted extensions: `.fastq`,
`.fastq.gz`, `.fq`, `.fq.gz`, `.fasta`, `.fasta.gz`, `.fa`, `.fa.gz`.

### Ground truth (optional — only needed for step 4's accuracy analysis)

Drop one ground-truth file per sample into `data/short-read/truth/` and `data/long-read/truth/`,
named `t{tissue}_s{sample}_true_counts.tsv` (short-read) or
`t{tissue}_s{sample}_transcript_info.tsv` (long-read, shared by ONT and PacBio). See each
pipeline's `analysis/<pipeline>.yaml` (`ground_truth:` section) for the expected column names.

---

## 2. Build indices

```bash
cd quantification_benchmarking
bash build_indices.sh                      # base (full-annotation) indices — always run this
bash build_reduced_annotation_indices.sh   # only if you're running a *_reduced_annotation
                                            # pipeline — see "Reduced-annotation pipelines" below,
                                            # you'll need to generate the reduced annotations first
```

Every step skips work that's already done, so it's safe to re-run these at any time.

---

## 3. Run a benchmarking pipeline

Every pipeline must be run from **inside its own `<pipeline>/benchmarking/` directory**:

```bash
cd quantification_benchmarking/short-read/benchmarking   # or ont/, pacbio/, or any *_reduced_annotation
bash pipeline.sh
```

`pipeline.sh` auto-discovers every FASTQ/FASTA file in the platform's input folder and processes
each one. Wall-clock time and peak memory for every tool invocation are logged to
`benchmarking/benchmark_log.tsv`.

### `short-read/`

Compares HISAT2, STAR, and Bowtie2 alignment (direct or via Bramble projection), quantified with
Salmon (quasi-mapping and alignment-based). Bowtie2 runs in two modes against the same
full-annotation transcriptome index — default and high-sensitivity (`-k 76`) — each with its own
Salmon quantification.

Outputs land in `benchmarking/results/<sample_id>/` — genomic and transcriptomic BAMs, tool logs,
and `salmon/` quant output.

### `ont/` and `pacbio/`

Compares minimap2 alignment (transcriptome, two presets) against Bramble's long-read-mode
genomic-to-transcriptomic projection, quantified with Oarfish and TranSigner. `pacbio/` runs the
same steps as `ont/` — only the minimap2 presets differ.

Outputs land in `benchmarking/results/<sample_id>/` — minimap2 BAMs, the Bramble-projected BAM,
and `oarfish/`/`transigner/` quant output.

### Reduced-annotation pipelines

`short-read_reduced_annotation/`, `ont_reduced_annotation/`, and `pacbio_reduced_annotation/`
repeat their non-reduced counterpart's benchmark across a grid of **reduced annotations** — GTFs
with only a percentage of transcripts kept, at several independent bootstrap replicates per
percentage — to measure how projection/quantification accuracy degrades as annotation
completeness drops. The grid (default: `80 60 40 20` percent kept, replicates `1`–`5`) is set in
`config.sh` (`percentages`, `bootstraps`).

**Step 0 — generate the reduced annotations** (once, shared by all three pipelines):

```bash
cd build_reduced_annotations
bash generate_annotations.sh    # writes GTFs (e.g. t_80_1.gtf) + a manifest TSV to data/ref/reduced_annotations/
bash generate_transcripts.sh    # extracts a matching transcript FASTA per GTF, via gffread
```

Run either with `--help` to see all options. Sampling is seeded (`seed` in `config.sh`, default
`42`), so re-running with an unchanged `config.sh` always reproduces the same set of reduced
annotations — don't override `--seed` unless you want a different random draw.

**Step 1 — build the reduced-annotation indices, then run the pipeline:**

```bash
cd quantification_benchmarking
bash build_reduced_annotation_indices.sh   # skips (with a warning) any combination Step 0 hasn't produced yet

cd short-read_reduced_annotation/benchmarking   # or ont_reduced_annotation/, pacbio_reduced_annotation/
bash pipeline.sh
```

Outputs land in `benchmarking/results/<sample_id>/pct<percentage>/rep<bootstrap>/` — e.g.
`results/t0_s0/pct80/rep1/`.

Note: HISAT2's alignment doesn't use the annotation, so it only runs once per sample (its output
is shared across every percentage/bootstrap combination); STAR's index bakes in splice junctions,
so it's rebuilt per combination. Both feed into their own per-combination Bramble projection.

---

## 4. Analyze the results

Two independent, optional analyses read the results from step 3.

### Quantification accuracy — `quantanalysis`

[`quantanalysis`](https://github.com/zrudnick/quantanalysis) (a separate tool, install it
separately) compares each pipeline's quant output against the ground truth in
`data/short-read/truth/`/`data/long-read/truth/` and produces accuracy plots. Each pipeline has
its own config at `quantification_benchmarking/<pipeline>/analysis/<pipeline>.yaml`, already
pointed at that pipeline's own results and the ground-truth files.

**Run every `quantanalysis` command from inside that pipeline's own `analysis/` directory** — its
config paths are relative to that location:

```bash
cd quantification_benchmarking/short-read/analysis   # or ont/analysis, pacbio/analysis
quantanalysis analysis --config short-read.yaml       # run this first
quantanalysis swarm --config short-read.yaml
quantanalysis scatter --config short-read.yaml
quantanalysis detect-errors --config short-read.yaml
```

For the three `*_reduced_annotation/` pipelines, use the grouped-mode commands instead:

```bash
cd quantification_benchmarking/short-read_reduced_annotation/analysis   # or ont_/pacbio_reduced_annotation/analysis
quantanalysis analysis --config short-read_reduced_annotation.yaml      # run this first
quantanalysis group-scatter --config short-read_reduced_annotation.yaml
quantanalysis group-box-plot --config short-read_reduced_annotation.yaml
quantanalysis group-detect-errors --config short-read_reduced_annotation.yaml
```

`analysis` must run first — it's the only command that reads the raw quant/ground-truth files.
Every other command reads back what `analysis` wrote to `results/analysis/` and renders plots
into `results/figures/` (both next to the config).

### Alignment precision/recall — `alignment_metrics/`

For each simulated read, does the aligner's BAM assign it back to the correct transcript?
Computed per-read, pooled across samples, and plotted by isoform-count bin. Run the platform's
pipeline (step 3) first — this only reads the BAMs already produced, it doesn't align anything
itself. Only base (non-reduced-annotation) pipelines are supported here.

Each platform has a working config already. Open `config/<platform>.yaml` and edit its `tools:`
list to whichever tools (see `tools.yaml` for the full list) you want plotted:

```yaml
# config/short_read.yaml
tools:
  - bramble_hisat2
  - bramble_star
  - bowtie2
  - bowtie2_k76
  - star_transcriptome
```

Then run it:

```bash
cd alignment_metrics
python3 run.py --config config/short_read.yaml   # or config/ont.yaml, config/pacbio.yaml
```

This writes `plots/<name>/{precision,recall,f1}_isoforms.png` (`--metrics` picks a different set:
also available are `hit_rate` and `avg_n_alignments`).
