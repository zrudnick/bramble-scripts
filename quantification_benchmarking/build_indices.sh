#!/usr/bin/env bash
# build_indices.sh: build every BASE (full-annotation / annotation-blind)
# index used across all six benchmarking pipelines — HISAT2, STAR, Bowtie2,
# and Salmon (short-read), plus minimap2 genome + transcriptome (ont/pacbio).
# Run once from the repo root before running any pipeline's pipeline.sh; each
# step below is idempotent (skips if its output already exists), so building
# indices for tools a particular run won't use is harmless.
#
# Run this before build_reduced_annotation_indices.sh — the reduced-annotation
# HISAT2 and minimap2-genome steps reuse the shared indices built here.
set -euo pipefail
source ./config.sh

echo "============================================================"
echo " Building base indices"
echo "============================================================"

mkdir -p "$index_dir"

bash scripts/HISAT2_index.sh
bash scripts/STAR_index.sh
bash scripts/bowtie2_index.sh
bash scripts/salmon_index.sh
bash scripts/minimap2_genome_index.sh
bash scripts/minimap2_transcriptome_index.sh

echo "============================================================"
echo " All base indices built successfully!"
echo "============================================================"
