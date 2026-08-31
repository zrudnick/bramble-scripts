#!/usr/bin/env bash
# build_reduced_annotation_indices.sh: build every per-percentage/bootstrap
# reduced-annotation index used across the three *_reduced_annotation
# pipelines — STAR and Bowtie2/Salmon indices (short-read_reduced_annotation),
# plus minimap2 transcriptome (ont_reduced_annotation/pacbio_reduced_annotation).
# HISAT2 has no reduced-annotation index of its own (see scripts/HISAT2.sh —
# it's annotation-blind and reuses the shared genome index as-is). Run once
# from the repo root before running any *_reduced_annotation pipeline's
# pipeline.sh.
#
# Requires build_indices.sh to have already been run — the shared HISAT2
# index and shared minimap2 genome index are reused as-is here, never
# rebuilt per combination. Also requires
# ../build_reduced_annotations/generate_annotations.sh and
# ../build_reduced_annotations/generate_transcripts.sh to have produced the
# per-combination GTFs/transcript FASTAs first (each step below skips, with
# a warning, any combination that isn't ready yet).
set -euo pipefail
source ./config.sh

echo "============================================================"
echo " Building reduced-annotation indices (${#percentages[@]} percentages x ${#bootstraps[@]} bootstraps)"
echo "============================================================"

bash scripts/STAR_reduced_annotation_index.sh
bash scripts/bowtie2_reduced_annotation_index.sh
bash scripts/salmon_reduced_annotation_index.sh
bash scripts/minimap2_transcriptome_reduced_annotation_index.sh

echo "============================================================"
echo " All reduced-annotation indices built successfully!"
echo "============================================================"
