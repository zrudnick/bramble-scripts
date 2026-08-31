#!/usr/bin/env bash
# scripts/bowtie2_index.sh: build the one shared Bowtie2 transcriptome index
# from the full annotation, used by short-read in base mode.
set -euo pipefail
source ./config.sh

cmd_bowtie2_build=($(cmd_in_env env_bowtie2 "${bowtie2_bin}-build"))

if [[ -f "$transcript_fa" && ! -f "${bowtie2_index}.1.bt2" ]]; then
    echo "Building Bowtie2 transcriptome index..."
    mkdir -p "$(dirname "$bowtie2_index")"
    "${cmd_bowtie2_build[@]}" -f --threads "$threads" \
        "$transcript_fa" \
        "$bowtie2_index"
fi
