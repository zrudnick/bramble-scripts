#!/usr/bin/env bash
# scripts/STAR_index.sh: build the one shared STAR genome index from the full
# annotation, used by short-read in base mode.
set -euo pipefail
source ./config.sh

cmd_STAR=($(cmd_in_env env_STAR "$STAR_bin"))

if [[ ! -d "$STAR_index" && -f "$genome_fa" && -f "$gtf_file" ]]; then
    echo "Building STAR genome index..."
    mkdir -p "$STAR_index"
    "${cmd_STAR[@]}" --runThreadN "$threads" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_index" \
        --genomeFastaFiles "$genome_fa" \
        --sjdbGTFfile "$gtf_file" \
        --sjdbOverhang 100
fi
