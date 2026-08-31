#!/usr/bin/env bash
# scripts/HISAT2_index.sh: build the one shared HISAT2 genome index.
# Annotation-blind — no splice-sites file is extracted or used, so this same
# index is reused as-is by short-read_reduced_annotation too.
set -euo pipefail
source ./config.sh

cmd_HISAT2_build=($(cmd_in_env env_HISAT2 "${HISAT2_bin}-build"))

if [[ ! -f "${HISAT2_index}.1.ht2" && -f "$genome_fa" ]]; then
    echo "Building HISAT2 index..."
    mkdir -p "$(dirname "$HISAT2_index")"
    "${cmd_HISAT2_build[@]}" -p "$threads" "$genome_fa" "$HISAT2_index"
fi
