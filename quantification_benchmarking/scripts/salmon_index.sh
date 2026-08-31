#!/usr/bin/env bash
# scripts/salmon_index.sh: build the shared Salmon standard index (from the
# full transcriptome) and decoy-aware index (transcriptome + genome), used by
# short-read in base mode.
set -euo pipefail
source ./config.sh

cmd_salmon=($(cmd_in_env env_salmon "$salmon_bin"))

if [[ -f "$transcript_fa" && ! -d "$salmon_std_index" ]]; then
    echo "Building Salmon standard index..."
    mkdir -p "$salmon_std_index"
    "${cmd_salmon[@]}" index --threads "$threads" \
        -t "$transcript_fa" \
        -i "$salmon_std_index"
fi

if [[ -f "$transcript_fa" && -f "$genome_fa" && ! -d "$salmon_decoy_index" ]]; then
    echo "Generating decoys & gentrome for Salmon decoy index..."
    mkdir -p "$salmon_decoy_index" "$(dirname "$decoys_txt")" "$(dirname "$gentrome_fa")"

    # Extract decoy chromosome/contig names from genome fasta
    grep "^>" "$genome_fa" | cut -d " " -f1 | sed 's/>//' > "$decoys_txt"

    # Combine transcript and genome sequences into gentrome
    cat "$transcript_fa" "$genome_fa" > "$gentrome_fa"

    "${cmd_salmon[@]}" index --threads "$threads" \
        -t "$gentrome_fa" \
        -d "$decoys_txt" \
        -i "$salmon_decoy_index"
fi
