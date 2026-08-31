#!/usr/bin/env bash
# scripts/salmon_reduced_annotation_index.sh: build one Salmon decoy-aware
# index per percentage/bootstrap combination (no separate standard index in
# reduced-annotation mode, matching each pipeline's existing behavior).
set -euo pipefail
source ./config.sh

cmd_salmon=($(cmd_in_env env_salmon "$salmon_bin"))

for pct in "${percentages[@]}"; do
    for bootstrap in "${bootstraps[@]}"; do
        this_transcript_fa=$(render_annotation_template "$reduced_annotation_transcript_fa_template" "$pct" "$bootstrap")
        this_gentrome_fa=$(render_annotation_template "$reduced_annotation_gentrome_template" "$pct" "$bootstrap")
        this_salmon_decoy_index=$(render_annotation_template "$reduced_annotation_salmon_decoy_index_template" "$pct" "$bootstrap")

        if [[ ! -f "$this_transcript_fa" ]]; then
            echo "   !! Missing reduced-annotation transcript FASTA: ${this_transcript_fa} — skipping"
            continue
        fi

        if [[ ! -d "$this_salmon_decoy_index" ]]; then
            echo ">> percentage ${pct}%, bootstrap ${bootstrap}: building Salmon decoy (gentrome) index..."
            mkdir -p "$(dirname "$this_gentrome_fa")" "$this_salmon_decoy_index"
            cat "$this_transcript_fa" "$genome_fa" > "$this_gentrome_fa"
            "${cmd_salmon[@]}" index --threads "$threads" \
                -t "$this_gentrome_fa" \
                -d "$decoys_txt" \
                -i "$this_salmon_decoy_index"
        fi
    done
done
