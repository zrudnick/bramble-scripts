#!/usr/bin/env bash
# scripts/minimap2_transcriptome_reduced_annotation_index.sh: build one
# minimap2 transcriptome index per percentage/bootstrap combination, from
# that combination's transcript FASTA.
set -euo pipefail
source ./config.sh

cmd_minimap2=($(cmd_in_env env_minimap2 "$minimap2_bin"))

for pct in "${percentages[@]}"; do
    for bootstrap in "${bootstraps[@]}"; do
        this_transcript_fa=$(render_annotation_template "$reduced_annotation_transcript_fa_template" "$pct" "$bootstrap")
        this_tx_mmi=$(render_annotation_template "$reduced_annotation_tx_mmi_template" "$pct" "$bootstrap")

        if [[ ! -f "$this_transcript_fa" ]]; then
            echo "   !! Missing reduced-annotation transcript FASTA: ${this_transcript_fa} — skipping"
            continue
        fi

        if [[ ! -f "$this_tx_mmi" ]]; then
            echo ">> percentage ${pct}%, bootstrap ${bootstrap}: building minimap2 transcript index..."
            mkdir -p "$(dirname "$this_tx_mmi")"
            "${cmd_minimap2[@]}" -d "$this_tx_mmi" "$this_transcript_fa"
        fi
    done
done
