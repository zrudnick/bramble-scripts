#!/usr/bin/env bash
# scripts/bowtie2_reduced_annotation_index.sh: build one Bowtie2
# transcriptome index per percentage/bootstrap combination, from that
# combination's transcript FASTA.
set -euo pipefail
source ./config.sh

cmd_bowtie2_build=($(cmd_in_env env_bowtie2 "${bowtie2_bin}-build"))

for pct in "${percentages[@]}"; do
    for bootstrap in "${bootstraps[@]}"; do
        this_transcript_fa=$(render_annotation_template "$reduced_annotation_transcript_fa_template" "$pct" "$bootstrap")
        this_bowtie2_index=$(render_annotation_template "$reduced_annotation_bowtie2_index_template" "$pct" "$bootstrap")

        if [[ ! -f "$this_transcript_fa" ]]; then
            echo "   !! Missing reduced-annotation transcript FASTA: ${this_transcript_fa} — skipping"
            continue
        fi

        if [[ ! -f "${this_bowtie2_index}.1.bt2" ]]; then
            echo ">> percentage ${pct}%, bootstrap ${bootstrap}: building Bowtie2 index..."
            mkdir -p "$(dirname "$this_bowtie2_index")"
            "${cmd_bowtie2_build[@]}" -f --threads "$threads" "$this_transcript_fa" "$this_bowtie2_index"
        fi
    done
done
