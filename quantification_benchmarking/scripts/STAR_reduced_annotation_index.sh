#!/usr/bin/env bash
# scripts/STAR_reduced_annotation_index.sh: build one full STAR genome index
# per percentage/bootstrap combination. STAR bakes splice junctions directly
# into its index, so unlike HISAT2 it cannot share one index across
# annotations — a shared index built from the full annotation would leak
# full-annotation knowledge into every combination's alignment.
set -euo pipefail
source ./config.sh

cmd_STAR=($(cmd_in_env env_STAR "$STAR_bin"))

for pct in "${percentages[@]}"; do
    for bootstrap in "${bootstraps[@]}"; do
        reduced_annotation_gtf=$(render_annotation_template "$reduced_annotation_gtf_template" "$pct" "$bootstrap")
        this_STAR_index=$(render_annotation_template "$reduced_annotation_STAR_index_template" "$pct" "$bootstrap")

        if [[ ! -f "$reduced_annotation_gtf" ]]; then
            echo "   !! Missing reduced-annotation GTF: ${reduced_annotation_gtf} — skipping"
            continue
        fi

        if [[ ! -d "$this_STAR_index" ]]; then
            echo ">> percentage ${pct}%, bootstrap ${bootstrap}: building STAR genome index..."
            mkdir -p "$this_STAR_index"
            "${cmd_STAR[@]}" --runThreadN "$threads" \
                --runMode genomeGenerate \
                --genomeDir "$this_STAR_index" \
                --genomeFastaFiles "$genome_fa" \
                --sjdbGTFfile "$reduced_annotation_gtf" \
                --sjdbOverhang 100
        fi
    done
done
