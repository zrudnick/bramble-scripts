#!/usr/bin/env bash
# scripts/minimap2_genome.sh: long-read mapping to the shared genome index.
# minimap2's spliced alignment takes no annotation input, so this step is
# always annotation-blind and always shared — used identically by ont/,
# ont_reduced_annotation/, pacbio/, and pacbio_reduced_annotation/ (never
# takes pct/bootstrap). The preset (mm2_gn_preset) comes from whichever
# pipeline's own config.sh sourced this script.
set -euo pipefail
source ./config.sh

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <input_reads> <sample_id> <tissue_id> <sample_num>"
    exit 1
fi

reads_in="$1"
sample_id="$2"
tissue="$3"
sample="$4"

sample_dir="${out_dir}/${sample_id}"
mkdir -p "$sample_dir"

cmd_minimap2=($(cmd_in_env env_minimap2 "$minimap2_bin"))
cmd_samtools=($(cmd_in_env env_samtools "$samtools_bin"))

echo "[minimap2] Mapping ${sample_id} to genome (shared index)..."
genome_bam="${sample_dir}/$(minimap2_genome_bam_name "$sample_id")"
genome_sam_tmp="${sample_dir}/${sample_id}_minimap2_genome.tmp.sam"

/usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tminimap2_genome\t%e\t%M" \
    -o "$timing_log" -a \
    "${cmd_minimap2[@]}" -a -x "$mm2_gn_preset" -t "$threads" \
        "$genome_mmi" "$reads_in" -o "$genome_sam_tmp"

"${cmd_samtools[@]}" sort -@ "$threads" -n "$genome_sam_tmp" -o "$genome_bam"
rm -f "$genome_sam_tmp"
