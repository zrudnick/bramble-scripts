#!/usr/bin/env bash
# scripts/bramble_short.sh: project a short-read genomic alignment (HISAT2 or
# STAR) to transcriptomic space via Bramble's short-read mode. Base mode uses
# the full annotation; pass <pct> <bootstrap> to project against that
# reduced-annotation combination's GTF instead.
set -euo pipefail
source ./config.sh

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <input_genomic_bam> <aligner_name> <sample_id> <tissue> <sample> [pct] [bootstrap]"
    exit 1
fi

genomic_bam_in="$1"
aligner_name="$2" # e.g. "HISAT2" or "STAR"
sample_id="$3"
tissue="$4"
sample="$5"

if [[ $# -ge 7 ]]; then
    pct="$6"
    bootstrap="$7"
    sample_dir="${out_dir}/${sample_id}/pct${pct}/rep${bootstrap}"
    this_gtf_file=$(render_annotation_template "$reduced_annotation_gtf_template" "$pct" "$bootstrap")
    timing_suffix="_pct${pct}_rep${bootstrap}"
else
    sample_dir="${out_dir}/${sample_id}"
    this_gtf_file="$gtf_file"
    timing_suffix=""
fi
mkdir -p "$sample_dir"

tx_bam_out="${sample_dir}/$(bramble_short_transcriptome_bam_name "$sample_id" "$aligner_name")"
unsorted_tmp="${sample_dir}/${sample_id}_bramble_${aligner_name}.tmp.bam"
tool_log="${sample_dir}/bramble_${aligner_name}.log"

echo "[Bramble] Projecting ${aligner_name} genomic alignments for ${sample_id}..."

cmd_bramble=($(cmd_in_env env_bramble "$bramble_bin"))
cmd_samtools=($(cmd_in_env env_samtools "$samtools_bin"))

/usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tBramble_${aligner_name}${timing_suffix}\t%e\t%M" \
    -o "$timing_log" -a \
    "${cmd_bramble[@]}" -p "$threads" \
        -G "$this_gtf_file" \
        -o "$unsorted_tmp" \
        "$genomic_bam_in" 2> "$tool_log"

echo "[Bramble] Sorting projected transcriptomic BAM by read name..."
"${cmd_samtools[@]}" sort -@ "$threads" -n "$unsorted_tmp" -o "$tx_bam_out"

rm -f "$unsorted_tmp"
