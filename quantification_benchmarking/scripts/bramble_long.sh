#!/usr/bin/env bash
# scripts/bramble_long.sh: long-mode projection of genomic alignments to
# transcriptomic space. Base mode uses the full annotation; pass <pct>
# <bootstrap> to project against that reduced-annotation combination's GTF
# instead. Output is the raw Bramble BAM plus a name-sorted copy — Oarfish
# and TranSigner both consume the sorted copy directly (no header rewrite
# needed).
set -euo pipefail
source ./config.sh

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <input_genome_bam> <sample_id> <tissue_id> <sample_num> [pct] [bootstrap]"
    exit 1
fi

genome_bam_in="$1"
sample_id="$2"
tissue="$3"
sample="$4"

if [[ $# -ge 6 ]]; then
    pct="$5"
    bootstrap="$6"
    sample_dir="${out_dir}/${sample_id}/pct${pct}/rep${bootstrap}"
    this_gtf_file=$(render_annotation_template "$reduced_annotation_gtf_template" "$pct" "$bootstrap")
    label_suffix=" (${pct}% kept, bootstrap ${bootstrap})"
    timing_suffix="_pct${pct}_rep${bootstrap}"
else
    sample_dir="${out_dir}/${sample_id}"
    this_gtf_file="$gtf_file"
    label_suffix=""
    timing_suffix=""
fi
mkdir -p "$sample_dir"

br_raw_bam="${sample_dir}/${sample_id}_bramble_transcriptome.bam"
br_sorted_bam="${sample_dir}/$(bramble_long_sorted_bam_name "$sample_id")"

echo "[Bramble] Projecting genomic alignments for ${sample_id}${label_suffix} (long mode)..."

cmd_bramble=($(cmd_in_env env_bramble "$bramble_bin"))
cmd_samtools=($(cmd_in_env env_samtools "$samtools_bin"))

/usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tBramble${timing_suffix}\t%e\t%M" \
    -o "$timing_log" -a \
    "${cmd_bramble[@]}" -p "$threads" \
        -S "$genome_fa" \
        --lr \
        -G "$this_gtf_file" \
        -o "$br_raw_bam" \
        "$genome_bam_in"

echo "[Bramble] Sorting projected BAM by read name..."
"${cmd_samtools[@]}" sort -@ "$threads" -n "$br_raw_bam" -o "$br_sorted_bam"
