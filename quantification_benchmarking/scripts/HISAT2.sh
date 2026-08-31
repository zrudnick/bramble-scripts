#!/usr/bin/env bash
# scripts/HISAT2.sh: genomic alignment using HISAT2 against the one shared
# genome index. Annotation-blind — no known-splicesite-infile is passed, so
# this step is always annotation-blind and always shared — used identically
# by short-read/ and short-read_reduced_annotation/ (never takes
# pct/bootstrap; see scripts/minimap2_genome.sh for the same pattern on the
# long-read side).
set -euo pipefail
source ./config.sh

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <input_fastq> <sample_id> <tissue_id> <sample_num>"
    exit 1
fi

fastq_in="$1"
sample_id="$2"
tissue="$3"
sample="$4"

sample_dir="${out_dir}/${sample_id}"
mkdir -p "$sample_dir"

bam_out="${sample_dir}/$(hisat2_genome_bam_name "$sample_id")"
sam_tmp="${sample_dir}/${sample_id}_HISAT2.tmp.sam"
tool_log="${sample_dir}/HISAT2.log"

echo "[HISAT2] Aligning ${sample_id}..."

cmd_HISAT2=($(cmd_in_env env_HISAT2 "$HISAT2_bin"))
cmd_samtools=($(cmd_in_env env_samtools "$samtools_bin"))

/usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tHISAT2\t%e\t%M" \
    -o "$timing_log" -a \
    "${cmd_HISAT2[@]}" -p "$threads" \
        -x "$HISAT2_index" \
        -U "$fastq_in" \
        -S "$sam_tmp" 2> "$tool_log"

echo "[HISAT2] Sorting BAM by read name..."
"${cmd_samtools[@]}" sort -@ "$threads" -n "$sam_tmp" -o "$bam_out"
rm -f "$sam_tmp"
