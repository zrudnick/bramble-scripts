#!/usr/bin/env bash
# scripts/STAR.sh: genomic alignment using STAR. STAR bakes splice junctions
# into its genome index (and re-inserts them again at alignment time), so
# unlike HISAT2 it cannot share one index across annotations — the base index
# is built from the full annotation; pass <pct> <bootstrap> to align against
# that reduced-annotation combination's own index instead (see
# scripts/STAR_index.sh / scripts/STAR_reduced_annotation_index.sh).
#
# In base mode this also emits STAR's native transcriptome BAM
# (--quantMode TranscriptomeSAM), which salmon.sh's "star" alignment-based
# quant mode consumes; the reduced-annotation pipeline never used that mode,
# so it's skipped when pct/bootstrap are given.
#
# STAR doesn't decompress gzipped FASTQ on its own — a .gz input needs
# --readFilesCommand zcat, so this script adds it automatically when needed.
set -euo pipefail
source ./config.sh

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <input_fastq> <sample_id> <tissue_id> <sample_num> [pct] [bootstrap]"
    exit 1
fi

fastq_in="$1"
sample_id="$2"
tissue="$3"
sample="$4"

reduced_mode=0
if [[ $# -ge 6 ]]; then
    reduced_mode=1
    pct="$5"
    bootstrap="$6"
    sample_dir="${out_dir}/${sample_id}/pct${pct}/rep${bootstrap}"
    this_STAR_index=$(render_annotation_template "$reduced_annotation_STAR_index_template" "$pct" "$bootstrap")
    this_gtf_file=$(render_annotation_template "$reduced_annotation_gtf_template" "$pct" "$bootstrap")
    timing_label="STAR_pct${pct}_rep${bootstrap}"
else
    sample_dir="${out_dir}/${sample_id}"
    this_STAR_index="$STAR_index"
    this_gtf_file="$gtf_file"
    timing_label="STAR"
fi
mkdir -p "$sample_dir"

prefix="${sample_dir}/"
genome_bam="${sample_dir}/$(star_genome_bam_name "$sample_id")"
tx_bam="${sample_dir}/$(star_transcriptome_bam_name "$sample_id")"
tool_log="${sample_dir}/STAR.log"

echo "[STAR] Aligning ${sample_id}..."

cmd_STAR=($(cmd_in_env env_STAR "$STAR_bin"))
cmd_samtools=($(cmd_in_env env_samtools "$samtools_bin"))

quant_mode_flags=()
if [[ "$reduced_mode" -eq 0 ]]; then
    quant_mode_flags=(--quantMode TranscriptomeSAM)
fi

readfiles_command_flags=()
if [[ "$fastq_in" == *.gz ]]; then
    readfiles_command_flags=(--readFilesCommand zcat)
fi

/usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\t${timing_label}\t%e\t%M" \
    -o "$timing_log" -a \
    "${cmd_STAR[@]}" --runThreadN "$threads" \
        --genomeDir "$this_STAR_index" \
        --readFilesIn "$fastq_in" \
        "${readfiles_command_flags[@]}" \
        --outSAMtype SAM \
        --sjdbGTFfile "$this_gtf_file" \
        --sjdbOverhang 100 \
        "${quant_mode_flags[@]}" \
        --outFileNamePrefix "$prefix" 2> "$tool_log"

echo "[STAR] Sorting genomic alignment BAM by read name..."
"${cmd_samtools[@]}" sort -@ "$threads" -n "${prefix}Aligned.out.sam" -o "$genome_bam"
rm -f "${prefix}Aligned.out.sam"

if [[ -f "${prefix}Aligned.toTranscriptome.out.bam" ]]; then
    echo "[STAR] Sorting native transcriptome BAM by read name..."
    "${cmd_samtools[@]}" sort -@ "$threads" -n "${prefix}Aligned.toTranscriptome.out.bam" -o "$tx_bam"
    rm -f "${prefix}Aligned.toTranscriptome.out.bam"
fi
