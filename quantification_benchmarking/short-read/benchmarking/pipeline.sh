#!/usr/bin/env bash
# pipeline.sh: short-read benchmark orchestration
set -euo pipefail
source ./config.sh

# Shared FASTQ data for all short-read pipelines (this one and
# short-read_reduced_annotation) — files are named with t{tissue}_s{sample}.
fastq_dir="$data_short_dir"

mkdir -p "$out_dir" "$(dirname "$timing_log")"

if [[ ! -f "$timing_log" ]]; then
    printf "Tissue\tSample\tSample_ID\tTool\tWall_Time_Sec\tMax_RSS_KB\n" > "$timing_log"
fi

# Discover all FASTQ files under fastq_dir matching fastq_extensions
find_args=()
for ext in "${fastq_extensions[@]}"; do
    find_args+=(-o -name "*.${ext}")
done
mapfile -t fastq_files < <(find "$fastq_dir" -type f \( "${find_args[@]:1}" \) | sort)

if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "Error: No FASTQ files found in directory '${fastq_dir}'"
    exit 1
fi

echo "============================================================"
echo " Found ${#fastq_files[@]} FASTQ file(s) for benchmarking"
echo "============================================================"

for fq in "${fastq_files[@]}"; do
    filename=$(basename "$fq")

    # Extract tissue (t) and sample (s) from the filename, e.g. t0_s1
    if [[ "$filename" =~ t([0-9]+)_s([0-9]+) ]]; then
        tissue="${BASH_REMATCH[1]}"
        sample="${BASH_REMATCH[2]}"
        sample_id="t${tissue}_s${sample}"
    else
        echo "Error: filename '${filename}' has no t{tissue}_s{sample} pattern"
        exit 1
    fi

    echo ">> Processing Sample: ${sample_id} (Tissue: ${tissue}, Sample: ${sample})"
    echo ">> Input FASTQ: ${fq}"

    # Step 1: Genomic alignments
    bash ../../scripts/HISAT2.sh "$fq" "$sample_id" "$tissue" "$sample"
    bash ../../scripts/STAR.sh "$fq" "$sample_id" "$tissue" "$sample"

    # Step 2: Direct transcriptome alignment
    bash ../../scripts/bowtie2.sh "$fq" "$sample_id" "$tissue" "$sample"

    # Step 3: Bramble transcriptomic projections
    bash ../../scripts/bramble_short.sh "${out_dir}/${sample_id}/$(hisat2_genome_bam_name "$sample_id")" "HISAT2" "$sample_id" "$tissue" "$sample"
    bash ../../scripts/bramble_short.sh "${out_dir}/${sample_id}/$(star_genome_bam_name "$sample_id")" "STAR" "$sample_id" "$tissue" "$sample"

    # Step 4: Salmon quantification
    bash ../../scripts/salmon.sh "$fq" "$sample_id" "$tissue" "$sample"

    echo ">> Completed ${sample_id}"
    echo "------------------------------------------------------------"
done

echo "============================================================"
echo " Benchmarking pipeline complete!"
echo " Logged to: ${timing_log}"
echo "============================================================"
