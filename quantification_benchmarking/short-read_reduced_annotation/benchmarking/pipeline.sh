#!/usr/bin/env bash
# pipeline.sh: reduced-annotation short-read benchmark orchestration
set -euo pipefail
source ./config.sh

# Shared FASTQ data for all short-read pipelines (this one and short-read) —
# files are named with t{tissue}_s{sample}; reads don't change across
# reduced-annotation percentages.
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
echo " Reduced-annotation short-read benchmarking: ${#fastq_files[@]} FASTQ file(s)"
echo "============================================================"

for fq in "${fastq_files[@]}"; do
    filename=$(basename "$fq")

    # Extract tissue (t) and sample (s) from the filename, e.g. t0_s1
    if [[ "$filename" =~ t([0-9]+)_s([0-9]+) ]]; then
        t="${BASH_REMATCH[1]}"
        s="${BASH_REMATCH[2]}"
        sample_id="t${t}_s${s}"
    else
        echo "Error: filename '${filename}' has no t{tissue}_s{sample} pattern"
        exit 1
    fi

    echo ">> Processing Sample: ${sample_id} (Tissue: ${t}, Sample: ${s})"
    echo ">> Input FASTQ: ${fq}"

    # Step 1: HISAT2 genomic alignment (shared plain index, annotation-blind
    # — no splice-sites file, so it produces the same alignment regardless of
    # pct/bootstrap; run once per sample, outside the loop, like
    # ont_reduced_annotation/pacbio_reduced_annotation's minimap2 genome step)
    bash ../../scripts/HISAT2.sh "$fq" "$sample_id" "$t" "$s"

    for pct in "${percentages[@]}"; do
        for bootstrap in "${bootstraps[@]}"; do
            echo ">>   ${pct}% kept, bootstrap ${bootstrap}"

            # Step 2: STAR genomic alignment (reduced-annotation index — STAR
            # bakes splice junctions into its index, so this must be redone
            # per pct/bootstrap to stay blind to the full annotation)
            bash ../../scripts/STAR.sh "$fq" "$sample_id" "$t" "$s" "$pct" "$bootstrap"

            # Step 3: Direct transcriptome alignment (reduced-annotation index)
            bash ../../scripts/bowtie2.sh "$fq" "$sample_id" "$t" "$s" "$pct" "$bootstrap"

            # Step 4: Bramble transcriptomic projections (reduced-annotation GTF)
            # — HISAT2's input BAM is the single shared one from Step 1;
            # bramble_short.sh still gets pct/bootstrap so its own projection
            # and output land in that combination's directory.
            bash ../../scripts/bramble_short.sh "${out_dir}/${sample_id}/$(hisat2_genome_bam_name "$sample_id")" "HISAT2" "$sample_id" "$t" "$s" "$pct" "$bootstrap"
            bash ../../scripts/bramble_short.sh "${out_dir}/${sample_id}/pct${pct}/rep${bootstrap}/$(star_genome_bam_name "$sample_id")" "STAR" "$sample_id" "$t" "$s" "$pct" "$bootstrap"

            # Step 5: Salmon quantification (reduced-annotation transcript reference)
            bash ../../scripts/salmon.sh "$fq" "$sample_id" "$t" "$s" "$pct" "$bootstrap"

            echo ">>   Done ${pct}% kept, bootstrap ${bootstrap}"
        done
    done

    echo ">> Completed ${sample_id}"
    echo "------------------------------------------------------------"
done

echo "============================================================"
echo " Reduced-annotation benchmarking pipeline complete!"
echo " Logged to: ${timing_log}"
echo "============================================================"
