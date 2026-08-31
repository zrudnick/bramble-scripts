#!/usr/bin/env bash
# pipeline.sh: reduced-annotation ONT long-read benchmark orchestration
set -euo pipefail
source ./config.sh

# Shared FASTQ data for both ONT pipelines (this one and ont) — files are
# named with t{tissue}_s{sample}; reads don't change across
# reduced-annotation percentages. Kept separate from data/long-read/FASTQ-pacbio/ — ONT and
# PacBio are distinct simulated datasets, not the same reads run through
# two aligners.
fastq_dir="$data_ont_dir"

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
echo " Reduced-annotation ONT benchmarking: ${#fastq_files[@]} FASTQ file(s)"
echo "============================================================"

for reads in "${fastq_files[@]}"; do
    filename=$(basename "$reads")

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
    echo ">> Input reads: ${reads}"

    # Step 1: minimap2 genome mapping (annotation-independent, run once per sample)
    bash ../../scripts/minimap2_genome.sh "$reads" "$sample_id" "$t" "$s"

    for pct in "${percentages[@]}"; do
        for bootstrap in "${bootstraps[@]}"; do
            echo ">>   ${pct}% kept, bootstrap ${bootstrap}"

            # Step 2: minimap2 transcriptome mapping (reduced-annotation index)
            bash ../../scripts/minimap2_transcriptome.sh "$reads" "$sample_id" "$t" "$s" "$pct" "$bootstrap"

            # Step 3: Bramble long-mode projection (reduced-annotation GTF)
            bash ../../scripts/bramble_long.sh "${out_dir}/${sample_id}/$(minimap2_genome_bam_name "$sample_id")" "$sample_id" "$t" "$s" "$pct" "$bootstrap"

            # Step 4: Quantification
            bash ../../scripts/oarfish.sh "$sample_id" "$t" "$s" "$pct" "$bootstrap"
            bash ../../scripts/transigner.sh "$sample_id" "$t" "$s" "$pct" "$bootstrap"

            echo ">>   Done ${pct}% kept, bootstrap ${bootstrap}"
        done
    done

    echo ">> Completed ${sample_id}"
    echo "------------------------------------------------------------"
done

echo "============================================================"
echo " Reduced-annotation ONT benchmarking pipeline complete!"
echo " Logged to: ${timing_log}"
echo "============================================================"
