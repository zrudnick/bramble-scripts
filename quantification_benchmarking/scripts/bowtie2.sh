#!/usr/bin/env bash
# scripts/bowtie2.sh: direct transcriptome mapping using Bowtie2. In base mode
# runs both default and high-sensitivity (-k 76) against the shared index; in
# reduced-annotation mode (pass <pct> <bootstrap>) runs only -k 76 against
# that combination's own index — matches each pipeline's existing behavior.
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

if [[ $# -ge 6 ]]; then
    pct="$5"
    bootstrap="$6"
    sample_dir="${out_dir}/${sample_id}/pct${pct}/rep${bootstrap}"
    this_bowtie2_index=$(render_annotation_template "$reduced_annotation_bowtie2_index_template" "$pct" "$bootstrap")
    timing_suffix="_pct${pct}_rep${bootstrap}"
else
    sample_dir="${out_dir}/${sample_id}"
    this_bowtie2_index="$bowtie2_index"
    timing_suffix=""
fi
mkdir -p "$sample_dir"

cmd_bowtie2=($(cmd_in_env env_bowtie2 "$bowtie2_bin"))
cmd_samtools=($(cmd_in_env env_samtools "$samtools_bin"))

run_bowtie2() {
    local tag="$1"          # "" for default, "_k76" for high-sensitivity
    local timing_label="$2" # e.g. "Bowtie2" or "Bowtie2_k76"
    local extra_flags="$3"
    local bam_name_fn="$4"  # which shared naming function names this run's BAM

    local bam_out="${sample_dir}/$("$bam_name_fn" "$sample_id")"
    local sam_tmp="${sample_dir}/${sample_id}_bowtie2${tag}.tmp.sam"
    local tool_log="${sample_dir}/bowtie2${tag}.log"

    echo "[Bowtie2] Mapping ${sample_id} to transcriptome (${timing_label})..."

    /usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\t${timing_label}\t%e\t%M" \
        -o "$timing_log" -a \
        "${cmd_bowtie2[@]}" $extra_flags -p "$threads" \
            -x "$this_bowtie2_index" \
            -U "$fastq_in" \
            -S "$sam_tmp" 2> "$tool_log"

    echo "[Bowtie2] Sorting BAM by read name (${timing_label})..."
    "${cmd_samtools[@]}" sort -@ "$threads" -n "$sam_tmp" -o "$bam_out"
    rm -f "$sam_tmp"
}

if [[ -n "$timing_suffix" ]]; then
    # Reduced-annotation mode only ever runs the -k 76 mapping, but names its
    # output the same as the "default" mode's BAM (there's only one to disambiguate)
    run_bowtie2 "" "Bowtie2${timing_suffix}" "-k 76" bowtie2_transcriptome_bam_name
else
    run_bowtie2 "" "Bowtie2" "" bowtie2_transcriptome_bam_name
    run_bowtie2 "_k76" "Bowtie2_k76" "-k 76" bowtie2_k76_transcriptome_bam_name
fi
