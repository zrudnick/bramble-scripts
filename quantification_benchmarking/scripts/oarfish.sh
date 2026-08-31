#!/usr/bin/env bash
# scripts/oarfish.sh: Oarfish alignment-based quantification (minimap2 &
# Bramble). Base mode reads from the flat sample directory; pass <pct>
# <bootstrap> to read that reduced-annotation combination's BAMs instead.
# Every call uses the standard thread count uniformly, including the
# Bramble-derived BAM.
set -euo pipefail
source ./config.sh

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <sample_id> <tissue_id> <sample_num> [pct] [bootstrap]"
    exit 1
fi

sample_id="$1"
tissue="$2"
sample="$3"

if [[ $# -ge 5 ]]; then
    pct="$4"
    bootstrap="$5"
    sample_dir="${out_dir}/${sample_id}/pct${pct}/rep${bootstrap}"
    timing_suffix="_pct${pct}_rep${bootstrap}"
else
    sample_dir="${out_dir}/${sample_id}"
    timing_suffix=""
fi
quant_dir="${sample_dir}/oarfish"
mkdir -p "$quant_dir"

cmd_oarfish=($(cmd_in_env env_oarfish "$oarfish_bin"))

run_oarfish() {
    local label="$1"
    local bam_file="$2"

    if [[ -f "$bam_file" ]]; then
        echo "[Oarfish] Quantifying (${label})..."
        /usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tOarfish_${label}${timing_suffix}\t%e\t%M" \
            -o "$timing_log" -a \
            "${cmd_oarfish[@]}" \
            -a "$bam_file" \
            -o "${quant_dir}/${label}" \
            -j ${threads} \
            --filter-group no-filters \
            --model-coverage
    fi
}

run_oarfish "mm2_N115" "${sample_dir}/$(minimap2_n115_transcriptome_bam_name "$sample_id")"
run_oarfish "mm2"      "${sample_dir}/$(minimap2_transcriptome_bam_name "$sample_id")"
run_oarfish "bramble"  "${sample_dir}/$(bramble_long_sorted_bam_name "$sample_id")"
