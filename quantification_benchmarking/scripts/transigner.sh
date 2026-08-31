#!/usr/bin/env bash
# scripts/transigner.sh: TranSigner alignment-based quantification (minimap2
# & Bramble). Base mode reads from the flat sample directory; pass <pct>
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
quant_dir="${sample_dir}/transigner"
mkdir -p "$quant_dir"

cmd_transigner=($(cmd_in_env env_transigner "$transigner_bin"))

run_transigner() {
    local label="$1"
    local bam_file="$2"
    local run_out_dir="${quant_dir}/${label}"

    if [[ -f "$bam_file" ]]; then
        echo "[TranSigner] Quantifying (${label})..."
        mkdir -p "$run_out_dir"
        /usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tTranSigner_${label}${timing_suffix}\t%e\t%M" \
            -o "$timing_log" -a \
            "${cmd_transigner[@]}" \
            "$bam_file" \
            -o "$run_out_dir" \
            -p "$threads" \
            --keep-low \
            --use-psw
    fi
}

run_transigner "mm2_N115" "${sample_dir}/$(minimap2_n115_transcriptome_bam_name "$sample_id")"
run_transigner "mm2"      "${sample_dir}/$(minimap2_transcriptome_bam_name "$sample_id")"
run_transigner "bramble"  "${sample_dir}/$(bramble_long_sorted_bam_name "$sample_id")"
