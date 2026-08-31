#!/usr/bin/env bash
# scripts/minimap2_transcriptome.sh: long-read mapping to the transcriptome,
# in both default and high-sensitivity (-N 115) presets. Base mode uses the
# shared transcriptome index; pass <pct> <bootstrap> to map against that
# reduced-annotation combination's own index instead. Presets
# (mm2_tx_preset, mm2_n115_flag) come from whichever pipeline's own
# config.sh sourced this script — they differ between ONT and PacBio.
set -euo pipefail
source ./config.sh

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <input_reads> <sample_id> <tissue_id> <sample_num> [pct] [bootstrap]"
    exit 1
fi

reads_in="$1"
sample_id="$2"
tissue="$3"
sample="$4"

if [[ $# -ge 6 ]]; then
    pct="$5"
    bootstrap="$6"
    sample_dir="${out_dir}/${sample_id}/pct${pct}/rep${bootstrap}"
    this_tx_mmi=$(render_annotation_template "$reduced_annotation_tx_mmi_template" "$pct" "$bootstrap")
    label_suffix=" (${pct}% kept, bootstrap ${bootstrap})"
    timing_suffix="_pct${pct}_rep${bootstrap}"
else
    sample_dir="${out_dir}/${sample_id}"
    this_tx_mmi="$tx_mmi"
    label_suffix=""
    timing_suffix=""
fi
mkdir -p "$sample_dir"

cmd_minimap2=($(cmd_in_env env_minimap2 "$minimap2_bin"))
cmd_samtools=($(cmd_in_env env_samtools "$samtools_bin"))

echo "[minimap2] Mapping ${sample_id} to transcriptome${label_suffix} (-N 115)..."
tx_n115_bam="${sample_dir}/$(minimap2_n115_transcriptome_bam_name "$sample_id")"
tx_n115_sam_tmp="${sample_dir}/${sample_id}_minimap2_N115.tmp.sam"

/usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tminimap2_N115${timing_suffix}\t%e\t%M" \
    -o "$timing_log" -a \
    "${cmd_minimap2[@]}" $mm2_n115_flag -a -x "$mm2_tx_preset" -t "$threads" \
        "$this_tx_mmi" "$reads_in" > "$tx_n115_sam_tmp"

"${cmd_samtools[@]}" sort -@ "$threads" -n "$tx_n115_sam_tmp" -o "$tx_n115_bam"
rm -f "$tx_n115_sam_tmp"

echo "[minimap2] Mapping ${sample_id} to transcriptome${label_suffix} (default)..."
tx_bam="${sample_dir}/$(minimap2_transcriptome_bam_name "$sample_id")"
tx_sam_tmp="${sample_dir}/${sample_id}_minimap2.tmp.sam"

/usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tminimap2${timing_suffix}\t%e\t%M" \
    -o "$timing_log" -a \
    "${cmd_minimap2[@]}" -a -x "$mm2_tx_preset" -t "$threads" \
        "$this_tx_mmi" "$reads_in" > "$tx_sam_tmp"

"${cmd_samtools[@]}" sort -@ "$threads" -n "$tx_sam_tmp" -o "$tx_bam"
rm -f "$tx_sam_tmp"
