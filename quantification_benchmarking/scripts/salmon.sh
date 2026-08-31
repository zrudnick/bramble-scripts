#!/usr/bin/env bash
# scripts/salmon.sh: Salmon quantification (quasi-mapping & alignment-based).
# In base mode runs quasi-mapping against both the standard and decoy indices,
# plus alignment-based quant against every alignment-based BAM this pipeline
# produces (Bowtie2 default+k76, native STAR, and both Bramble projections).
# In reduced-annotation mode (pass <pct> <bootstrap>) runs quasi-mapping
# against only that combination's decoy index, plus alignment-based quant
# against the 3 BAMs the reduced-annotation pipeline actually produces
# (Bowtie2 k76-only, and both Bramble projections) — matches each pipeline's
# existing behavior.
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
    this_transcript_fa=$(render_annotation_template "$reduced_annotation_transcript_fa_template" "$pct" "$bootstrap")
    this_salmon_decoy_index=$(render_annotation_template "$reduced_annotation_salmon_decoy_index_template" "$pct" "$bootstrap")
    timing_suffix="_pct${pct}_rep${bootstrap}"
else
    sample_dir="${out_dir}/${sample_id}"
    this_transcript_fa="$transcript_fa"
    this_salmon_decoy_index="$salmon_decoy_index"
    timing_suffix=""
fi
salmon_out="${sample_dir}/salmon"
mkdir -p "$salmon_out"

cmd_salmon=($(cmd_in_env env_salmon "$salmon_bin"))

# 1. Quasi-mapping mode
if [[ "$reduced_mode" -eq 0 ]]; then
    echo "[Salmon] Running quasi-mapping (standard index)..."
    /usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tSalmon_Quasi_Std\t%e\t%M" \
        -o "$timing_log" -a \
        "${cmd_salmon[@]}" quant --threads "$threads" \
            -i "$salmon_std_index" \
            -l SF \
            -r "$fastq_in" \
            -o "${salmon_out}/quasi_std" 2> "${salmon_out}/quasi_std.log"
fi

if [[ -d "$this_salmon_decoy_index" ]]; then
    echo "[Salmon] Running quasi-mapping (decoy index)..."
    /usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tSalmon_Quasi_Decoy${timing_suffix}\t%e\t%M" \
        -o "$timing_log" -a \
        "${cmd_salmon[@]}" quant --threads "$threads" \
            -i "$this_salmon_decoy_index" \
            -l SF \
            -r "$fastq_in" \
            -o "${salmon_out}/quasi_decoy" 2> "${salmon_out}/quasi_decoy.log"
fi

# 2. Alignment-based mode (consuming projected/transcriptomic BAMs)
run_align_quant() {
    local mode_name="$1"
    local bam_file="$2"

    if [[ -f "$bam_file" ]]; then
        echo "[Salmon] Running alignment quantification (${mode_name})..."
        /usr/bin/time -f "${tissue}\t${sample}\t${sample_id}\tSalmon_Align_${mode_name}${timing_suffix}\t%e\t%M" \
            -o "$timing_log" -a \
            "${cmd_salmon[@]}" quant --threads "$threads" \
                -t "$this_transcript_fa" \
                -l SF \
                -a "$bam_file" \
                -o "${salmon_out}/align_${mode_name}" 2> "${salmon_out}/align_${mode_name}.log"
    fi
}

if [[ "$reduced_mode" -eq 0 ]]; then
    run_align_quant "bowtie2"      "${sample_dir}/$(bowtie2_transcriptome_bam_name "$sample_id")"
    run_align_quant "bowtie2_k76"  "${sample_dir}/$(bowtie2_k76_transcriptome_bam_name "$sample_id")"
    run_align_quant "STAR"         "${sample_dir}/$(star_transcriptome_bam_name "$sample_id")"
else
    run_align_quant "bowtie2"      "${sample_dir}/$(bowtie2_transcriptome_bam_name "$sample_id")"
fi
run_align_quant "bramble_HISAT2" "${sample_dir}/$(bramble_short_transcriptome_bam_name "$sample_id" "HISAT2")"
run_align_quant "bramble_STAR"   "${sample_dir}/$(bramble_short_transcriptome_bam_name "$sample_id" "STAR")"
