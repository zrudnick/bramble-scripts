#!/usr/bin/env bash
# config.sh: shared configuration for the quantification_benchmarking pipelines
# (short-read, ont, pacbio, and their *_reduced_annotation variants),
# build_reduced_annotations/'s GTF-reduction scripts, and
# alignment_metrics/'s precision/recall analysis

_config_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ==========================================
# Conda Environments
# Leave as "" if installed in base/active environment
# ==========================================
env_bramble=""
env_samtools=""
env_gffread=""
# short-read
env_HISAT2=""
env_STAR=""
env_bowtie2=""
env_salmon=""
# long-read (ont/pacbio)
env_minimap2=""
env_oarfish=""
env_transigner=""

# ==========================================
# Environment manager & environments
# Set conda_mamba_bin to "conda", "mamba", etc.
# ==========================================
conda_mamba_bin="mamba"

# ==========================================
# Tool binaries (ensure these are in PATH or set absolute paths)
# ==========================================
export bramble_bin=""
export samtools_bin=""
export gffread_bin=""
# short-read
export HISAT2_bin=""
export STAR_bin=""
export bowtie2_bin=""
export salmon_bin=""
# long-read (ont/pacbio)
export minimap2_bin=""
export oarfish_bin=""
export transigner_bin=""

# ==========================================
# Resource allocation
# ==========================================
export threads=16
export timing_log="benchmark_log.tsv"

# ==========================================
# Directory structures
# ==========================================
export ref_dir="${_config_dir}/data/ref"
export index_dir="${_config_dir}/data/index"
export data_short_dir="${_config_dir}/data/short-read/FASTQ"
export data_ont_dir="${_config_dir}/data/long-read/FASTQ-ont"
export data_pacbio_dir="${_config_dir}/data/long-read/FASTQ-pacbio"
export out_dir="results"

fastq_extensions=("fastq" "fastq.gz" "fq" "fq.gz" "fasta" "fasta.gz" "fa" "fa.gz")

# ==========================================
# Full-annotation reference files
# Used directly by every pipeline running in "base" (non-reduced-annotation)
# mode, as the source for the shared/base indices below, and as the GTF
# build_reduced_annotations/generate_annotations.sh reduces.
# ==========================================
export genome_fa="$ref_dir/hg38_p12_ucsc.no_alts.no_fixs.fa"
export gtf_file="$ref_dir/chess2.2_ALL.gtf"
export transcript_fa="$ref_dir/chess2.2_ALL.no_alts.no_fixs.fa"
export gentrome_fa="$ref_dir/chess2.2_ALL.gentrome.no_alts.no_fixs.fa"
export decoys_txt="$ref_dir/decoys.no_alts.no_fixs.txt"

# ==========================================
# Shared/base index paths (short-read: HISAT2/STAR/Bowtie2/Salmon;
# long-read: minimap2). Built once by build_indices.sh, reused by every
# pipeline in "base" mode and by HISAT2/minimap2-genome in reduced-annotation
# mode too (see reduced-annotation templates below for what isn't shared).
# ==========================================
export HISAT2_index="${index_dir}/HISAT2/gn"
export STAR_index="${index_dir}/STAR_gn"
export bowtie2_index="${index_dir}/bowtie2/tx"
export salmon_std_index="${index_dir}/salmon_std"
export salmon_decoy_index="${index_dir}/salmon_decoy"
export tx_mmi="${index_dir}/minimap2_tx/transcripts_chess2.2_ALL_no_alts.no_fixs.mmi"
export genome_mmi="${index_dir}/minimap2_gn/genome_chess2.2_ALL_no_alts.no_fixs.mmi"

# ==========================================
# Reduced-annotation percentage/bootstrap grid
# Every *_reduced_annotation pipeline.sh, build_reduced_annotation_indices.sh,
# and generate_annotations.sh loop/pass this same percentage grid.
# ==========================================
percentages=("80" "60" "40" "20")
bootstraps=("1" "2" "3" "4" "5")
export bootstrap_count=5

# ==========================================
# Reduced-annotation generation
# (build_reduced_annotations/generate_annotations.sh + generate_transcripts.sh)
# and the inputs those scripts produce, consumed downstream by the
# *_reduced_annotation pipelines. {pct} (percentage of transcripts kept) and
# {bootstrap} (replicate index) in the templates below are substituted by
# render_annotation_template().
# ==========================================
export seed=42
export mode="transcript"
mode_prefix="t"

export reduced_annotations_dir="${ref_dir}/reduced_annotations"
export reduced_annotation_transcripts_dir="${ref_dir}/reduced_annotation_transcripts"
export reduced_annotation_gtf_template="${reduced_annotations_dir}/${mode_prefix}_{pct}_{bootstrap}.gtf"
export reduced_annotation_transcript_fa_template="${reduced_annotation_transcripts_dir}/${mode_prefix}_{pct}_{bootstrap}.fa"

# Per-percentage/bootstrap index & reference derivatives built by
# build_reduced_annotation_indices.sh
export reduced_annotation_STAR_index_template="${index_dir}/STAR_reduced_annotation/gn_{pct}_{bootstrap}"
export reduced_annotation_bowtie2_index_template="${index_dir}/bowtie2_reduced_annotation/tx_{pct}_{bootstrap}"
export reduced_annotation_gentrome_template="${ref_dir}/reduced_annotation_gentromes/chess2.2_ALL.gentrome.{pct}_{bootstrap}.fa"
export reduced_annotation_salmon_decoy_index_template="${index_dir}/salmon_decoy_reduced_annotation/tx_{pct}_{bootstrap}"
export reduced_annotation_tx_mmi_template="${index_dir}/minimap2_tx_reduced_annotation/transcripts_chess2.2_ALL_no_alts.no_fixs.{pct}_{bootstrap}.mmi"

# Substitute {pct}/{bootstrap} placeholders in a path template
render_annotation_template() {
    local template="$1" pct="$2" bootstrap="$3"
    template="${template//\{pct\}/$pct}"
    template="${template//\{bootstrap\}/$bootstrap}"
    echo "$template"
}

# ==========================================
# BAM output filename conventions
# ==========================================
hisat2_genome_bam_name() { echo "${1}_HISAT2_genome.bam"; }
star_genome_bam_name() { echo "${1}_STAR_genome.bam"; }
star_transcriptome_bam_name() { echo "${1}_STAR_transcriptome.bam"; }
bowtie2_transcriptome_bam_name() { echo "${1}_bowtie2_transcriptome.bam"; }
bowtie2_k76_transcriptome_bam_name() { echo "${1}_bowtie2_k76_transcriptome.bam"; }
bramble_short_transcriptome_bam_name() { echo "${1}_bramble_${2}_transcriptome.bam"; }  # $2=aligner_name
minimap2_genome_bam_name() { echo "${1}_minimap2_genome.bam"; }
minimap2_transcriptome_bam_name() { echo "${1}_minimap2_transcriptome.bam"; }
minimap2_n115_transcriptome_bam_name() { echo "${1}_minimap2_N115_transcriptome.bam"; }
bramble_long_sorted_bam_name() { echo "${1}_bramble_transcriptome.sorted_name.bam"; }

# ==========================================
# Helper function for command building
# ==========================================
cmd_in_env() {
    local env_var="$1"
    shift
    local env_name="${!env_var:-}"

    if [[ -n "$env_name" ]]; then
        # Force MAMBA_USE_LOCKFILES=false to prevent process lock collisions
        echo "env MAMBA_USE_LOCKFILES=false ${conda_mamba_bin} run -n ${env_name} $*"
    else
        echo "$*"
    fi
}
