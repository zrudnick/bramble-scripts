#!/usr/bin/env bash
# generate_transcripts.sh: extracts transcript FASTAs via gffread (see README for usage).
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/../config.sh"

gtf_dir="$reduced_annotations_dir"
genome="$genome_fa"
output_dir="$reduced_annotation_transcripts_dir"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Run gffread on every GTF found under GTF_DIR, extracting transcript
sequences as FASTA into a mirrored directory structure under OUTPUT_DIR.
Intended to be run after generate_annotations.sh.

Options:
  -g, --gtf-dir DIR    Directory to search for *.gtf files (recursive)
                       (default: ${gtf_dir})
  -f, --genome FILE    Genome FASTA to extract sequences from
                       (default: ${genome})
  -o, --output DIR     Output directory for transcript FASTAs
                       (default: ${output_dir})
  -h, --help           Show this help and exit

Output: OUTPUT_DIR/<same relative path as the GTF, .gtf -> .fa>
        (e.g. .../t_80_1.gtf -> OUTPUT_DIR/t_80_1.fa)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--gtf-dir) gtf_dir="$2"; shift 2 ;;
        -f|--genome) genome="$2"; shift 2 ;;
        -o|--output) output_dir="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# ---- sanity checks --------------------------------------------------------
if [[ ! -d "$gtf_dir" ]]; then
    echo "ERROR: GTF directory not found: $gtf_dir" >&2
    echo "       (run generate_annotations.sh first)" >&2
    exit 1
fi

if [[ ! -f "$genome" ]]; then
    echo "ERROR: Genome FASTA not found: $genome" >&2
    exit 1
fi

# ---- collect GTF files ----------------------------------------------------
mapfile -t gtf_files < <(find "$gtf_dir" -name "*.gtf" | sort)

if [[ ${#gtf_files[@]} -eq 0 ]]; then
    echo "ERROR: No .gtf files found under $gtf_dir" >&2
    exit 1
fi

echo "============================================================"
echo " Extracting transcript FASTAs"
echo "   Found ${#gtf_files[@]} GTF file(s) under: ${gtf_dir}"
echo "   Genome:  ${genome}"
echo "   Output:  ${output_dir}"
echo "============================================================"

mkdir -p "$output_dir"

total=${#gtf_files[@]}
idx=0
for gtf in "${gtf_files[@]}"; do
    idx=$((idx + 1))

    # Mirror directory structure under output_dir
    rel="${gtf#"$gtf_dir"/}"          # e.g. t_80_1.gtf
    out_fa="$output_dir/${rel%.gtf}.fa"  # e.g. t_80_1.fa
    mkdir -p "$(dirname "$out_fa")"

    "$gffread_bin" "$gtf" -g "$genome" -w "$out_fa"

    n_seqs=$(grep -c "^>" "$out_fa" 2>/dev/null || echo 0)
    echo "[$idx/$total] $out_fa  ($n_seqs transcripts)"
done

echo ""
echo "Done. Transcript FASTAs written to: ${output_dir}"
