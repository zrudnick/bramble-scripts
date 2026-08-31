#!/usr/bin/env bash
# generate_annotations.sh: wraps reduce_annotations.py (see README for usage).
set -euo pipefail

# Resolve our own location so this can be run from anywhere, not just this directory
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
generator="${script_dir}/reduce_annotations.py"
source "${script_dir}/../config.sh"

input_gtf="$gtf_file"
output_dir="$reduced_annotations_dir"
bootstrap="$bootstrap_count"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Generate reduced-annotation GTFs by randomly keeping a percentage of
transcripts (or genes) from a full annotation, bootstrapped across
multiple independent replicates per percentage level.

Options:
  -i, --input FILE          Full annotation GTF to reduce
                             (default: ${input_gtf})
  -o, --output DIR          Output directory
                             (default: ${output_dir})
  -p, --percentages "P..."  Space-separated percentages to keep
                             (default: "${percentages[*]}")
  -b, --bootstrap N         Number of independent replicates per percentage
                             (default: ${bootstrap})
  -s, --seed N              Base random seed (default: ${seed})
      --mode gene|transcript
                             Reduce at the gene or transcript level
                             (default: ${mode})
  -h, --help                Show this help and exit

Output: OUTPUT/t_{percentage}_{replicate}.gtf
        (e.g. reduced_annotations/t_80_1.gtf)
        plus a manifest TSV recording every keep/remove decision.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input) input_gtf="$2"; shift 2 ;;
        -o|--output) output_dir="$2"; shift 2 ;;
        -p|--percentages) read -r -a percentages <<< "$2"; shift 2 ;;
        -b|--bootstrap) bootstrap="$2"; shift 2 ;;
        -s|--seed) seed="$2"; shift 2 ;;
        --mode) mode="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ ! -f "$input_gtf" ]]; then
    echo "ERROR: Input GTF not found: $input_gtf" >&2
    exit 1
fi

if [[ "$mode" != "gene" && "$mode" != "transcript" ]]; then
    echo "ERROR: --mode must be 'gene' or 'transcript', got: $mode" >&2
    exit 1
fi

mkdir -p "$output_dir"
manifest="${output_dir}/${mode}_manifest.tsv"

echo "============================================================"
echo " Generating reduced annotations (${mode}-level)"
echo "   Input:        ${input_gtf}"
echo "   Output:       ${output_dir}/"
echo "   Percentages:  ${percentages[*]}"
echo "   Bootstrap:    ${bootstrap} replicate(s) per percentage"
echo "   Manifest:     ${manifest}"
echo "============================================================"

python3 "$generator" \
    -i "$input_gtf" \
    -o "$output_dir" \
    -m "$mode" \
    -p "${percentages[@]}" \
    --bootstrap "$bootstrap" \
    --seed "$seed" \
    --manifest "$manifest" \
    --batch \
    --verbose

echo ""
echo "Done. Reduced annotations written to: ${output_dir}/"
echo "Manifest: ${manifest}"
