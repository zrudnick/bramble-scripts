#!/usr/bin/env bash
# scripts/minimap2_genome_index.sh: build the one shared minimap2 genome
# index. No reduced-annotation counterpart — minimap2's spliced alignment
# takes no annotation input, so this index is always shared, for every
# long-read pipeline.
set -euo pipefail
source ./config.sh

cmd_minimap2=($(cmd_in_env env_minimap2 "$minimap2_bin"))

if [[ -f "$genome_fa" && ! -f "$genome_mmi" ]]; then
    echo "Building minimap2 genome index..."
    mkdir -p "$(dirname "$genome_mmi")"
    "${cmd_minimap2[@]}" -d "$genome_mmi" "$genome_fa"
fi
