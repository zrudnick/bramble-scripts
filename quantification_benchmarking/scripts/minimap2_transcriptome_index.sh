#!/usr/bin/env bash
# scripts/minimap2_transcriptome_index.sh: build the one shared minimap2
# transcriptome index from the full annotation, used by ont/pacbio in base
# mode.
set -euo pipefail
source ./config.sh

cmd_minimap2=($(cmd_in_env env_minimap2 "$minimap2_bin"))

if [[ -f "$transcript_fa" && ! -f "$tx_mmi" ]]; then
    echo "Building minimap2 transcriptome index..."
    mkdir -p "$(dirname "$tx_mmi")"
    "${cmd_minimap2[@]}" -d "$tx_mmi" "$transcript_fa"
fi
