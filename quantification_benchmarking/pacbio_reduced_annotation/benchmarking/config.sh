#!/usr/bin/env bash
# config.sh: PacBio reduced-annotation pipeline configuration — pulls in the
# shared config, then sets the platform-specific overrides PacBio needs (same
# presets as pacbio/, differ from ONT).
source ../../../config.sh

# minimap2 mapping presets — differ between ONT and PacBio, and are only used
# at alignment time (never index-build time), so they live here rather than
# in the shared config.
export mm2_tx_preset="lr:hq"
export mm2_gn_preset="splice:hq"
export mm2_n115_flag="-N 115"
