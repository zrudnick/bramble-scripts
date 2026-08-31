#!/usr/bin/env bash
# config.sh: ONT reduced-annotation pipeline configuration — pulls in the
# shared config, then sets the platform-specific overrides ONT needs (same
# presets as ont/, differ from PacBio).
source ../../../config.sh

# minimap2 mapping presets — differ between ONT and PacBio, and are only used
# at alignment time (never index-build time), so they live here rather than
# in the shared config.
export mm2_tx_preset="map-ont"
export mm2_gn_preset="splice"
export mm2_n115_flag="-N 115"
