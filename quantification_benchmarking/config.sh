#!/usr/bin/env bash
# config.sh: pulls in the shared root config. Exists so that `scripts/*.sh`
# (which do `source ./config.sh`, resolved relative to the caller's cwd) keep
# working whether they're invoked from here (by build_indices.sh /
# build_reduced_annotation_indices.sh) or from a pipeline's benchmarking/
# directory (which has its own thin config.sh doing the same thing, one
# level deeper). No overrides needed here — index-building never touches the
# platform-specific values (minimap2 presets) that pipeline.sh runs do.
source ../config.sh
