#!/bin/bash
# Use the Slurm scripts' JULIA_CPU_TARGET so Julia shares their precompiled caches
# (overrides any value in ~/.bashrc, which other projects may need).
#   ./julia_hpc.sh [julia args]   run Julia in this project, e.g. ./julia_hpc.sh -p 4 script.jl
#   source julia_hpc.sh           only export the target in the current shell
here=$(dirname "${BASH_SOURCE[0]}")
eval "$(grep '^export JULIA_CPU_TARGET=' "$here"/precompile_slurm.sh)"
(return 0 2>/dev/null) && return 0   # sourced: stop here
exec julia --project="$here" "$@"
