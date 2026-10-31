#!/bin/bash
# Start Julia in this project with the same JULIA_CPU_TARGET as the Slurm scripts, so the REPL
# shares their precompiled caches (overrides any value exported in ~/.bashrc, which other
# projects may need). Usage: ./julia_hpc.sh [julia args]   e.g. ./julia_hpc.sh -p 4 script.jl
here=$(dirname "$0")
eval "$(grep '^export JULIA_CPU_TARGET=' "$here"/precompile_slurm.sh)"
exec julia --project="$here" "$@"
