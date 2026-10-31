#!/bin/bash
# Run Julia in this project with the Slurm JULIA_CPU_TARGET:  ./julia_hpc.sh [julia args]
# (for plain `julia` in your shell, `source hpc_env.sh` instead)
source "$(dirname "${BASH_SOURCE[0]}")"/hpc_env.sh
exec julia "$@"
