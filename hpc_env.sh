# Source this (don't run it) before using julia interactively on the HPC:  source hpc_env.sh
# Sets JULIA_CPU_TARGET to the Slurm scripts' value (shared precompiled caches) and
# JULIA_PROJECT to this directory, so plain `julia` works. Only sets variables.
_here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
eval "$(grep '^export JULIA_CPU_TARGET=' "$_here"/precompile_slurm.sh)"
export JULIA_PROJECT=$_here
unset _here
