#!/bin/bash
# Precompile the project (packages and EkmanCommon) and its plotting environment, plot/.
# The Julia settings, including the CPU target shared by all x86_64 nodes, come from
# ~/.config/julia/arch_env.sh. Run on the head node after instantiating a new machine,
# changing packages, or editing EkmanCommon:
#   ./precompile.sh            # precompile what is stale
#   ./precompile.sh --check    # then check that a compute node loads without recompiling
# Run one precompile at a time, and not while Jupyter kernels are compiling.
# SRUN_ARGS (e.g. "-p mypartition -A myaccount") is passed to the check's srun.

set -euo pipefail
repo=$(cd "$(dirname "$0")" && pwd)
source ~/.config/julia/arch_env.sh

check=false
for arg in "$@"; do
    case $arg in
        --check) check=true ;;
        *) echo "unknown option: $arg" >&2; exit 1 ;;
    esac
done

julia --project="$repo" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project="$repo/plot" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

if $check; then
    # should print ok with no "Precompiling" messages
    srun ${SRUN_ARGS:-} -N1 -n1 julia --project="$repo" -e \
        'using EkmanCommon, SlurmClusterManager; println("compute node ", Sys.CPU_NAME, ": ok")'
fi
