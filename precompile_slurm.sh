#!/bin/bash
# Precompile the project (packages and EkmanCommon) for every node type on the HPC:
# the head node (znver2) and the Sandy Bridge compute nodes share ~/.julia, so the
# caches must hold code for both. Run on the head node, in or next to the project, after
# instantiating a new machine, changing packages, or editing EkmanCommon:
#   ./precompile_slurm.sh            # precompile what is stale
#   ./precompile_slurm.sh --clean    # first delete this Julia's compiled caches; needed after
#                                    # changing JULIA_CPU_TARGET, which Pkg's check doesn't see
#   ./precompile_slurm.sh --check    # afterwards, check that a compute node loads without recompiling
# Options combine. SRUN_ARGS (e.g. "-p mypartition -A myaccount") is passed to the check's srun.

set -euo pipefail
cd "$(dirname "$0")"

# compile for every node type; run_slurm.sh must use the same target
export JULIA_CPU_TARGET="generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1);x86-64-v4,-rdrnd,base(1)"

clean=false; check=false
for arg in "$@"; do
    case $arg in
        --clean) clean=true ;;
        --check) check=true ;;
        *) echo "unknown option: $arg" >&2; exit 1 ;;
    esac
done

if $clean; then
    # e.g. ~/.julia/compiled/v1.13; only caches, rebuilt below (other environments rebuild on next use)
    compiled=$(julia -e 'print(joinpath(DEPOT_PATH[1], "compiled", "v$(VERSION.major).$(VERSION.minor)"))')
    echo "removing $compiled"
    rm -rf "$compiled"
fi

julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

if $check; then
    # should print ok with no "Precompiling" messages
    srun ${SRUN_ARGS:-} -N1 -n1 julia --project -e \
        'using EkmanCommon, SlurmClusterManager; println("compute node ", Sys.CPU_NAME, ": ok")'
fi
