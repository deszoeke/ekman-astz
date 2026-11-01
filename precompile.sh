#!/bin/bash
# Precompile the project (packages and EkmanCommon) and plot/ for every node type on the HPC:
# the head node (znver2) and the Sandy Bridge compute nodes share ~/.julia, so the
# caches must hold code for both. Run on the head node, in or next to the project, after
# instantiating a new machine, changing packages, or editing EkmanCommon:
#   ./precompile.sh            # precompile what is stale
#   ./precompile.sh --clean    # first delete this Julia's compiled caches; needed after
#                                    # changing JULIA_CPU_TARGET, which Pkg's check doesn't see
#   ./precompile.sh --check    # afterwards, check that a compute node loads without recompiling
# Options combine. SRUN_ARGS (e.g. "-p mypartition -A myaccount") is passed to the check's srun.

set -euo pipefail
cd "$(dirname "$0")"

# JULIA_CPU_TARGET and JULIA_PROJECT come from env.sh
source ./env.sh

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

julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
# the plotting environment, which the Jupyter kernel (kernel.sh) uses
(cd plot && JULIA_PROJECT=. julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()')

if $check; then
    # should print ok with no "Precompiling" messages
    srun ${SRUN_ARGS:-} -N1 -n1 julia -e \
        'using EkmanCommon, SlurmClusterManager; println("compute node ", Sys.CPU_NAME, ": ok")'
fi
