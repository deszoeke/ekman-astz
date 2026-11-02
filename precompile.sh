#!/bin/bash
# Precompile the project (packages and EkmanCommon) and its plotting environment, plot/.
# Runs only as a Slurm job: the head node's per-process limits (8 GB virtual memory, 300 s CPU
# time) make compiles segfault or get killed
# (e.g. ZMQ for IJulia, ProcessSignaled(11)). From the repository directory:
#   sbatch precompile.sh            # precompile what is stale
#   sbatch precompile.sh --check    # then check that a compute node loads without recompiling
# Output goes to precompile-<jobid>.out. Run after instantiating a new machine, changing
# packages, or editing EkmanCommon; one precompile at a time, and not while Jupyter kernels
# are compiling. The Julia settings, including the CPU target shared by all x86_64 nodes,
# come from ~/.config/julia/arch_env.sh.
#SBATCH --job-name=precompile
#SBATCH --partition=ceoas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.out

set -euo pipefail
if [[ -z ${SLURM_JOB_ID:-} ]]; then
    echo "run as a Slurm job, from the repository directory: sbatch precompile.sh" >&2
    exit 1
fi
repo=$SLURM_SUBMIT_DIR
source ~/.config/julia/arch_env.sh
# one parallel compile per allocated CPU; Julia would otherwise count all the node's CPUs
export JULIA_NUM_PRECOMPILE_TASKS=${SLURM_CPUS_PER_TASK:-4}

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
    srun -N1 -n1 julia --project="$repo" -e \
        'using EkmanCommon, SlurmClusterManager; println("compute node ", Sys.CPU_NAME, ": ok")'
fi
