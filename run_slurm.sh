#!/bin/bash
# Run a compositing script on Slurm, one Julia worker per Slurm task, on any of the job's nodes.
#   sbatch run_slurm.sh ekman_mse.jl                                     # 25 pentads: 27 tasks
#   sbatch --nodes=1 --ntasks-per-node=8 run_slurm.sh bsiso_phase.jl     # 8 phases:   8 tasks
# Options given to sbatch override the #SBATCH defaults below. Each worker takes whole
# pentads (or phases), so tasks beyond their number sit idle.
#SBATCH --job-name=ekman
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=9
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
##SBATCH --partition=...
##SBATCH --account=...
# --mem is per node: ~1 GB per worker, plus ~6 GB on the first node for the main
# process, which collects and stacks all the composites.

set -euo pipefail
script=${1:-ekman_mse.jl}
cd "$SLURM_SUBMIT_DIR"

# the same compile target as precompile_slurm.sh (head and compute nodes differ), taken from
# its export line so the two can't drift apart; a different target would rebuild the caches
eval "$(grep '^export JULIA_CPU_TARGET=' precompile_slurm.sh)"

# precompile ahead of time with ./precompile_slurm.sh on the head node; this only loads the
# cache (or builds anything stale, once, before any worker starts)
julia --project -e 'using EkmanCommon, SlurmClusterManager'

# the main process runs here; the script's addprocs(SlurmManager()) srun's the workers
julia --project "$script"
