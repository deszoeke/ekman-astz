#!/bin/bash
# Run a compositing script on Slurm, one Julia worker per Slurm task, on any of the job's nodes.
#   sbatch run_slurm.sh ekman_mse.jl                                     # 25 pentads: 27 tasks
#   sbatch --nodes=1 --ntasks-per-node=8 run_slurm.sh bsiso_phase.jl     # 8 phases:   8 tasks
#   sbatch run_slurm.sh ekman_mse.jl oaflux                              # OAFlux 1985-2009 instead of ERA5
# Options given to sbatch override the #SBATCH defaults below. Each worker takes whole
# pentads (or phases), so tasks beyond their number sit idle.
#SBATCH --job-name=ekman
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=9
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=ceoas
##SBATCH --account=...
# --mem=32G is the allowance per node. Estimated need for ekman_mse.jl (25 pentads):
# ~2.2 GB per worker (~1 GB, plus ~1.2 GB of netCDF chunk cache, EkmanCommon.chunk_cache),
# plus ~9 GB for the main process on the first node, which collects and stacks all the
# composites: ~29 GB on the first node, ~20 GB on the others.

set -euo pipefail
script=${1:-ekman_mse.jl}
cd "$SLURM_SUBMIT_DIR"

# Julia settings for this node's architecture (CPU target, depot, juliaup), the same as
# precompile.sh used
source ~/.config/julia/arch_env.sh

# precompile ahead of time with sbatch precompile.sh; this only loads the
# cache (or builds anything stale, once, before any worker starts)
julia --project=. -e 'using EkmanCommon, SlurmClusterManager'

# the main process runs here; the script's addprocs(SlurmManager()) srun's the workers,
# which get this project from SlurmManager (JULIA_PROJECT = the active project)
julia --project=. "$script" "${@:2}" # further arguments go to the script, e.g. oaflux
