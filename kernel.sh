#!/bin/bash
# Jupyter kernel for this project: Julia with env.sh's settings. jupyter_kernel.sh installs a
# kernel spec that runs this script, so a kernel started by JupyterHub (which doesn't read
# ~/.bashrc or env.sh) compiles for the same JULIA_CPU_TARGET as the shells and Slurm jobs.
# A different target overwrites the shared caches in ~/.julia/compiled, because their file names
# don't include the target.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/env.sh"
unset JULIA_PROJECT     # --project=@. finds the notebook's Project.toml (plot/)
exec "${JULIA:-$HOME/.juliaup/bin/julia}" -i --color=yes --project=@. \
    -e 'import IJulia; IJulia.run_kernel()' "$@"
