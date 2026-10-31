# The only place this project sets environment variables. Source it (don't run it):
#   source env.sh
# Every entry point (precompile.sh, run_slurm.sh, your interactive shell) sources this, so
# the precompiled caches and the runs always use the same settings.
#
# JULIA_CPU_TARGET: the head node (znver2) and compute nodes (Sandy Bridge) share ~/.julia,
#   so compile for both. Changing this needs ./precompile.sh --clean.
# JULIA_PROJECT: this directory, so plain `julia` uses the project without --project.
export JULIA_CPU_TARGET="generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1);x86-64-v4,-rdrnd,base(1)"
export JULIA_PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
