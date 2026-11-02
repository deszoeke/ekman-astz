#!/bin/bash
# Read-only check of the Julia setup: which target, depot, julia and kernels this process
# would use, and which packages of the project and plot/ would recompile on this node.
# Changes nothing. Run on any node:
#   ./check_julia.sh
# From a notebook, to see the kernel's own environment:
#   run(`bash /home/ceoas/deszoeks/projects/ekman-astz/check_julia.sh`)
set -u
repo=$(cd "$(dirname "$0")" && pwd)

echo "== node"
echo "host $(hostname), arch $(uname -m), Slurm job ${SLURM_JOB_ID:-none}"
echo "ulimit -v $(ulimit -v) KB, ulimit -t $(ulimit -t) s"
if [[ -z ${SLURM_JOB_ID:-} && $(ulimit -v) != unlimited ]]; then
    echo "WARNING: memory-limited and not in a Slurm job (head node?): compile only with sbatch precompile.sh"
fi

echo "== environment"
for v in JULIA_CPU_TARGET JULIA_DEPOT_PATH JULIAUP_DEPOT_PATH JULIA_PROJECT JULIA_LOAD_PATH IJULIA_NODEFAULTKERNEL; do
    echo "$v=${!v-<unset>}"
done
if [[ $(uname -m) == x86_64 && -z ${JULIA_CPU_TARGET:-} ]]; then
    echo "WARNING: JULIA_CPU_TARGET unset on x86_64: caches would be built as \"native\"; source ~/.config/julia/arch_env.sh"
fi
echo "julia on PATH: $(command -v julia || echo none)"

echo "== Jupyter kernel specs"
kernels=${JUPYTER_DATA_DIR:-$HOME/.local/share/jupyter}/kernels
for k in "$kernels"/*/kernel.json; do
    [[ -e $k ]] || { echo "none in $kernels"; break; }
    echo "-- $k"
    cat "$k"; echo
    case $k in */julia-[0-9]*) echo "WARNING: IJulia's default kernel; it sets no target (native)";; esac
done

command -v julia >/dev/null || exit 0
echo "== precompile status on this node (packages that would recompile)"
for p in "$repo" "$repo/plot"; do
    julia --startup-file=no --project="$p" -e '
        using Pkg
        println("-- ", Base.active_project(), ": CPU ", Sys.CPU_NAME,
                ", target ", get(ENV, "JULIA_CPU_TARGET", "native"))
        deps = try
            Pkg.dependencies()
        catch err
            println("cannot read the manifest: ", sprint(showerror, err)); exit()
        end
        isstd(d) = startswith(d.source, Sys.STDLIB)
        stale(f) = sort([d.name for (u, d) in deps if f(d) && !Base.isprecompiled(Base.PkgId(u, d.name))])
        pk, sl = stale(!isstd), stale(isstd)
        npk = count(!isstd, values(deps))
        println(npk - length(pk), " of ", npk, " packages compiled for this node",
                isempty(pk) ? "" : "; would recompile: " * join(pk, ", "))
        isempty(sl) || println("standard libraries that would compile into the depot: ", join(sl, ", "))'
done
