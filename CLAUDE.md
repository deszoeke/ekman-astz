# ekman-astz

`bsiso_phase.jl`: composites of ERA5 surface fields and air-sea Ekman terms, grouped by BSISO1 phase (1–8), 2012–2025.

`ekman_mse.jl`: April–July pentad climatology of the same fields and terms, 2012–2025.

`EkmanCommon/`: local package holding the code both scripts share: the physics, `daily_terms`, `composite` and `save_composite`. It's precompiled, with a PrecompileTools workload on a tiny synthetic grid, so workers load compiled code instead of recompiling it. Scripts load it with `@everywhere using EkmanCommon`. Editing its source triggers a recompile on the next load. Each script sets its own data root, `ceoasdir`, and passes `dir=era5dir(ceoasdir)`.

Each script defines its helpers and groups of days (phases or pentads) as functions, and `main(; kw...)` calls `composite(groups)`, saves the file and returns `(comp, nobs, lon, lat, terms)`. `main()` runs only when the file is run as a script, not on `include`; `derived_terms(comp, lon, lat)` computes the enthalpy/MSE terms from a composite. Run locally with `julia --project -p N <script>.jl` on Julia 1.13. On the HPC use Slurm: `sbatch run_slurm.sh <script>.jl`. Inside a job (`SLURM_JOB_ID` set), each script calls `addprocs(SlurmManager())` to start one worker per Slurm task across the job's nodes. Only the main process loads SlurmClusterManager.

### Julia settings on the HPC
The project sets no environment variables. Julia settings are user-level, in `~/.config/julia/arch_env.sh` (one copy). `~/.bashrc`, `precompile.sh`, `run_slurm.sh` and the Jupyter kernel launcher source it. For each architecture it sets the depot, the juliaup directory and the CPU target:
- x86_64 (Sandy Bridge compute nodes, znver2 head node, znver4 chen): depot `~/.julia`, target `sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(0);x86-64-v4,-rdrnd,base(0)`, one string for all x86_64 nodes.
- aarch64 (GH200): depot `~/.julia-aarch64`, target unset (native; one CPU type).

How targets and caches work (Julia 1.13 `base/loading.jl`):
- A cache file's name hashes the target string with the project path, binary and flags (`compilecache_path`, ~3353). So different target strings, or different projects, never share a file.
- An unset target hashes as "native", the same name on every CPU. Native builds on different nodes overwrite each other.
- At load, Julia accepts any cache whose clones run on this CPU, whatever its target string. So a cache built with another target can satisfy `Pkg.precompile` on one node and fail on another.
- Rule: never leave the target unset on x86_64, and use one string.

monsoon-breeze still uses its own `helpers/arch_env.sh` (`haswell,-rdrnd`) for its jobs. Its caches are separate (different project path).

| To | Run |
|---|---|
| use the REPL | `julia --project` (in the repo) |
| run a script locally | `julia --project -p N script.jl` |
| test in the REPL | `julia --project`, `include("ekman_mse.jl")`; for workers `addprocs(4); @everywhere using EkmanCommon`; then e.g. `main(years=2012:2012, pentads=19:20)` |
| rebuild caches after package or EkmanCommon changes (head node) | `./precompile.sh` (project and `plot/`) |
| run on Slurm | `sbatch run_slurm.sh script.jl` |

Workers get the project from Distributed/SlurmManager, which pass `JULIA_PROJECT` = the active project.

Run one precompile at a time, never while Jupyter kernels are compiling. Never `rm -rf ~/.julia/compiled/...` while any Julia is running: it deletes every project's caches, and open files leave `.nfs*` placeholders.

Jupyter on the HPC (JupyterHub on chen): kernel "Julia", `~/.local/share/jupyter/kernels/julia/` (user-level, not in this repo). Its `kernel.json` runs `kernel.sh` in the same directory, which sources `~/.config/julia/arch_env.sh` and starts IJulia with `--project=@.`. IJulia's own `julia-1.x` kernels set no target (native); don't use them.

`pmap` gives each worker process whole groups, so there are no shared accumulators and no locks. Use processes, not threads, because NetCDF-C/HDF5 are not thread-safe.

Setup on a new machine: `julia --project -e 'using Pkg; Pkg.instantiate()'`. `Manifest.toml` is gitignored.

`plot/`: all plotting, with its own `plot/Project.toml` (PythonPlot, NCDatasets, Revise, IJulia, and EkmanCommon developed from `../EkmanCommon`) and CondaPkg in `plot/.CondaPkg`. Plots are made interactively in Jupyter notebooks, using the "Julia" kernel (`~/Library/Jupyter/kernels/julia` on the Mac; on the HPC see above). It runs the juliaup launcher, so it follows Julia updates, with `--project=@.`, so it finds `plot/Project.toml`. Don't pin Julia versions in kernels or manifests; let Julia update, recompile and re-resolve manifests itself. Each notebook loads a plotting library with `using Revise; includet(...)`. Libraries hold plot defaults and reusable functions, and no top-level plotting. Notebooks keep the data being plotted in variables, and adjust appearance there. Plot code takes constants and analysis functions from EkmanCommon, not local copies; short formulas (e.g. m_a = c_pa T + Lv q) may be written out in library functions for readability.
- `ekman_V_ma_mo_plot.jl` (library), `ekman_plots.ipynb` (notebook): maps of group k of a composite NetCDF. M, m_a, m_o and their gradients come from the composite means (exact, since these are linear). The advection, s∇⋅M and flux-form terms come from the saved `adv_*`/`sdiv_*` composites of daily products. Each figure is a list of `panel`s passed to `plot_panels`.
- `plot_bsiso_phase.jl` plots the phases in `BSISO_phase.txt`.

Keep PythonPlot out of the project and out of the compositing scripts. Otherwise every worker process loads PythonPlot and CondaPkg, and they recompile and collide on CondaPkg's lock files.

## Who reads code
The main session runs on Opus (`.claude/settings.json`).
- **File and line already known:** read them directly with `sed` or `grep`; keep it short.
- **Anything else** (finding code, tracing calls, summarizing `.jl` files or variables, reading large files): send it to the `code-reader` agent (Sonnet, read-only, `.claude/agents/code-reader.md`). It reports facts with file:line references.
- **Judgment calls, planning, editing and writing code** stay with the main session. Don't give them to `code-reader`.

## Status
- Tested only on synthetic ERA5-format files; not yet run on real data.
- `rdadir` (NCAR RDA d633000 path for the iews/inss wind stress) is a placeholder relative path; point it at the real copy.
- SST, t2m and d2m come from `/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC/{sst,t2m,d2m}`.
- `Ro_rad = 250e3` m (equatorial curtailment of the stress) is a placeholder, not a derived value.

## Conventions
- Ocean Ekman mass transport is M = (τy, −τx)/f; the atmosphere's is −M. Signs in `scalars` encode this.
- Nonlinear terms (`adv_s`, `sdiv_s` for s in sst, t2, q) are computed from DAILY means inside the loop, then composited. The flux-form tendency is `adv_s − sdiv_s`.
- q is computed from HOURLY dewpoint (Buck 1981, p = 1000 hPa), then averaged to daily.
- m_a and h_o use temperatures in K, i.e. measured from 0 K.
- BSISO composites use only active days, where BSISO1 amplitude √(PC1² + PC2²) > `ampmin` = 1. Index days with fill values (−999.9, e.g. 2015-06-21) are missing: they get phase 0 in `BSISO_phase.txt` and are never composited.
