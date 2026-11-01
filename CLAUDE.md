# ekman-astz

`bsiso_phase.jl`: composites of ERA5 surface fields and air-sea Ekman terms, grouped by BSISO1 phase (1–8), 2012–2025.

`ekman_mse.jl`: April–July pentad climatology of the same fields and terms, 2012–2025.

`EkmanCommon/`: local package holding the code both scripts share: the physics, `daily_terms`, `composite` and `save_composite`. It's precompiled, with a PrecompileTools workload on a tiny synthetic grid, so workers load compiled code instead of recompiling it. Scripts load it with `@everywhere using EkmanCommon`. Editing its source triggers a recompile on the next load. Each script sets its own data root, `ceoasdir`, and passes `dir=era5dir(ceoasdir)`.

Each script defines its helpers and groups of days (phases or pentads) as functions, and `main(; kw...)` calls `composite(groups)`, saves the file and returns `(comp, nobs, lon, lat, terms)`. `main()` runs only when the file is run as a script, not on `include`; `derived_terms(comp, lon, lat)` computes the enthalpy/MSE terms from a composite. Run locally with `julia --project -p N <script>.jl` on Julia 1.13. On the HPC use Slurm: `sbatch run_slurm.sh <script>.jl`. Inside a job (`SLURM_JOB_ID` set), each script calls `addprocs(SlurmManager())` to start one worker per Slurm task across the job's nodes. Only the main process loads SlurmClusterManager.

HPC nodes differ in CPU: the head node is znver2, the compute nodes Sandy Bridge, and they share `~/.julia`. So everything is precompiled for both with one `JULIA_CPU_TARGET`.

`env.sh` is the only place the project sets environment variables (`JULIA_CPU_TARGET`, `JULIA_PROJECT`). `precompile.sh` and `run_slurm.sh` source it; in your own shell, `source env.sh` (only this file is sourced; run the others). A different target, e.g. `haswell,-rdrnd` from `~/.bashrc`, builds a second cache and crashes on the Sandy Bridge nodes.

| To | Run |
|---|---|
| use the REPL | `source env.sh`, then `julia` |
| run a script locally | `source env.sh`, then `julia -p N script.jl` |
| test in the REPL | `source env.sh`, `julia`, `include("ekman_mse.jl")`; for workers `addprocs(4); @everywhere using EkmanCommon`; then e.g. `main(years=2012:2012, pentads=19:20)` |
| rebuild caches after package or EkmanCommon changes (head node) | `./precompile.sh` |
| run on Slurm | `sbatch run_slurm.sh script.jl` |

Use `./precompile.sh --clean` after changing the target in `env.sh`: Pkg's up-to-date check ignores the CPU target and keeps stale caches, which gives "Precompiled image … not available with flags".

`pmap` gives each worker process whole groups, so there are no shared accumulators and no locks. Use processes, not threads, because NetCDF-C/HDF5 are not thread-safe.

Setup on a new machine: `julia --project -e 'using Pkg; Pkg.instantiate()'`. `Manifest.toml` is gitignored.

`plot/`: all plotting, with its own `plot/Project.toml` (PythonPlot, NCDatasets) and CondaPkg in `plot/.CondaPkg`. Run from `plot/`, e.g. `julia --project=. ekman_V_ma_mo_plot.jl ../file.nc k`.
- `ekman_V_ma_mo_plot.jl` maps group k of a composite NetCDF: M, m_a, m_o and their gradients, from the composite means (exact, since these are linear); and the advection, s∇⋅M and flux-form terms from the saved `adv_*`/`sdiv_*` composites of daily products.
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
