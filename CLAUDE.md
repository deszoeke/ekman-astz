# ekman-astz

`bsiso_phase.jl`: composites of ERA5 surface fields and air-sea Ekman terms, grouped by BSISO1 phase (1–8), 2012–2025, with OISST SST.

`ekman_mse.jl`: April–July pentad climatology of the same fields and terms, 2012–2025.

`EkmanCommon/`: local package holding the code both scripts share: the physics, `daily_terms`, `composite` and `save_composite`. It's precompiled, with a PrecompileTools workload on a tiny synthetic grid, so workers load compiled code instead of recompiling it. Scripts load it with `@everywhere using EkmanCommon`. Editing its source triggers a recompile on the next load. Each script sets its own data root, `ceoasdir`, and passes `dir=era5dir(ceoasdir)`.

Each script defines its helpers and groups of days (phases or pentads) as functions, and `main(; kw...)` calls `composite(groups)`, saves the file and returns `(comp, nobs, lon, lat)`. `main()` runs only when the file is run as a script, not on `include`. Enthalpy terms (h_o, h_a) are computed where they are plotted (`plot/`), from the saved composites. Run locally with `julia --project -p N <script>.jl` on Julia 1.13. On the HPC use Slurm: `sbatch run_slurm.sh <script>.jl`. Inside a job (`SLURM_JOB_ID` set), each script calls `addprocs(SlurmManager())` to start one worker per Slurm task across the job's nodes. Only the main process loads SlurmClusterManager.

### Julia on the HPC
Julia settings (depot, juliaup, CPU target per architecture) are user-level, in `~/.config/julia/arch_env.sh`, sourced by `~/.bashrc`, `precompile.sh`, `run_slurm.sh` and the Jupyter kernel launcher. The project sets no environment variables. Targets, caches, kernels and the don'ts are in the user-level `hpc-julia` skill.

| To | Run |
|---|---|
| use the REPL | `julia --project` (in the repo) |
| run a script locally | `julia --project -p N script.jl` |
| test in the REPL | `julia --project`, `include("ekman_mse.jl")`; for workers `addprocs(4); @everywhere using EkmanCommon`; then e.g. `main(years=2012:2012, pentads=19:20)` |
| rebuild caches after package or EkmanCommon changes | `sbatch precompile.sh` (project and `plot/`; never on the head node) |
| run on Slurm | `sbatch run_slurm.sh script.jl` |

Jupyter on the HPC: kernel "Julia" (`~/.local/share/jupyter/kernels/julia/`, user-level).

When something recompiles or fails to load: `./check_julia.sh` (read-only, any node).

`pmap` gives each worker process whole groups, so there are no shared accumulators and no locks. Use processes, not threads, because NetCDF-C/HDF5 are not thread-safe.

Setup on a new machine: `julia --project -e 'using Pkg; Pkg.instantiate()'`. `Manifest.toml` is gitignored.

`plot/`: all plotting, in Jupyter notebooks with their own `plot/Project.toml`. **For any plot, figure or notebook task, read `plot/CLAUDE.md` first** (environment, the `composites.jl` and `maps.jl` libraries, how figures are written, and the notebooks).

Keep PythonPlot out of the project and out of the compositing scripts. Otherwise every worker process loads PythonPlot and CondaPkg, and they recompile and collide on CondaPkg's lock files.

## Who reads code
The main session runs on Opus (`.claude/settings.json`).
- **File and line already known:** read them directly with `sed` or `grep`; keep it short.
- **Anything else** (finding code, tracing calls, summarizing `.jl` files or variables, reading large files): send it to the `code-reader` agent (Sonnet, read-only, `.claude/agents/code-reader.md`). It reports facts with file:line references.
- **Judgment calls, planning, editing and writing code** stay with the main session. Don't give them to `code-reader`.

## Status
- Tested only on synthetic ERA5-format files; not yet run on real data.
- `rdadir` (NCAR RDA d633000 path for the iews/inss wind stress) is a placeholder relative path; point it at the real copy.
- Surface flux data on the HPC: OAFlux daily (WHOI, 1°) in `data/oaflux/daily/` (`oafluxdir`), with turbulent fluxes and evaporation 1985–2022 and radiation (ISCCP) and qnet 1985–2009; CERES SYN1deg-Day Ed4.2 adjusted all-sky surface and TOA fluxes in `data/ceres/adjusted/` (`ceresdir`; a symlink to `/storage/ceoas-datasets/datasets/ERA5/staging/simon/ceres/adjusted`), March 2000–August 2026 (download list `data/ceres/urls.txt`, `wget -c -i urls.txt`). Data sets (`main(data=...)` or the command-line word): `oaflux` (all variables, 1985–2009), `oafluxturb` (turbulent fields) and `ceres`, matched, 2000-04-01 to 2022-12-31. EkmanCommon `composite_data` picks the source. Each data set is composited in its own run and file (`<prefix>_pentad_clim_<first_day>-<last_day>.nc`, `<prefix>_bsiso1_phase_<first_day>-<last_day>.nc`, dates yyyymmdd, so different epochs go to different files). Turbulent (OAFlux) and radiative (CERES) terms over the same days are combined after compositing, e.g. the net surface heat flux; the file attributes `first_day` and `last_day` give the period of each composite, so files of matching epochs can be checked before their terms are added. Not yet run on the real files.
- SST, t2m and d2m come from `/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC/{sst,t2m,d2m}`.
- t2m and d2m are extended through 2026 in `/storage/ceoas-datasets/datasets/ERA5/staging/simon/{t2,d2}` (point the `data/t2`, `data/d2` symlinks there). Both scripts take any of these command-line words, in any order (`script_args`): a data set, `era5` (default), `oaflux`, `oafluxturb` or `ceres`; an SST product for `era5`, `oisst` (default) or `e5sst`; a first and last year; and `peryear` (one file per year), e.g. `ekman_mse.jl e5sst 2012 2020` or `ekman_mse.jl 2023 2025 peryear`. `oisst` is daily NOAA OISST v2.1 (`data/oisst/YYYY/oisst-avhrr-v02r01.YYYYMMDD.nc`, `oisstdir`, 1981–2026) in place of ERA5 SST, which is then not read: each ERA5 grid point, the corner of four OISST cells offset by 1/8°, gets the mean of those that are ocean (`oisst_sst`). Defaults: `oisst` 2012–2025; `e5sst` 2012–2020 (Andrea's SST stops in 2020). Files are named `<data>_<grouping>_<sst>_<period>.nc`, e.g. `ekman_pentad_clim_oisst_20120401-20250803.nc`, `ekman_bsiso1_phase_e5sst_20120101-20201231.nc`; data sets that use no SST product (`oaflux`, `oafluxturb`, `ceres`) have no `<sst>_`.
- Compositing uses the uncurtailed Ekman transport. The equatorial curtailment r = `curtail(lat, Ro_rad)` = 1 − exp(−y²/Ro_rad²) depends only on latitude, so it is applied in postprocessing (`ekman_fields(d; Ro_rad)`; the notebook sets `Ro_rad = 250e3` m, an ad hoc plotting choice). Results should not depend on it; it matters only within ~5° of the equator.

## Conventions
- Ocean Ekman mass transport is M = (τy, −τx)/f; the atmosphere's is −M. Signs in `scalars` encode this.
- Nonlinear terms (`adv_s`, `sdiv_s`, `s_taux`, `s_tauy` for s in sst, t2, q) are computed from DAILY means inside the loop, then composited. The flux-form tendency is `adv_s − sdiv_s` = −∇⋅(±M s), with M s = (s_tauy, −s_taux)/f.
- q is computed from HOURLY dewpoint (Buck 1981, p = 1000 hPa), then averaged to daily.
- h_a = c_pa T2 + Lv q (surface-air moist enthalpy, i.e. surface MSE) and h_o = c_po SST use temperatures in K, i.e. measured from 0 K.
- ASTZ: −M_o⋅∇(h_o − h_a) = −M_o⋅∇(h_o − h_*) − M_o⋅∇h'_a. `adv_star` = −M_o⋅∇h_* is linearized (c_pa + Lv q_sat′ at the composite-mean SST times the composite of daily −M_o⋅∇SST; −M⋅∇q_sat(SST) is not composited); `adv_ha_prime` is the remainder, so the split sums exactly. d(h_o − h_*)/dSST = c_po − c_pa − Lv q_sat′ changes sign near 25 °C.
- BSISO composites use only active days, where BSISO1 amplitude √(PC1² + PC2²) > `ampmin` = 1. Index days with fill values (−999.9, e.g. 2015-06-21) are missing: they get phase 0 in `BSISO_phase.txt` and are never composited.
