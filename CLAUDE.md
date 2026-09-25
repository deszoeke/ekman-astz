# ekman-astz

`bsiso_phase.jl`: composites of ERA5 surface fields and air-sea Ekman terms, grouped by BSISO1 phase (1–8), 2012–2025.

`ekman_mse.jl`: April–July pentad climatology of the same fields and terms, 2012–2025.

`EkmanCommon/`: local package holding the code both scripts share: the physics, `daily_terms`, `composite` and `save_composite`. It's precompiled, with a PrecompileTools workload on a tiny synthetic grid, so workers load compiled code instead of recompiling it. Scripts load it with `@everywhere using EkmanCommon`. Editing its source triggers a recompile on the next load. Each script sets its own data root, `ceoasdir`, and passes `dir=era5dir(ceoasdir)`.

Each script defines its groups of days (phases or pentads) and calls `composite(groups)`. Run locally with `julia --project -p N <script>.jl` on Julia 1.13. On the HPC use Slurm: `sbatch run_slurm.sh <script>.jl`. Inside a job (`SLURM_JOB_ID` set), each script calls `addprocs(SlurmManager())` to start one worker per Slurm task across the job's nodes. Only the main process loads SlurmClusterManager. `pmap` gives each worker process whole groups, so there are no shared accumulators and no locks. Use processes, not threads, because NetCDF-C/HDF5 are not thread-safe.

Setup on a new machine: `julia --project -e 'using Pkg; Pkg.instantiate()'`. `Manifest.toml` is gitignored.

`plots.jl`: all plotting. It reads outputs (e.g. `BSISO_phase.txt`) and runs in the default environment, which has PythonPlot: `julia plots.jl`, without `--project`. Keep PythonPlot out of the project and out of the compositing scripts. Otherwise every worker process loads PythonPlot and CondaPkg, and they recompile and collide on CondaPkg's lock files.

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
