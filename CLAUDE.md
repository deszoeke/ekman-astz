# ekman-astz

`bsiso_phase.jl`: composites of ERA5 surface fields and air-sea Ekman terms, grouped by BSISO1 phase (1–8), 2012–2026.

`ekman_mse.jl`: April–July pentad climatology of the same fields and terms. Run with `julia --project -p N ekman_mse.jl`: `pmap` gives each worker process whole pentads, so there are no shared accumulators and no locks. Use processes, not threads, because NetCDF-C/HDF5 are not thread-safe.

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
