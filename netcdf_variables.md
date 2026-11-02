# NetCDF output variables

Both scripts write through `save_composite` in `EkmanCommon/src/EkmanCommon.jl`, so the variable set is the same. Only the third dimension (`dim`) and the file attributes differ.

Dimension order is (longitude, latitude, dim) as written from Julia. Tools that read the C order, such as `ncdump`, show it reversed.

Only `comp` and `nobs` are saved. Derived quantities (M, m_a, m_o, gradients, enthalpy fluxes) are computed from them in `plot/`.

## Common to both files

File attributes: `history` (`"<now()> <basename(PROGRAM_FILE)>"`; the script name is empty if run from the REPL), plus the script-specific ones below.

Every composite field is Float32 with `_FillValue` = NaN. It is the mean of the finite daily values at each grid point. It is NaN where the count is 0, for example SST over land.

Every `nobs_<k>` is Int16, with `long_name` "number of daily values in the mean of <k>".

| variable | dimensions | description |
|---|---|---|
| longitude | (longitude) | degrees_east, read from the ERA5 stress file (Float64) |
| latitude | (latitude) | degrees_north, read from the ERA5 stress file (Float64). Order is whatever the file has; the code says ERA5 is decreasing |
| *dim* (`phase` or `pentad`) | (*dim*) | Int32 group coordinate; see each script below |
| taux | (longitude, latitude, dim) | Zonal wind stress, N m-2. ERA5 `iews`, hourly averaged to daily. Sign convention is ERA5's; the code does not change it |
| tauy | (longitude, latitude, dim) | Meridional wind stress, N m-2 (ERA5 `inss`) |
| sst | (longitude, latitude, dim) | Sea surface temperature, K |
| t2 | (longitude, latitude, dim) | 2 m air temperature, K |
| d2 | (longitude, latitude, dim) | 2 m dewpoint temperature, K |
| q | (longitude, latitude, dim) | Specific humidity, kg kg-1. Computed from the hourly dewpoint (Buck 1981, p = 1000 hPa), then averaged to daily |
| adv_sst | (longitude, latitude, dim) | K kg m-2 s-1. +1 × (−M·∇sst), the advection of SST by the ocean Ekman mass transport M. Computed from daily means, then composited |
| adv_t2 | (longitude, latitude, dim) | K kg m-2 s-1. −1 × (−M·∇t2), advection of t2 by the atmospheric Ekman transport −M |
| adv_q | (longitude, latitude, dim) | kg m-2 s-1. −1 × (−M·∇q), advection of q by the atmospheric Ekman transport −M |
| sdiv_sst | (longitude, latitude, dim) | K kg m-2 s-1. +1 × sst·∇·M |
| sdiv_t2 | (longitude, latitude, dim) | K kg m-2 s-1. −1 × t2·∇·M |
| sdiv_q | (longitude, latitude, dim) | kg m-2 s-1. −1 × q·∇·M |
| sst_taux, sst_tauy | (longitude, latitude, dim) | K N m-2. sst × taux, sst × tauy (no sign), from daily means. The flux sst·M = (sst_tauy, −sst_taux)/f |
| t2_taux, t2_tauy | (longitude, latitude, dim) | K N m-2. t2 × taux, t2 × tauy |
| q_taux, q_tauy | (longitude, latitude, dim) | N m-2. q × taux, q × tauy |
| nobs_taux, nobs_tauy, nobs_sst, nobs_t2, nobs_d2, nobs_q, nobs_adv_sst, nobs_adv_t2, nobs_adv_q, nobs_sdiv_sst, nobs_sdiv_t2, nobs_sdiv_q, nobs_sst_taux, … nobs_q_tauy | (longitude, latitude, dim) | Int16 count of finite daily values in the mean of the matching field, at each grid point |

Notes on M and the signs:

- M = (Mx, My) = (tauy/f, −taux/f), the ocean Ekman mass transport in kg m-1 s-1, NOT curtailed at the equator. At f = 0 it is infinite, so the adv and sdiv terms there (and sdiv on the neighboring rows) are not finite and drop out of the means (nobs = 0, value NaN).
- Curtailment r = 1 − exp(−y²/Ro_rad²) depends only on latitude, so it is applied after compositing (`ekman_fields(d; Ro_rad)` in `plot/`): r·adv, −∇⋅(r M s) from the s·τ composites, and s∇⋅(rM) = r·adv + ∇⋅(r M s).
- The code comment says the flux-form tendency is `adv_s − sdiv_s`.
- The `units` strings are hard-coded in `save_composite`. For adv and sdiv they follow from M in kg m-1 s-1 times the units of s.
- Gradients and divergence are centered differences on the sphere. They are NaN at the first and last latitude. They are also NaN across longitude unless the grid is global (`isglobal`).

## bsiso_phase.jl

File name: `ekman_bsiso1_phase_<firstyear>-<lastyear>.nc`, e.g. `ekman_bsiso1_phase_2012-2020.nc`. It is written to the current directory and uses the years of `first(days)` and `last(days)`.

The group coordinate is `phase`, with values 1:8.

| variable | dimensions | description |
|---|---|---|
| phase | (phase) | Int32 1..8, `long_name` "BSISO1 phase". Only active days are composited, i.e. BSISO1 amplitude > `ampmin` and not missing |

File attributes: `title` = "ERA5 composites of surface fields and Ekman terms by BSISO1 phase", `days` = "<first> to <last>", `BSISO1_amplitude_min` = `ampmin` (default 1.0). All other variables are those in the common table, on (longitude, latitude, phase).

A side output, `BSISO_phase.txt`, is a text file and not NetCDF. It has the columns YEAR DOY MONTH DATE BSISO1_phase BSISO2_phase.

## ekman_mse.jl

File name: `ekman_pentad_clim_<firstyear>-<lastyear>.nc`, e.g. `ekman_pentad_clim_2012-2020.nc` for the default years.

The group coordinate is `pentad`, default 19:43. The values come from `pentad(d)`, which uses the 365-day calendar with Feb 29 in pentad 12.

| variable | dimensions | description |
|---|---|---|
| pentad | (pentad) | Int32 pentad of year (1-73). `long_name` is "pentad of year (1-73), 365-day calendar, Feb 29 in pentad 12" |

File attributes: `title` = "ERA5 pentad climatology of surface fields and Ekman terms", `years` = "<first>-<last>". All other variables are those in the common table, on (longitude, latitude, pentad). Each mean is over every day of that pentad in all the years.

## Other writers

`ekman_mse_test.jl` also calls `save_composite` (file `ekman_pentad<p>_<year>_test.nc`). I did not inspect its attributes or its coordinate. It uses the same variables.
