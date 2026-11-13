"""
    EkmanCommon

Shared code for ekman_mse.jl and bsiso_phase.jl: ERA5 file access, thermodynamics,
Ekman transport, finite differences on the sphere, and compositing daily fields over
groups of days on worker processes, for ERA5 (`composite`, optionally with OISST SST), OAFlux (`oaflux_composite`) or CERES (`ceres_composite`); `composite_data` picks one. Load it on every process:
    @everywhere using EkmanCommon
"""
module EkmanCommon

using Dates
using Distributed
using NCDatasets
using PrecompileTools

export era5dir, era5code, era5file, month_files, month_blocks, chunk_cache, daily_terms, composite_mean, composite, save_composite,
       fieldkeys, nlkeys, keys_all, scalars, c_pa, c_po, rho_o, Lv, esat_buck, qsat_dew,
       coriolis, curtail, ekman_transport, ekman_transport_xy, gradient, divergence, advection,
       dailymean, accumulate!,
       oaflux_files, oaflux_turbulent, oafluxfile, oaflux_mean, oaflux_composite,
       ceres_vars, ceresfiles, ceres_mean, ceres_composite, pmap_groups, composite_data,
       oisstfile, oisst_sst, script_args

# ERA5 hourly surface analysis from NCAR RDA d633000, one file per variable per month:
#   <dir>/201201/e5.oper.an.sfc.128_229_iews.ll025sc.2012010100_2012013123.nc
# rdadir  = "d633000/e5.oper.an.sfc" # not used
# ECMWF parameter codes: iews, inss (N/m^2), sstk, 2t, 2d (K)
const era5code = (taux="128_229_iews", tauy="128_230_inss", sst="128_034_sstk",
                  t2="128_167_2t", d2="128_168_2d")

"""
    era5dir(root="./data")

Directories holding each variable's files under `root`, e.g.
"/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC": YYYYMM subdirectories
for stress, t2 and d2; SST files are grouped by month directly in root/sst/.
"""
era5dir(root="./data") = (taux=joinpath(root, "stress"), tauy=joinpath(root, "stress"),
                          sst=joinpath(root, "sst"),
                          t2=joinpath(root, "t2"), d2=joinpath(root, "d2"))

"ERA5 RDA file in directory dir for parameter code in the month containing date d"
function era5file(dir, code, d)
    m0 = firstdayofmonth(d)
    m1 = lastdayofmonth(d)
    ym = Dates.format(m0, "yyyymm")
    span = Dates.format(m0, "yyyymmdd") * "00_" * Dates.format(m1, "yyyymmdd") * "23"
    if code=="128_034_sstk"
	    joinpath(dir, "ERA5_SFC_sst_$(ym)_r1440x721_hr.nc")
    else
    	joinpath(dir, ym, "e5.oper.an.sfc.$(code).ll025sc.$(span).nc")
    end
end

# netCDF chunk cache (bytes) per file: one time-span of chunks, so that consecutive days
# reuse decompressed chunks instead of reading them again. Chunks (time, lat, lon): sst
# (124, 121, 240), ~520 MB per 124 h; RDA stress, t2, d2 (27, 139, 277), ~150 MB per 27 h.
const chunk_cache = (taux=160_000_000, tauy=160_000_000, sst=550_000_000,
                     t2=160_000_000, d2=160_000_000)

"the month's ERA5 files (NamedTuple like `dir`) for the month containing date d"
month_files(d; dir=era5dir(), code=era5code) = map((dk, ck) -> era5file(dk, ck, d), dir, code)

"sorted `days` split into vectors of consecutive days in the same month"
function month_blocks(days)
    ds = sort(collect(days))
    [filter(d -> firstdayofmonth(d) == m, ds) for m in unique(firstdayofmonth.(ds))]
end

"the data variable in an RDA file (the only variable with 3 dimensions)"
datavar(ds) = first(k for k in keys(ds) if ndims(ds[k]) == 3)

"the time coordinate variable name in a file (sst uses valid_time, others use time)"
timevar(ds) = haskey(ds, "valid_time") ? "valid_time" : "time"

"sorted indices as a UnitRange when contiguous (fast NetCDF hyperslab read), else left as-is"
asrange(idx) = !isempty(idx) && idx == first(idx):last(idx) ? (first(idx):last(idx)) : idx

"""
    withdatasets(f, files; cache=nothing)

Open each file in the NamedTuple `files` (fill values read as NaN), call `f` on the
NamedTuple of datasets, and close them all afterwards, like `NCDataset(file) do ds`
for several files. Files opened before a failure are closed too. `cache`, a NamedTuple
with the keys of `files`, gives each file's netCDF chunk cache size in bytes.
"""
function withdatasets(f, files::NamedTuple; cache=nothing)
    opened = NCDataset[]
    default = NCDatasets.get_chunk_cache()
    try
        for (k, file) in pairs(files)
            isnothing(cache) || NCDatasets.set_chunk_cache(size=cache[k]) # for files opened next
            push!(opened, NCDataset(file; maskingvalue=NaN))
        end
        f(NamedTuple{keys(files)}(Tuple(opened)))
    finally
        isnothing(cache) || NCDatasets.set_chunk_cache(size=default[1])
        foreach(close, opened)
    end
end

"average hourly data h (lon x lat x hour) to daily (lon x lat), applying f to each hour"
function dailymean(h, f=identity)
    x = zeros(size(h, 1), size(h, 2))
    for k in axes(h, 3), j in axes(h, 2), i in axes(h, 1)
        x[i,j] += f(h[i,j,k])
    end
    x ./= size(h, 3)
end

# thermodynamic constants
const c_pa = 1004.6  # J/kg/K specific heat of dry air at constant pressure
const c_po = 3985.0  # J/kg/K specific heat of seawater
const rho_o = 1025.0 # kg/m^3 reference density of seawater
const Lv   = 2.501e6 # J/kg latent heat of vaporization
const Rd   = 287.04  # J/kg/K
const Rv   = 461.5   # J/kg/K
const ε    = Rd/Rv
const p0   = 1000.0  # hPa, assumed surface pressure

"Buck (1981) saturation vapor pressure (hPa) over water, with enhancement factor, T in K, p in hPa"
function esat_buck(T, p=p0)
    Tc = T - 273.15
    f = 1.0007 + 3.46e-6*p
    f * 6.1121 * exp(17.502*Tc / (240.97 + Tc))
end

"specific humidity (kg/kg) from dewpoint Td (K) at pressure p (hPa)"
function qsat_dew(Td, p=p0)
    e = esat_buck(Td, p)
    ε*e / (p - (1 - ε)*e)
end

const Ω = 7.292e-5 # rad/s Earth rotation rate
const a = 6.371e6  # m Earth radius

"Coriolis parameter (1/s) at latitude lat (degrees)"
coriolis(lat) = 2Ω*sind(lat)

"""
    curtail(lat, Ro_rad)

Equatorial curtailment factor r = 1 - exp(-y^2/Ro_rad^2) at latitude lat (degrees), with
y the distance (m) from the equator: 0 at the equator, > 0.99 poleward of ~2.2 Ro_rad.
Depends only on latitude, so it can be applied to composites after compositing.
"""
curtail(lat, Ro_rad) = -expm1(-(a*deg2rad(lat)/Ro_rad)^2)

"""
    ekman_transport(taux, tauy, lat; Ro_rad=nothing)

Ekman mass transport (Mx, My) = (tauy, -taux)/f (kg/m/s) from wind stress (N/m^2)
at latitude lat (degrees). If Ro_rad (m) is given, curtail the effective stress
near the equator by `curtail(lat, Ro_rad)` to account for the pressure gradient balancing the wind stress and keep
the transport finite (-> 0 at the equator). Without Ro_rad the transport is
infinite (NaN for zero stress) at the equator.
"""
function ekman_transport(taux, tauy, lat; Ro_rad=nothing)
    f = coriolis(lat)
    if isnothing(Ro_rad)
        r = 1.0
    else
        f == 0 && return (0.0*taux, 0.0*tauy) # limit of r/f ~ y/(beta*Ro_rad^2) -> 0; keeps NaN stress NaN
        r = curtail(lat, Ro_rad)
    end
    (r*tauy/f, -r*taux/f)
end

"Ekman transport component arrays (Mx, My) from stress arrays (lon x lat x ...)"
function ekman_transport_xy(taux, tauy, lat; Ro_rad=nothing)
    M = ekman_transport.(taux, tauy, lat'; Ro_rad=Ro_rad)
    first.(M), last.(M)
end

"true if longitudes lon (degrees, uniform) span the whole globe"
isglobal(lon) = isapprox(mod(lon[end] - lon[1] + (lon[2] - lon[1]), 360), 0, atol=1e-6) ||
                isapprox(mod(lon[end] - lon[1] + (lon[2] - lon[1]), 360), 360, atol=1e-6)

"""
    gradient(φ, lon, lat; periodic=isglobal(lon))

Centered-difference gradient (∂φ/∂x, ∂φ/∂y) on the sphere of a scalar field φ
(lon x lat x ...), with lon, lat in degrees (lat may be decreasing, as in ERA5).
∂/∂x = ∂/∂λ / (a cos(lat)), ∂/∂y = ∂/∂lat / a. Longitude wraps if periodic;
otherwise, and at the first/last latitude and the poles, the gradient is NaN.
"""
function gradient(φ, lon, lat; periodic=isglobal(lon))
    nx, ny = size(φ, 1), size(φ, 2)
    dφdx = fill(NaN, size(φ))
    dφdy = fill(NaN, size(φ))
    # precompute per-i, per-j quantities that don't depend on the other index or on k
    dλ = [deg2rad(mod(lon[mod1(i+1, nx)] - lon[mod1(i-1, nx)], 360)) for i in 1:nx]
    c  = cosd.(lat)
    ady = [1 < j < ny ? a*deg2rad(lat[j+1] - lat[j-1]) : NaN for j in 1:ny]
    for k in CartesianIndices(size(φ)[3:end]), j in 1:ny, i in 1:nx
        # zonal
        if (periodic || 1 < i < nx) && c[j] > 1e-10
            ip, im = mod1(i+1, nx), mod1(i-1, nx)
            dφdx[i,j,k] = (φ[ip,j,k] - φ[im,j,k]) / (a*c[j]*dλ[i])
        end
        # meridional
        if 1 < j < ny
            dφdy[i,j,k] = (φ[i,j+1,k] - φ[i,j-1,k]) / ady[j]
        end
    end
    dφdx, dφdy
end

"""
    divergence(Fx, Fy, lon, lat; periodic=isglobal(lon))

Centered-difference divergence on the sphere of a vector field (Fx, Fy) (lon x lat x ...):
∇⋅F = ∂Fx/∂x + ∂(Fy cos(lat))/∂y / cos(lat). NaN where the gradient is NaN.
"""
function divergence(Fx, Fy, lon, lat; periodic=isglobal(lon))
    c = cosd.(lat')
    dFxdx, _ = gradient(Fx, lon, lat; periodic=periodic)
    _, dFydy = gradient(Fy .* c, lon, lat; periodic=periodic)
    @. dFxdx + dFydy / c
end

"advective tendency -(Fx ∂φ/∂x + Fy ∂φ/∂y) of scalar φ by transport (Fx, Fy)"
function advection(Fx, Fy, φ, lon, lat; periodic=isglobal(lon))
    dφdx, dφdy = gradient(φ, lon, lat; periodic=periodic)
    @. -(Fx*dφdx + Fy*dφdy)
end


# scalars s and the sign of the Ekman transport of their medium:
# ocean (sst) is carried by M, atmosphere (t2, q) by the equal and opposite -M
const scalars = (sst=1, t2=-1, q=-1)
# nonlinear terms, computed from daily means, then composited:
#   adv_s  = -(±M)⋅∇s   advection by Ekman transport of the medium
#   sdiv_s = s ∇⋅(±M)   s times Ekman transport divergence of the medium
#   s_taux, s_tauy      s times the stress (no sign), so that s M = (s_tauy, -s_taux)/f
# so the flux-form tendency is -∇⋅(±M s) = adv_s - sdiv_s. M is NOT curtailed at the
# equator here (Inf where f = 0, so those points drop out of the composites); curtailment
# by a function of latitude is applied afterwards, e.g. ekman_transport(...; Ro_rad).
const nlkeys = (Symbol.(:adv_,  keys(scalars))..., Symbol.(:sdiv_, keys(scalars))...,
                Symbol.(keys(scalars), :_taux)..., Symbol.(keys(scalars), :_tauy)...)
const fieldkeys = (keys(era5code)..., :q) # taux, tauy, sst, t2, d2, q
const keys_all  = (fieldkeys..., nlkeys...)

"add the finite values of x to the sums S, and count them in N"
function accumulate!(S, N, x)
    for i in eachindex(S, N, x)
        if isfinite(x[i])
            S[i] += x[i]
            N[i] += 1
        end
    end
end

"""
    daily_terms(ds, times, d, lon, lat; scalars=scalars, sst=nothing)

Daily means (lon x lat) on date d of the ERA5 fields, q from hourly dewpoint, and the
nonlinear Ekman terms computed from them (uncurtailed M), as a NamedTuple with keys `keys_all`.
`ds` are the month's open datasets (see `withdatasets`) and `times` their dates, one per hour.
`sst` (lon x lat, K), e.g. OISST, replaces the daily mean ERA5 SST, which `ds` then need not hold.
Returns nothing, with a warning, if hours for d are missing.
"""
function daily_terms(ds, times, d, lon, lat; scalars=scalars, sst=nothing)
    it = map(t -> asrange(findall(==(d), t)), times)
    if any(isempty, it)
        @warn "missing hours on $d"
        return nothing
    end
    # average hourly data to daily; q from hourly dewpoint, then daily mean
    hours(k) = ds[k][datavar(ds[k])][:, :, it[k]]
    hd2 = hours(:d2)
    x = (taux=dailymean(hours(:taux)), tauy=dailymean(hours(:tauy)),
         sst=isnothing(sst) ? dailymean(hours(:sst)) : sst,
         t2=dailymean(hours(:t2)), d2=dailymean(hd2), q=dailymean(hd2, qsat_dew))
    # daily Ekman transport (uncurtailed; Inf at f = 0) and its divergence
    Mx, My = ekman_transport_xy(x.taux, x.tauy, lat)
    divM = divergence(Mx, My, lon, lat)
    adv  = Tuple(sgn .* advection(Mx, My, x[s], lon, lat) for (s, sgn) in pairs(scalars))
    sdiv = Tuple(sgn .* divM .* x[s] for (s, sgn) in pairs(scalars))
    staux = Tuple(x[s] .* x.taux for s in keys(scalars))
    stauy = Tuple(x[s] .* x.tauy for s in keys(scalars))
    merge(x, NamedTuple{Symbol.(:adv_,  keys(scalars))}(adv),
             NamedTuple{Symbol.(:sdiv_, keys(scalars))}(sdiv),
             NamedTuple{Symbol.(keys(scalars), :_taux)}(staux),
             NamedTuple{Symbol.(keys(scalars), :_tauy)}(stauy))
end

"""
    composite_mean(days, lon, lat; progress=nothing, dir=era5dir(), code=era5code, scalars=scalars, oisst=nothing)

Mean over `days` of the finite daily values (see `daily_terms`). The days are read a month
at a time: each month's files are opened once, with a chunk cache (`chunk_cache`), and its
days read from them in order. Returns (mean, nobs), NamedTuples of lon x lat arrays; the
mean is NaN where nobs is 0. If `progress` is a channel, put 1 on it as each day finishes.
With `oisst` (the OISST root directory), SST is the daily OISST on the ERA5 grid (`oisst_sst`)
instead of ERA5's, whose files are then not read; a day without an OISST file is skipped.
"""
function composite_mean(days, lon, lat; progress=nothing, dir=era5dir(), code=era5code, scalars=scalars,
                        oisst=nothing)
    S = NamedTuple{keys_all}(Tuple(zeros(length(lon), length(lat)) for k in keys_all))
    N = NamedTuple{keys_all}(Tuple(zeros(Int16, length(lon), length(lat)) for k in keys_all))
    done(n) = isnothing(progress) || foreach(_ -> put!(progress, 1), 1:n)
    for mdays in month_blocks(days)
        files = month_files(first(mdays); dir, code)
        isnothing(oisst) || (files = Base.structdiff(files, NamedTuple{(:sst,)}))
        if !all(isfile, files)
            @warn "missing ERA5 files for $(Dates.format(first(mdays), "yyyy-mm")) ($(length(mdays)) days)" filter(!isfile, collect(files))
            done(length(mdays))
            continue
        end
        withdatasets(files; cache=chunk_cache) do ds
            # every source must be on the stress grid (sst and t2/d2 come from a different archive)
            for (k, f) in pairs(ds)
                f["latitude"][:] ≈ lat && f["longitude"][:] ≈ lon ||
                    error("$k grid in $(files[k]) differs from the stress grid")
            end
            times = map(f -> Date.(f[timevar(f)][:]), ds) # sst: valid_time
            for d in mdays
                sst = nothing
                if !isnothing(oisst)
                    f = oisstfile(oisst, d)
                    if !isfile(f)
                        @warn "missing OISST file $f"
                        done(1)
                        continue
                    end
                    sst = oisst_sst(f, lon, lat)
                end
                x = daily_terms(ds, times, d, lon, lat; scalars, sst)
                isnothing(x) || foreach(accumulate!, S, N, x)
                done(1)
            end
        end
    end
    map((s, n) -> s ./ n, S, N), N
end

"""
    composite(groups; dir=era5dir(), code=era5code, scalars=scalars, oisst=nothing)

Composite means (lon x lat x group) of the ERA5 fields and the nonlinear Ekman terms,
one for each vector of dates in `groups`. Means are taken over the finite daily values
at each grid point, so a missing hour or a land point drops out of that point's mean
instead of making it NaN. The groups are distributed over the worker processes, each
worker summing whole groups, so no two processes share an accumulator. Progress is a
counter of finished days, redrawn in place. With `oisst` (root directory), SST is the daily
OISST (see `composite_mean`). Returns (comp, nobs, lon, lat).
"""
function composite(groups; dir=era5dir(), code=era5code, scalars=scalars, oisst=nothing)
    lon, lat = NCDataset(era5file(dir.taux, code.taux, first(Iterators.flatten(groups)))) do ds
        Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:])
    end
    # workers report each finished day to the main process, which counts and prints
    progress = RemoteChannel(() -> Channel{Int}(256))
    total = sum(length, groups)
    printer = @async for (n, _) in enumerate(progress)
        print("\r$n/$total days"); flush(stdout)
    end
    r = pmap(days -> composite_mean(days, lon, lat; progress=progress,
                                    dir=dir, code=code, scalars=scalars, oisst=oisst),
             groups)
    close(progress); wait(printer); println()
    # stack the groups (lon x lat x group)
    comp = map(k -> stack(ri[1][k] for ri in r), NamedTuple{keys_all}(keys_all))
    nobs = map(k -> stack(ri[2][k] for ri in r), NamedTuple{keys_all}(keys_all))
    comp, nobs, lon, lat
end

# NOAA OISST v2.1 (AVHRR-only) daily SST, 0.25 degrees, one file per day:
#   <root>/YYYY/oisst-avhrr-v02r01.YYYYMMDD.nc, variable sst (lon, lat, zlev, time) in Celsius,
# on cell centres 0.125:0.25:359.875 E, -89.875:0.25:89.875 N, offset by 1/8 degree from ERA5's
# 0:0.25:359.75 E, 90:-0.25:-90 N.

"OISST v2.1 file for date d under root"
oisstfile(root, d) = joinpath(root, Dates.format(d, "yyyy"), "oisst-avhrr-v02r01.$(Dates.format(d, "yyyymmdd")).nc")

"""
    oisst_sst(file, lon, lat)

Daily OISST SST (K) in `file` on the ERA5 grid `lon`, `lat`: each ERA5 point is the corner of
four OISST cells, and gets the mean of those that are ocean (finite), so coastal points keep
the available ocean SST; NaN where none is. The grids are checked.
"""
function oisst_sst(file, lon, lat)
    ongrid(x, r) = length(x) == length(r) && x ≈ r
    ongrid(lon, 0:0.25:359.75) && ongrid(lat, 90:-0.25:-90) || error("oisst_sst: lon, lat are not the ERA5 0.25 degree grid")
    so, kelvin = NCDataset(file; maskingvalue=NaN) do ds
        ongrid(ds["lon"][:], 0.125:0.25:359.875) && ongrid(ds["lat"][:], -89.875:0.25:89.875) ||
            error("$file is not on the OISST 0.25 degree grid")
        v = ds["sst"]
        Float64.(v[:, :, 1, 1]), occursin(r"^(K|kelvin)$"i, get(v.attrib, "units", ""))
    end
    kelvin || (so .+= 273.15) # Celsius to K
    nx, ny = size(so)
    out = fill(NaN, length(lon), length(lat))
    for (j, φ) in enumerate(lat), (i, λ) in enumerate(lon)
        i1 = mod1(round(Int, λ/0.25), nx)    # OISST cells west and east of λ (wrapping at 0 E)
        k1 = round(Int, (φ + 90)/0.25)        # OISST cells south (k1) and north (k1 + 1) of φ
        s, n = 0.0, 0
        for ii in (i1, mod1(i1 + 1, nx)), kk in (k1, k1 + 1)
            if 1 <= kk <= ny && isfinite(so[ii, kk])
                s += so[ii, kk]; n += 1
            end
        end
        n > 0 && (out[i, j] = s/n)
    end
    out
end

# OAFlux daily 1-degree fields (WHOI), one file per variable per year under `root`, e.g.
#   <root>/turbulence/lh_oaflux_1985.nc, <root>/radiation/sw_isccp_1985.nc, <root>/netheat/qnet_1985.nc
# turbulence and evaporation 1985-2022; radiation (ISCCP) and qnet 1985-2009.
const oaflux_files = (lh="turbulence/lh_oaflux", sh="turbulence/sh_oaflux",
                      evapr="evaporation/evapr_oaflux",
                      qa="turbulence/qa_oaflux", ta="turbulence/ta_oaflux",
                      ts="turbulence/ts_oaflux", ws="turbulence/ws_oaflux",
                      lw="radiation/lw_isccp", sw="radiation/sw_isccp", qnet="netheat/qnet")
const oaflux_turbulent = (:lh, :sh, :evapr, :qa, :ta, :ts, :ws) # the variables through 2022

# CERES SYN1deg-Day Ed4.2 adjusted all-sky radiative fluxes (W m-2), 1-degree like OAFlux,
# in files of ~2.5 years each (CERES_SYN1deg-Day_*.nc under `ceresroot`), March 2000-August 2026.
# Short name => variable; lw down at the TOA (zero) is left out.
const ceres_vars = (sw_dn_toa="adj_atmos_sw_down_all_toa_daily", sw_up_toa="adj_atmos_sw_up_all_toa_daily",
                    lw_up_toa="adj_atmos_lw_up_all_toa_daily",
                    sw_dn_sfc="adj_atmos_sw_down_all_surface_daily", sw_up_sfc="adj_atmos_sw_up_all_surface_daily",
                    lw_dn_sfc="adj_atmos_lw_down_all_surface_daily", lw_up_sfc="adj_atmos_lw_up_all_surface_daily")

"OAFlux file of variable k (a key of oaflux_files) for year y under root"
oafluxfile(root, k, y) = joinpath(root, "$(oaflux_files[k])_$y.nc")

"CERES SYN1deg daily files under root, in time order"
ceresfiles(root) = sort(filter(f -> startswith(basename(f), "CERES_SYN1deg-Day") && endswith(f, ".nc"),
                               readdir(root; join=true)))

"the data variable in an OAFlux file (3 dimensions lon x lat x time; not an err field)"
oafluxvar(ds) = first(k for k in keys(ds) if ndims(ds[k]) == 3 && !startswith(lowercase(k), "err"))

"longitude and latitude of a lon x lat x time variable v, from its first two dimensions"
datagrid(ds, v) = (Float64.(ds[dimnames(v)[1]][:]), Float64.(ds[dimnames(v)[2]][:]))

"lon x lat zero sums and counts for each of vars"
zerosums(vars, lon, lat) = (NamedTuple{vars}(Tuple(zeros(length(lon), length(lat)) for k in vars)),
                            NamedTuple{vars}(Tuple(zeros(Int16, length(lon), length(lat)) for k in vars)))

"""
    oaflux_accumulate!(S, N, days, lon, lat; root, vars, done=n -> nothing)

Add the OAFlux daily fields `vars` on `days` to the sums S and counts N. Each file holds one
field per day of its year, in order. A year with any missing file is skipped with a warning.
`done(n)` is called with the number of days of each year as it finishes.
"""
function oaflux_accumulate!(S, N, days, lon, lat; root, vars, done=n -> nothing)
    for y in unique(year.(days))
        ydays = filter(d -> year(d) == y, days)
        files = map(k -> oafluxfile(root, k, y), NamedTuple{vars}(vars))
        if !all(isfile, files)
            @warn "missing OAFlux files for $y ($(length(ydays)) days)" filter(!isfile, collect(files))
            done(length(ydays))
            continue
        end
        for k in vars
            NCDataset(files[k]; maskingvalue=NaN) do ds
                v = ds[oafluxvar(ds)]
                datagrid(ds, v) == (lon, lat) || error("grid of $(files[k]) differs from the first file's")
                size(v, 3) == daysinyear(y) || error("$(files[k]) has $(size(v, 3)) days, not $(daysinyear(y))")
                for d in ydays
                    accumulate!(S[k], N[k], Float64.(v[:, :, dayofyear(d)]))
                end
            end
        end
        done(length(ydays))
    end
end

"""
    ceres_accumulate!(S, N, days, lon, lat; files, vars)

Add the CERES daily fields `vars` (keys of ceres_vars) on `days` to the sums S and counts N,
finding the days by the time coordinate of each of `files`. Warns of days in no file.
`done(n)` is called with the number of days found in each file, and of days not found.
"""
function ceres_accumulate!(S, N, days, lon, lat; files, vars, done=n -> nothing)
    want = Set(days)
    found = 0
    for file in files
        NCDataset(file; maskingvalue=NaN) do ds
            it = findall(in(want), Date.(ds["time"][:]))
            isempty(it) && return
            for k in vars
                v = ds[ceres_vars[k]]
                datagrid(ds, v) == (lon, lat) || error("grid of $(ceres_vars[k]) in $file differs from the first file's")
                for i in it
                    accumulate!(S[k], N[k], Float64.(v[:, :, i]))
                end
            end
            found += length(it)
            done(length(it))
        end
    end
    if found < length(days)
        @warn "$(length(days) - found) of $(length(days)) days are not in the CERES files"
        done(length(days) - found)
    end
end

"""
    oaflux_mean(days, lon, lat; progress=nothing, root, vars)
    ceres_mean(days, lon, lat; progress=nothing, files, vars)

Mean over `days` of the OAFlux or CERES daily fields `vars`, like `composite_mean`:
returns (mean, nobs), NamedTuples of lon x lat arrays, NaN where nobs is 0.
"""
function oaflux_mean(days, lon, lat; progress=nothing, root, vars)
    S, N = zerosums(vars, lon, lat)
    oaflux_accumulate!(S, N, days, lon, lat; root, vars, done=progress_done(progress))
    map((s, n) -> s ./ n, S, N), N
end
function ceres_mean(days, lon, lat; progress=nothing, files, vars)
    S, N = zerosums(vars, lon, lat)
    ceres_accumulate!(S, N, days, lon, lat; files, vars, done=progress_done(progress))
    map((s, n) -> s ./ n, S, N), N
end

"a function that puts n ones on the channel `progress` (nothing: does nothing)"
progress_done(progress) = n -> isnothing(progress) || foreach(_ -> put!(progress, 1), 1:n)

"""
    pmap_groups(f, groups)

Call f(days, progress) on the workers for each group of days (each worker takes whole groups),
printing a counter of days as f puts 1 on the channel `progress` for each. f returns (mean, nobs),
NamedTuples of lon x lat arrays; these are stacked to lon x lat x group and returned.
"""
function pmap_groups(f, groups)
    progress = RemoteChannel(() -> Channel{Int}(256))
    total = sum(length, groups)
    printer = @async for (n, _) in enumerate(progress)
        print("\r$n/$total days"); flush(stdout)
    end
    r = pmap(days -> f(days, progress), groups)
    close(progress); wait(printer); println()
    ks = keys(first(r)[1])
    map(k -> stack(ri[1][k] for ri in r), NamedTuple{ks}(ks)),
    map(k -> stack(ri[2][k] for ri in r), NamedTuple{ks}(ks))
end

"""
    oaflux_composite(groups; root, vars=keys(oaflux_files))
    ceres_composite(groups; root, vars=keys(ceres_vars))

Composite OAFlux or CERES SYN1deg daily fields `vars` over each group of days, like `composite`
for ERA5. Returns (comp, nobs, lon, lat, units); comp and nobs are NamedTuples of
lon x lat x group arrays, units a NamedTuple of each variable's units attribute.
"""
function oaflux_composite(groups; root, vars=keys(oaflux_files))
    y0 = year(first(Iterators.flatten(groups)))
    lon, lat = NCDataset(ds -> datagrid(ds, ds[oafluxvar(ds)]), oafluxfile(root, first(vars), y0))
    units = map(NamedTuple{vars}(vars)) do k
        NCDataset(ds -> get(ds[oafluxvar(ds)].attrib, "units", ""), oafluxfile(root, k, y0))
    end
    comp, nobs = pmap_groups((days, progress) -> oaflux_mean(days, lon, lat; progress, root, vars), groups)
    comp, nobs, lon, lat, units
end
function ceres_composite(groups; root, vars=keys(ceres_vars))
    files = ceresfiles(root)
    isempty(files) && error("no CERES_SYN1deg-Day files in $root")
    lon, lat, units = NCDataset(first(files)) do ds
        (datagrid(ds, ds[ceres_vars[first(vars)]])...,
         map(k -> get(ds[ceres_vars[k]].attrib, "units", ""), NamedTuple{vars}(vars)))
    end
    comp, nobs = pmap_groups((days, progress) -> ceres_mean(days, lon, lat; progress, files, vars), groups)
    comp, nobs, lon, lat, units
end

"""
    composite_data(data, groups; sst=:oisst, era5root, oafluxroot, ceresroot, oisstroot=nothing)

Composite each group of days of the data set `data`, each in its own run and file:
`:era5`: ERA5 fields and Ekman terms (`composite`), from era5dir(era5root), with the SST
product `sst`: `:oisst` (default), daily OISST v2.1 from oisstroot on the ERA5 grid, or
`:e5sst`, ERA5's own SST (Andrea's copy stops in 2020);
`:oaflux`: all OAFlux fields, with ISCCP radiation and qnet (through 2009);
`:oafluxturb`: the OAFlux turbulent fields only (through 2022);
`:ceres`: CERES SYN1deg adjusted all-sky surface and TOA radiative fluxes (from March 2000).
OAFlux and CERES composites of the same days can be combined after compositing (same grid).
Returns (comp, nobs, lon, lat, units, prefix, sstname, title): prefix names the data set
and sstname the SST product ("oisst", "e5sst"; "" for data sets that use neither) in the
output file name, <prefix>_<grouping>_<sstname>_<period>.nc, and title in its title.
"""
function composite_data(data, groups; sst=:oisst, era5root, oafluxroot, ceresroot, oisstroot=nothing)
    if data == :era5 && sst == :oisst
        isnothing(oisstroot) && error("SST :oisst needs oisstroot, the OISST directory")
        (composite(groups; dir=era5dir(era5root), oisst=oisstroot)..., units,
         "ekman", "oisst", "ERA5 surface fields with OISST v2.1 SST, and Ekman terms")
    elseif data == :era5 && sst == :e5sst
        (composite(groups; dir=era5dir(era5root))..., units,
         "ekman", "e5sst", "ERA5 surface fields with ERA5 SST, and Ekman terms")
    elseif data == :era5
        error("sst must be :oisst or :e5sst, not :$sst")
    elseif data == :oaflux
        (oaflux_composite(groups; root=oafluxroot)..., "oaflux", "", "OAFlux daily surface fluxes and meteorology")
    elseif data == :oafluxturb
        (oaflux_composite(groups; root=oafluxroot, vars=oaflux_turbulent)...,
         "oafluxturb", "", "OAFlux daily turbulent fluxes and meteorology")
    elseif data == :ceres
        (ceres_composite(groups; root=ceresroot)...,
         "ceres", "", "CERES SYN1deg daily adjusted all-sky surface and TOA radiative fluxes")
    else
        error("data must be :era5, :oaflux, :oafluxturb or :ceres, not :$data")
    end
end

"""
    script_args(args=ARGS)

Options of a compositing script from its command-line words, in any order: a data set
(era5, oaflux, oafluxturb or ceres; default era5), an SST product for era5 (oisst or
e5sst; default oisst), a first and last year (two integers), and peryear (one file per
year). Returns (; data, sst, years, peryear); years is nothing unless given.
"""
function script_args(args=ARGS)
    data, sst, yrs, peryear = :era5, :oisst, Int[], false
    for a in args
        if a in ("era5", "oaflux", "oafluxturb", "ceres")
            data = Symbol(a)
        elseif a in ("oisst", "e5sst")
            sst = Symbol(a)
        elseif a == "peryear"
            peryear = true
        elseif !isempty(a) && all(isdigit, a)
            push!(yrs, parse(Int, a))
        else
            error("unknown argument \"$a\": expected era5, oaflux, oafluxturb, ceres, oisst, e5sst, peryear or a year")
        end
    end
    length(yrs) in (0, 2) || error("give a first and a last year, or no years")
    (; data, sst, years=isempty(yrs) ? nothing : yrs[1]:yrs[2], peryear)
end

# units of the saved fields; M in kg/m/s, so M⋅∇s and s∇⋅M are in kg/m^2/s times the units of s
const units = (taux="N m-2", tauy="N m-2", sst="K", t2="K", d2="K", q="kg kg-1",
         adv_sst="K kg m-2 s-1", adv_t2="K kg m-2 s-1", adv_q="kg m-2 s-1",
         sdiv_sst="K kg m-2 s-1", sdiv_t2="K kg m-2 s-1", sdiv_q="kg m-2 s-1",
         sst_taux="K N m-2", t2_taux="K N m-2", q_taux="N m-2",
         sst_tauy="K N m-2", t2_tauy="K N m-2", q_tauy="N m-2")

"""
    save_composite(file, comp, nobs, lon, lat, dim, dimvals; dimattrib=[], attrib=[], fieldunits=units)

Write composites to NetCDF `file` (overwritten): each field k of `comp` as Float32
(NaN fill), and its count of daily values as Int16 `nobs_k`, on (longitude, latitude, dim),
with coordinate `dim` = `dimvals`. `dimattrib` and `attrib` are the attributes of the
coordinate and of the file; the history is added to the latter. `fieldunits` gives each
field's units (default: the ERA5 fields and Ekman terms; OAFlux units come from `oaflux_composite`).
"""
function save_composite(file, comp, nobs, lon, lat, dim, dimvals; dimattrib=[], attrib=[], fieldunits=units)
    NCDataset(file, "c", attrib=[attrib...,
                                 "history" => "$(now()) $(basename(PROGRAM_FILE))"]) do ds
        defVar(ds, "longitude", lon, ("longitude",), attrib=["units" => "degrees_east"])
        defVar(ds, "latitude",  lat, ("latitude",),  attrib=["units" => "degrees_north"])
        defVar(ds, dim, collect(Int32, dimvals), (dim,), attrib=dimattrib)
        dims = ("longitude", "latitude", dim)
        for k in keys(comp)
            defVar(ds, String(k), Float32.(comp[k]), dims; fillvalue=NaN32, deflatelevel=4,
                   attrib=["units" => get(fieldunits, k, "")])
            defVar(ds, "nobs_$k", nobs[k], dims; deflatelevel=4,
                   attrib=["long_name" => "number of daily values in the mean of $k"])
        end
    end
end

# compile the analysis code on a tiny synthetic grid when the package is precompiled
@setup_workload begin
    lon = collect(0.0:30.0:330.0)
    lat = collect(90.0:-30.0:-90.0)
    @compile_workload begin
        mktempdir() do root
            dir = era5dir(root)
            d = Date(2000, 1, 1)
            t = collect(DateTime(d):Hour(1):DateTime(d) + Hour(23))
            for k in keys(era5code)
                f = era5file(dir[k], era5code[k], d)
                mkpath(dirname(f))
                tn = k == :sst ? "valid_time" : "time"
                NCDataset(f, "c") do ds
                    defVar(ds, "longitude", lon, ("longitude",))
                    defVar(ds, "latitude", lat, ("latitude",))
                    defVar(ds, tn, t, (tn,), attrib=["units" => "hours since 1900-01-01"])
                    defVar(ds, uppercase(String(k)), fill(k == :d2 ? 290f0 : 0.1f0, length(lon), length(lat), 24),
                           ("longitude", "latitude", tn); fillvalue=NaN32)
                end
            end
            composite_mean([d], lon, lat; dir=dir)
        end
    end
end

end # module
