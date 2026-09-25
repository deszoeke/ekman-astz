"""
    EkmanCommon

Shared code for ekman_mse.jl and bsiso_phase.jl: ERA5 file access, thermodynamics,
Ekman transport, finite differences on the sphere, and compositing daily fields over
groups of days on worker processes. Load it on every process:
    @everywhere using EkmanCommon
"""
module EkmanCommon

using Dates
using Distributed
using NCDatasets
using PrecompileTools

export era5dir, era5code, era5file, daily_terms, composite_mean, composite, save_composite,
       fieldkeys, nlkeys, keys_all, scalars, Ro_rad, c_pa, c_po, Lv, esat_buck, qsat_dew,
       coriolis, ekman_transport, ekman_transport_xy, gradient, divergence, advection,
       dailymean, accumulate!

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

"the data variable in an RDA file (the only variable with 3 dimensions)"
datavar(ds) = first(k for k in keys(ds) if ndims(ds[k]) == 3)

"the time coordinate variable name in a file (sst uses valid_time, others use time)"
timevar(ds) = haskey(ds, "valid_time") ? "valid_time" : "time"

"sorted indices as a UnitRange when contiguous (fast NetCDF hyperslab read), else left as-is"
asrange(idx) = !isempty(idx) && idx == first(idx):last(idx) ? (first(idx):last(idx)) : idx

"""
    withdatasets(f, files)

Open each file in the NamedTuple `files` (fill values read as NaN), call `f` on the
NamedTuple of datasets, and close them all afterwards, like `NCDataset(file) do ds`
for several files. Files opened before a failure are closed too.
"""
function withdatasets(f, files::NamedTuple)
    opened = NCDataset[]
    try
        for file in files
            push!(opened, NCDataset(file; maskingvalue=NaN))
        end
        f(NamedTuple{keys(files)}(Tuple(opened)))
    finally
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
    ekman_transport(taux, tauy, lat; Ro_rad=nothing)

Ekman mass transport (Mx, My) = (tauy, -taux)/f (kg/m/s) from wind stress (N/m^2)
at latitude lat (degrees). If Ro_rad (m) is given, curtail the effective stress
near the equator by 1 - exp(-y^2/Ro_rad^2), with y the distance (m) from the
equator, to account for the pressure gradient balancing the wind stress and keep
the transport finite (-> 0 at the equator). Without Ro_rad the transport is
infinite (NaN for zero stress) at the equator.
"""
function ekman_transport(taux, tauy, lat; Ro_rad=nothing)
    f = coriolis(lat)
    if isnothing(Ro_rad)
        r = 1.0
    else
        y = a*deg2rad(lat)
        f == 0 && return (0.0*taux, 0.0*tauy) # limit of r/f ~ y/(beta*Ro_rad^2) -> 0; keeps NaN stress NaN
        r = -expm1(-(y/Ro_rad)^2) # 1 - exp(-y^2/Ro_rad^2)
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


const Ro_rad = 250e3 # m, equatorial Rossby radius at which pressure balances tau

# scalars s and the sign of the Ekman transport of their medium:
# ocean (sst) is carried by M, atmosphere (t2, q) by the equal and opposite -M
const scalars = (sst=1, t2=-1, q=-1)
# nonlinear terms, computed from daily means, then composited:
#   adv_s  = -(±M)⋅∇s   advection by Ekman transport of the medium
#   sdiv_s = s ∇⋅(±M)   s times Ekman transport divergence of the medium
# so the flux-form tendency is -∇⋅(±M s) = adv_s - sdiv_s
const nlkeys = (Symbol.(:adv_,  keys(scalars))..., Symbol.(:sdiv_, keys(scalars))...)
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
    daily_terms(d, lon, lat; dir=era5dir(), code=era5code, Ro_rad=Ro_rad, scalars=scalars)

Daily means (lon x lat) on date d of the ERA5 fields, q from hourly dewpoint, and the
nonlinear Ekman terms computed from them, as a NamedTuple with keys `keys_all`.
Returns nothing, with a warning, if the files or hours for d are missing.
"""
function daily_terms(d, lon, lat; dir=era5dir(), code=era5code, Ro_rad=Ro_rad, scalars=scalars)
    files = map((dk, ck) -> era5file(dk, ck, d), dir, code)
    if !all(isfile, files)
        @warn "missing ERA5 files for $d" filter(!isfile, collect(files))
        return nothing
    end
    x = withdatasets(files) do ds
        # every source must be on the stress grid (sst and t2/d2 come from a different archive)
        for (k, f) in pairs(ds)
            f["latitude"][:] ≈ lat && f["longitude"][:] ≈ lon ||
                error("$k grid in $(files[k]) differs from the stress grid")
        end
        it = map(f -> asrange(findall(==(d), Date.(f[timevar(f)][:]))), ds) # sst: valid_time
        if any(isempty, it)
            @warn "missing hours on $d"
            return nothing
        end
        # average hourly data to daily; q from hourly dewpoint, then daily mean
        hours(k) = ds[k][datavar(ds[k])][:, :, it[k]]
        hd2 = hours(:d2)
        (taux=dailymean(hours(:taux)), tauy=dailymean(hours(:tauy)), sst=dailymean(hours(:sst)),
         t2=dailymean(hours(:t2)), d2=dailymean(hd2), q=dailymean(hd2, qsat_dew))
    end
    isnothing(x) && return nothing
    # daily Ekman transport and its divergence
    Mx, My = ekman_transport_xy(x.taux, x.tauy, lat; Ro_rad=Ro_rad)
    divM = divergence(Mx, My, lon, lat)
    adv  = Tuple(sgn .* advection(Mx, My, x[s], lon, lat) for (s, sgn) in pairs(scalars))
    sdiv = Tuple(sgn .* divM .* x[s] for (s, sgn) in pairs(scalars))
    merge(x, NamedTuple{Symbol.(:adv_,  keys(scalars))}(adv),
             NamedTuple{Symbol.(:sdiv_, keys(scalars))}(sdiv))
end

"""
    composite_mean(days, lon, lat; progress=nothing, kw...)

Mean over `days` of the finite daily values (see `daily_terms`, which takes `kw`).
Returns (mean, nobs), NamedTuples of lon x lat arrays; the mean is NaN where nobs is 0.
If `progress` is a channel, put 1 on it as each day finishes.
"""
function composite_mean(days, lon, lat; progress=nothing, kw...)
    S = NamedTuple{keys_all}(Tuple(zeros(length(lon), length(lat)) for k in keys_all))
    N = NamedTuple{keys_all}(Tuple(zeros(Int16, length(lon), length(lat)) for k in keys_all))
    for d in days
        x = daily_terms(d, lon, lat; kw...)
        isnothing(x) || foreach(accumulate!, S, N, x)
        isnothing(progress) || put!(progress, 1)
    end
    map((s, n) -> s ./ n, S, N), N
end

"""
    composite(groups; dir=era5dir(), code=era5code, Ro_rad=Ro_rad, scalars=scalars)

Composite means (lon x lat x group) of the ERA5 fields and the nonlinear Ekman terms,
one for each vector of dates in `groups`. Means are taken over the finite daily values
at each grid point, so a missing hour or a land point drops out of that point's mean
instead of making it NaN. The groups are distributed over the worker processes, each
worker summing whole groups, so no two processes share an accumulator. Progress is a
counter of finished days, redrawn in place. Returns (comp, nobs, lon, lat).
"""
function composite(groups; dir=era5dir(), code=era5code, Ro_rad=Ro_rad, scalars=scalars)
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
                                    dir=dir, code=code, Ro_rad=Ro_rad, scalars=scalars),
             groups)
    close(progress); wait(printer); println()
    # stack the groups (lon x lat x group)
    comp = map(k -> stack(ri[1][k] for ri in r), NamedTuple{keys_all}(keys_all))
    nobs = map(k -> stack(ri[2][k] for ri in r), NamedTuple{keys_all}(keys_all))
    comp, nobs, lon, lat
end

# units of the saved fields; M in kg/m/s, so M⋅∇s and s∇⋅M are in kg/m^2/s times the units of s
const units = (taux="N m-2", tauy="N m-2", sst="K", t2="K", d2="K", q="kg kg-1",
         adv_sst="K kg m-2 s-1", adv_t2="K kg m-2 s-1", adv_q="kg m-2 s-1",
         sdiv_sst="K kg m-2 s-1", sdiv_t2="K kg m-2 s-1", sdiv_q="kg m-2 s-1")

"""
    save_composite(file, comp, nobs, lon, lat, dim, dimvals; dimattrib=[], attrib=[])

Write composites to NetCDF `file` (overwritten): each field k of `comp` as Float32
(NaN fill), and its count of daily values as Int16 `nobs_k`, on (longitude, latitude, dim),
with coordinate `dim` = `dimvals`. `dimattrib` and `attrib` are the attributes of the
coordinate and of the file; Ro_rad and the history are added to the latter.
"""
function save_composite(file, comp, nobs, lon, lat, dim, dimvals; dimattrib=[], attrib=[])
    NCDataset(file, "c", attrib=[attrib..., "Ro_rad_m" => Ro_rad,
                                 "history" => "$(now()) $(basename(PROGRAM_FILE))"]) do ds
        defVar(ds, "longitude", lon, ("longitude",), attrib=["units" => "degrees_east"])
        defVar(ds, "latitude",  lat, ("latitude",),  attrib=["units" => "degrees_north"])
        defVar(ds, dim, collect(Int32, dimvals), (dim,), attrib=dimattrib)
        dims = ("longitude", "latitude", dim)
        for k in keys(comp)
            defVar(ds, String(k), Float32.(comp[k]), dims; fillvalue=NaN32, deflatelevel=4,
                   attrib=["units" => get(units, k, "")])
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
