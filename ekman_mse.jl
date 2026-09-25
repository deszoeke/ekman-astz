using PythonPlot
using Dates
using NCDatasets
using Statistics

# Pentad (5-day) climatologies of ERA5 surface fields and air-sea Ekman terms for April-July.
# loop over years, and the months containing the selected pentads, opening each month's files once
#   average hourly data to daily
#   compute daily Ekman transport and the nonlinear terms
#   add finite values to the pentad sums, counting them at each grid point
# climatological mean = sum / count at each grid point (NaN where nobs is 0, e.g. sst over land)

# ERA5 hourly surface analysis from NCAR RDA d633000, one file per variable per month:
# adjust to local directory structure
#   <dir>/201201/e5.oper.an.sfc.128_229_iews.ll025sc.2012010100_2012013123.nc
# rdadir  = "d633000/e5.oper.an.sfc" # not used
# ceoasdir = "/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC"
ceoasdir = "./data"
# ECMWF parameter codes: iews, inss (N/m^2), sstk, 2t, 2d (K)
era5code = (taux="128_229_iews", tauy="128_230_inss", sst="128_034_sstk",
            t2="128_167_2t", d2="128_168_2d")
# directory holding the YYYYMM subdirectories for each variable
era5dir = (taux=joinpath(ceoasdir, "stress"), tauy=joinpath(ceoasdir, "stress"), 
           sst=joinpath(ceoasdir, "sst"),
           t2=joinpath(ceoasdir, "t2"), d2=joinpath(ceoasdir, "d2"))
# SST files are grouped by month directly in the directory data/sst/

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


Ro_rad = 250e3 # m, equatorial Rossby radius at which pressure balances tau

# scalars s and the sign of the Ekman transport of their medium:
# ocean (sst) is carried by M, atmosphere (t2, q) by the equal and opposite -M
scalars = (sst=1, t2=-1, q=-1)
# nonlinear terms, computed from daily means, then composited:
#   adv_s  = -(±M)⋅∇s   advection by Ekman transport of the medium
#   sdiv_s = s ∇⋅(±M)   s times Ekman transport divergence of the medium
# so the flux-form tendency is -∇⋅(±M s) = adv_s - sdiv_s
nlkeys = (Symbol.(:adv_,  keys(scalars))..., Symbol.(:sdiv_, keys(scalars))...)
fieldkeys = (keys(era5code)..., :q) # taux, tauy, sst, t2, d2, q

"pentad (1-73) of date d, on the 365-day calendar (Feb 29 joins pentad 12)"
pentad(d) = (dayofyear(d) - (isleapyear(d) && month(d) > 2) - 1) ÷ 5 + 1

"days of year y in the pentads"
pentaddays(y, pentads) = filter(d -> pentad(d) in pentads, Date(y,1,1):Day(1):Date(y,12,31))

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
    pentad_climatology(years, pentads; dir=era5dir, code=era5code, Ro_rad=Ro_rad, scalars=scalars)

Climatological mean (lon x lat x pentad) over `years` of the ERA5 fields and the
nonlinear Ekman terms for each pentad in the range `pentads`. Means are taken over
the finite daily values at each grid point, so a missing hour or a land point
drops out of that point's mean instead of making it NaN. Returns (comp, nobs, lon, lat).
"""
function pentad_climatology(years, pentads::UnitRange; dir=era5dir, code=era5code,
                            Ro_rad=Ro_rad, scalars=scalars)
    lon, lat = NCDataset(era5file(dir.taux, code.taux, first(pentaddays(first(years), pentads)))) do ds
        Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:])
    end
    nx, ny, np = length(lon), length(lat), length(pentads)
    keys_all = (fieldkeys..., nlkeys...)
    comp  = Dict(k => zeros(nx, ny, np) for k in keys_all)
    nobs  = Dict(k => zeros(Int16, nx, ny, np) for k in keys_all)

    for y in years
        days = pentaddays(y, pentads)
        for m in unique(firstdayofmonth.(days))
            files = map((dk, ck) -> era5file(dk, ck, m), dir, code)
            if !all(isfile, files)
                @warn "missing ERA5 files for $(Dates.format(m, "yyyy-mm"))" filter(!isfile, collect(files))
                continue
            end
            withdatasets(files) do ds
                # every source must be on the stress grid (sst and t2/d2 come from a different archive)
                for (k, d) in pairs(ds)
                    d["latitude"][:] ≈ lat && d["longitude"][:] ≈ lon ||
                        error("$k grid in $(files[k]) differs from the stress grid")
                end
                vars = map(d -> d[datavar(d)], ds)
                daystamp = map(d -> Date.(d[timevar(d)][:]), ds) # sst: valid_time
                for d in filter(d -> firstdayofmonth(d) == m, days)
                    ip = pentad(d) - first(pentads) + 1
                    it = map(t -> asrange(findall(==(d), t)), daystamp)
                    if any(isempty, it)
                        @warn "missing hours on $d"
                        continue
                    end
                    # average hourly data to daily, reading each variable's hours once;
                    # q from hourly dewpoint, then daily mean
                    x = Dict{Symbol,Matrix{Float64}}()
                    for k in keys(vars)
                        h = vars[k][:, :, it[k]]
                        x[k] = dailymean(h)
                        k == :d2 && (x[:q] = dailymean(h, qsat_dew))
                    end
                    # daily Ekman transport and its divergence
                    Mx, My = ekman_transport_xy(x[:taux], x[:tauy], lat; Ro_rad=Ro_rad)
                    divM = divergence(Mx, My, lon, lat)
                    for (s, sgn) in pairs(scalars)
                        x[Symbol(:adv_, s)]  = sgn .* advection(Mx, My, x[s], lon, lat)
                        x[Symbol(:sdiv_, s)] = sgn .* divM .* x[s]
                    end
                    for k in keys_all
                        accumulate!(view(comp[k], :, :, ip), view(nobs[k], :, :, ip), x[k])
                    end
                end
            end
        end
    end
    # climatological mean; NaN where there are no finite values
    for k in keys_all
        comp[k] ./= nobs[k]
    end
    comp, nobs, lon, lat
end

# units of the saved fields; M in kg/m/s, so M⋅∇s and s∇⋅M are in kg/m^2/s times the units of s
units = (taux="N m-2", tauy="N m-2", sst="K", t2="K", d2="K", q="kg kg-1",
         adv_sst="K kg m-2 s-1", adv_t2="K kg m-2 s-1", adv_q="kg m-2 s-1",
         sdiv_sst="K kg m-2 s-1", sdiv_t2="K kg m-2 s-1", sdiv_q="kg m-2 s-1")

"""
    save_climatology(file, comp, nobs, lon, lat, years, pentads)

Write the pentad climatology to NetCDF `file` (overwritten): each field k of `comp` as
Float32 (NaN fill), and its count of daily values as Int16 `nobs_k`, on (longitude, latitude, pentad).
"""
function save_climatology(file, comp, nobs, lon, lat, years, pentads)
    NCDataset(file, "c", attrib=["title" => "ERA5 pentad climatology of surface fields and Ekman terms",
                                 "years" => "$(first(years))-$(last(years))",
                                 "Ro_rad_m" => Ro_rad,
                                 "history" => "$(now()) ekman_mse.jl"]) do ds
        defVar(ds, "longitude", lon, ("longitude",), attrib=["units" => "degrees_east"])
        defVar(ds, "latitude",  lat, ("latitude",),  attrib=["units" => "degrees_north"])
        defVar(ds, "pentad", collect(Int32, pentads), ("pentad",),
               attrib=["long_name" => "pentad of year (1-73), 365-day calendar, Feb 29 in pentad 12"])
        dims = ("longitude", "latitude", "pentad")
        for k in keys(comp)
            defVar(ds, String(k), Float32.(comp[k]), dims; fillvalue=NaN32, deflatelevel=4,
                   attrib=["units" => get(units, k, "")])
            defVar(ds, "nobs_$k", nobs[k], dims; deflatelevel=4,
                   attrib=["long_name" => "number of daily values in the mean of $k"])
        end
    end
end

# April-July: pentads 19 (Apr 1-5) through 43 (Jul 30-Aug 3)
years   = 2012:2026
pentads = 19:43
comp, nobs, lon, lat = pentad_climatology(years, pentads)
# comp, nobs, lon, lat = pentad_climatology(2012:2012, 19:20) # short test run
save_climatology("ekman_pentad_clim_$(first(years))-$(last(years)).nc", comp, nobs, lon, lat, years, pentads)

taux, tauy, sst, t2, d2, q = (comp[k] for k in fieldkeys)
adv_sst, adv_t2, adv_q, sdiv_sst, sdiv_t2, sdiv_q = (comp[k] for k in nlkeys)

# compute full advection and scalar divergence of 
adv_h_o = c_po * adv_sst
adv_h_a = c_pa * adv_t2
adv_m_a = adv_h_a + L * adv_q
sdiv_h_o = c_po * sdiv_sst
sdiv_h_a = c_pa * sdiv_t2
sdiv_m_a = adv_h_a * sdiv_t2 + L*sdiv_q
# NOTE the signs of these are defined as if 
# they are on opposite sides of the equation!
# adv  = -M⋅∇
# sdiv = +∇⋅M

# Compute advection of the mean by the mean Ekman mass transport

# mean specific moist static energy of surface air, m_a = c_pa*Ta + Lv*q (J/kg; T in K, z=2 m term neglected)
m_a = @. c_pa*t2 + Lv*q
# specific enthalpy of ocean surface water, h_o = c_po*SST (J/kg; SST in K)
h_o = @. c_po*sst

# mean Ekman mass transport of climatological mean stress (lon x lat x pentad), curtailed within ~Ro_rad of the equator
Mx, My = ekman_transport_xy(taux, tauy, lat; Ro_rad=Ro_rad)
divM = divergence(Mx, My, lon, lat) # mean divergence of Ekman mass transport
# gradients of ocean enthalpy and surface air moist static energy (J/kg/m)
dhodx, dhody = gradient(h_o, lon, lat)
dmadx, dmady = gradient(m_a, lon, lat)
# advection of climatological means by climatological Ekman mass transport (W/m^2): -M⋅∇
# ocean Ekman transport is M; atmospheric Ekman transport is equal and opposite, -M
adv_ho = @. -(Mx*dhodx + My*dhody)
adv_ma = @.  (Mx*dmadx + My*dmady)
