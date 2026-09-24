using DelimitedFiles
using PythonPlot
using Dates
using NCDatasets
using Statistics

"convert PC1, PC2 components of BSISO into canonical phase index (1-8)"
phase(b1, b2) = mod(floor(8/(2*pi) * atan(b1, -b2)), 8) + 1


ds = readdlm("data/BSISO.INDEX.NORM.LY.data", header=true)
B = ds[1]
header = ds[2]
dt = @. Date(B[:,1], 1, 1) + Day(B[:,2] - 1)

plot( phase.(B[:,3], B[:,4]) ) # increases
plot( phase.(B[:,5], B[:,6]) )

# save BSISO phases
# YEAR DOY MONTH DATE BSISO1_phase BSISO2_phase
ph1 = Int.(phase.(B[:,3], B[:,4]))
ph2 = Int.(phase.(B[:,5], B[:,6]))
writedlm("BSISO_phase.txt",
         [["YEAR" "DOY" "MONTH" "DATE" "BSISO1_phase" "BSISO2_phase"];
          year.(dt) dayofyear.(dt) month.(dt) day.(dt) ph1 ph2])

# preallocate stress, sst, t2, d2 0.25x0.25 degee arrays for 8 BSISO phases
# preallocate counter array(s)
# loop over days for 2012-2026 
#   load ERA5 wind stress, sst, t2, and d2 each day
#   open and read file(s); average hourly data to daily
#   look up BSISO phase index by day
#   increment stress, sst, t2, d2 arrays for appropriate phase
#   increment counter for phase
# get composite mean by dividing each array by the count for each phase


# look up BSISO1 phase by day
phaselookup = Dict(zip(dt, ph1))

# ERA5 hourly surface analysis from NCAR RDA d633000, one file per variable per month:
#   <dir>/201201/e5.oper.an.sfc.128_229_iews.ll025sc.2012010100_2012013123.nc
rdadir  = "d633000/e5.oper.an.sfc" # adjust to local copy
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

"average hourly data to daily (lon x lat) for time indices it, applying f to each hour; missing -> NaN"
function dailymean(v, it, f=identity)
    x = f.(Float64.(coalesce.(v[:,:,it], NaN)))
    dropdims(mean(x, dims=3), dims=3)
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
        f == 0 && return (zero(taux/1), zero(tauy/1)) # limit of r/f ~ y/(beta*Ro_rad^2) -> 0
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
    for k in CartesianIndices(size(φ)[3:end]), j in 1:ny, i in 1:nx
        # zonal
        if periodic || 1 < i < nx
            ip, im = mod1(i+1, nx), mod1(i-1, nx)
            dλ = deg2rad(mod(lon[ip] - lon[im], 360))
            c = cosd(lat[j])
            if c > 1e-10
                dφdx[i,j,k] = (φ[ip,j,k] - φ[im,j,k]) / (a*c*dλ)
            end
        end
        # meridional
        if 1 < j < ny
            dφdy[i,j,k] = (φ[i,j+1,k] - φ[i,j-1,k]) / (a*deg2rad(lat[j+1] - lat[j-1]))
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

# preallocate stress, sst, t2, d2 0.25x0.25 degree arrays for 8 BSISO phases
months = Date(2012,1,1):Month(1):Date(2026,12,1)
lon, lat = NCDataset(era5file(era5dir.taux, era5code.taux, first(months))) do ds
    ds["longitude"][:], ds["latitude"][:]
end
nx, ny = length(lon), length(lat)
fieldkeys = (keys(era5code)..., :q) # taux, tauy, sst, t2, d2, q
comp = Dict(k => zeros(nx, ny, 8) for k in (fieldkeys..., nlkeys...))
# preallocate counter array(s)
nday = zeros(Int, 8)

# loop over days for 2012-2026, opening each month's files once
# for m in months
for m in months[1:3] # short test run
    files = Dict(k => era5file(era5dir[k], era5code[k], m) for k in keys(era5code))
    if !all(isfile, values(files))
        @warn "missing ERA5 files for $(Dates.format(m, "yyyy-mm"))"
        continue
    end
    # open and read files
    ds = Dict(k => NCDataset(f) for (k, f) in files)
    vars = Dict(k => d[datavar(d)] for (k, d) in ds)
    daystamp = Dict(k => Date.(d["time"][:]) for (k, d) in ds) # per file: sources differ
    for d in m:Day(1):lastdayofmonth(m)
        haskey(phaselookup, d) || continue # no BSISO index for this day
        # look up BSISO phase index by day
        p = phaselookup[d]
        it = Dict(k => findall(==(d), daystamp[k]) for k in keys(era5code))
        if any(isempty, values(it))
            @warn "missing hours on $d"
            continue
        end
        # average hourly data to daily
        x = Dict(k => dailymean(vars[k], it[k]) for k in keys(era5code))
        # specific humidity from hourly dewpoint, then daily mean
        x[:q] = dailymean(vars[:d2], it[:d2], qsat_dew)
        # daily Ekman transport and its divergence
        Mx, My = ekman_transport_xy(x[:taux], x[:tauy], lat; Ro_rad=Ro_rad)
        divM = divergence(Mx, My, lon, lat)
        # increment stress, sst, t2, d2, q arrays for appropriate phase
        for k in fieldkeys
            comp[k][:,:,p] .+= x[k] # sst NaN over land
        end
        # increment nonlinear terms from daily means
        for (s, sgn) in pairs(scalars)
            comp[Symbol(:adv_, s)][:,:,p]  .+= sgn .* advection(Mx, My, x[s], lon, lat)
            comp[Symbol(:sdiv_, s)][:,:,p] .+= sgn .* divM .* x[s]
        end
        # increment counter for phase
        nday[p] += 1
    end
    foreach(close, values(ds))
end

# get composite mean by dividing each array by the count for each phase
for A in values(comp), p in 1:8
    A[:,:,p] ./= nday[p]
end
taux, tauy, sst, t2, d2, q = (comp[k] for k in fieldkeys)
adv_sst, adv_t2, adv_q, sdiv_sst, sdiv_t2, sdiv_q = (comp[k] for k in nlkeys)

# specific moist static energy of surface air, m_a = c_pa*Ta + Lv*q (J/kg; T in K, z=2 m term neglected)
m_a = @. c_pa*t2 + Lv*q
# specific enthalpy of ocean surface water, h_o = c_po*SST (J/kg; SST in K)
h_o = @. c_po*sst

# Ekman mass transport of composite mean stress (lon x lat x phase), curtailed within ~Ro_rad of the equator
Mx, My = ekman_transport_xy(taux, tauy, lat; Ro_rad=Ro_rad)
divM = divergence(Mx, My, lon, lat)

# gradients of ocean enthalpy and surface air moist static energy (J/kg/m)
dhodx, dhody = gradient(h_o, lon, lat)
dmadx, dmady = gradient(m_a, lon, lat)

# advection of composite means by composite Ekman mass transport (W/m^2): -M⋅∇
# ocean Ekman transport is M; atmospheric Ekman transport is equal and opposite, -M
adv_ho = @. -(Mx*dhodx + My*dhody)
adv_ma = @.  (Mx*dmadx + My*dmady)
