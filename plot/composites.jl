# Composites: reading the composite NetCDF files written by save_composite and computing the
# fields to plot. No plotting here (that is maps.jl). Sections, from general to specific:
#   General: any composite file   read_group, group_index, read_mean, group_days, anom, box_mean, ekman_fields
#   Pentad climatology            (ekman_mse.jl output) pentads for a period of dates
#   BSISO phases                  (bsiso_phase.jl output) phases
# The pentad and BSISO sections use only the general section, not each other.
#
# Fields: the Ekman mass transport (Mx, My), surface-air moist enthalpy h_a = c_pa T2 + Lv q
# (moist static energy at the surface) and ocean-surface enthalpy h_o = c_po SST; the gradients of h_o and h_a and
# their Ekman advection -M⋅∇s; and s ∇⋅M and the flux-form tendency -∇⋅(M s) = -M⋅∇s - s ∇⋅M.
# The ocean's transport is M_o = M = (tauy, -taux)/f; the atmosphere's is M_a = -M_o.
# M, h_a, h_o and the gradients are linear in the composited fields, so computing them from
# the composite means is exact. The nonlinear terms are not: -M⋅∇s is the saved adv_*
# composite of daily products, and the flux M s = (s τy, -s τx)/f comes from the saved
# s_taux, s_tauy composites (EkmanCommon.daily_terms), so they include the transient (eddy)
# transport. The saved adv_* signs already include M_a = -M_o.
# -M⋅∇s > 0 is advection from higher s, a positive tendency (as adv = -M⋅∇ in ekman_mse.jl).
#
# The composites are of the uncurtailed M. The equatorial curtailment r = curtail(lat, Ro_rad)
# depends only on latitude, so it is applied here: M -> r M, -M⋅∇s -> -r M⋅∇s, and
# -∇⋅(M s) -> -∇⋅(r M s); then s ∇⋅(r M) = -r M⋅∇s + ∇⋅(r M s).
#
# Notebooks load it with Revise, so edits take effect without restarting the kernel:
#   using Revise; includet("composites.jl")
using NCDatasets, Dates, Statistics
using EkmanCommon: c_pa, c_po, rho_o, Lv, qsat_dew, curtail, ekman_transport_xy, gradient, divergence

## General: any composite file (group dimension: pentad or phase)

const composite_vars = ("taux", "tauy", "sst", "t2", "q", "adv_sst", "adv_t2", "adv_q",
                        "sst_taux", "sst_tauy", "t2_taux", "t2_tauy", "q_taux", "q_tauy")

"""
    read_group(file, k; vars=composite_vars, lonlim=nothing, latlim=nothing)

Composite fields at index k along the group (pentad/phase) dimension of a composite NetCDF
file, as a NamedTuple: lon, lat, each of `vars` (lon x lat, Float64, NaN where missing),
dim and dimval (the group dimension's name and value at k), file and k. With a range or
vector k, e.g. 1:8, each variable is lon x lat x group, indexed d.taux[:, :, q]. `lonlim`,
`latlim` (degrees, e.g. (40, 110), (-5, 30)) read only that region from the file.
"""
function read_group(file, k; vars=composite_vars, lonlim=nothing, latlim=nothing)
    NCDataset(file; maskingvalue=NaN) do ds
        dim = group_dim(ds)
        inlim(x, lim) = isnothing(lim) ? (1:length(x)) : (r = findall(v -> lim[1] <= v <= lim[2], x); first(r):last(r))
        i = inlim(ds["longitude"][:], lonlim)
        j = inlim(ds["latitude"][:], latlim)
        (; lon=ds["longitude"][i], lat=ds["latitude"][j],
           (Symbol(v) => Float64.(ds[v][i, j, k]) for v in vars)...,
           dim, dimval=ds[dim][k], file, k)
    end
end

"name of the group dimension (pentad or phase) of open composite file ds: the last dimension of its 3-D variables"
group_dim(ds) = last(dimnames(ds[first(k for k in keys(ds) if ndims(ds[k]) == 3)]))

"""
    composite_file(stem, years; dir="../out")

Path of the composite file in `dir` named <stem>_<first>-<last>.nc whose period runs from
year first(years) to year last(years), e.g. composite_file("ekman_pentad_clim_oisst", 2012:2025)
for the climatology over 2012-2025, or composite_file("ekman_pentad_clim_oisst", 2023:2023)
for the single year 2023. The stem is <data>_<grouping>_<sst> (no <sst> for data sets
without one); first and last are dates, yyyymmdd, or years, yyyy. Errors unless exactly one
file matches.
"""
function composite_file(stem, years; dir="../out")
    fs = filter(readdir(dir)) do f
        startswith(f, stem * "_") || return false
        # the period after the stem: yyyy[mmdd]-yyyy[mmdd].nc
        m = match(r"^(\d{4})(?:\d{4})?-(\d{4})(?:\d{4})?\.nc$", f[ncodeunits(stem)+2:end])
        !isnothing(m) && parse(Int, m[1]) == first(years) && parse(Int, m[2]) == last(years)
    end
    length(fs) == 1 ||
        error("$(length(fs)) files $(stem)_<first>-<last>.nc in $dir for years $(first(years))-$(last(years)): $fs")
    joinpath(dir, only(fs))
end

"name and values of the group dimension of a composite file, e.g. (\"phase\", [1, ..., 8])"
group_values(file) = NCDataset(file) do ds
    dim = group_dim(ds)
    dim, ds[dim][:]
end

"""
    group_index(file, values)

Indices along the group dimension of the group `values` (e.g. pentads 25:30, phases 1:8).
Errors, saying what the file has, if any value is not in the file.
"""
function group_index(file, values)
    dim, have = group_values(file)
    missingv = setdiff(values, have)
    isempty(missingv) || error("$dim $missingv not in $file (it has $(extrema(have)))")
    [findfirst(==(v), have) for v in values]
end

"""
    read_mean(file, k; vars=composite_vars, lonlim=nothing, latlim=nothing)

Mean of the composite fields over the groups at indices k, weighting each group's mean by
its count of daily values (nobs_v), so it is the mean over all the days. Returns the same
NamedTuple as read_group (lon x lat fields), with dimval the group values averaged.
"""
function read_mean(file, k; vars=composite_vars, lonlim=nothing, latlim=nothing)
    d = read_group(file, k; vars, lonlim, latlim)
    w = read_group(file, k; vars=map(v -> "nobs_$v", vars), lonlim, latlim)
    function wmean(x, n)
        ok = @. isfinite(x) & (n > 0)
        s = dropdims(sum(@. ifelse(ok, n*x, 0.0); dims=3); dims=3)
        c = dropdims(sum(@. ifelse(ok, n, 0.0); dims=3); dims=3)
        @. ifelse(c > 0, s/c, NaN)
    end
    (; lon=d.lon, lat=d.lat,
       (Symbol(v) => wmean(d[Symbol(v)], w[Symbol("nobs_$v")]) for v in vars)...,
       dim=d.dim, dimval=collect(d.dimval), file, k)
end

"number of days composited in group index k (the maximum of nobs_taux over the region)"
group_days(file, k; lonlim=nothing, latlim=nothing) =
    Int(maximum(read_group(file, k; vars=("nobs_taux",), lonlim, latlim).nobs_taux))

"anomaly of x from its mean over dimension `dims` (e.g. phase minus the 8-phase mean)"
anom(x; dims=3) = x .- mean(x; dims)

"""
    box_mean(x, lon, lat, lonlim, latlim; land=nothing)

Area-weighted (cos latitude) mean of x (lon x lat, or lon x lat x group) over the box
`lonlim`, `latlim` (degrees), over the points where x is finite and `land` (lon x lat) is
false: a number for lon x lat x, else a vector over the groups (a time series of pentads).
"""
function box_mean(x, lon, lat, lonlim, latlim; land=nothing)
    i = findall(v -> lonlim[1] <= v <= lonlim[2], lon)
    j = findall(v -> latlim[1] <= v <= latlim[2], lat)
    w = cosd.(lat[j])' # 1 x lat
    ocean = isnothing(land) ? trues(length(i), length(j)) : .!land[i, j]
    m = map(eachslice(view(x, i, j, :); dims=3)) do xs
        ok = isfinite.(xs) .& ocean
        sum(ifelse.(ok, w .* xs, 0.0)) / sum(ifelse.(ok, w, 0.0))
    end
    ndims(x) == 2 ? only(m) : m
end

"dq_sat/dT (kg/kg/K) at T (K): centered difference of the Buck saturation humidity qsat_dew, as used for q"
dqsat_dT(T; dT=0.01) = (qsat_dew(T + dT) - qsat_dew(T - dT)) / 2dT

"""
    ekman_fields(d; Ro_rad)

Fields to plot from composites `d` (`read_group`, `read_mean`, `read_period`, `read_phase`),
with the Ekman transport curtailed within ~`Ro_rad` (m) of the equator by
r = curtail(lat, Ro_rad), as a NamedTuple (lon x lat arrays, or lon x lat x group when d is):
SST sst, 2 m air temperature t2 (K), specific humidity q (kg/kg); wind stress taux, tauy (N/m^2, not curtailed); Mx, My (kg/m/s); their divergence divM
(kg/m^2/s) and Ekman pumping w_ek = divM/rho_o = curl(r τ/f)/rho_o (m/s, upward > 0); h_a, h_o (J/kg); their
gradients dhodx, dhody, dhadx, dhady (J/kg/m); Ekman advection of temperature advT_o (SST),
advT_a (t2) in K kg/m^2/s, and, in W/m^2, Ekman advection adv_o, adv_a, s times transport
divergence sdiv_o, sdiv_a, and flux-form tendency flux_o, flux_a; for the air-sea transition
zone (ASTZ), the air moist enthalpy saturated at SST h_star = h_* = c_pa SST + Lv q_sat(SST) and
ha_prime = h'_a = h_* - h_a = c_pa (SST - t2) + Lv (q_sat(SST) - q) (J/kg), and in W/m^2 adv_star =
-M_o⋅∇h_* and adv_ha_prime = -M_o⋅∇h'_a, so the ASTZ sum adv_o + adv_a = -M_o⋅∇(h_o - h_a) =
(adv_o - adv_star) + adv_ha_prime exactly; r; land (lon x lat, true
where SST is NaN, in the first group); lon, lat; and dim, dimval from d.
"""
function ekman_fields(d; Ro_rad)
    r  = curtail.(d.lat', Ro_rad)
    rx(x) = @. ifelse(r == 0, 0.0, r*x) # r x; 0 at the equator, where uncurtailed M is Inf
    Mx, My = ekman_transport_xy(d.taux, d.tauy, d.lat; Ro_rad)
    # Ekman pumping: ∇⋅M = curl(r τ/f), exact from the composite-mean stress (M is linear in τ)
    divM = divergence(Mx, My, d.lon, d.lat) # kg/m^2/s
    w_ek = divM ./ rho_o                     # m/s, upward > 0
    h_a = @. c_pa*d.t2 + Lv*d.q # J/kg, T in K
    h_o = @. c_po*d.sst         # J/kg, SST in K
    # gradients of the composite means
    dhodx, dhody = gradient(h_o, d.lon, d.lat)
    dhadx, dhady = gradient(h_a, d.lon, d.lat)
    # -r M⋅∇s from the saved composites of daily terms; the atmospheric ones already use M_a = -M
    adv_o = rx(c_po*d.adv_sst)              # -M_o⋅∇h_o
    adv_a = rx(c_pa*d.adv_t2 + Lv*d.adv_q)  # -M_a⋅∇h_a
    advT_o = rx(d.adv_sst)                  # -M_o⋅∇SST
    advT_a = rx(d.adv_t2)                   # -M_a⋅∇t2
    # flux form -∇⋅(r M s), with r M s = ekman_transport of the composite s τ
    fluxdiv(staux, stauy) = divergence(ekman_transport_xy(staux, stauy, d.lat; Ro_rad)..., d.lon, d.lat)
    flux_o = -c_po*fluxdiv(d.sst_taux, d.sst_tauy)                                # -∇⋅(M_o h_o)
    flux_a = c_pa*fluxdiv(d.t2_taux, d.t2_tauy) + Lv*fluxdiv(d.q_taux, d.q_tauy)  # -∇⋅(M_a h_a), M_a = -M_o
    sdiv_o = adv_o - flux_o                                                       # h_o ∇⋅M_o
    sdiv_a = adv_a - flux_a                                                       # h_a ∇⋅M_a
    # ASTZ: air moist enthalpy saturated at SST, h_* = c_pa SST + Lv q_sat(SST), and h_a = h_* - h'_a
    h_star = @. c_pa*d.sst + Lv*qsat_dew(d.sst) # J/kg
    ha_prime = h_star .- h_a                     # J/kg, h'_a = c_pa (SST - t2) + Lv (q_sat(SST) - q)
    # -M_o⋅∇h_* = (c_pa + Lv q_sat'(SST)) (-M_o⋅∇SST): linearized, q_sat' at the composite-mean SST
    # times the composite of daily -M_o⋅∇SST (-M⋅∇q_sat(SST) was not composited)
    adv_star = rx(@. (c_pa + Lv*dqsat_dT(d.sst)) * d.adv_sst) # -M_o⋅∇h_*
    adv_ha_prime = adv_star + adv_a # -M_o⋅∇h'_a = -M_o⋅∇h_* + M_o⋅∇h_a (M_a = -M_o), the remainder of the ASTZ sum
    land = isnan.(selectdim(d.sst, 3, 1)) # lon x lat for 2-D and 3-D d
    (; d.lon, d.lat, d.sst, d.t2, d.q, d.taux, d.tauy, Mx, My, divM, w_ek, h_a, h_o, dhodx, dhody, dhadx, dhady,
       advT_o, advT_a, adv_o, adv_a, sdiv_o, sdiv_a, flux_o, flux_a,
       h_star, ha_prime, adv_star, adv_ha_prime, r, land, d.dim, d.dimval)
end

## Pentad climatology (ekman_mse.jl output: group dimension "pentad")

"pentad (1-73) of date d, on the 365-day calendar (Feb 29 joins pentad 12), as in ekman_mse.jl"
pentad(d::Date) = (dayofyear(d) - (isleapyear(d) && month(d) > 2) - 1) ÷ 5 + 1

"first and last day of pentad p, on the 365-day calendar of 2001 (not a leap year)"
pentad_dates(p) = (Date(2001, 1, 1) + Day(5p - 5), Date(2001, 1, 1) + Day(5p - 1))

"""
    pentads_for(first_day, last_day; n=nothing)

Pentads (a range) covering the dates first_day:last_day. With `n`, the n consecutive
pentads that overlap those dates by the most days; ties go to the pentads whose middle is
nearest the middle of the dates. E.g. the 6 pentads most closely aligned with May:
pentads_for(Date(2001, 5, 1), Date(2001, 5, 31); n=6) == 25:30 (May 1-30).
"""
function pentads_for(first_day::Date, last_day::Date; n=nothing)
    ps = pentad.(first_day:Day(1):last_day)
    isnothing(n) && return minimum(ps):maximum(ps)
    mid = (pentad(first_day) + pentad(last_day)) / 2 # in pentads
    score(p0) = (count(in(p0:p0+n-1), ps), -abs(p0 + (n - 1)/2 - mid))
    p0 = argmax(score, 1:74-n)
    p0:p0+n-1
end

"the days covered by consecutive `pentads`, e.g. \"May 1-30\" or \"May 26-Jun 4\""
function period_label(pentads)
    a, b = first(pentad_dates(first(pentads))), last(pentad_dates(last(pentads)))
    Dates.format(a, "u d") * "-" * Dates.format(b, month(a) == month(b) ? "d" : "u d")
end

"`pentads` for file names, e.g. \"p25-30\""
period_tag(pentads) = "p$(first(pentads))-$(last(pentads))"

"""
    read_period(file, first_day, last_day; n=nothing, kw...)

Mean over all the days of the pentads `pentads_for(first_day, last_day; n)` (year ignored)
of a pentad composite file (`read_mean`); `kw`: vars, lonlim, latlim. dimval = the pentads.
"""
read_period(file, first_day::Date, last_day::Date; n=nothing, kw...) =
    read_mean(file, group_index(file, pentads_for(first_day, last_day; n)); kw...)

"""
    ekman_fields(file, first_day::Date, last_day::Date; n=nothing, lonlim=nothing, latlim=nothing, Ro_rad)

Fields of the mean over a period of a pentad composite file:
`ekman_fields(read_period(file, first_day, last_day; n, lonlim, latlim); Ro_rad)`.
E.g. May 1-30: `ekman_fields(file, Date(2001, 5, 1), Date(2001, 5, 31); n=6, lonlim, latlim, Ro_rad)`.
"""
ekman_fields(file, first_day::Date, last_day::Date; n=nothing, lonlim=nothing, latlim=nothing, Ro_rad) =
    ekman_fields(read_period(file, first_day, last_day; n, lonlim, latlim); Ro_rad)

"""
    flux_fields(turbfile, radfile, first_day::Date, last_day::Date; n=nothing, lonlim=nothing, latlim=nothing)

Surface heat fluxes, positive from ocean to atmosphere (W/m^2), averaged over the pentads
`pentads_for(first_day, last_day; n)` of an OAFlux turbulent composite file `turbfile`
(lh, sh, positive upward; ts, SST in °C) and a CERES composite file `radfile` (surface
sw_dn, sw_up, lw_dn, lw_up, positive magnitudes) on the same grid and over the same days
(errors if the files' first_day/last_day differ). Returns a NamedTuple: lh, sh, turb = lh + sh;
sw = sw_dn - sw_up (net shortwave into the ocean), lw = lw_up - lw_dn (net longwave out of the
ocean), rad = lw - sw; net Q = turb + rad; ts; land (true where lh is NaN; OAFlux is
ocean only); lon, lat, dim, dimval.
"""
function flux_fields(turbfile, radfile, first_day::Date, last_day::Date; n=nothing, lonlim=nothing, latlim=nothing)
    period(f) = NCDataset(ds -> (ds.attrib["first_day"], ds.attrib["last_day"]), f)
    period(turbfile) == period(radfile) ||
        error("different periods: $turbfile $(period(turbfile)), $radfile $(period(radfile))")
    t = read_period(turbfile, first_day, last_day; n, vars=("lh", "sh", "ts"), lonlim, latlim)
    r = read_period(radfile, first_day, last_day; n, vars=("sw_dn_sfc", "sw_up_sfc", "lw_dn_sfc", "lw_up_sfc"), lonlim, latlim)
    t.lon == r.lon && t.lat == r.lat || error("$turbfile and $radfile are on different grids")
    turb = t.lh .+ t.sh
    sw = r.sw_dn_sfc .- r.sw_up_sfc
    lw = r.lw_up_sfc .- r.lw_dn_sfc
    rad = lw .- sw
    (; t.lon, t.lat, t.lh, t.sh, turb, sw, lw, rad, Q=turb .+ rad, t.ts, land=isnan.(t.lh), t.dim, t.dimval)
end

## BSISO phases (bsiso_phase.jl output: group dimension "phase", active days only)

"""
    read_phase(file, p; kw...)

Composites of BSISO phase p (1-8) of a phase composite file: lon x lat fields for one
phase, lon x lat x phase for several (e.g. 1:8, indexed d.taux[:, :, q] for p[q]).
`kw`: vars, lonlim, latlim (see `read_group`).
"""
function read_phase(file, p; kw...)
    dim, _ = group_values(file)
    dim == "phase" || error("$file has group dimension $dim, not phase")
    k = group_index(file, p isa Integer ? [p] : p)
    read_group(file, p isa Integer ? only(k) : k; kw...)
end

"phase p and its number of active days, e.g. \"phase 4, 37 days\" (BSISO1)"
phase_label(file, p::Integer; lonlim=nothing, latlim=nothing) =
    "phase $p, $(group_days(file, only(group_index(file, [p])); lonlim, latlim)) days"

"phases for file names, e.g. \"phase4\", \"phase1-8\""
phase_tag(p) = p isa Integer ? "phase$p" : "phase$(first(p))-$(last(p))"

"""
    ekman_fields(file, p::Integer | phases; lonlim=nothing, latlim=nothing, Ro_rad)

Fields of BSISO phase p, or of several phases (lon x lat x phase), of a phase composite file:
`ekman_fields(read_phase(file, p; lonlim, latlim); Ro_rad)`.
"""
ekman_fields(file, p::Union{Integer, AbstractVector{<:Integer}}; lonlim=nothing, latlim=nothing, Ro_rad) =
    ekman_fields(read_phase(file, p; lonlim, latlim); Ro_rad)
