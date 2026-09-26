# Maps of the Ekman mass transport (Mx, My), surface-air moist static energy m_a and
# ocean-surface enthalpy m_o (h_o in the analysis code), from a composite NetCDF file
# written by save_composite (ekman_mse.jl, ekman_mse_test.jl, bsiso_phase.jl); and of
# the gradients of m_o and m_a and their Ekman advection -M⋅∇s. The ocean's transport is
# M_o = M = (tauy, -taux)/f; the atmosphere's is equal and opposite, M_a = -M_o. Both
# come from the composite mean stress: advection of the mean by the mean.
# -M⋅∇s > 0 is advection from higher s, a positive tendency (as adv = -M⋅∇ in ekman_mse.jl).
# Runs in the plot/ environment, which has PythonPlot and NCDatasets, keeping them (and
# CondaPkg's plot/.CondaPkg) out of the compositing project:
#   julia --project=. ekman_V_ma_mo_plot.jl [file.nc] [index along the pentad/phase dimension, default 1]
# Saves ekman_V_ma_mo_<file>_<index>.png and ekman_grad_adv_<file>_<index>.png next to this script.
using NCDatasets
using PythonPlot

# sans-serif text in the first of these fonts that is installed (not all systems have all),
# and math ($...$: symbols, super/subscripts) in the same font; the few symbols it lacks,
# e.g. ∇, come from STIX sans, bundled with matplotlib, rather than DejaVu Sans
let rc = PythonPlot.matplotlib.rcParams
    rc["font.family"] = "sans-serif"
    rc["font.sans-serif"] = PythonPlot.PythonCall.pylist(["Liberation Sans", "Nimbus Sans", "Arial", "Helvetica", "Verdana"])
    rc["mathtext.fontset"] = "custom"
    rc["mathtext.rm"], rc["mathtext.it"], rc["mathtext.bf"] = "sans", "sans:italic", "sans:bold"
    rc["mathtext.fallback"] = "stixsans"
end

file = get(ARGS, 1, "../ekman_pentad31_2020_test.nc")
k    = parse(Int, get(ARGS, 2, "1"))

# constants and Ekman transport as in EkmanCommon
const c_pa = 1004.6  # J/kg/K specific heat of dry air at constant pressure
const c_po = 3985.0  # J/kg/K specific heat of seawater
const Lv   = 2.501e6 # J/kg latent heat of vaporization
const Ω    = 7.292e-5 # rad/s Earth rotation rate
const a    = 6.371e6  # m Earth radius

"Ekman mass transport (Mx, My) = r (tauy, -taux)/f (kg/m/s), curtailed by r = 1 - exp(-y^2/Ro_rad^2) near the equator"
function ekman_transport(taux, tauy, lat, Ro_rad)
    f = 2Ω*sind(lat)
    f == 0 && return (0.0*taux, 0.0*tauy)
    r = -expm1(-(a*deg2rad(lat)/Ro_rad)^2)
    (r*tauy/f, -r*taux/f)
end

"centered-difference gradient (∂φ/∂x, ∂φ/∂y) on the sphere of φ (lon x lat, global lon); NaN at the first/last latitude"
function gradient(φ, lon, lat)
    nx, ny = size(φ)
    dφdx, dφdy = fill(NaN, nx, ny), fill(NaN, nx, ny)
    for j in 1:ny, i in 1:nx
        ip, im = mod1(i+1, nx), mod1(i-1, nx)
        c = cosd(lat[j])
        c > 1e-10 && (dφdx[i,j] = (φ[ip,j] - φ[im,j]) / (a*c*deg2rad(mod(lon[ip] - lon[im], 360))))
        1 < j < ny && (dφdy[i,j] = (φ[i,j+1] - φ[i,j-1]) / (a*deg2rad(lat[j+1] - lat[j-1])))
    end
    dφdx, dφdy
end

lon, lat, taux, tauy, sst, t2, q, Ro_rad, dim, dimval = NCDataset(file; maskingvalue=NaN) do ds
    dim = last(dimnames(ds["taux"]))
    (ds["longitude"][:], ds["latitude"][:],
     (Float64.(ds[v][:, :, k]) for v in ("taux", "tauy", "sst", "t2", "q"))...,
     ds.attrib["Ro_rad_m"], dim, ds[dim][k])
end

M  = ekman_transport.(taux, tauy, lat', Ro_rad)
Mx, My = first.(M), last.(M)
m_a = @. c_pa*t2 + Lv*q # J/kg, T in K
m_o = @. c_po*sst       # J/kg, SST in K

# gradients (J/kg/m) and Ekman advection (W/m^2); ocean M_o = M, atmosphere M_a = -M
dmodx, dmody = gradient(m_o, lon, lat)
dmadx, dmady = gradient(m_a, lon, lat)
adv_o = @. -(Mx*dmodx + My*dmody) # -M_o⋅∇m_o
adv_a = @.  (Mx*dmadx + My*dmady) # -M_a⋅∇m_a, M_a = -M

# land mask from the NaN SST, applied to every panel
land = isnan.(sst)
mask(x) = ifelse.(land, NaN, x)

"symmetric color limit: the 98th percentile of |x| over finite points"
function symlim(x)
    v = sort(abs.(filter(isfinite, vec(x))))
    isempty(v) ? 1.0 : v[max(1, round(Int, 0.98*length(v)))]
end

fig, axs = subplots(2, 2, figsize=(14, 7), sharex=true, sharey=true, layout="constrained")
panels = [(mask(Mx), raw"Ekman transport $M_x$", raw"kg m$^{-1}$ s$^{-1}$", "RdBu_r", true),
          (mask(My), raw"Ekman transport $M_y$", raw"kg m$^{-1}$ s$^{-1}$", "RdBu_r", true),
          (mask(m_a ./ 1e3), raw"surface air MSE $m_a = c_{pa} T_2 + L_v q$", raw"kJ kg$^{-1}$", "viridis", false),
          (mask(m_o ./ 1e3), raw"ocean surface enthalpy $m_o = c_{po}$ SST", raw"kJ kg$^{-1}$", "viridis", false)]
# axs is a Python array of axes, indexed from 0 as axs[row, col]
for (n, (x, title, units, cmap, diverging)) in enumerate(panels)
    ax = axs[(n-1) ÷ 2, (n-1) % 2]
    lims = diverging ? (-symlim(x), symlim(x)) : extrema(filter(isfinite, x))
    # x is lon x lat; pcolormesh takes (rows=lat, cols=lon)
    pc = ax.pcolormesh(lon, lat, permutedims(x), cmap=cmap, vmin=lims[1], vmax=lims[2], shading="nearest")
    colorbar(pc, ax=ax, label=units)
    ax.set_title(title)
    (n-1) % 2 == 0 && ax.set_ylabel("latitude")
    (n-1) ÷ 2 == 1 && ax.set_xlabel("longitude")
end
fig.suptitle("$(basename(file)): $dim $dimval (Ro_rad = $(Ro_rad/1e3) km)")

out = joinpath(@__DIR__, "ekman_V_ma_mo_$(splitext(basename(file))[1])_$k.png")
savefig(out, dpi=150)
println("saved $out")

# gradients and Ekman advection; columns: ocean, atmosphere; rows: ∂/∂x, ∂/∂y, -M⋅∇; gradients per km.
# color limits ±lim, or ±symlim(x) where lim is nothing; both advection panels share ±180 W/m^2
fig, axs = subplots(3, 2, figsize=(14, 10), sharex=true, sharey=true, layout="constrained")
panels = [(mask(1e3dmodx), raw"$\partial m_o/\partial x$", raw"J kg$^{-1}$ km$^{-1}$", nothing), (mask(1e3dmadx), raw"$\partial m_a/\partial x$", raw"J kg$^{-1}$ km$^{-1}$", nothing),
          (mask(1e3dmody), raw"$\partial m_o/\partial y$", raw"J kg$^{-1}$ km$^{-1}$", nothing), (mask(1e3dmady), raw"$\partial m_a/\partial y$", raw"J kg$^{-1}$ km$^{-1}$", nothing),
          (mask(adv_o), raw"ocean Ekman advection $-M_o \cdot \nabla m_o$", raw"W m$^{-2}$", 180),
          (mask(adv_a), raw"atmosphere Ekman advection $-M_a \cdot \nabla m_a$, $M_a = -M_o$", raw"W m$^{-2}$", 180)]
for (n, (x, title, units, lim)) in enumerate(panels)
    row, col = (n-1) ÷ 2, (n-1) % 2
    ax = axs[row, col]
    l = something(lim, symlim(x))
    pc = ax.pcolormesh(lon, lat, permutedims(x), cmap="RdBu_r", vmin=-l, vmax=l, shading="nearest")
    colorbar(pc, ax=ax, label=units)
    ax.set_title(title)
    col == 0 && ax.set_ylabel("latitude")
    row == 2 && ax.set_xlabel("longitude")
end
fig.suptitle("$(basename(file)): $dim $dimval (Ro_rad = $(Ro_rad/1e3) km)")

out = joinpath(@__DIR__, "ekman_grad_adv_$(splitext(basename(file))[1])_$k.png")
savefig(out, dpi=150)
println("saved $out")
