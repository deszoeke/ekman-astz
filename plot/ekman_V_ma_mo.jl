# Maps from a composite NetCDF file written by save_composite (ekman_mse.jl,
# ekman_mse_test.jl, bsiso_phase.jl):
#  1. Ekman mass transport (Mx, My), surface-air MSE m_a and ocean-surface enthalpy h_o
#  2. gradients of h_o and m_a, and their Ekman advection M⋅∇s. The ocean's transport is
#     M_o = M = (tauy, -taux)/f; the atmosphere's is equal and opposite, M_a = -M_o. Both
#     come from the composite mean stress: advection of the mean by the mean.
#     M⋅∇s > 0 is transport toward higher s (the tendency -M⋅∇s is negative).
# Runs in the plot/ environment:
#   julia --project=plot plot/ekman_V_ma_mo.jl [file.nc] [index along the pentad/phase dimension, default 1]
# Saves plot/ekman_V_ma_mo_<file>_<index>.png and plot/ekman_grad_adv_<file>_<index>.png.
using NCDatasets
using PythonPlot

file = get(ARGS, 1, "ekman_pentad31_2020_test.nc")
k    = parse(Int, get(ARGS, 2, "1"))
name = "$(splitext(basename(file))[1])_$k"

# constants, Ekman transport and gradient as in EkmanCommon
const c_pa = 1004.6   # J/kg/K specific heat of dry air at constant pressure
const c_po = 3985.0   # J/kg/K specific heat of seawater
const Lv   = 2.501e6  # J/kg latent heat of vaporization
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
Mx, My = first.(M), last.(M) # ocean, M_o; the atmosphere's is -M
m_a = @. c_pa*t2 + Lv*q      # J/kg, T in K
h_o = @. c_po*sst            # J/kg, SST in K

dhodx, dhody = gradient(h_o, lon, lat) # J/kg/m
dmadx, dmady = gradient(m_a, lon, lat)
adv_o = @.  Mx*dhodx + My*dhody        # M_o⋅∇h_o, W/m^2
adv_a = @. -Mx*dmadx - My*dmady        # M_a⋅∇m_a = -M⋅∇m_a, W/m^2

"symmetric color limit: the 98th percentile of |x| over finite points"
function symlim(x)
    v = sort(abs.(filter(isfinite, vec(x))))
    isempty(v) ? 1.0 : v[max(1, round(Int, 0.98*length(v)))]
end

"a figure of pcolormesh panels (x, title, units, diverging), 2 per row; saved as plot/<prefix>_<name>.png"
function mapfigure(prefix, panels)
    nrow = cld(length(panels), 2)
    fig, axs = subplots(nrow, 2, figsize=(14, 3.5nrow), sharex=true, sharey=true, layout="constrained")
    for (n, (x, title, units, diverging)) in enumerate(panels)
        row, col = (n-1) ÷ 2, (n-1) % 2
        ax = axs[row, col] # Python array of axes, indexed from 0
        lims = diverging ? (-symlim(x), symlim(x)) : extrema(filter(isfinite, x))
        # x is lon x lat; pcolormesh takes (rows=lat, cols=lon)
        pc = ax.pcolormesh(lon, lat, permutedims(x), cmap=diverging ? "RdBu_r" : "viridis",
                           vmin=lims[1], vmax=lims[2], shading="nearest")
        colorbar(pc, ax=ax, label=units)
        ax.set_title(title)
        col == 0 && ax.set_ylabel("latitude")
        row == nrow-1 && ax.set_xlabel("longitude")
    end
    fig.suptitle("$(basename(file)): $dim $dimval (Ro_rad = $(Ro_rad/1e3) km)")
    out = joinpath(@__DIR__, "$(prefix)_$name.png")
    savefig(out, dpi=150)
    println("saved $out")
end

mapfigure("ekman_V_ma_mo",
          [(Mx, "Ekman transport Mx", "kg m⁻¹ s⁻¹", true),
           (My, "Ekman transport My", "kg m⁻¹ s⁻¹", true),
           (m_a ./ 1e3, "surface air MSE m_a = c_pa T₂ + L_v q", "kJ kg⁻¹", false),
           (h_o ./ 1e3, "ocean surface enthalpy h_o = c_po SST", "kJ kg⁻¹", false)])

# columns: ocean, atmosphere; rows: ∂/∂x, ∂/∂y, M⋅∇; gradients per km
mapfigure("ekman_grad_adv",
          [(1e3dhodx, "∂h_o/∂x", "J kg⁻¹ km⁻¹", true), (1e3dmadx, "∂m_a/∂x", "J kg⁻¹ km⁻¹", true),
           (1e3dhody, "∂h_o/∂y", "J kg⁻¹ km⁻¹", true), (1e3dmady, "∂m_a/∂y", "J kg⁻¹ km⁻¹", true),
           (adv_o, "ocean Ekman advection M_o⋅∇h_o", "W m⁻²", true),
           (adv_a, "atmosphere Ekman advection M_a⋅∇m_a, M_a = −M_o", "W m⁻²", true)])
