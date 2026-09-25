# Maps of the Ekman mass transport (Mx, My), surface-air moist static energy m_a and
# ocean-surface enthalpy m_o (h_o in the analysis code), from a composite NetCDF file
# written by save_composite (ekman_mse.jl, ekman_mse_test.jl, bsiso_phase.jl).
# Runs in the plot/ environment, which has PythonPlot and NCDatasets, keeping them (and
# CondaPkg's plot/.CondaPkg) out of the compositing project:
#   julia --project=plot plot/ekman_V_ma_mo.jl [file.nc] [index along the pentad/phase dimension, default 1]
# Saves plot/ekman_V_ma_mo_<file>_<index>.png next to this script.
using NCDatasets
using PythonPlot

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

"symmetric color limit: the 98th percentile of |x| over finite points"
function symlim(x)
    v = sort(abs.(filter(isfinite, vec(x))))
    isempty(v) ? 1.0 : v[max(1, round(Int, 0.98*length(v)))]
end

fig, axs = subplots(2, 2, figsize=(14, 7), sharex=true, sharey=true, layout="constrained")
panels = [(Mx, "Ekman transport Mx", "kg m⁻¹ s⁻¹", "RdBu_r", true),
          (My, "Ekman transport My", "kg m⁻¹ s⁻¹", "RdBu_r", true),
          (m_a ./ 1e3, "surface air MSE m_a = c_pa T₂ + L_v q", "kJ kg⁻¹", "viridis", false),
          (m_o ./ 1e3, "ocean surface enthalpy m_o = c_po SST", "kJ kg⁻¹", "viridis", false)]
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
