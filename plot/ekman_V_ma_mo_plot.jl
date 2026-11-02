# Plotting library for composite NetCDF files written by save_composite (ekman_mse.jl,
# ekman_mse_test.jl, bsiso_phase.jl): maps of the Ekman mass transport (Mx, My),
# surface-air moist static energy m_a and ocean-surface enthalpy m_o (h_o in the analysis
# code); the gradients of m_o and m_a and their Ekman advection -M⋅∇s; and s ∇⋅M and the
# flux-form tendency -∇⋅(M s) = -M⋅∇s - s ∇⋅M. The ocean's transport is
# M_o = M = (tauy, -taux)/f; the atmosphere's is equal and opposite, M_a = -M_o.
# M, m_a, m_o and the gradients are linear in the composited fields, so computing them from
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
# Plots are made in ekman_plots.ipynb, which loads this file with Revise:
#   using Revise; includet("ekman_V_ma_mo_plot.jl")
# so edits here take effect in the notebook without restarting the kernel. Runs in the plot/
# environment (PythonPlot, NCDatasets, Revise, IJulia, and EkmanCommon developed from
# ../EkmanCommon), keeping PythonPlot (and CondaPkg's plot/.CondaPkg) out of the
# compositing project.
using NCDatasets
using PythonPlot
using EkmanCommon: c_pa, c_po, Lv, curtail, ekman_transport_xy, gradient, divergence

"""
    plot_defaults!(; fonts, dpi=100)

Matplotlib defaults: sans-serif text in the first of `fonts` that is installed (not all
systems have all), and math (\$...\$: symbols, super/subscripts) in the same font; the few
symbols it lacks, e.g. ∇, come from STIX sans, bundled with matplotlib, rather than DejaVu Sans.
"""
function plot_defaults!(; fonts=["Liberation Sans", "Nimbus Sans", "Arial", "Helvetica", "Verdana"], dpi=100)
    rc = PythonPlot.matplotlib.rcParams
    rc["font.family"] = "sans-serif"
    rc["font.sans-serif"] = PythonPlot.PythonCall.pylist(fonts)
    rc["mathtext.fontset"] = "custom"
    rc["mathtext.rm"], rc["mathtext.it"], rc["mathtext.bf"] = "sans", "sans:italic", "sans:bold"
    rc["mathtext.fallback"] = "stixsans"
    rc["figure.dpi"] = dpi
    nothing
end
plot_defaults!()

## data

const composite_vars = ("taux", "tauy", "sst", "t2", "q", "adv_sst", "adv_t2", "adv_q",
                        "sst_taux", "sst_tauy", "t2_taux", "t2_tauy", "q_taux", "q_tauy")

"""
    read_group(file, k; vars=composite_vars)

Composite fields at index k along the group (pentad/phase) dimension of a composite NetCDF
file, as a NamedTuple: lon, lat, each of `vars` (lon x lat, Float64, NaN where missing),
dim and dimval (the group dimension's name and value at k), file and k.
"""
function read_group(file, k; vars=composite_vars)
    NCDataset(file; maskingvalue=NaN) do ds
        dim = last(dimnames(ds["taux"]))
        (; lon=ds["longitude"][:], lat=ds["latitude"][:],
           (Symbol(v) => Float64.(ds[v][:, :, k]) for v in vars)...,
           dim, dimval=ds[dim][k], file, k)
    end
end

"""
    ekman_fields(d; Ro_rad)

Fields to plot from `d = read_group(file, k)`, with the Ekman transport curtailed within
~`Ro_rad` (m) of the equator by r = curtail(lat, Ro_rad), as a NamedTuple (lon x lat arrays):
Mx, My (kg/m/s), m_a, m_o (J/kg), their gradients dmodx, dmody, dmadx, dmady (J/kg/m),
and, in W/m^2, Ekman advection adv_o, adv_a, s times transport divergence sdiv_o, sdiv_a,
and flux-form tendency flux_o, flux_a; r; land (true where SST is NaN); and lon, lat.
"""
function ekman_fields(d; Ro_rad)
    r  = curtail.(d.lat', Ro_rad)
    rx(x) = @. ifelse(r == 0, 0.0, r*x) # r x; 0 at the equator, where uncurtailed M is Inf
    Mx, My = ekman_transport_xy(d.taux, d.tauy, d.lat; Ro_rad)
    m_a = @. c_pa*d.t2 + Lv*d.q # J/kg, T in K
    m_o = @. c_po*d.sst         # J/kg, SST in K
    # gradients of the composite means
    dmodx, dmody = gradient(m_o, d.lon, d.lat)
    dmadx, dmady = gradient(m_a, d.lon, d.lat)
    # -r M⋅∇s from the saved composites of daily terms; the atmospheric ones already use M_a = -M
    adv_o = rx(c_po*d.adv_sst)              # -M_o⋅∇m_o
    adv_a = rx(c_pa*d.adv_t2 + Lv*d.adv_q)  # -M_a⋅∇m_a
    # flux form -∇⋅(r M s), with r M s = ekman_transport of the composite s τ
    fluxdiv(staux, stauy) = divergence(ekman_transport_xy(staux, stauy, d.lat; Ro_rad)..., d.lon, d.lat)
    flux_o = -c_po*fluxdiv(d.sst_taux, d.sst_tauy)                                # -∇⋅(M_o m_o)
    flux_a = c_pa*fluxdiv(d.t2_taux, d.t2_tauy) + Lv*fluxdiv(d.q_taux, d.q_tauy)  # -∇⋅(M_a m_a), M_a = -M_o
    sdiv_o = adv_o - flux_o                                                       # m_o ∇⋅M_o
    sdiv_a = adv_a - flux_a                                                       # m_a ∇⋅M_a
    land = isnan.(d.sst)
    (; d.lon, d.lat, Mx, My, m_a, m_o, dmodx, dmody, dmadx, dmady,
       adv_o, adv_a, sdiv_o, sdiv_a, flux_o, flux_a, r, land)
end

"x with NaN where land is true"
mask(x, land) = ifelse.(land, NaN, x)

"symmetric color limit: the p quantile of |x| over finite points"
function symlim(x; p=0.98)
    v = sort(abs.(filter(isfinite, vec(x))))
    isempty(v) ? 1.0 : v[max(1, round(Int, p*length(v)))]
end

## panels and figures

"""
    panel(x, title, units; cmap="RdBu_r", kw...)

A map panel of x (lon x lat). `kw` are passed to pcolormesh, e.g. `vmin=-180, vmax=180`;
with neither vmin nor vmax the color limits are ±symlim of the plotted (land-masked) x, and
`vmin=nothing, vmax=nothing` lets matplotlib scale to the data. Change one in a list with
e.g. `panels[3] = panel(f.m_a, "m_a", "J/kg"; cmap="magma")`.
"""
panel(x, title, units; cmap="RdBu_r", kw...) = (; x, title, units, cmap, kw=NamedTuple(kw))

"""
    plot_panels(lon, lat, panels; ncols=2, land=nothing, suptitle="", figsize=(7ncols, 3.5nrows))

Map each panel in a grid with ncols columns, sharing axes, with land (if given) masked.
Returns (fig, axs); axs is a Python array of axes indexed from 0 as axs[row, col].
"""
function plot_panels(lon, lat, panels; ncols=2, land=nothing, suptitle="",
                     figsize=(7ncols, 3.5cld(length(panels), ncols)))
    nrows = cld(length(panels), ncols)
    fig, axs = subplots(nrows, ncols; figsize, sharex=true, sharey=true, layout="constrained", squeeze=false)
    for (n, p) in enumerate(panels)
        row, col = divrem(n-1, ncols)
        ax = axs[row, col]
        x = isnothing(land) ? p.x : mask(p.x, land)
        lims = haskey(p.kw, :vmin) || haskey(p.kw, :vmax) ? (;) : (l = symlim(x); (vmin=-l, vmax=l))
        # x is lon x lat; pcolormesh takes (rows=lat, cols=lon)
        pc = ax.pcolormesh(lon, lat, permutedims(x); cmap=p.cmap, shading="nearest", lims..., p.kw...)
        colorbar(pc; ax, label=p.units)
        ax.set_title(p.title)
        col == 0 && ax.set_ylabel("latitude")
        row == nrows-1 && ax.set_xlabel("longitude")
    end
    for n in length(panels)+1:nrows*ncols # hide unused axes
        row, col = divrem(n-1, ncols)
        axs[row, col].set_visible(false)
    end
    fig.suptitle(suptitle)
    fig, axs
end

# unit labels for panels
const flux_units = raw"W m$^{-2}$"
const grad_units = raw"J kg$^{-1}$ km$^{-1}$"
