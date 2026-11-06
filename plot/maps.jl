# Maps: plotting lon x lat fields with PythonPlot. No reading or analysis here (that is
# composites.jl). Sections, from general to specific:
#   Settings                      plot_defaults!
#   General maps                  panel, arrows, map_panel!, plot_panels (any lon x lat fields)
#   Figures of ekman_fields       fig_transport_mse, fig_adv, fig_flux (any period or phase)
#
# Figures are shown at the end of the notebook cell that makes them (PythonPlot does this on
# the IJulia kernel), so there is no display call; adjust a figure in the cell that makes it.
#
# Notebooks load it with Revise, so edits take effect without restarting the kernel:
#   using Revise; includet("maps.jl")
# It runs in the plot/ environment, keeping PythonPlot (and CondaPkg's plot/.CondaPkg) out
# of the compositing project.
redirect_stdio(stdout=devnull, stderr=devnull) do
    @eval using PythonPlot
end

## Settings

"""
    plot_defaults!(; fonts, dpi=100)

Matplotlib defaults: sans-serif text in those of `fonts` that matplotlib already has in its
font list (not all systems have all), else its default, DejaVu Sans; math (\$...\$: symbols,
super/subscripts) in the same font; the few symbols it lacks, e.g. ∇, come from STIX sans,
bundled with matplotlib. Fonts not in the list are never requested: a missing font can make
matplotlib search all of the system's fonts, which hangs for minutes on the HPC.
`dpi` sets the size of inline figures on screen (figsize in inches times dpi, in pixels).
Returns the fonts used.
"""
function plot_defaults!(; fonts=["Liberation Sans", "Nimbus Sans", "Arial", "Helvetica", "Verdana"], dpi=100)
    rc = PythonPlot.matplotlib.rcParams
    known = Set(PythonPlot.PythonCall.pyconvert(String, f.name) for f in PythonPlot.matplotlib.font_manager.fontManager.ttflist)
    use = filter(in(known), fonts)
    if !isempty(use)
        rc["font.family"] = "sans-serif"
        rc["font.sans-serif"] = PythonPlot.PythonCall.pylist(use)
        rc["mathtext.fontset"] = "custom"
        rc["mathtext.rm"], rc["mathtext.it"], rc["mathtext.bf"] = "sans", "sans:italic", "sans:bold"
        rc["mathtext.fallback"] = "stixsans"
    end
    rc["figure.dpi"] = dpi
    isempty(use) ? ["DejaVu Sans (default)"] : use
end
plot_defaults!()

## General maps

# unit labels for panels
const flux_units = raw"W m$^{-2}$"
const grad_units = raw"J kg$^{-1}$ km$^{-1}$"
const advT_units = raw"K kg m$^{-2}$ s$^{-1}$"
const tau_units  = raw"N m$^{-2}$"

"x with NaN where land is true"
mask(x, land) = ifelse.(land, NaN, x)

"symmetric color limit: the p quantile of |x| over finite points"
function symlim(x; p=0.98)
    v = sort(abs.(filter(isfinite, vec(x))))
    isempty(v) ? 1.0 : v[max(1, round(Int, p*length(v)))]
end

"""
    panel(x, title, units; cmap="RdBu_r", vec=nothing, kw...)

A map panel of x (lon x lat). `vec`: arrows drawn on top (see `arrows`); `con`: contour
lines drawn on top (see `contours`). `kw` are passed
to pcolormesh, e.g. `vmin=-180, vmax=180`; with neither vmin nor vmax the color limits are
±symlim of the plotted (land-masked) x, and `vmin=nothing, vmax=nothing` lets matplotlib
scale to the data. Change one in a list with
e.g. `panels[3] = panel(f.m_a, "m_a", "J/kg"; cmap="magma")`.
"""
panel(x, title, units; cmap="RdBu_r", vec=nothing, con=nothing, kw...) = (; x, title, units, cmap, vec, con, kw=NamedTuple(kw))

"""
    contours(x; levels=nothing, colors="k", linewidths=0.7, label=true, fmt="%g", kw...)

Contour lines of x (lon x lat) to draw on a panel with `panel(...; con=contours(...))`, at
`levels` (default: matplotlib's), labelled with `fmt` if `label`. `kw` go to contour,
e.g. `linestyles="-"`. E.g. SST in °C: `contours(f.sst .- 273.15; levels=24:31)`.
"""
contours(x; levels=nothing, colors="k", linewidths=0.7, label=true, fmt="%g", kw...) =
    (; x, levels, label, fmt, kw=(; colors, linewidths, kw...))

"""
    arrows(u, v; every=12, scale=0.25, key=nothing, keylabel="", keypos=(0.9, 0.05), kw...)

Vectors (u, v) (lon x lat) to draw on a panel with `panel(...; vec=arrows(...))`: every
`every`-th grid point in each direction, `scale` data units per inch of arrow (the same on
every panel and figure), and a reference arrow of length `key` labelled `keylabel` at
`keypos` (axes fraction; default inside the lower right corner, label above) if `key` is
given. `kw` go to quiver, e.g. `color="w", width=0.003`.
E.g. wind stress: `arrows(f.taux, f.tauy; key=0.1, keylabel=raw"0.1 N m\$^{-2}\$")`.
"""
arrows(u, v; every=12, scale=0.25, key=nothing, keylabel="", keypos=(0.9, 0.05), kw...) =
    (; u, v, every, scale, key, keylabel, keypos, kw=NamedTuple(kw))

"""
    map_panel!(ax, lon, lat, p; land=nothing)

Draw panel `p` (see `panel`) on matplotlib axis `ax`: a rasterized pcolormesh map of p.x
(lon x lat), land (if given) masked, colorbar with p.units, title p.title, and the arrows
p.vec and contours p.con if given. Color limits from p.kw, else ±symlim of the masked x.
Returns (pc, cb, q, c), the mesh, the colorbar, the quiver and the contour set (or nothing),
for further changes,
e.g. `cb.ax.tick_params(labelsize=8)`. Use it in a notebook's own loop over its own axes:
    fig, axs = subplots(4, 2; figsize=(12, 14), sharex=true, sharey=true, layout="compressed", squeeze=false)
    for (p, ax) in zip(panels, axs.flat); map_panel!(ax, lon, lat, p; land); end
"""
function map_panel!(ax, lon, lat, p; land=nothing)
    x = isnothing(land) ? p.x : mask(p.x, land)
    lims = haskey(p.kw, :vmin) || haskey(p.kw, :vmax) ? (;) : (l = symlim(x); (vmin=-l, vmax=l))
    # x is lon x lat; pcolormesh takes (rows=lat, cols=lon). rasterized: in SVG/PDF the mesh is
    # one embedded image, not a path per grid cell (tens of MB, minutes to render); text stays vector.
    pc = ax.pcolormesh(lon, lat, permutedims(x); cmap=p.cmap, shading="nearest", rasterized=true, lims..., p.kw...)
    cb = colorbar(pc; ax, label=p.units)
    ax.set_title(p.title)
    c = isnothing(p.con) ? nothing : contour_panel!(ax, lon, lat, p.con; land)
    q = isnothing(p.vec) ? nothing : quiver_panel!(ax, lon, lat, p.vec; land)
    pc, cb, q, c
end

"draw `contours` c on axis ax, land masked; returns the contour set"
function contour_panel!(ax, lon, lat, c; land=nothing)
    x = permutedims(isnothing(land) ? c.x : mask(c.x, land))
    cs = isnothing(c.levels) ? ax.contour(lon, lat, x; c.kw...) : ax.contour(lon, lat, x, collect(c.levels); c.kw...)
    c.label && ax.clabel(cs; fmt=c.fmt, fontsize="small")
    cs
end

"draw `arrows` a on axis ax, land masked; returns the quiver"
function quiver_panel!(ax, lon, lat, a; land=nothing)
    i, j = 1:a.every:length(lon), 1:a.every:length(lat)
    m(x) = permutedims((isnothing(land) ? x : mask(x, land))[i, j])
    q = ax.quiver(lon[i], lat[j], m(a.u), m(a.v); scale=a.scale, scale_units="inches", a.kw...)
    isnothing(a.key) || ax.quiverkey(q, a.keypos..., a.key, a.keylabel; labelpos="N", coordinates="axes")
    q
end

"""
    plot_panels(lon, lat, panels; ncols=2, land=nothing, suptitle="", xlim=nothing, ylim=nothing,
                aspect="equal", figsize=(7ncols, 3.5nrows))
    plot_panels(f, panels; kw...)

Map each panel with `map_panel!` in a grid with ncols columns, sharing axes, with land (if
given) masked, latitude labels on the first column and longitude labels on the last row.
`xlim`, `ylim` (degrees) and `aspect` apply to all the (shared) axes; aspect "equal" makes a
degree of latitude and of longitude the same length. The second form takes lon, lat and land
from the fields `f` (e.g. of `ekman_fields`). Returns (fig, axs, pcs, cbs, qs, cs); axs is a
Python array of axes indexed from 0 as axs[row, col], pcs, cbs and qs the meshes, colorbars
quivers and contour sets in panel order (`fig, axs = plot_panels(...)` also works). The "compressed"
layout keeps colorbars as tall as their maps. For other layouts, make the axes yourself and
call `map_panel!` (see there).
"""
function plot_panels(lon, lat, panels; ncols=2, land=nothing, suptitle="", xlim=nothing,
                     ylim=nothing, aspect="equal", figsize=(7ncols, 3.5cld(length(panels), ncols)))
    nrows = cld(length(panels), ncols)
    fig, axs = subplots(nrows, ncols; figsize, sharex=true, sharey=true, layout="compressed", squeeze=false)
    pcs, cbs, qs, cs = [], [], [], []
    for (n, p) in enumerate(panels)
        row, col = divrem(n-1, ncols)
        ax = axs[row, col]
        pc, cb, q, c = map_panel!(ax, lon, lat, p; land)
        push!(pcs, pc); push!(cbs, cb); push!(qs, q); push!(cs, c)
        col == 0 && ax.set_ylabel("latitude")
        row == nrows-1 && ax.set_xlabel("longitude")
    end
    for n in length(panels)+1:nrows*ncols # hide unused axes
        row, col = divrem(n-1, ncols)
        axs[row, col].set_visible(false)
    end
    for ax in axs.flat # axes are shared, but aspect is per axis
        isnothing(xlim) || ax.set_xlim(xlim...)
        isnothing(ylim) || ax.set_ylim(ylim...)
        isnothing(aspect) || ax.set_aspect(aspect)
    end
    fig.suptitle(suptitle)
    fig, axs, pcs, cbs, qs, cs
end
plot_panels(f, panels; kw...) = plot_panels(f.lon, f.lat, panels; land=f.land, kw...)

## Figures of ekman_fields: f from ekman_fields (one period or phase: lon x lat fields)

"""
    stress_arrows(f; key=0.1, scale=0.25, kw...)

Wind stress arrows (`arrows`) of f.taux, f.tauy: `scale` N m^-2 per inch (smaller is longer),
with a `key` N m^-2 arrow. `kw` go to `arrows`, e.g. `keypos=(0.92, 0.88)`, `every=8`.
"""
stress_arrows(f; key=0.1, scale=0.25, kw...) =
    arrows(f.taux, f.tauy; scale, key, keylabel="$key" * raw" N m$^{-2}$", kw...)

"""
    sst_contours(f; levels=24:0.5:31, colors="0.25", kw...)

SST contours (`contours`) of f.sst in °C. `kw` go to `contours`.
"""
sst_contours(f; levels=24:0.5:31, colors="0.25", kw...) = contours(f.sst .- 273.15; levels, colors, kw...)

"""
    fig_transport_mse(f; title="", mo_lims=(1200, 1212), ma_lims=(345, 365), advlim=180, tau, kw...)

Ekman transport Mx, My, ocean surface enthalpy m_o and surface-air MSE m_a, with wind stress
arrows `tau` (default: every 12th point, 0.25 N m^-2 per inch, 0.1 N m^-2 key), and ocean and
atmosphere Ekman advection. Color limits `mo_lims`, `ma_lims` (kJ/kg), ±`advlim` (W/m^2).
`kw` go to `plot_panels`, e.g. `ylim=(0, 25)`. Returns (fig, axs, pcs, cbs, qs, cs).
"""
function fig_transport_mse(f; title="", mo_lims=(1200, 1212), ma_lims=(345, 365), advlim=180,
        tau=stress_arrows(f), kw...)
    M_units = raw"kg m$^{-1}$ s$^{-1}$"
    plot_panels(f, [
        panel(f.Mx, raw"Ekman transport $M_x$ (ocean), wind stress", M_units; vec=tau),
        panel(f.My, raw"Ekman transport $M_y$ (ocean), wind stress", M_units; vec=tau),
        panel(f.m_o ./ 1e3, raw"ocean surface enthalpy $m_o = c_{po}$ SST", raw"kJ kg$^{-1}$";
              vec=tau, cmap="viridis", vmin=first(mo_lims), vmax=last(mo_lims)),
        panel(f.m_a ./ 1e3, raw"surface air MSE $m_a = c_{pa} T_2 + L_v q$", raw"kJ kg$^{-1}$";
              vec=tau, cmap="viridis", vmin=first(ma_lims), vmax=last(ma_lims)),
        panel(f.adv_o, raw"ocean Ekman advection $-M_o \cdot \nabla m_o$", flux_units; vmin=-advlim, vmax=advlim),
        panel(f.adv_a, raw"atmosphere Ekman advection $-M_a \cdot \nabla m_a$, $M_a = -M_o$", flux_units; vmin=-advlim, vmax=advlim)];
        suptitle=title, kw...)
end

"""
    fig_adv(f; title="", Tlim=0.05, Hlim=60, kw...)

Ocean and atmosphere Ekman advection, their sum and difference (rows). Columns: temperature
(±`Tlim`, K kg/m^2/s) and (moist) enthalpy (±`Hlim`, W/m^2). `kw` go to `plot_panels`.
"""
fig_adv(f; title="", Tlim=0.05, Hlim=60, kw...) = plot_panels(f, [
    panel(f.advT_o, raw"ocean advection $-M_o \cdot \nabla$SST", advT_units; vmin=-Tlim, vmax=Tlim),
    panel(f.adv_o, raw"ocean enthalpy advection $-M_o \cdot \nabla m_o$, $m_o = c_{po}$SST", flux_units; vmin=-Hlim, vmax=Hlim),
    panel(f.advT_a, raw"atmosphere advection $-M_a \cdot \nabla T_2$", advT_units; vmin=-Tlim, vmax=Tlim),
    panel(f.adv_a, raw"atmosphere moist enthalpy advection $-M_a \cdot \nabla m_a$, $m_a = c_{pa}T_2 + L_v q$", flux_units; vmin=-Hlim, vmax=Hlim),
    panel(f.advT_o .+ f.advT_a, "oc+atm T advection sum", advT_units; vmin=-Tlim, vmax=Tlim),
    panel(f.adv_o .+ f.adv_a, "oc+atm advection sum", flux_units; vmin=-Hlim, vmax=Hlim),
    panel(f.advT_o .- f.advT_a, "oc-atm T advection difference", advT_units; vmin=-Tlim, vmax=Tlim),
    panel(f.adv_o .- f.adv_a, "oc-atm advection difference", flux_units; vmin=-Hlim, vmax=Hlim)];
    suptitle=title, kw...)

"""
    fig_flux(f; title="", kw...)

Ocean and atmosphere m ∇⋅M, flux form -∇⋅(M m), and their sum (color limits ±symlim).
`kw` go to `plot_panels`.
"""
fig_flux(f; title="", kw...) = plot_panels(f, [
    panel(f.sdiv_o, raw"ocean $m_o \nabla \cdot M_o$", flux_units),
    panel(f.sdiv_a, raw"atmosphere $m_a \nabla \cdot M_a$", flux_units),
    panel(f.flux_o, raw"ocean flux form $-\nabla \cdot (M_o m_o)$", flux_units),
    panel(f.flux_a, raw"atmosphere flux form $-\nabla \cdot (M_a m_a)$", flux_units),
    panel(f.flux_o .+ f.sdiv_o, raw"ocean sum $-\nabla \cdot (M_o m_o) + m_o \nabla \cdot M_o$", flux_units),
    panel(f.flux_a .+ f.sdiv_a, raw"atmosphere sum $-\nabla \cdot (M_a m_a) + m_a \nabla \cdot M_a$", flux_units)];
    suptitle=title, kw...)

"""
    fig_ocean_advection(fs, titles; advlim=100, sst=sst_contours.(fs), tau=stress_arrows.(fs), kw...)

One panel per f in `fs` (e.g. periods to compare), titled by `titles`: colors ocean Ekman
advection -M_o⋅∇m_o (±`advlim` W/m^2), contours SST `sst`, quivers wind stress `tau`
(one `contours` and one `arrows` per panel). `kw` go to `plot_panels`, e.g. `ncols=1`.
"""
fig_ocean_advection(fs, titles; advlim=100, sst=sst_contours.(fs), tau=stress_arrows.(fs), kw...) =
    plot_panels(first(fs), [panel(f.adv_o, t * raw": ocean Ekman advection $-M_o \cdot \nabla m_o$", flux_units;
                                  vmin=-advlim, vmax=advlim, con=c, vec=a)
                            for (f, t, c, a) in zip(fs, titles, sst, tau)]; kw...)
