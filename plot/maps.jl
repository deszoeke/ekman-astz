# Maps: plotting lon x lat fields with PythonPlot. No reading or analysis here (that is
# composites.jl). Sections, from general to specific:
#   Settings                      plot_defaults!
#   General maps                  panel, arrows, contours, plot_panels, draw_*! (any lon x lat fields)
#   Figures of ekman_fields       stress_arrows, sst_contours, ocean_advection_panel,
#                                 fig_transport_mse, fig_adv, fig_flux (any period or phase)
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
#
# A figure is a list of panel specs, drawn by plot_panels. Each spec is a NamedTuple made by
# panel (colors), with optional overlays made by arrows and contours. A value keeps its name
# all the way down: panel(...; arrows, contours) -> panel.arrows, panel.contours ->
# draw_arrows!(ax, ..., arrows), draw_contours!(ax, ..., contours) -> ax.quiver, ax.contour.

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
    panel(x, title, units; cmap="RdBu_r", arrows=nothing, contours=nothing, kw...)

A map panel of x (lon x lat) in colors, with `arrows` (see `arrows`) and `contours` (see
`contours`) drawn on top if given. `kw` go to pcolormesh, e.g. `vmin=-180, vmax=180`; with
neither vmin nor vmax the color limits are ±symlim of the plotted (land-masked) x, and
`vmin=nothing, vmax=nothing` lets matplotlib scale to the data.
"""
panel(x, title, units; cmap="RdBu_r", arrows=nothing, contours=nothing, kw...) =
    (; x, title, units, cmap, arrows, contours, kw=NamedTuple(kw))

"""
    arrows(u, v; every=12, scale=0.25, key=nothing, keylabel="", keypos=(0.9, 0.08), kw...)

Vectors (u, v) (lon x lat) for `panel(...; arrows)`: every `every`-th grid point in each
direction, `scale` data units per inch of arrow (smaller is longer; the same on every panel
and figure), and a key arrow of length `key` labelled `keylabel` at `keypos` (axes fraction,
(1, 1) is the upper right corner; label above the arrow) if `key` is given. `kw` go to
quiver, e.g. `color="w", width=0.003`. Drawn by `draw_arrows!`.
"""
arrows(u, v; every=12, scale=0.25, key=nothing, keylabel="", keypos=(0.9, 0.08), kw...) =
    (; u, v, every, scale, key, keylabel, keypos, kw=NamedTuple(kw))

"""
    contours(x; levels=nothing, colors="k", linewidths=0.7, label=true, fmt="%g", kw...)

Contour lines of x (lon x lat) for `panel(...; contours)`, at `levels` (default:
matplotlib's), labelled with `fmt` if `label`. `kw` go to contour, e.g. `linestyles="-"`.
Drawn by `draw_contours!`.
"""
contours(x; levels=nothing, colors="k", linewidths=0.7, label=true, fmt="%g", kw...) =
    (; x, levels, label, fmt, kw=(; colors, linewidths, kw...))

"""
    plot_panels(f, panels; ncols=2, suptitle="", xlim=nothing, ylim=nothing, aspect="equal",
                figsize=(7ncols, 3.5nrows))

Draw `panels` (see `panel`) on the grid f.lon, f.lat (e.g. `f` from `ekman_fields`) with
f.land masked, in a grid of `ncols` columns sharing axes, latitude labels on the first column
and longitude labels on the last row. `xlim`, `ylim` (degrees) and `aspect` apply to all
axes; aspect "equal" makes a degree of latitude and of longitude the same length.

Returns `(; fig, axs, meshes, cbars, quivers, contoursets)`, which can also be unpacked in
order, `fig, axs = plot_panels(...)`: axs is a Python array of axes indexed from 0 as
axs[row, col]; the others are lists in panel order of the matplotlib pcolormeshes,
colorbars, quivers and contour sets (nothing for a panel without arrows or contours).
"""
function plot_panels(f, panels; ncols=2, suptitle="", xlim=nothing, ylim=nothing, aspect="equal",
                     figsize=(7ncols, 3.5cld(length(panels), ncols)))
    nrows = cld(length(panels), ncols)
    fig, axs = subplots(nrows, ncols; figsize, sharex=true, sharey=true, layout="compressed", squeeze=false)
    meshes, cbars, quivers, contoursets = [], [], [], []
    for (n, panel) in enumerate(panels)
        row, col = divrem(n-1, ncols)
        ax = axs[row, col]
        mesh, cbar, quiver, contourset = draw_panel!(ax, f.lon, f.lat, panel; f.land)
        push!(meshes, mesh); push!(cbars, cbar); push!(quivers, quiver); push!(contoursets, contourset)
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
    (; fig, axs, meshes, cbars, quivers, contoursets)
end

"""
    draw_panel!(ax, lon, lat, panel; land=nothing)

Draw `panel` (see `panel`) on matplotlib axis `ax`: a rasterized pcolormesh of panel.x
(lon x lat), land masked, a colorbar labelled panel.units, title panel.title, and
panel.arrows and panel.contours if given. Returns (mesh, cbar, quiver, contourset).
"""
function draw_panel!(ax, lon, lat, panel; land=nothing)
    x = isnothing(land) ? panel.x : mask(panel.x, land)
    lims = haskey(panel.kw, :vmin) || haskey(panel.kw, :vmax) ? (;) : (l = symlim(x); (vmin=-l, vmax=l))
    # x is lon x lat; pcolormesh takes (rows=lat, cols=lon). rasterized: in SVG/PDF the mesh is
    # one embedded image, not a path per grid cell (tens of MB, minutes to render); text stays vector.
    mesh = ax.pcolormesh(lon, lat, permutedims(x); panel.cmap, shading="nearest", rasterized=true, lims..., panel.kw...)
    cbar = colorbar(mesh; ax, label=panel.units)
    ax.set_title(panel.title)
    contourset = isnothing(panel.contours) ? nothing : draw_contours!(ax, lon, lat, panel.contours; land)
    quiver = isnothing(panel.arrows) ? nothing : draw_arrows!(ax, lon, lat, panel.arrows; land)
    mesh, cbar, quiver, contourset
end

"draw `contours` (see `contours`) on axis ax, land masked; returns the contour set"
function draw_contours!(ax, lon, lat, contours; land=nothing)
    x = permutedims(isnothing(land) ? contours.x : mask(contours.x, land))
    levels = isnothing(contours.levels) ? () : (collect(contours.levels),)
    contourset = ax.contour(lon, lat, x, levels...; contours.kw...)
    contours.label && ax.clabel(contourset; contours.fmt, fontsize="small")
    contourset
end

"draw `arrows` (see `arrows`) on axis ax, land masked, with the key at arrows.keypos; returns the quiver"
function draw_arrows!(ax, lon, lat, arrows; land=nothing)
    i, j = 1:arrows.every:length(lon), 1:arrows.every:length(lat)
    thin(x) = permutedims((isnothing(land) ? x : mask(x, land))[i, j])
    quiver = ax.quiver(lon[i], lat[j], thin(arrows.u), thin(arrows.v); arrows.scale, scale_units="inches", arrows.kw...)
    isnothing(arrows.key) || ax.quiverkey(quiver, arrows.keypos..., arrows.key, arrows.keylabel; labelpos="N", coordinates="axes")
    quiver
end

## Figures of ekman_fields: f from ekman_fields (one period or phase: lon x lat fields)

"""
    stress_arrows(f; scale=0.25, key=0.1, kw...)

Wind stress `arrows` of f.taux, f.tauy: `scale` N m^-2 per inch (smaller is longer), and a
key arrow of `key` N m^-2, labelled. `kw` go to `arrows`, e.g. `keypos=(0.9, 0.85)`, `every=8`.
"""
stress_arrows(f; scale=0.25, key=0.1, kw...) =
    arrows(f.taux, f.tauy; scale, key, keylabel="$key " * tau_units, kw...)

"""
    sst_contours(f; levels=24:0.5:31, colors="0.25", kw...)

SST `contours` of f.sst in °C. `kw` go to `contours`.
"""
sst_contours(f; levels=24:0.5:31, colors="0.25", kw...) = contours(f.sst .- 273.15; levels, colors, kw...)

"""
    ocean_advection_panel(f, title; advlim=100, contours=sst_contours(f), arrows=stress_arrows(f))

A `panel` of ocean Ekman advection -M_o⋅∇m_o in colors (±`advlim` W m^-2), titled `title`,
with SST `contours` and wind stress `arrows`.
"""
ocean_advection_panel(f, title; advlim=100, contours=sst_contours(f), arrows=stress_arrows(f)) =
    panel(f.adv_o, title * raw": ocean Ekman advection $-M_o \cdot \nabla m_o$", flux_units;
          vmin=-advlim, vmax=advlim, contours, arrows)

"""
    fig_transport_mse(f; title="", mo_lims=(1200, 1212), ma_lims=(345, 365), advlim=180,
                      arrows=stress_arrows(f), kw...)

Ekman transport Mx, My, ocean surface enthalpy m_o and surface-air MSE m_a, with wind stress
`arrows`, and ocean and atmosphere Ekman advection. Color limits `mo_lims`, `ma_lims`
(kJ/kg), ±`advlim` (W/m^2). `kw` go to `plot_panels`, e.g. `ylim=(0, 25)`.
"""
function fig_transport_mse(f; title="", mo_lims=(1200, 1212), ma_lims=(345, 365), advlim=180,
                           arrows=stress_arrows(f), kw...)
    M_units = raw"kg m$^{-1}$ s$^{-1}$"
    plot_panels(f, [
        panel(f.Mx, raw"Ekman transport $M_x$ (ocean), wind stress", M_units; arrows),
        panel(f.My, raw"Ekman transport $M_y$ (ocean), wind stress", M_units; arrows),
        panel(f.m_o ./ 1e3, raw"ocean surface enthalpy $m_o = c_{po}$ SST", raw"kJ kg$^{-1}$";
              arrows, cmap="viridis", vmin=first(mo_lims), vmax=last(mo_lims)),
        panel(f.m_a ./ 1e3, raw"surface air MSE $m_a = c_{pa} T_2 + L_v q$", raw"kJ kg$^{-1}$";
              arrows, cmap="viridis", vmin=first(ma_lims), vmax=last(ma_lims)),
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
