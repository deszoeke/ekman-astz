# plot/

All plotting for ekman-astz. The physics, conventions and data notes are in the root `CLAUDE.md`; this file covers how plots are made.

## Environment
`plot/` has its own `plot/Project.toml` (PythonPlot, NCDatasets, Revise, IJulia, and EkmanCommon developed from `../EkmanCommon`) and CondaPkg in `plot/.CondaPkg`. Plots are made interactively in Jupyter notebooks, using the "Julia" kernel (`~/Library/Jupyter/kernels/julia` on the Mac; on the HPC see the root `CLAUDE.md`). It runs the juliaup launcher, so it follows Julia updates, with `--project=@.`, so it finds `plot/Project.toml`. Don't pin Julia versions in kernels or manifests; let Julia update, recompile and re-resolve manifests itself.

## Libraries
Two libraries, loaded by notebooks with `using Revise; includet(...)`, each in sections from general to specific:
- `composites.jl` (no plotting): **general**, any composite file (`composite_file(stem, years)` finds the file in `../out` for a data set and period, e.g. `composite_file("ekman_pentad_clim_oisst", 2012:2025)`, or `2023:2023` for one year; `read_group`, `group_index`, `read_mean`, `group_days`, `anom`, `ekman_fields(d; Ro_rad)`); **pentad climatology** (`pentads_for`, `period_label`, `period_tag`, `read_period`, `ekman_fields(file, first_day::Date, last_day::Date; n, lonlim, latlim, Ro_rad)`); **BSISO phases** (`read_phase`, `phase_label`, `phase_tag`, `ekman_fields(file, p; ...)` for a phase or phases). The pentad and BSISO sections use only the general section, never each other.
- `maps.jl` (PythonPlot only): **settings** (`plot_defaults!`); **general maps**: a figure is a list of panel specs, `panel(x, title, units; arrows=arrows(u, v; scale, key, keypos), contours=contours(y; levels))`, drawn by `plot_panels(f, panels; xlim, ylim, aspect, background)` (background: light gray where NaN, e.g. land), which returns `(; fig, axs, meshes, cbars, quivers, contoursets)` via `draw_panel!`, `draw_arrows!`, `draw_contours!`; **overlays of `ekman_fields`** (`stress_arrows(f)`, `sst_contours(f)`). Unit labels: `flux_units`, `M_units`, `h_units`, `advT_units`, `tau_units`. Maps are rasterized, so SVG and PDF output stays small.

## Writing figures
- There are no figure functions: each notebook cell writes its figure out as `plot_panels(f, [panel(...), ...]; suptitle, ylim)`, with color limits set at the top of the cell, so titles, fields, units and limits are visible where the figure is made.
- Notebooks choose the data with one call, e.g. `f = ekman_fields(file, Date(2001,5,1), Date(2001,5,31); n=6, lonlim, latlim, Ro_rad)`, make figures with the shared functions, and adjust appearance after the call in the same cell.
- A value keeps one name as it is passed down (e.g. `arrows`: keyword, `panel.arrows`, `draw_arrows!` argument); in notebooks the wind stress arrows are `tau` (`arrows` is the library function).
- Figures are shown at the end of the cell by the IJulia kernel ("Julia"); there is no `display(fig)`. Display settings that depend on the platform (dpi, fonts, `PythonPlot.svg`) are in each notebook's first cell.
- Plot code takes constants and analysis functions from EkmanCommon, not local copies; short formulas (e.g. h_a = c_pa T + Lv q) may be written out in library functions for readability.
- Enthalpy terms (h_o, h_a) are computed here, from the saved composites.

## Notebooks and scripts
- `ekman_plots.ipynb`: pentad climatology averaged over a period (May set up; June as an example of another period; ASTZ split of the Ekman advection sum with h_* = c_pa SST + Lv q_sat(SST) and h'_a = h_* − h_a; a May-June comparison of ocean Ekman advection with SST contours and wind stress). M, h_a, h_o and their gradients come from the composite means (exact, since these are linear). The advection, s∇⋅M and flux-form terms come from the saved `adv_*`/`sdiv_*` composites of daily products.
- `ekman_years.ipynb`: May Ekman advection of single years 2023–2025 (per-year files, `composite_file(stem, y:y)`), and their anomalies from the 2012–2025 climatology.
- `ekman_flux_form.ipynb`: why the advection form is used: in the flux form, h ∇⋅M and −∇⋅(M h) nearly cancel (the flux-form figure), May.
- `bsiso_plots.ipynb`: BSISO1 phase composites: one phase `p`, and all 8 phases (`fp = ekman_fields(file, 1:8; ...)`, lon x lat x phase).
- `plot_bsiso_phase.jl` plots the phases in `BSISO_phase.txt`.
