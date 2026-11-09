using Distributed, Dates
using EkmanCommon

# Pentad (5-day) climatologies of ERA5 surface fields and air-sea Ekman terms for April-July.
#   script:  julia --project -p 8 ekman_mse.jl    or    sbatch run_slurm.sh ekman_mse.jl
#            other data (oaflux, oafluxturb, ceres): julia --project -p 8 ekman_mse.jl ceres   or   sbatch run_slurm.sh ekman_mse.jl ceres
#   REPL:    include("ekman_mse.jl") defines the functions only; then, for workers,
#              addprocs(4); @everywhere using EkmanCommon
#            and e.g. main(years=2012:2012, pentads=19:20)
# for each pentad, loop over its days in every year (each worker takes whole pentads)
#   average hourly data to daily
#   compute daily Ekman transport and the nonlinear terms
#   add finite values to the pentad sums, counting them at each grid point
# climatological mean = sum / count at each grid point (NaN where nobs is 0, e.g. sst over land)

# ERA5 data root, e.g. "/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC"
ceoasdir = "./data"
# OAFlux daily data root (turbulence/, evaporation/, radiation/, netheat/), on the HPC in the repo
oafluxdir = "data/oaflux/daily"
# CERES SYN1deg-Day files (CERES_SYN1deg-Day_*.nc), on the HPC in the repo
ceresdir = "data/ceres"
# years of each data set: ERA5 (Andrea's data stop at 2020); all OAFlux variables (radiation
# and qnet stop at 2009); OAFlux turbulent fields and CERES radiation over the same years, to be
# combined after compositing (CERES starts in March 2000, OAFlux turbulence stops in 2022)
default_years = (era5=2012:2020, oaflux=1985:2009, oafluxturb=2000:2022, ceres=2000:2022) # 2000: from April 1 (pentad 19)

"pentad (1-73) of date d, on the 365-day calendar (Feb 29 joins pentad 12)"
pentad(d) = (dayofyear(d) - (isleapyear(d) && month(d) > 2) - 1) ÷ 5 + 1

"days of year y in the pentads"
pentaddays(y, pentads) = filter(d -> pentad(d) in pentads, Date(y,1,1):Day(1):Date(y,12,31))

"""
    pentad_climatology(years, pentads; kw...)

Climatological mean (lon x lat x pentad) over `years` of the ERA5 fields and the
nonlinear Ekman terms for each pentad in the range `pentads` (see `composite`,
which takes `kw`). Returns (comp, nobs, lon, lat).
"""
pentad_climatology(years, pentads::UnitRange; kw...) =
    composite([[d for y in years for d in pentaddays(y, p:p)] for p in pentads]; kw...)

"""
    main(; data=:era5, years, pentads, ceoasdir, oafluxdir, ceresdir)

Pentad climatology of the data set `data` (see `composite_data`): save the NetCDF file
<prefix>_pentad_clim_<first_day>-<last_day>.nc (yyyymmdd) and return (comp, nobs, lon, lat). Pentads default to
April-July, 19 (Apr 1-5) through 43 (Jul 30-Aug 3); years to `default_years[data]`.
"""
function main(; data=:era5, years=default_years[data], pentads=19:43,
                ceoasdir=ceoasdir, oafluxdir=oafluxdir, ceresdir=ceresdir)
    groups = [[d for y in years for d in pentaddays(y, p:p)] for p in pentads]
    comp, nobs, lon, lat, units, prefix, title =
        composite_data(data, groups; era5root=ceoasdir, oafluxroot=oafluxdir, ceresroot=ceresdir)
    yrs = "$(first(years))-$(last(years))"
    day0, day1 = minimum(minimum, groups), maximum(maximum, groups) # period of the composite
    ymd(d) = Dates.format(d, "yyyymmdd")
    save_composite("$(prefix)_pentad_clim_$(ymd(day0))-$(ymd(day1)).nc", comp, nobs, lon, lat, "pentad", pentads;
                   dimattrib=["long_name" => "pentad of year (1-73), 365-day calendar, Feb 29 in pentad 12"],
                   fieldunits=units,
                   attrib=["title" => "$title: pentad climatology", "years" => yrs,
                           # period of the composite, to match files before combining their terms
                           "first_day" => string(day0), "last_day" => string(day1)])
    return (; comp, nobs, lon, lat)
end

# Run only when started as a script (`julia -p N script.jl`, or sbatch run_slurm.sh), not on include
if abspath(PROGRAM_FILE) == @__FILE__
    # inside a Slurm job, start one worker per Slurm task, on any of its nodes
    if haskey(ENV, "SLURM_JOB_ID")
        using SlurmClusterManager
        addprocs(SlurmManager(); exeflags="--project=$(Base.active_project())")
    end
    @everywhere using EkmanCommon   # after the workers exist
    main(; data=Symbol(get(ARGS, 1, "era5"))) # e.g. julia --project -p 8 ekman_mse.jl oaflux
end
