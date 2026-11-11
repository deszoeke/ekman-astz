using Distributed, Dates
using EkmanCommon

# Pentad (5-day) climatologies of ERA5 surface fields and air-sea Ekman terms for April-July.
#   script:  julia --project -p 8 ekman_mse.jl    or    sbatch run_slurm.sh ekman_mse.jl
#            followed by any of these words, in any order (see script_args):
#              data set: era5 (default), oaflux, oafluxturb or ceres
#              SST product for era5: oisst (default) or e5sst
#              first and last year, and peryear for one file per year
#            e.g. sbatch run_slurm.sh ekman_mse.jl e5sst 2012 2020
#                 sbatch run_slurm.sh ekman_mse.jl 2023 2025 peryear
#                 sbatch run_slurm.sh ekman_mse.jl ceres
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
ceresdir = "data/ceres/adjusted" # symlink to /storage/ceoas-datasets/datasets/ERA5/staging/simon/ceres/adjusted
# NOAA OISST v2.1 daily SST files (YYYY/oisst-avhrr-v02r01.YYYYMMDD.nc), 1981-present, on the HPC in the repo
oisstdir = "data/oisst"
# composite NetCDF files are written here
outdir = "out"
# default years of each data set: ERA5 with OISST SST (t2, d2 extended through 2026) or with
# ERA5 SST (Andrea's stops at 2020); all OAFlux variables (radiation and qnet stop at 2009);
# OAFlux turbulent fields and CERES radiation over the same years, to be combined after
# compositing (CERES starts in March 2000, OAFlux turbulence stops in 2022)
default_years(data, sst=:oisst) = data == :era5 ? (sst == :e5sst ? (2012:2020) : (2012:2025)) :
    (oaflux=1985:2009, oafluxturb=2000:2022, ceres=2000:2022)[data]

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
    main(; data=:era5, sst=:oisst, years, pentads, ceoasdir, oafluxdir, ceresdir, oisstdir)

Pentad climatology of the data set `data` with SST product `sst` (see `composite_data`):
save the NetCDF file <prefix>_pentad_clim_<sstname>_<first_day>-<last_day>.nc (yyyymmdd; no
<sstname>_ for data sets without an SST product) and return
(comp, nobs, lon, lat). Pentads default to April-July, 19 (Apr 1-5) through 43 (Jul 30-Aug 3);
years to `default_years(data, sst)`.
"""
function main(; data=:era5, sst=:oisst, years=default_years(data, sst), pentads=19:43,
                ceoasdir=ceoasdir, oafluxdir=oafluxdir, ceresdir=ceresdir, oisstdir=oisstdir)
    groups = [[d for y in years for d in pentaddays(y, p:p)] for p in pentads]
    comp, nobs, lon, lat, units, prefix, sstname, title =
        composite_data(data, groups; sst, era5root=ceoasdir, oafluxroot=oafluxdir, ceresroot=ceresdir,
                       oisstroot=oisstdir)
    yrs = "$(first(years))-$(last(years))"
    day0, day1 = minimum(minimum, groups), maximum(maximum, groups) # period of the composite
    ymd(d) = Dates.format(d, "yyyymmdd")
    mkpath(outdir)
    sstpart = isempty(sstname) ? "" : "$(sstname)_"
    save_composite(joinpath(outdir, "$(prefix)_pentad_clim_$(sstpart)$(ymd(day0))-$(ymd(day1)).nc"), comp, nobs, lon, lat, "pentad", pentads;
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
    # words in any order: data set, SST product, first and last year, peryear (see script_args)
    (; data, sst, years, peryear) = script_args(ARGS)
    years = something(years, default_years(data, sst))
    if peryear # a climatology of each year, each in its own file
        foreach(y -> main(; data, sst, years=y:y), years)
    else
        main(; data, sst, years)
    end
end
