using Distributed, Dates
using EkmanCommon

# Pentad (5-day) climatologies of ERA5 surface fields and air-sea Ekman terms for April-July.
#   script:  julia --project -p 8 ekman_mse.jl    or    sbatch run_slurm.sh ekman_mse.jl
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
    main(; years, pentads, ceoasdir)

Pentad climatology: save the NetCDF file and return (comp, nobs, lon, lat). Defaults:
April-July, pentads 19 (Apr 1-5) through 43 (Jul 30-Aug 3), years 2012-2020 (Andrea's data stops at 2020; 2026 is incomplete).
"""
function main(; years=2012:2020, pentads=19:43, ceoasdir=ceoasdir)
    comp, nobs, lon, lat = pentad_climatology(years, pentads; dir=era5dir(ceoasdir))
    save_composite("ekman_pentad_clim_$(first(years))-$(last(years)).nc", comp, nobs, lon, lat, "pentad", pentads;
                   dimattrib=["long_name" => "pentad of year (1-73), 365-day calendar, Feb 29 in pentad 12"],
                   attrib=["title" => "ERA5 pentad climatology of surface fields and Ekman terms",
                           "years" => "$(first(years))-$(last(years))"])
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
    main()
end
