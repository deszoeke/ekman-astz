using Distributed, Dates, DelimitedFiles
using EkmanCommon

# Composites of ERA5 surface fields and air-sea Ekman terms by BSISO1 phase (1-8).
#   script:  julia --project -p 8 bsiso_phase.jl    or    sbatch --nodes=1 --ntasks-per-node=8 run_slurm.sh bsiso_phase.jl
#            other data: add the argument oaflux, oafluxturb or ceres, e.g. sbatch run_slurm.sh bsiso_phase.jl ceres
#   REPL:    include("bsiso_phase.jl") defines the functions only; then, for workers,
#              addprocs(4); @everywhere using EkmanCommon
#            and e.g. main(days=Date(2012,1,1):Day(1):Date(2012,3,31))
# for each phase, loop over its days (each worker takes whole phases)
#   average hourly data to daily
#   compute daily Ekman transport and the nonlinear terms
#   add finite values to the phase sums, counting them at each grid point
# composite mean = sum / count at each grid point (NaN where nobs is 0, e.g. sst over land)

# ERA5 data root, e.g. "/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC"
ceoasdir = "./data"
# OAFlux daily data root (turbulence/, evaporation/, radiation/, netheat/), on the HPC in the repo
oafluxdir = "data/oaflux/daily"
# CERES SYN1deg-Day files (CERES_SYN1deg-Day_*.nc), on the HPC in the repo
ceresdir = "data/ceres"
# days of each data set: ERA5 (Andrea's data stop at 2020); all OAFlux variables (radiation
# and qnet stop at 2009); OAFlux turbulent fields and CERES radiation over the same days, to be
# combined after compositing (CERES starts on 2000-03-01; from 2000-04-01, as the pentad climatology; OAFlux turbulence stops in 2022)
default_days = (era5=Date(2012,1,1):Day(1):Date(2020,12,31), oaflux=Date(1985,1,1):Day(1):Date(2009,12,31),
                oafluxturb=Date(2000,4,1):Day(1):Date(2022,12,31), ceres=Date(2000,4,1):Day(1):Date(2022,12,31))

"convert PC1, PC2 components of BSISO into canonical phase index (1-8)"
phase(b1, b2) = mod(floor(8/(2*pi) * atan(b1, -b2)), 8) + 1

"""
    bsiso_phases(indexfile; ampmin=1.0, outfile="BSISO_phase.txt")

Read the BSISO index file and write the phases of every day to `outfile` (0 for missing days;
columns YEAR DOY MONTH DATE BSISO1_phase BSISO2_phase). Returns a Dict from date to BSISO1
phase for active days only (amplitude > `ampmin`, not missing).
"""
function bsiso_phases(indexfile="data/BSISO.INDEX.NORM.LY.data"; ampmin=1.0, outfile="BSISO_phase.txt")
    B, header = readdlm(indexfile, header=true)
    dt = @. Date(B[:,1], 1, 1) + Day(B[:,2] - 1)

    # missing days are -999.9 in the index file (e.g. 2015-06-21)
    valid = vec(all(B[:,3:6] .> -999, dims=2))

    ph1 = Int.(phase.(B[:,3], B[:,4])) .* valid
    ph2 = Int.(phase.(B[:,5], B[:,6])) .* valid
    writedlm(outfile,
             [["YEAR" "DOY" "MONTH" "DATE" "BSISO1_phase" "BSISO2_phase"];
              year.(dt) dayofyear.(dt) month.(dt) day.(dt) ph1 ph2])

    # BSISO1 amplitude sqrt(PC1^2 + PC2^2) (= column BSISO1 of the index file); only active days
    amp1 = hypot.(B[:,3], B[:,4])
    active = valid .& (amp1 .> ampmin)
    return Dict(zip(dt[active], ph1[active]))
end

"days of each BSISO1 phase 1:8 among `days`; weak or missing days are left out"
phase_days(phaselookup, days) = [filter(d -> get(phaselookup, d, 0) == p, days) for p in 1:8]

"""
    main(; data=:era5, days, ampmin, ceoasdir, oafluxdir, ceresdir, indexfile)

Composite the data set `data` (see `composite_data`) by BSISO1 phase over `days` (default
`default_days[data]`), save the NetCDF file <prefix>_bsiso1_phase_<first_day>-<last_day>.nc (yyyymmdd), and return
(comp, nobs, lon, lat).
"""
function main(; data=:era5, days=default_days[data], ampmin=1.0,
                ceoasdir=ceoasdir, oafluxdir=oafluxdir, ceresdir=ceresdir,
                indexfile="data/BSISO.INDEX.NORM.LY.data")
    phaselookup = bsiso_phases(indexfile; ampmin)
    comp, nobs, lon, lat, units, prefix, title =
        composite_data(data, phase_days(phaselookup, days);
                       era5root=ceoasdir, oafluxroot=oafluxdir, ceresroot=ceresdir)
    ymd(d) = Dates.format(d, "yyyymmdd")
    save_composite("$(prefix)_bsiso1_phase_$(ymd(first(days)))-$(ymd(last(days))).nc",
                   comp, nobs, lon, lat, "phase", 1:8;
                   dimattrib=["long_name" => "BSISO1 phase"], fieldunits=units,
                   attrib=["title" => "$title: composites by BSISO1 phase",
                           # period of the composite, to match files before combining their terms
                           "first_day" => string(first(days)), "last_day" => string(last(days)),
                           "BSISO1_amplitude_min" => ampmin])
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
    main(; data=Symbol(get(ARGS, 1, "era5"))) # e.g. julia --project -p 8 bsiso_phase.jl oaflux
end
