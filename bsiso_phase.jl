using Distributed, Dates, DelimitedFiles
using EkmanCommon

# Composites of ERA5 surface fields and air-sea Ekman terms by BSISO1 phase (1-8).
#   script:  julia --project -p 8 bsiso_phase.jl    or    sbatch --nodes=1 --ntasks-per-node=8 run_slurm.sh bsiso_phase.jl
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
    main(; days, ampmin, ceoasdir, indexfile)

Composite by BSISO1 phase over `days` (default 2012-2025; 2026 is incomplete), save the NetCDF
file, and return (comp, nobs, lon, lat).
"""
function main(; days=Date(2012,1,1):Day(1):Date(2025,12,31), ampmin=1.0,
                ceoasdir=ceoasdir, indexfile="data/BSISO.INDEX.NORM.LY.data")
    phaselookup = bsiso_phases(indexfile; ampmin)
    comp, nobs, lon, lat = composite(phase_days(phaselookup, days); dir=era5dir(ceoasdir))
    save_composite("ekman_bsiso1_phase_$(year(first(days)))-$(year(last(days))).nc", comp, nobs, lon, lat,
                   "phase", 1:8;
                   dimattrib=["long_name" => "BSISO1 phase"],
                   attrib=["title" => "ERA5 composites of surface fields and Ekman terms by BSISO1 phase",
                           "days" => "$(first(days)) to $(last(days))",
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
    main()
end
