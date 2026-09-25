using Distributed
# inside a Slurm job (see run_slurm.sh), start one worker per Slurm task, on any of its nodes;
# otherwise use `julia --project -p N` for local workers, or none
if haskey(ENV, "SLURM_JOB_ID")
    using SlurmClusterManager
    addprocs(SlurmManager(); exeflags="--project=$(Base.active_project())")
end
using DelimitedFiles

# Composites of ERA5 surface fields and air-sea Ekman terms by BSISO1 phase (1-8).
# Run with worker processes, e.g. `julia --project -p 8 bsiso_phase.jl`, or on Slurm
# `sbatch --nodes=1 --ntasks-per-node=8 run_slurm.sh bsiso_phase.jl`; each worker composites whole phases.
# for each phase, loop over its days
#   average hourly data to daily
#   compute daily Ekman transport and the nonlinear terms
#   add finite values to the phase sums, counting them at each grid point
# composite mean = sum / count at each grid point (NaN where nobs is 0, e.g. sst over land)

using Dates
# shared code (package EkmanCommon/), loaded on all worker processes too
@everywhere using EkmanCommon

# ERA5 data root, e.g. "/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC"
ceoasdir = "./data"

"convert PC1, PC2 components of BSISO into canonical phase index (1-8)"
phase(b1, b2) = mod(floor(8/(2*pi) * atan(b1, -b2)), 8) + 1

ds = readdlm("data/BSISO.INDEX.NORM.LY.data", header=true)
B = ds[1]
header = ds[2]
dt = @. Date(B[:,1], 1, 1) + Day(B[:,2] - 1)

# missing days are -999.9 in the index file (e.g. 2015-06-21)
valid = vec(all(B[:,3:6] .> -999, dims=2))

# save BSISO phases, 0 for missing days
# YEAR DOY MONTH DATE BSISO1_phase BSISO2_phase
ph1 = Int.(phase.(B[:,3], B[:,4])) .* valid
ph2 = Int.(phase.(B[:,5], B[:,6])) .* valid
writedlm("BSISO_phase.txt",
         [["YEAR" "DOY" "MONTH" "DATE" "BSISO1_phase" "BSISO2_phase"];
          year.(dt) dayofyear.(dt) month.(dt) day.(dt) ph1 ph2])

# BSISO1 amplitude sqrt(PC1^2 + PC2^2) (= column BSISO1 of the index file)
amp1 = hypot.(B[:,3], B[:,4])
ampmin = 1.0 # composite only active days, amplitude > ampmin
active = valid .& (amp1 .> ampmin)

# look up BSISO1 phase by day, for active days only
phaselookup = Dict(zip(dt[active], ph1[active]))

# days of each BSISO1 phase, 2012-2025 (2026 is incomplete); weak or missing days are left out
days = Date(2012,1,1):Day(1):Date(2025,12,31)
# days = Date(2012,1,1):Day(1):Date(2012,3,31) # short test run
phasedays = [filter(d -> get(phaselookup, d, 0) == p, days) for p in 1:8]
comp, nobs, lon, lat = composite(phasedays; dir=era5dir(ceoasdir))
save_composite("ekman_bsiso1_phase_$(year(first(days)))-$(year(last(days))).nc", comp, nobs, lon, lat,
               "phase", 1:8;
               dimattrib=["long_name" => "BSISO1 phase"],
               attrib=["title" => "ERA5 composites of surface fields and Ekman terms by BSISO1 phase",
                       "days" => "$(first(days)) to $(last(days))",
                       "BSISO1_amplitude_min" => ampmin])

taux, tauy, sst, t2, d2, q = (comp[k] for k in fieldkeys)
adv_sst, adv_t2, adv_q, sdiv_sst, sdiv_t2, sdiv_q = (comp[k] for k in nlkeys)

# composite advection and scalar divergence of ocean enthalpy and air enthalpy and moist static energy
adv_h_o = c_po * adv_sst
adv_h_a = c_pa * adv_t2
adv_m_a = adv_h_a + Lv * adv_q
sdiv_h_o = c_po * sdiv_sst
sdiv_h_a = c_pa * sdiv_t2
sdiv_m_a = sdiv_h_a + Lv * sdiv_q
# NOTE the signs of these are defined as if
# they are on opposite sides of the equation!
# adv  = -M⋅∇
# sdiv = +∇⋅M

# specific moist static energy of surface air, m_a = c_pa*Ta + Lv*q (J/kg; T in K, z=2 m term neglected)
m_a = @. c_pa*t2 + Lv*q
# specific enthalpy of ocean surface water, h_o = c_po*SST (J/kg; SST in K)
h_o = @. c_po*sst

# Ekman mass transport of composite mean stress (lon x lat x phase), curtailed within ~Ro_rad of the equator
Mx, My = ekman_transport_xy(taux, tauy, lat; Ro_rad=Ro_rad)
divM = divergence(Mx, My, lon, lat)

# gradients of ocean enthalpy and surface air moist static energy (J/kg/m)
dhodx, dhody = gradient(h_o, lon, lat)
dmadx, dmady = gradient(m_a, lon, lat)

# advection of composite means by composite Ekman mass transport (W/m^2): -M⋅∇
# ocean Ekman transport is M; atmospheric Ekman transport is equal and opposite, -M
adv_ho = @. -(Mx*dhodx + My*dhody)
adv_ma = @.  (Mx*dmadx + My*dmady)
