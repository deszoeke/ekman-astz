using Distributed
# inside a Slurm job (see run_slurm.sh), start one worker per Slurm task, on any of its nodes;
# otherwise use `julia --project -p N` for local workers, or none
if haskey(ENV, "SLURM_JOB_ID")
    using SlurmClusterManager
    addprocs(SlurmManager(); exeflags="--project=$(Base.active_project())")
end

# Pentad (5-day) climatologies of ERA5 surface fields and air-sea Ekman terms for April-July.
# Run with worker processes, e.g. `julia --project -p 8 ekman_mse.jl`, or on Slurm `sbatch run_slurm.sh ekman_mse.jl`;
# each worker composites whole pentads.
# for each pentad, loop over its days in every year
#   average hourly data to daily
#   compute daily Ekman transport and the nonlinear terms
#   add finite values to the pentad sums, counting them at each grid point
# climatological mean = sum / count at each grid point (NaN where nobs is 0, e.g. sst over land)

using Dates
# shared code (package EkmanCommon/), loaded on all worker processes too
@everywhere using EkmanCommon

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

# April-July: pentads 19 (Apr 1-5) through 43 (Jul 30-Aug 3)
years   = 2012:2020 # Andrea's data stops at 2020; 2026 is incomplete
pentads = 19:43
comp, nobs, lon, lat = pentad_climatology(years, pentads; dir=era5dir(ceoasdir))
# comp, nobs, lon, lat = pentad_climatology(2012:2012, 19:20; dir=era5dir(ceoasdir)) # short test run
save_composite("ekman_pentad_clim_$(first(years))-$(last(years)).nc", comp, nobs, lon, lat, "pentad", pentads;
               dimattrib=["long_name" => "pentad of year (1-73), 365-day calendar, Feb 29 in pentad 12"],
               attrib=["title" => "ERA5 pentad climatology of surface fields and Ekman terms",
                       "years" => "$(first(years))-$(last(years))"])

taux, tauy, sst, t2, d2, q = (comp[k] for k in fieldkeys)
adv_sst, adv_t2, adv_q, sdiv_sst, sdiv_t2, sdiv_q = (comp[k] for k in nlkeys)

# compute full advection and scalar divergence of 
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

# Compute advection of the mean by the mean Ekman mass transport

# mean specific moist static energy of surface air, m_a = c_pa*Ta + Lv*q (J/kg; T in K, z=2 m term neglected)
m_a = @. c_pa*t2 + Lv*q
# specific enthalpy of ocean surface water, h_o = c_po*SST (J/kg; SST in K)
h_o = @. c_po*sst

# mean Ekman mass transport of climatological mean stress (lon x lat x pentad), curtailed within ~Ro_rad of the equator
Mx, My = ekman_transport_xy(taux, tauy, lat; Ro_rad=Ro_rad)
divM = divergence(Mx, My, lon, lat) # mean divergence of Ekman mass transport
# gradients of ocean enthalpy and surface air moist static energy (J/kg/m)
dhodx, dhody = gradient(h_o, lon, lat)
dmadx, dmady = gradient(m_a, lon, lat)
# advection of climatological means by climatological Ekman mass transport (W/m^2): -M⋅∇
# ocean Ekman transport is M; atmospheric Ekman transport is equal and opposite, -M
adv_ho = @. -(Mx*dhodx + My*dhody)
adv_ma = @.  (Mx*dmadx + My*dmady)
