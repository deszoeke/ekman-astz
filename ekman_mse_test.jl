# Tester for ekman_mse.jl: the mean of a single pentad in a single year, the pentad
# containing 2023-06-04 (pentad 31, May 31 - Jun 4, which spans two monthly files).
# One group, so no workers are needed:  julia --project ekman_mse_test.jl
using Dates
using EkmanCommon

# ERA5 data root, e.g. "/storage/ceoas-datasets/datasets/ERA5/staging/andrea/1hr/SFC"
ceoasdir = "./data"

"pentad (1-73) of date d, on the 365-day calendar (Feb 29 joins pentad 12)"
pentad(d) = (dayofyear(d) - (isleapyear(d) && month(d) > 2) - 1) ÷ 5 + 1

d0 = Date(2023, 6, 4)
p  = pentad(d0)
days = filter(d -> pentad(d) == p, Date(year(d0),1,1):Day(1):Date(year(d0),12,31))
println("pentad $p of $(year(d0)): $(first(days)) to $(last(days))")

t = @elapsed (comp, nobs, lon, lat) = composite([days]; dir=era5dir(ceoasdir))
println("composited $(length(days)) days on a $(length(lon))x$(length(lat)) grid in $(round(t, digits=1)) s")

# summary of each field: finite points, range and mean, and the range of its daily counts
println(rpad("field", 10), lpad("finite", 9), lpad("min", 12), lpad("max", 12), lpad("mean", 12), "   nobs")
for k in keys(comp)
    x = filter(isfinite, comp[k])
    lo, hi = isempty(x) ? (NaN, NaN) : extrema(x)
    m = isempty(x) ? NaN : sum(x) / length(x)
    println(rpad(k, 10), lpad(length(x), 9), lpad(round(lo, sigdigits=4), 12), lpad(round(hi, sigdigits=4), 12),
            lpad(round(m, sigdigits=4), 12), "   ", join(extrema(nobs[k]), "-"))
end

save_composite("ekman_pentad$(p)_$(year(d0))_test.nc", comp, nobs, lon, lat, "pentad", p:p;
               dimattrib=["long_name" => "pentad of year (1-73), 365-day calendar, Feb 29 in pentad 12"],
               attrib=["title" => "TEST: ERA5 single-pentad mean of surface fields and Ekman terms",
                       "days" => "$(first(days)) to $(last(days))"])
println("saved ekman_pentad$(p)_$(year(d0))_test.nc")
