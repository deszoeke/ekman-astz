# Tester for ekman_mse.jl: the mean of a single pentad in a single year, the pentad
# containing 2020-06-04 (pentad 31, May 31 - Jun 4, which spans two monthly files).
# One group, so no workers are needed:  julia --project ekman_mse_test.jl
include("ekman_mse.jl") # pentad, pentaddays, pentad_climatology, ceoasdir; main() does not run

d0 = Date(2020, 6, 4)
p  = pentad(d0)
days = pentaddays(year(d0), p:p)
println("pentad $p of $(year(d0)): $(first(days)) to $(last(days))")

t = @elapsed (comp, nobs, lon, lat) = pentad_climatology(year(d0):year(d0), p:p; dir=era5dir(ceoasdir))
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
