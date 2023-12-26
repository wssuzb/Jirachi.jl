
nanmean(arr) = mean(filter(!isnan, arr))
nanmedian(arr) = median(filter(!isnan, arr))

function Base.show(io::IO, lc::lightcurve) # ::MIME"text/plain", 
    
    flux_err_mag = flux2mag(lc.flux, lc.err)

    println(io, "Loading $lightcurve, band used: ", lc.band, "-band")
    println(" ")
    println(io, )
    println(io, "\t Time: ")
    println(io, "\t\t min time: ", round.(minimum(lc.time), digits=2))
    println(io, "\t\t max time: ", round.(maximum(lc.time), digits=2))
    println(io, "\t\t mean cadence: ", round.(nanmean(diff(lc.time)), digits=2))
    println(io, "\t\t median cadence: ", round.(nanmedian(diff(lc.time)), digits=2))
    println(" ")
    println(io, "\t Flux: ")
    println(io, "\t\t mean flux: ", round.(nanmean(lc.flux), digits=2))
    println(io, "\t\t median flux: ", round.(nanmedian(lc.flux), digits=2))
    println(io, "\t\t mean mag flux: ", round.(nanmean(flux_err_mag[1]), digits=2))
    println(io, "\t\t std mag flux: ", round.(std(flux_err_mag[1]), digits=2))
    println(" ")
    println(io, "\t Error: ")
    println(io, "\t\t mean error: ", nanmean(lc.err))
    println(io, "\t\t median error: ", nanmedian(lc.err))
    println(io, "\t\t mean mag err: ", round.(nanmean(flux_err_mag[2]), digits=2))
    println(io, "\t\t std mag err: ", round.(std(flux_err_mag[2]), digits=2))
end
