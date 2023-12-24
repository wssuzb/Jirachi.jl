
function Base.show(io::IO, lc::lightcurve) # ::MIME"text/plain", 
    println(io, "Loading $lightcurve, band used: ", lc.band, "-band")
    println(" ")
    println(io, )
    println(io, "\t Time: ")
    println(io, "\t\t min time: ", minimum(lc.time))
    println(io, "\t\t max time: ", maximum(lc.time))
    println(io, "\t\t mean cadence: ", mean(lc.time))
    println(io, "\t\t median cadence: ", median(lc.time))
    println(" ")
    println(io, "\t Flux: ")
    println(io, "\t\t mean flux: ", mean(lc.flux))
    println(io, "\t\t median flux: ", median(lc.flux))
    println(" ")
    println(io, "\t Error: ")
    println(io, "\t\t mean error: ", mean(lc.err))
    println(io, "\t\t median error: ", median(lc.err))
end
