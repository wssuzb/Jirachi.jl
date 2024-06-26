export fractional_variability
struct mag end
struct flux end
function fractional_variability(flux, err)
    """
    reference: Vaughan+2003.
    """
    N = length(flux)
    x̄ = mean(flux)
    
    S² = 1 / (N-1) * sum((flux .- x̄) .^ 2)
    err² = sum(err .^ 2) / N
    
    println("S²: ", S²)
    println("err²: ", err²)

    σ_xs2 = S² - err² # σ²_rms
    
    σ_nxs2 = σ_xs2 / x̄ ^ 2


    fvar = sqrt((S² - (err²)) / x̄ ^ 2)

    err_nxs2 = sqrt((sqrt(2 / N) * err² / x̄ ^ 2) ^ 2 + ( sqrt(err² / N) * 2 * fvar / x̄ ) ^ 2) # equation 11

    # println("equation 11, first terms: ", sqrt(2 / N) * (sum(err .^ 2) / N) / mean(flux) ^ 2)
    # println("equation 11, second terms: ", sqrt(sum(err .^ 2) / N / N) * 2 * fvar / mean(flux))


    err_fvar = sqrt(
        (sqrt(1/2/N) * err² / ( x̄ ^ 2 * fvar)) ^ 2 + (sqrt( err² / N ) * ( 1/ x̄ )) ^ 2
        ) # equatuon B2
    

    σ²_rms = σ_xs2

    err_σ²_rms =  (S² / err² >= 100) ? err_σ²_rms = sqrt((sqrt(2 / N) * err² ) ^ 2 + (sqrt( err² / N ) * 2 * fvar * x̄ ) ^ 2) : sqrt(2 / N) * (sum(err .^ 2) / N)

    return fvar, err_fvar, σ_nxs2, err_nxs2, σ²_rms, err_σ²_rms
end

function fractional_variability(flux, err, ::mag)
    """
    reference: Vaughan+2003.
    """
    N = length(flux)
    x̄ = mean(flux)
    
    S² = 1 / (N-1) * sum((flux .- x̄) .^ 2)
    err² = sum(err .^ 2) / N
    
    println("S²: ", S²)
    println("err²: ", err²)

    σ_xs2 = S² - err² # σ²_rms

    fvar = sqrt((S² - (err²)) / x̄ ^ 2)

    σ²_rms = σ_xs2

    err_σ²_rms =  (S² / err² >= 100) ? err_σ²_rms = sqrt((sqrt(2 / N) * err² ) ^ 2 + (sqrt( err² / N ) * 2 * fvar * x̄ ) ^ 2) : sqrt(2 / N) * (sum(err .^ 2) / N)

    return σ²_rms, err_σ²_rms
end