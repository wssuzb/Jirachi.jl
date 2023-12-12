export fractional_variability
# using Unitful

function fractional_variability(flux, err)
    """
    reference: Vaughan+2003.
    """
    N = length(flux)
    S2 = 1 / (N-1) * sum((flux .- mean(flux)) .^ 2)
    
    sigma_xs2 = S2 - sum(err .^ 2) / N
    sigma_nxs2 = sigma_xs2 / mean(flux) ^ 2

    fvar = sqrt((S2 - (sum(err .^ 2) / N)) / mean(flux) ^ 2)

    err_nxs = sqrt((sqrt(2 / N) * (sum(err .^ 2) / N) / mean(flux) ^ 2) ^ 2 + (sqrt(sum(err .^ 2) / N / N) * 2 * fvar / mean(flux)) ^ 2)
	
    err_fvar = sqrt(
        (sqrt(1/2/N) * (sum(err .^ 2) / N) / (mean(flux) ^ 2 * fvar)) ^ 2 + (sqrt((sum(err .^ 2) / N)/N) * (1/mean(flux))) ^ 2
        )
    return fvar, err_fvar #sigma_nxs2, err_nxs,
end