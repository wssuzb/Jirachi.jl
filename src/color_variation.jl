export color_variation, binned_color_variation, flux2mag


isEqual(lc1::lightcurve, lc2::lightcurve) = diff(lc1.time) == diff(lc2.time) ? println(" ") : error("The Light curves should be quasi-simultaneou observed!!!")

function flux2mag(flux, err)
    mag = @. -2.5 * log10(flux) - 21.175
    mag_err = @. 2.5 / log(10) * (err / flux) # approximately!!!
    return mag, mag_err
end

"""
    color_variation(lc1::lightcurve, lc2::lightcurve, nsigma=3, erron=true, mode="mag";debug=false)

calculate color variation, reference: Su et al., (2023)
"""

function color_variation(lc1::lightcurve, lc2::lightcurve, nsigma=3, erron=true, mode="mag"; showhist=false, err_method="sun2014")

    isEqual(lc1, lc2)
    
    t = lc1.time
    
    if mode == "mag"
        y1, err1 = flux2mag(lc1.flux, lc1.err)
        y2, err2 = flux2mag(lc2.flux, lc2.err)
    else
        y1, err1 = lc1.flux, lc1.err
        y2, err2 = lc2.flux, lc2.err
    end

    time_2d = repeat(t', length(t))
    dt_all = time_2d - transpose(time_2d)

    x_wmin_2d = repeat(y1', length(y1))
    dx_wmin_all = x_wmin_2d - transpose(x_wmin_2d)

    x_wmax_2d = repeat(y2', length(y2))
    dx_wmax_all = x_wmax_2d - transpose(x_wmax_2d)
    
    idx = dt_all .> 0

    dt_all = dt_all[idx]
    dx_wmin_all = dx_wmin_all[idx]
    dx_wmax_all = dx_wmax_all[idx]
    
    num_all = dt_all
    
    theta = dx_wmax_all ./ dx_wmin_all

    if erron
        mmd = @. sqrt(dx_wmin_all ^ 2 + dx_wmax_all ^ 2)

        # generate 2d matrix and transpose it.
        ex_wmin_1e = repeat(err1', length(err1))
        ex_wmax_1e = repeat(err2' ,length(err2))
        ex_wmin_2e, ex_wmax_2e = transpose(ex_wmin_1e), transpose(ex_wmax_1e)

        # select the index where `dt_all>0`.
        ex_wmin_1e, ex_wmax_1e = ex_wmin_1e[idx], ex_wmax_1e[idx]
        ex_wmin_2e, ex_wmax_2e = ex_wmin_2e[idx], ex_wmax_2e[idx]

        if err_method == "sun2014"
            # calculate tanθ.
            dx_ratio2 = @. (dx_wmax_all / dx_wmin_all) ^ 2        
            
            # calculate error ellipse
            # σ1	
            err_mmd1 = @. ex_wmin_1e * ex_wmax_1e * sqrt((1+dx_ratio2) / (ex_wmax_1e ^ 2 + ex_wmin_1e ^ 2 * dx_ratio2))
            # σ2
            err_mmd2 = @. ex_wmin_2e * ex_wmax_2e * sqrt((1+dx_ratio2) / (ex_wmax_2e ^ 2 + ex_wmin_2e ^ 2 * dx_ratio2))
            err_mmd = @. sqrt(err_mmd1 ^ 2 + err_mmd2 ^2)

        else
            err_mmd = @. sqrt(ex_wmin_1e ^ 2 + ex_wmax_1e ^2 + ex_wmin_2e ^ 2 + ex_wmax_2e ^ 2)
        end

        err_mmd_nsigma = nsigma * err_mmd

        criterion = mmd .>= err_mmd_nsigma

        dt_all = dt_all[criterion]
        theta = theta[criterion]
    end
    
    
    int0 = sortperm(dt_all)
    dt_all_sort = dt_all[int0]
    theta_sort = theta[int0]

    idx_pos = theta_sort .> 0
    num_cut = dt_all_sort
    num_pos = dt_all_sort[idx_pos]
    
    if showhist
        return (cv=cv(dt_all_sort[idx_pos], theta_sort[idx_pos]), num_all = num_all, num_cut = num_cut, num_pos = num_pos)
    else
        return cv(dt_all_sort[idx_pos], theta_sort[idx_pos])
    end
end




function binned_color_variation(data::cv, bin_edges::AbstractArray=1:0.2:5)
    tau = log10.(data.tau)
    color = data.color
    
    bin_edges = bin_edges

    bin_all = bin(bin_edges, tau, color)
    bin_yerr = zeros(length(bin_all))

    
    for i=1: lastindex(bin_all)
        isempty(bin_all[i]) ? (bin_yerr[i] = 0) : (bin_yerr[i] = err_bootstraped(bin_all[i]))
    end

    bin_value = bin(median, bin_edges, tau, color)
    _bin_edges = collect(bin_edges)
    bin_width = round.([(_bin_edges[i+1] - _bin_edges[i]) / 2 for i=1: lastindex(_bin_edges)-1], digits=2)
    bin_center = round.(_bin_edges[2:end] .- bin_width, digits=2)
    
    res = binned_result(bin_center, bin_width, bin_value, bin_yerr)

    return res
end

