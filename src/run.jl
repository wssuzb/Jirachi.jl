export  runall

"""
runall(
    lc1::lightcurve, lc2::lightcurve, lc_edges::AbstractArray; 
    tunits = 24 * 3600, 
    label::String,
    sf_bin_edges=1:0.1:5, 
    cv_bin_edges=1:0.2:5, 
    nsigma=3, 
    erron=true, 
    nsim=10, 
    fi_np::String="./run_all.h5", 
    lower_bounds = [0, 0, 0, 0.001], 
    upper_bounds = [10, 2e4, 2, 0.1], 
    p0=[], 
    mode="both", 
    t_fit = 10 .^ range(log10(1), log10(6e4), step=0.1, sf_noise_sigma=2)
    )


    The pipeline for run all results.
"""
function runall(
    lc1::lightcurve, lc2::lightcurve, lc_edges::AbstractArray; 
    tunits = 24 * 3600, label::String,
    sf_bin_edges=1:0.1:5, cv_bin_edges=1:0.2:5, nsigma=3, erron=true, nsim=10, fi_np::String="./run_all.h5", lower_bounds = [0, 0, 0, 0.001], upper_bounds = [10, 2e4, 2, 0.1], p0=[], mode="both", t_fit = 10 .^ range(log10(1), log10(6e4), step=0.1), sf_noise_sigma=2
    )
    
    # here the first time epoch of observation are the same

    # otherwise, 
    #   t_start = maximum([lc1.time[1], lc2.time[1]])
    #   lc1.time = lc1.time .- t_start
    #   lc2.time = lc2.time .- t_start

    lc1.time = lc1.time .- lc1.time[1]
    lc1.time = round.(lc1.time * tunits, digits=2)

    lc2.time = lc2.time .- lc2.time[1]
    lc2.time = round.(lc2.time * tunits, digits=2)

    lc1_bin_nan = bin_light_curve(lc1; lc_edges = lc_edges)
    
    lc2_bin_nan = bin_light_curve(lc2; lc_edges = lc_edges)

    lc1_bin = remove_lc_nan(lc1_bin_nan)
    lc2_bin = remove_lc_nan(lc2_bin_nan)

    flux_ratio = mean(lc1_bin.flux) / mean(lc2_bin.flux)

    println(lc1_bin)

    println(lc2_bin)

    fit_sf1 = fitsf_mcmc(lc1_bin; nsim=nsim, lb = lower_bounds , ub = upper_bounds, sf_bin_edges=sf_bin_edges, p0=p0, mode = mode)
    fit_sf2 = fitsf_mcmc(lc2_bin; nsim=nsim, lb = lower_bounds, ub = upper_bounds, sf_bin_edges=sf_bin_edges, p0=p0, mode = mode)
    
    binsf1, binsf2 = fit_sf1.binsf, fit_sf2.binsf
    
    par_1, par_1_err = fit_sf1.param, fit_sf1.param_err
    par_2, par_2_err = fit_sf2.param, fit_sf2.param_err

    chisq1, chisq2 = fit_sf1.chisquare, fit_sf2.chisquare

    t_break_1 = find_t_break(binsf1)
    t_break_2 = find_t_break(binsf2)
    
    itp1 = find_t_min(binsf1, par_1; sf_noise_sigma=sf_noise_sigma, t_fit=t_fit)
    itp2 = find_t_min(binsf2, par_2; sf_noise_sigma=sf_noise_sigma, t_fit=t_fit)

    t_min_1, sf_min_1 = itp1.t_min, itp1.sf_min
    t_min_2, sf_min_2 = itp2.t_min, itp2.sf_min

    t_fit_1, sf_fit_1 = itp1.t_fit, itp1.sf_fit
    t_fit_2, sf_fit_2 = itp2.t_fit, itp2.sf_fit

    t_used_min = maximum([t_min_1, t_min_2])
    t_used_max = minimum([t_break_1, t_break_2])


    lc1_bin_com, lc2_bin_com = get_common_lc(lc1_bin, lc2_bin)

    # save_data(lc1_bin_com, lc2_bin_com; fi_np = "./data/bin_lightcurve/lc_bin_" * string(step(lc_edges)) * "_" * label * "_nsim_" * string(nsim) * "_mode_" * string(mode) *  "_nsigma_" * string(nsigma) * "_band_" * lc1.band * "_" * lc2.band * ".txt")


    nsigma = nsigma
    erron = erron
    
    # cv in flux-Flux
    cv_flux_res = color_variation(lc1_bin_com, lc2_bin_com, nsigma, erron, "flux"; showhist=true)
    
    cv_flux = cv_flux_res.cv

    num_all = cv_flux_res.num_all
    num_cut = cv_flux_res.num_cut
    num_pos = cv_flux_res.num_pos

    bincv_flux = binned_color_variation(cv_flux, cv_bin_edges)

    # cv in mag-mag
    cv_mag_res = color_variation(lc1_bin_com, lc2_bin_com, nsigma, erron, "mag")
    cv_mag = cv_mag_res

    bincv_mag = binned_color_variation(cv_mag, cv_bin_edges)
    
    _sf = zeros(length(binsf1.x), 4, 2)
    _sf[:,:,1] = [binsf1.x binsf1.xerr binsf1.y binsf1.yerr]
    _sf[:, :, 2] = [binsf2.x  binsf2.xerr binsf2.y binsf2.yerr]
    
    _cv = zeros(length(bincv_flux.x), 4, 2)
    _cv[:,:,1] = [bincv_flux.x bincv_flux.xerr bincv_flux.y bincv_flux.yerr]
    _cv[:,:,2] = [bincv_mag.x bincv_mag.xerr bincv_mag.y bincv_mag.yerr]

    _fit = zeros(length(t_fit_1), 2, 2)
    _fit[:, :, 1] = [t_fit_1 sf_fit_1]
    _fit[:, :, 2] = [t_fit_2 sf_fit_2]

    _par = zeros(length(4), 8, 2)
    _par[:, :, 1] = [par_1[1] par_1_err[1] par_1[2] par_1_err[2] par_1[3] par_1_err[3] par_1[4] par_1_err[4]]

    _par[:, :, 2] = [par_2[1] par_2_err[1] par_2[2] par_2_err[2] par_2[3] par_2_err[3] par_2[4] par_2_err[4]]

    band_pair = [lc1.band lc2.band]
    

    # _fit_mcmc_t = zeros(nsim, length(sf_bin_edges)-1, 2)
    # _fit_mcmc_t[:,:,1] = fit_sf1.t
    # _fit_mcmc_t[:,:,2] = fit_sf2.t

    # _fit_mcmc_sf = zeros(nsim, length(sf_bin_edges)-1, 2)
    # _fit_mcmc_sf[:,:,1] = fit_sf1.sf
    # _fit_mcmc_sf[:,:,2] = fit_sf2.sf

    h5open(fi_np, "w") do file
        file["flux_ratio"] = [mean(flux_ratio), median(flux_ratio)]
        file["band"] = band_pair
        file["sf"] = _sf
        file["cv"] = _cv
        
        file["fit"] = _fit
        file["par"] = _par
        
        file["fit_chisqure"] = [chisq1, chisq2]

        file["sf_min"] = [sf_min_1, sf_min_2]
        file["t_min"] = [t_min_1, t_min_2]
        file["t_max"] = [t_break_1, t_break_2]

        file["num_all"] = num_all
        file["num_cut"] = num_cut
        file["num_pos"] = num_pos
    end

    result = (
        sf = _sf, cv = _cv, fit = _fit, par = _par, sf_min = [sf_min_1, sf_min_2], t_min = [t_min_1, t_min_2], t_max =  [t_break_1, t_break_2], num_all = num_all, num_cut = num_cut, num_pos = num_pos
    )
    
    return result
end
