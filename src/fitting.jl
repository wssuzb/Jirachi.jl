export jmodel, model, fitsf, fitsf_mcmc, find_t_min, find_t_break

# @pyimport scipy.optimize as so

@. jmodel(t, p) = sqrt(p[1] ^ 2 * (1-exp(-(t / p[2]) ^ p[3])) + 2 * p[4] ^ 2)

# _model_ = pyeval("""lambda f: lambda a, b, c, d, e: f(a, b, c, d, e)""")
# model = _model_((t, sf, tau, beta, sigma) -> sqrt.(sf ^ 2 * (1 .- exp.(-(t ./ tau) .^ beta)) .+ 2 * sigma .^ 2))

function find_t_min(data::binned_result, p::Vector{Float64}; noise_sigma=2, t_fit=10 .^ range(log10(1), log10(6e4), step=0.1))

    idx = all.(isfinite, data.y)
    t_new, sf_new = data.x[idx], data.y[idx]

    t_br = find_t_break(data)
    idx_tmax = t_new .<= t_br

    noise_cut = (noise_sigma * sqrt(2) * p[4])
    
    itp = linear_interpolation(sf_new[idx_tmax], t_new[idx_tmax], extrapolation_bc=Line())

    t_used_min = itp(noise_cut)
    # modeling
    sf_fit = jmodel(t_fit, p)
    # tmp_min = (noise_sigma * sqrt(2) * p[4])
    # noise_cut = sf_fit .>= tmp_min
    # t_used_min = t_fit[noise_cut][1]
    
    return (t_min = t_used_min, sf_min = noise_cut, t_fit = t_fit, sf_fit = sf_fit)
end

function find_t_break(data::binned_result; check_plot::Bool=false)

    idx = all.(isfinite, data.y)

    t_new, sf_new = data.x[idx], data.y[idx]

    pks, vals = findmaxima(sf_new; strict=false)

    pks, proms = peakproms(pks, sf_new; strict=false)
    _, proms = peakproms!(pks, sf_new; minprom=1e-3)
    pks, widths, leftedge, rightedge = peakwidths(pks, sf_new, proms; strict=false)
    _, widths = peakwidths(pks, sf_new, proms; minwidth=1e-3)
    
    idx_max_ = findmax(widths)[2]

    if check_plot
        println(" plotting")
        plt = plotpeaks(t_new, sf_new, peaks=[pks[idx_max_]], prominences=true, widths=true)
        display(plt)
    end

    return t_new[pks[idx_max_]]
end


function fitsf(data::binned_result, t_break::Float64, lb = [0, 0, 0, 0.001], ub = [10, 2e4, 2, 0.1])

    idx = all.(isfinite, data.y)

    t_new, sf_new = data.x[idx], data.y[idx]

    idx_tmax = t_new .<= t_break
    
    p0_bounds = (lb .+ ub) / 2

    fit_bounds = curve_fit(jmodel, 10 .^ t_new[idx_tmax], sf_new[idx_tmax], p0_bounds, lower=lb, upper=ub)

    p = fit_bounds.param
    # fit_bounds = so.curve_fit(model, 10 .^ t_new[idx_tmax], sf_new[idx_tmax], p0=collect(p0_bounds))
    # p = fit_bounds[1]    
    return p
end

check_bounds(p0_tau, t_max) = p0_tau < t_max ? println(" ") : error("The initial guess of tau parameters values in p0 should smaller than t max!!!")

function fitsf_mcmc(data::lc; nsim=1000, lb = [0, 0, 0, 0.001], ub = [10, 2e4, 2, 0.1], sf_bin_edges=1:0.1:5, p0 = [], mode = "both")



    _tmp_sf = zeros(nsim, length(sf_bin_edges)-1)
    _tmp_t = zeros(nsim, length(sf_bin_edges)-1)

    _sf = structure_function(data.time, data.flux)
        
    binsf = binned_structure_function(_sf, sf_bin_edges)
    
    if mode != "none"
        _nseed_ = ifelse(isempty(data.band), 1, Int(only(data.band)) - Int('a') + 1)

        Threads.@threads for i=1: nsim
            
            nseed = i + _nseed_ * nsim
            
            tmp = lc_bootstrapped(data; seed = nseed, mode = mode)
            sf_tmp = structure_function(tmp.time, tmp.flux)
            
            binned_sf_tmp = binned_structure_function(sf_tmp, sf_bin_edges)
        
            _tmp_sf[i, :] = binned_sf_tmp.y
            _tmp_t[i, :] = binned_sf_tmp.x

        end
        
        sf_err = zeros(length(sf_bin_edges) - 1, 1)

        for i=1: length(sf_bin_edges) - 1
            try
                # percent = percentile_16_50_84(filter(!isnan, _tmp_sf[:, i])) 
                sf_err[i, 1] = std(filter(!isnan, _tmp_sf[:, i]))
                # sf_err[i, 1] = percent.low
                # sf_err[i, 2] = percent.hig
            catch y
                # warn("Exception: ", y) # What to do on error.
                println("pass")
            end
        end

        binsf.yerr = sf_err[:, 1]

    end

    idx = all.(isfinite, binsf.y)

    t_new, sf_new, err_new = binsf.x[idx], binsf.y[idx],  binsf.yerr[idx]
    
    t_br = find_t_break(binsf)
    
    idx_tmax = t_new .<= t_br

    # cheak initial guess of p0 values is under t_max
    check_bounds(log10(p0[2]), t_br)

    p0_bounds = ifelse(isempty(p0), (lb .+ ub) / 2, p0)
    # p0_bounds = (lb .+ ub) / 3 # we have to start inside the bounds
    # p0_bounds[2] = 10 ^ t_br
    wt = 1  ./ err_new[idx_tmax] .^ 2
    fit_bounds = curve_fit(jmodel, 10 .^ t_new[idx_tmax], sf_new[idx_tmax], wt, p0_bounds, lower=lb, upper=ub)

    # fit_bounds = so.curve_fit(model, 10 .^ t_new[idx_tmax], sf_new[idx_tmax], sigma=collect(binsf.yerr[idx][idx_tmax]), p0=collect(p0_bounds))
    # p  = fit_bounds[1]
    
    # p_std_err = [sqrt(fit_bounds[2][i, i]) for i=1: 4]
    p = fit_bounds.param
    p_std_err = stderror(fit_bounds)

    # cf = coef(fit_bounds)
    # ci = confidence_interval(fit_bounds, 0.05)
    # @printf("%f %F %f %F", Inf, Inf, NaN, NaN)
    println(" ")
    @printf "sigma: %.2f +/- %.2f" p[4]/1e-2 p_std_err[4]/1e-2
    println(" ")
    @printf "beta: %.2f +/- %.2f" p[3] p_std_err[3]
    println(" ")
    @printf "tau: %.2f +/- %.2f" p[2]/1e3 p_std_err[2]/1e3
    println(" ")
    @printf "SF: %.2f +/- %.2f" p[1]/1e-2  p_std_err[1]/1e-2

    println(" ")
    println(" ")
    # str = @sprintf("Y = (%.5f +/- %.5f) + (%.5f +/- %.5f) + (%.5f +/- %.5f) + (%.5f +/- %.5f)",
        #   cf[1],diff([ci[1]...])[1]/2, cf[2], diff([ci[2]...])[1]/2, cf[3],diff([ci[3]...])[1]/2, cf[4], diff([ci[4]...])[1]/2)


    # p = (sf=cf[1],sf_err=diff([ci[1]...])[1]/2, tau=cf[2],tau_err= diff([ci[2]...])[1]/2, beta= cf[3], beta_err= diff([ci[3]...])[3]/2,sigma= cf[4],sigma_err = diff([ci[4]...])[1]/2)

    # p = fit_bounds.param


    #TODO
    #Chi-square r = ydata - jmodel(xdata, param), chisq = @. sum((r / sigma) ^ 2).

    return (binsf = binsf, param = p, param_err = p_std_err)

end