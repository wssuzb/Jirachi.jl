export structure_function, binned_structure_function, err_bootstraped

function structure_function(time, flux, mode=:flux)
    
    mode ∉ [:mag, :flux] && throw(DomainError(mode, "mode must be :flux or :mag!!!"))
    # method ∉ [:mean, :iqr] && throw(DomainError(mode, "SF method must be :mean or :iqr!!!"))
    
    mode ==:flux ? flux = -2.5 * log10.(flux) : flux = flux

    time_2d = repeat(time', length(time))
    dt_all = time_2d - transpose(time_2d)

    x_flux_2d = repeat(flux', length(flux))
    dx_flux_all = x_flux_2d - transpose(x_flux_2d)
    
    idx = dt_all .> 0

    dt_all = dt_all[idx]
    dx_flux_all = dx_flux_all[idx]

    int0 = sortperm(dt_all)
    dt_all_sort = dt_all[int0]

    flux_sort = dx_flux_all[int0]
    # method ==:mean ? flux_sort = sqrt.(pi / 2 * dx_flux_all[int0] .^ 2) : 
    
    return sf(dt_all_sort, flux_sort)
end

function err_bootstraped(x, num::Int64=500)
    res = []
    # num = 500
    for i in range(1, num, step=1)
        idx = collect(rand(1:length(x), length(x)))
        idx_unique = unique(idx)
        push!(res, median([x[d] for d in idx_unique]))
    end
    return std(res)
end

function binned_structure_function(data::sf, bin_edges::AbstractArray=1:0.1:5, method::Function=mean)
    
    # typeof(method) !<: Function && 
    
    method ∉ [mean, iqr, std] && throw(DomainError(method, "SF method must be mean OR iqr OR std!!!"))
    
    tau = log10.(data.tau)
    sf = abs.(data.sf)
    
    bin_edges = bin_edges

    bin_all = bin(bin_edges, tau, sf)
    bin_yerr = zeros(length(bin_all))
    
    for i=1: lastindex(bin_all)
        isempty(bin_all[i]) ? (bin_yerr[i] = 0) : (bin_yerr[i] = std(bin_all[i]) / sqrt(length(bin_all[i])))
    end

    bin_value = bin(method, bin_edges, tau, sf)
    _bin_edges = collect(bin_edges)
    bin_width = round.([(_bin_edges[i+1] - _bin_edges[i]) / 2 for i=1: lastindex(_bin_edges)-1], digits=2)
    bin_center = round.(_bin_edges[2:end] .- bin_width, digits=2)
    
    method == mean ? res = binned_result(bin_center, bin_width, sqrt.(pi/2 * bin_value .^ 2), bin_yerr) :
    method == iqr ? res = binned_result(bin_center, bin_width, 0.741 * bin_value, bin_yerr) :
    method == std ? res = binned_result(bin_center, bin_width, bin_value, bin_yerr) : throw(DomainError(method, "SF method must be :mean or :iqr!!!"))

    return res
end