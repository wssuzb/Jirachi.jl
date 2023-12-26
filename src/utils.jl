export lightcurve, cv, sf, binned_result,percentile_16_50_84, load_data, save_data, lc_bootstrapped, find_nearest, select_time, bin_light_curve, get_common_lc, bin_lc_edges, remove_lc_outlier, remove_lc_nan


@kwdef mutable struct parameters
    # for light curve bins
    t_start::Float64
    t_end::Float64
    t_binsize::Float64

    # for structure function
    sf_binsize::Float64
    sf_start::Float64
    sf_end::Float64

    # for structure function fitting
    mode::String # bootstrapped lc
    nsim::Int64
    lower_bounds::Vector{Any}
    upper_bounds::Vector{Any}
    p0::Vector{Any}

    # for color variation
    cv_binsize::Float64
    cv_start::Float64
    cv_end::Float64
    nsigma::Int64
    erroron::Bool
end


bin_lc_edges(binsize, t_start, t_end) = (t_start - (binsize / 2)):binsize:(t_end + (binsize / 2))

# range(start=0, stop=1, step=2)

"""
    lightcurve(time, flux, err)

structure for loading light curve.
"""
@kwdef mutable struct lightcurve
    time::Vector{Float64}
    flux::Vector{Float64}
    err::Vector{Float64}
    band
end
lightcurve(time, flux, err) = lightcurve(time, flux, err, [])


struct cv
    tau::Vector{Float64}
    color::Vector{Float64}
end

struct sf
    tau::Vector{Float64}
    sf::Vector{Float64}
end


@kwdef mutable struct binned_result
    x::Vector{Float64}
    xerr::Vector{Float64}
    y::Vector{Float64}
    yerr::Vector{Float64}
end


function percentile_16_50_84(x::T) where {T}
    med = percentile(x, 50)
    low = med - percentile(x, 16)
    hig = percentile(x, 84) - med
    return (med=med, low=low, hig=hig)
end

function uniquecount(data::T) where{T}
    unique_array = unique(data)
    counts = Dict(unique_array .=> 0)
    for (i, c) in enumerate(data)
       counts[c] += 1     
    end
    return (keys=collect(keys(counts)), vals=collect(values(counts)))
 end


"""
    lc_bootstrapped(data::lightcurve; seed=1, mode="both")
 
 
 return bootstrapped light curves, input data must be with format of ::lc.
 mode: 
    - 1, rss only
    - 2, fr only
    - 0, fr/rss
"""
function lc_bootstrapped(data::lightcurve; seed=1, mode="both")

    Random.seed!(seed)

    idx = collect(rand(1:length(data.time), length(data.time)))
    unique_array = uniquecount(idx)
    idx_unique, counts = unique_array.keys, unique_array.vals
    
    t_rss = data.time[idx_unique]
    y_rss = data.flux[idx_unique]
    e_rss = @. data.err[idx_unique] / sqrt(counts)

    y_fr_rss = [rand(Normal(y_rss[i], e_rss[i]), 1)[1] for i=1: lastindex(y_rss)]
    idx_sort = sortperm(t_rss)

    mode == "rss" && return lightcurve(t_rss[idx_sort], y_rss[idx_sort], e_rss[idx_sort], data.band)

    mode == "both" && return lightcurve(t_rss[idx_sort], y_fr_rss[idx_sort], e_rss[idx_sort], data.band)

    lightcurve(
        data.time,
        [rand(Normal(data.flux[i], data.err[i]), 1)[1] for i=1: lastindex(data.flux)],
        data.err,
        data.band
    )
    
end


"""
    load_data(fi_np::String, usecols=[1, 2, 3]; band=[])

Loading light curves from txt files.
"""
function load_data(fi_np::String, usecols=[1, 2, 3]; band=[])
    tmp = readdlm(fi_np)
    data = lightcurve(tmp[:, usecols[1]], tmp[:, usecols[2]], tmp[:, usecols[3]], band)
    return data
end


function hcatlc(lc::lightcurve)
    T = typeof(lc)
    names = fieldnames(T)
    return reduce(hcat, (getfield(lc, names[i]) for i=1: length(names)-1)) 
end


function save_data(lc::lightcurve...; fi_np="./test.txt")
    open(fi_np,"w") do io
        writedlm(io, reduce(hcat, (hcatlc(l) for l in lc)))
    end
end
# write(io, "# Light curves are saved as time, flux, err\n")

function save_data(arr::Vector{Float64}...; fi_np="./test.txt")
    open(fi_np,"w") do io
        writedlm(io, hcat(arr...))
    end
end


find_nearest(arr, val) = argmin(abs.(arr .- val))
    
function select_time(data1::lightcurve, data2::lightcurve, fi_np::String)

    time, time_, flux, err = [], [], [], []
    idx_list = []

    for i=1: lastindex(data1.time)
        idx = find_nearest(data2.time, data1.time[i])
        idx_list = push!(idx_list, idx)
        
        time = push!(time, round(data1.time[i], digits=2))
        time_ = push!(time_, round(data2.time[idx],digits=2))
        flux = push!(flux, data2.flux[idx])
        err = push!(err, data2.err[idx])
    end

    if length(unique(idx_list)) !=  length(idx_list)
        # print(data2.band)
        print(length(unique(idx_list)))
        print(length(idx_list))
        print("WARNING")
    end

    data = hcat(time, time_, flux, err)
    
    open(fi_np, "w") do io
        writedlm(io, data)
    end

    return data
end



lc_bin_err(err::Vector{T}) where T =  sqrt.(sum(err .^ 2)) / length(err)

function bin_light_curve(lc::lightcurve; lc_edges::AbstractArray)
    _bin_edges = collect(lc_edges)
    _bin_width = round.([(_bin_edges[i+1] - _bin_edges[i]) / 2 for i=1: lastindex(_bin_edges)-1], digits=2)
    bin_center = round.(_bin_edges[2:end] .- _bin_width, digits=2)

    flux_bin = bin(mean, lc_edges, lc.time, lc.flux)
    err_bin = bin(lc_bin_err, lc_edges, lc.time, lc.err)
    
    return lightcurve(bin_center, flux_bin, err_bin, lc.band)
end

function remove_lc_nan(lc::lightcurve)
    idx = all(!isnan, lc.flux; dims=2) |> vec
    return lightcurve(lc.time[idx], lc.flux[idx], lc.err[idx], lc.band)
end


# function get_common_lc(lc1::lightcurve, lc2::lightcurve)
#     # find the index of values that are finite.
#     idx = all(!isnan, hcat(lc1.flux, lc2.flux); dims=2) |> vec
    
#     return lightcurve(lc1.time[idx], lc1.flux[idx], lc1.err[idx], lc1.band), lightcurve(lc2.time[idx], lc2.flux[idx], lc2.err[idx], lc2.band)
# end


function get_common_lc(lc1::lightcurve, lc2::lightcurve)
    # find the common values of two given array: arr1, arr2
    common_t = intersect(lc1.time, lc2.time)
   
    # find the index of the common values in: arr1, arr2, respectively.
    ind1 = indexin(common_t, lc1.time)
    ind2 = indexin(common_t, lc2.time)

    return lightcurve(lc1.time[ind1], lc1.flux[ind1], lc1.err[ind1], lc1.band), lightcurve(lc2.time[ind2], lc2.flux[ind2], lc2.err[ind2], lc2.band)
end

function remove_lc_outlier(lc::lightcurve; cut=2, criteria="err")
    if criteria == "err"
        mean_err, std_err = mean(lc.err), std(lc.err)
        idx = @. abs(lc.err - mean_err) <= (cut * std_err)
    else
        mean_flux, std_flux = mean(lc.flux), std(lc.flux)
        idx = @. abs(lc.flux - mean_flux) <= (cut * std_flux)
    end
    
    return lightcurve(lc.time[idx], lc.flux[idx], lc.err[idx], lc.band)
end