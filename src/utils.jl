export par, lightcurve, cv, sf, binned_result,percentile_16_50_84, load_data, save_data, lc_bootstrapped, find_nearest, select_time, bin_light_curve, get_common_lc
# """
#     par(med, low, hig)
# structure for loading statistical results.
# """
# # Examples
# ```jldoctest
# julia> par(med, low, hig)
# ```
struct par
    med::Float64
    low::Float64
    hig::Float64
end


# # Examples
# ```jldoctest
# julia> lightcurve(time, flux, err)
# ```
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



# # Examples
# ```jldoctest
# julia> cv(tau, color)
# ```
"""
    cv(tau, color)

structure for loading color variation temp results.

"""
struct cv
    tau::Vector{Float64}
    color::Vector{Float64}
end

# # Examples
# ```jldoctest
# julia> sf(tau, sf)
# ```
"""
    sf(tau, sf)

structure for loading structure function results.

"""
struct sf
    tau::Vector{Float64}
    sf::Vector{Float64}
end


@kwdef mutable struct binned_result
    x::Vector{Float64}
    xerr::Vector{Float64}
    y::Vector{Float64}
    yerr::Vector{Float64}
    yerr_::Vector{Float64}
end
binned_result(x, xerr, y, yerr) = binned_result(x, xerr, y, yerr, [])

# function save_struct(data::Any, fi_np::String)
#     """
#     Convert struct to dictionary.
#     """
#     struct_to_dict(s) = Dict(key => getfield(s, key) for key in propertynames(s))

#     test = struct_to_dict(data)
#     open(fi_np, "w") do io
#         JSON3.write(io, test)
#     end

#     println("File have save to: " * "\"" * fi_np * "\"")
#     println(" ")
#     println("The data set: ")
#     dump(data)
#     println(" ")
#     return test
# end

# Examples
# ```jldoctest
# julia> percentile_16_50_84([1, 2, 3])
# ```
"""
    percentile_16_50_84(array)

return 16%, 50% and 84% values with given array.
"""
function percentile_16_50_84(x::T) where {T}
    med = percentile(x, 50)
    low = med - percentile(x, 16)
    hig = percentile(x, 84) - med
    return par(med, low, hig)
end



function uniquecount(data::T) where{T}
    unique_array = unique(data)
    counts = Dict(unique_array .=> 0)
    for (i, c) in enumerate(data)
       counts[c] += 1     
    end
    return (keys=collect(keys(counts)), vals=collect(values(counts)))
 end



#  # Examples
#  ```jldoctest
#  julia> lc_bootstrapped(data)
#  ```

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
    
    # if mode == "rss"
    
    #     return lightcurve(t_rss[idx_sort], y_rss[idx_sort], e_rss[idx_sort], data.band)
    
    # elseif mode == "both"
    
    #     return lightcurve(t_rss[idx_sort], y_fr_rss[idx_sort], e_rss[idx_sort], data.band)
    
    # else
    #     return lightcurve(
    #         data.time,
    #         [rand(Normal(data.flux[i], data.err[i]), 1)[1] for i=1: lastindex(data.flux)],
    #         data.err,
    #         data.band
    #     )
    # end

end

# """
#     idx = findmat(x->x > 0, dt_all)
# """

function findmat(f, A::AbstractMatrix)
    m,n = size(A)
    out = []
    for i in 1:m, j in 1:n
      f(A[i,j]) && push!(out,(i,j))
    end
    return out
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
# load_data(fi_np::String; usecols=[1, 2, 3]) = load_data(fi_np, usecols)


function save_data(arr1::Vector{Float64}, arr2::Vector{Float64}, arr3::Vector{Float64}; fi_np::String="./test_data/color_variability.txt")
    res = []
    for i=1: lastindex(arr1)
        push!(res, [arr1[i], arr2[i], arr3[i]])
    end

    open(fi_np, "w") do io
        writedlm(io, [map(x->x[1], res) map(x->x[2], res) map(x->x[3], res)], ' ')
    end
end

function save_data(lc1::lightcurve, lc2::lightcurve; fi_np::String="./test_data/color_variability.txt")
    res = []
    for i=1: lastindex(lc1.time)
        push!(res, [lc1.time[i], lc1.flux[i], lc1.err[i], lc2.time[i], lc2.flux[i], lc2.err[i]])
    end

    open(fi_np, "w") do io
        writedlm(io, [map(x->x[1], res) map(x->x[2], res) map(x->x[3], res) map(x->x[4], res) map(x->x[5], res) map(x->x[6], res)], ' ')
    end
end
# save_data(data::sf; fi_np::String="./test_data/color_variability.txt") = save_data(data.tau, data.sf; fi_np)
# save_data(data::cv; fi_np::String="./test_data/color_variability.txt") = save_data(data.tau, data.color; fi_np)


# function save_data(arr1::Vector{Float64}, arr2::Vector{Float64}, arr3::Vector{Float64}, arr4::Vector{Float64}; fi_np::String="./test_data/binned_color_variability.txt")
#     res = []
#     for i=1: lastindex(arr1)
#         push!(res, [arr1[i], arr2[i], arr3[i], arr4[i]])
#     end

#     open(fi_np, "w") do io
#         writedlm(io, [map(x->x[1], res) map(x->x[2], res) map(x->x[3], res) map(x->x[4], res)], ' ')
#     end
# end

# save_data(data::binned_result;  fi_np::String="./test_data/binned_color_variability.txt") = save_data(data.x, data.xerr, data.y, data.yerr; fi_np)



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

    data = []
    for i=1: lastindex(time)
        push!(data, [time[i], time_[i], flux[i], err[i]])
    end

    open(fi_np, "w") do io
        writedlm(io, [map(x->x[1], data) map(x->x[2], data) map(x->x[3], data) map(x->x[4], data)], ' ')
    end

    return time, time_, flux, err
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
# bin_light_curve(lc, lc_edges) = bin_light_curve(lc::lightcurve; lc_edges::AbstractArray)

function get_common_lc(lc1::lightcurve, lc2::lightcurve)
    # find the index of values that are finite.
    idx = all(!isnan, hcat(lc1.flux, lc2.flux); dims=2) |> vec
    
    lc1.time, lc1.flux, lc1.err = lc1.time[idx], lc1.flux[idx], lc1.err[idx]
    lc2.time, lc2.flux, lc2.err = lc2.time[idx], lc2.flux[idx], lc2.err[idx]
    
    # lc1.time[t1_idx] ∩ lc2.time[t2_idx]
    # filter(row -> all(x -> !(x isa Number && isnan(x)), row), t_bin)
    return lc1, lc2
end

# function get_common_lc(arr1::T, arr2::T) where{T}
#     # find the common values of two given array: arr1, arr2
#     common_t = intersect(arr1, arr2)
   
#     # find the index of the common values in: arr1, arr2, respectively.
#     ind1 = indexin(common_t, arr1)
#     ind2 = indexin(common_t, arr2)
#     return ind1, ind2
# end
