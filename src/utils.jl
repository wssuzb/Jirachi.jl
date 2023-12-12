export par, lc, cv, sf, binned_result,percentile_16_50_84, load_data, save_data, lc_bootstrapped, find_nearest, select_time
"""
par(med, low, hig)

structure for loading statistical results.

# Examples
```jldoctest
julia> par(med, low, hig)
julia> par.med, par.low, par.hig
```
"""
struct par
    med::Float64
    low::Float64
    hig::Float64
end

"""
lc(time, flux, err)

structure for loading light curve.

# Examples
```jldoctest
julia> res = lc(time, flux, err)
julia> res.time
julia> res.flux
julia> res.err
```
"""
@kwdef mutable struct lc
    time::Vector{Float64}
    flux::Vector{Float64}
    err::Vector{Float64}
    band
end
lc(time, flux, err) = lc(time, flux, err, [])


"""
cv(tau, color)

structure for loading color variation results.

# Examples
```jldoctest
julia> res = cv(tau, color)
julia> res.tau
julia> res.color
```
"""
struct cv
    tau::Vector{Float64}
    color::Vector{Float64}
end

"""
sf(tau, sf)

structure for loading structure function results.

# Examples
```jldoctest
julia> res = sf(tau, sf)
julia> res.tau
julia> res.sf
```
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

"""
percentile_16_50_84(array)

return 16%, 50% and 84% values with given array.

# Examples
```jldoctest
julia> x = [1, 2, 3]
julia> percentile_16_50_84(x)
```
"""
function percentile_16_50_84(x::T) where {T}
    med = percentile(x, 50)
    low = med - percentile(x, 16)
    hig = percentile(x, 84) - med
    return par(med, low, hig)
end



function uniquecount(data)
    unique_array = unique(data)
    counts = Dict(unique_array .=> 0)
    for (i, c) in enumerate(data)
       counts[c] += 1     
    end
    return (keys=collect(keys(counts)), vals=collect(values(counts)))
 end



 """
 lc_bootstrapped(data)
 
 mode: 
    - 1, rss only
    - 2, fr only
    - 0, fr/rss

 return bootstrapped light curves, input data must be with format of ::lc.
 
 # Examples
 ```jldoctest
 julia> _lc = lc_bootstrapped(data)
 ```
 """ 

function lc_bootstrapped(data::lc; seed=1, mode="both")

    Random.seed!(seed)

    idx = collect(rand(1:length(data.time), length(data.time)))
    unique_array = uniquecount(idx)
    idx_unique, counts = unique_array.keys, unique_array.vals
    
    t_rss = data.time[idx_unique]
    y_rss = data.flux[idx_unique]
    e_rss = @. data.err[idx_unique] / sqrt(counts)

    y_fr_rss = [rand(Normal(y_rss[i], e_rss[i]), 1)[1] for i=1: lastindex(y_rss)]
    idx_sort = sortperm(t_rss)

    if mode == "rss"
    
        return lc(t_rss[idx_sort], y_rss[idx_sort], e_rss[idx_sort], data.band)
    
    elseif mode == "both"
    
        return lc(t_rss[idx_sort], y_fr_rss[idx_sort], e_rss[idx_sort], data.band)
    
    else
        return lc(
            data.time,
            [rand(Normal(data.flux[i], data.err[i]), 1)[1] for i=1: lastindex(data.flux)],
            data.err,
            data.band
        )
    end

end

"""
    idx = findmat(x->x > 0, dt_all)
"""

function findmat(f, A::AbstractMatrix)
    m,n = size(A)
    out = []
    for i in 1:m, j in 1:n
      f(A[i,j]) && push!(out,(i,j))
    end
    return out
end 


"""
Loading light curves from txt files.

# Examples
```jldoctest
julia> load_data(fi_np, usecols)
julia> percentile_16_50_84(x)
```
"""
function load_data(fi_np::String, usecols=[1, 2, 3]; band=[])
    tmp = readdlm(fi_np)
    data = lc(tmp[:, usecols[1]], tmp[:, usecols[2]], tmp[:, usecols[3]], band)
    return data
end
# load_data(fi_np::String; usecols=[1, 2, 3]) = load_data(fi_np, usecols)


function save_data(arr1::Vector{Float64}, arr2::Vector{Float64}; fi_np::String="./test_data/color_variability.txt")
    res = []
    for i=1: lastindex(arr1)
        push!(res, [arr1[i], arr2[i]])
    end

    open(fi_np, "w") do io
        writedlm(io, [map(x->x[1], res) map(x->x[2], res)], ' ')
    end
end

# save_data(data::sf; fi_np::String="./test_data/color_variability.txt") = save_data(data.tau, data.sf; fi_np)
# save_data(data::cv; fi_np::String="./test_data/color_variability.txt") = save_data(data.tau, data.color; fi_np)


function save_data(arr1::Vector{Float64}, arr2::Vector{Float64}, arr3::Vector{Float64}, arr4::Vector{Float64}; fi_np::String="./test_data/binned_color_variability.txt")
    res = []
    for i=1: lastindex(arr1)
        push!(res, [arr1[i], arr2[i], arr3[i], arr4[i]])
    end

    open(fi_np, "w") do io
        writedlm(io, [map(x->x[1], res) map(x->x[2], res) map(x->x[3], res) map(x->x[4], res)], ' ')
    end
end

# save_data(data::binned_result;  fi_np::String="./test_data/binned_color_variability.txt") = save_data(data.x, data.xerr, data.y, data.yerr; fi_np)



function find_nearest(arr, val)
    idx = argmin(abs.(arr .- val))
    return idx
end

function select_time(data1::lc, data2::lc, fi_np::String)

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

