module Jirachi

# Write your package code here.

using Pkg
using HDF5
using Interpolations
using Random
using StatsBase
using UncertainData
using DelimitedFiles
using LaTeXStrings
using Peaks
using Printf
using LsqFit
using LinearAlgebra
# using PyCall

# ENV["PYTHON"] = "/usr/bin/python"
# Pkg.build("PyCall")
# ENV["PYTHON"] = ""
# Pkg.build("PyCall")
include("utils.jl")
include("structure_function.jl")
include("color_variation.jl")
include("fractional_variability.jl")
include("fitting.jl")
include("run.jl")

# const so = PyNULL()

# function __init__()
#     copy!(so, pyimport_conda("scipy.optimize", "scipy"))
# end


export par, lightcurve, cv, sf, binned_result,percentile_16_50_84, load_data, save_data, lc_bootstrapped, find_nearest, select_time

export color_variation, binned_color_variation, structure_function, binned_structure_function, err_bootstraped, flux2mag

export jmodel, model, fitsf, fitsf_mcmc, find_t_min, find_t_break

export runall


end
