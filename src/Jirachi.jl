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
using CairoMakie
using MakiePublication

include("utils.jl")
include("structure_function.jl")
include("color_variation.jl")
include("fractional_variability.jl")
include("fitting.jl")
include("showinfo.jl")
include("run.jl")
include("plotting.jl")
include("generateDRW.jl")

export lightcurve, cv, sf, binned_result, percentile_16_50_84, load_data, save_data, lc_bootstrapped, find_nearest, select_time, get_common_lc, bin_light_curve, bin_lc_edges, remove_lc_outlier, remove_lc_nan, hcatlc, mergelc

export color_variation, binned_color_variation, structure_function, binned_structure_function, err_bootstraped, flux2mag

export jmodel, model, fitsf, fitsf_mcmc, find_t_min, find_t_break

export runall

export stochastic_process

export IntegerTicks, theme_lc, theme_sf, plotlc, plotsf, plotcv


end
