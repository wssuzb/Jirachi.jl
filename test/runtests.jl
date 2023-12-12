# using Jirachi
# using Test

# @testset "Jirachi.jl" begin
    # Write your tests here.
# end
include("../src/Jirachi.jl")
include("../src/utils.jl")
include("../src/structure_function.jl")
include("../src/color_variation.jl")
include("../src/fractional_variability.jl")
include("../src/run.jl")
include("../src/fitting.jl")
using UncertainData, Random, Peaks, LaTeXStrings, DelimitedFiles, LsqFit, Interpolations, StatsBase, HDF5, Printf, LinearAlgebra, CairoMakie

band = ["g", "r", "i", "z"]

for b in band
    lc1 = load_data("test/data/montano22_n1_g.txt")
    lc1.time = lc1.time * 3600 * 24
    lc2 = load_data("test/data/montano22_n1_" * b * ".txt")
    lc2.time = lc2.time * 3600 * 24
    select_time(lc1, lc2, "test/data/montano22_n1_" * b * "_binned.txt")
end

lc1 = load_data("/Users/suzhenbo/Mylibrary/Projects/lib_julia_external/VariabilityTools/test_data/lc/montano22_n1_i_binned.txt", [1, 3, 4]; band = "i")
lc2 = load_data("/Users/suzhenbo/Mylibrary/Projects/lib_julia_external/VariabilityTools/test_data/lc/montano22_n1_z_binned.txt", [1, 3, 4]; band = "z")


lc1 = load_data("test/data/montano22_n1_i_binned.txt", [1, 3, 4]; band = "i")
# lc1.time =  lc1.time* 3600 * 24

lc2 = load_data("test/data/montano22_n1_z_binned.txt", [1, 3, 4]; band = "z")
# lc2.time =  lc2.time * 3600 * 24

# lc1 == lc1_


sf_bin_edges=0:0.1:5
cv_bin_edges=0:0.2:5
nsigma=3
erron=true
nsim=100
mode = "both"
fi_np::String="./test/run_all.h5"
lower_bounds = [0, 0, 0, 0.001]
upper_bounds = [10, 2e4, 2, 0.1]
p0 = [1, 1e3, 1, 0.05]


fit_sf1 = fitsf_mcmc(lc1; nsim=nsim, lb = lower_bounds , ub = upper_bounds, sf_bin_edges=sf_bin_edges, p0=p0, mode = mode)
fit_sf2 = fitsf_mcmc(lc2; nsim=nsim, lb = lower_bounds, ub = upper_bounds, sf_bin_edges=sf_bin_edges, p0=p0, mode = mode)

binsf1, binsf2 = fit_sf1.binsf, fit_sf2.binsf
    
par_1, par_1_err = fit_sf1.param, fit_sf1.param_err
par_2, par_2_err = fit_sf2.param, fit_sf2.param_err

t_break_1 = find_t_break(binsf1)
t_break_2 = find_t_break(binsf2)

itp1 = find_t_min(binsf1, par_1)
itp2 = find_t_min(binsf2, par_2)

t_min_1, sf_min_1 = itp1.t_min, itp1.sf_min
t_min_2, sf_min_2 = itp2.t_min, itp2.sf_min

t_fit_1, sf_fit_1 = itp1.t_fit, itp1.sf_fit
t_fit_2, sf_fit_2 = itp2.t_fit, itp2.sf_fit

# cv in flux-Flux
cv_flux_res = color_variation(lc1, lc2, nsigma, erron, "flux"; debug=true)

cv_flux = cv_flux_res.cv

num_all = cv_flux_res.num_all
num_cut = cv_flux_res.num_cut
num_pos = cv_flux_res.num_pos

bincv_flux = binned_color_variation(cv_flux, cv_bin_edges)

# cv in mag-mag
cv_mag_res = color_variation(lc1, lc2, nsigma, erron, "mag")
cv_mag = cv_mag_res

bincv_mag = binned_color_variation(cv_mag, cv_bin_edges)
  

fig = Figure()
function mysf(fig::Figure)
    # set_theme!(fonts = (; regular = "Times New Roman"))
    set_theme!(fonts = (; regular = "/Users/suzhenbo/opt/anaconda3/lib/python3.8/site-packages/smplotlib/ttf/AVHersheyComplexMedium.ttf"))

  
    ax = Axis(fig[1, 1], aspect = 1, xscale=log10, yscale=log10, 

    xminorticks = IntervalsBetween(10),
    yminorticks = IntervalsBetween(10),
    xticksmirrored = true, yticksmirrored = true, 
    xminorticksvisible = true,
    xminorgridvisible = false,      
    yminorticksvisible = true,
    yminorgridvisible = false,  
    xminortickalign = 1,
    yminortickalign = 1,
    xtickalign = 1,
    ytickalign = 1,
    xgridvisible = false,
    ygridvisible = false)

    l1 = lines!(ax, t_fit_1, sf_fit_1, )
    l2 = lines!(ax, t_fit_2, sf_fit_2)
    s1 = scatter!(ax, 10 .^ binsf1.x, binsf1.y)
    s2 = scatter!(ax, 10 .^ binsf2.x, binsf2.y)

    e1 = errorbars!(ax, 10 .^ binsf1.x, binsf1.y, binsf1.yerr, binsf1.yerr, whiskerwidth = 0.2)
    e2 = errorbars!(ax, 10 .^ binsf2.x, binsf2.y, binsf2.yerr, binsf2.yerr, whiskerwidth = 0.2; transparency=true)
    # vspan!([10 ^ t_min_1], [10 ^ t_break_1], color = (:blue))

    band!(ax, 10 ^ t_min_1:10^t_break_1, 0.006, 0.1, color= (:red, 0.1) )
    # hspan!(-1.1, -0.9, color = (:blue, 0.5))

    xlims!(80, 3e4)
    ylims!(0.006, 0.1)

    # resize_to_layout!(fig)
    # fig[1, 1] = ax
    return fig
end

function myhist!(fig::Figure)
    # set_theme!(fonts = (; regular = "Times New Roman"))
    set_theme!(fonts = (; regular = "/Users/suzhenbo/opt/anaconda3/lib/python3.8/site-packages/smplotlib/ttf/AVHersheyComplexMedium.ttf"))

  
    ax = Axis(fig[1, 2], aspect = 1, xscale=log10, yscale=log10, 

    xminorticks = IntervalsBetween(10),
    yminorticks = IntervalsBetween(10),
    xticksmirrored = true, yticksmirrored = true, 
    xminorticksvisible = true,
    xminorgridvisible = false,      
    yminorticksvisible = true,
    yminorgridvisible = false,  
    xminortickalign = 1,
    yminortickalign = 1,
    xtickalign = 1,
    ytickalign = 1,
    xgridvisible = false,
    ygridvisible = false)

    l1 = lines!(ax, t_fit_1, sf_fit_1, )
    l2 = lines!(ax, t_fit_2, sf_fit_2)
    s1 = scatter!(ax, 10 .^ binsf1.x, binsf1.y)
    s2 = scatter!(ax, 10 .^ binsf2.x, binsf2.y)

    e1 = errorbars!(ax, 10 .^ binsf1.x, binsf1.y, binsf1.yerr, binsf1.yerr, whiskerwidth = 0.2)
    e2 = errorbars!(ax, 10 .^ binsf2.x, binsf2.y, binsf2.yerr, binsf2.yerr, whiskerwidth = 0.2; transparency=true)
    # vspan!([10 ^ t_min_1], [10 ^ t_break_1], color = (:blue))

    band!(ax, 10 ^ t_min_1:10^t_break_1, 0.006, 0.1, color= (:red, 0.1) )
    # hspan!(-1.1, -0.9, color = (:blue, 0.5))

    xlims!(80, 3e4)
    ylims!(0.006, 0.1)

    # resize_to_layout!(fig)
    # fig[1, 1] = ax
    return fig
end
function mycv_flux!(f::Figure)

    ax2 = Axis(f[2, 1], aspect = 1, xscale=log10, 
    
    xminorticks = IntervalsBetween(10),
    yminorticks = IntervalsBetween(10),
    xticksmirrored = true, yticksmirrored = true, 
    xminorticksvisible = true,
    xminorgridvisible = false,      
    yminorticksvisible = true,
    yminorgridvisible = false,  
    xminortickalign = 1,
    yminortickalign = 1,
    xtickalign = 1,
    ytickalign = 1,
    xgridvisible = false,
    ygridvisible = false)

    s1 = scatter!(ax2, 10 .^ bincv_flux.x, bincv_flux.y)

    e1 = errorbars!(ax2, 10 .^ bincv_flux.x, bincv_flux.y, bincv_flux.yerr, bincv_flux.yerr, whiskerwidth = 0.2)
    e2 = errorbars!(ax2, 10 .^ bincv_flux.x, bincv_flux.y, bincv_flux.xerr, bincv_flux.xerr, whiskerwidth = 0.2, direction = :x)
    band!(ax2, 10 ^ t_min_1:10^t_break_1, 0, 2, color= (:red, 0.1) )
    # hspan!(-1.1, -0.9, color = (:blue, 0.5))

    xlims!(80, 3e4)
    ylims!(0.5, 1.1)
 
    return f
end

function mycv_mag!(f::Figure)

    ax3 = Axis(f[2, 2], aspect = 1, xscale=log10, 
    
    xminorticks = IntervalsBetween(10),
    yminorticks = IntervalsBetween(10),
    xticksmirrored = true, yticksmirrored = true, 
    xminorticksvisible = true,
    xminorgridvisible = false,      
    yminorticksvisible = true,
    yminorgridvisible = false,  
    xminortickalign = 1,
    yminortickalign = 1,
    xtickalign = 1,
    ytickalign = 1,
    xgridvisible = false,
    ygridvisible = false)

    s1 = scatter!(ax3, 10 .^ bincv_mag.x, bincv_mag.y)

    e1 = errorbars!(ax3, 10 .^ bincv_mag.x, bincv_mag.y, bincv_mag.yerr, bincv_mag.yerr, whiskerwidth = 0.2)
    e2 = errorbars!(ax3, 10 .^ bincv_mag.x, bincv_mag.y, bincv_mag.xerr, bincv_mag.xerr, whiskerwidth = 0.2, direction = :x)
    band!(ax3, 10 ^ t_min_1:10^t_break_1, 0., 2, color= (:red, 0.1) )
    # hspan!(-1.1, -0.9, color = (:blue, 0.5))

    xlims!(80, 3e4)
    ylims!(0.7, 1.3)
 
    return f
end

fig = mysf(Figure(size=(500, 500)))
myhist!(fig)
mycv_flux!(fig)
mycv_mag!(fig)
resize_to_layout!(fig)