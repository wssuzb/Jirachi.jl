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
using UncertainData, Random, Peaks, LaTeXStrings, DelimitedFiles, LsqFit, Interpolations, StatsBase, HDF5, Printf, LinearAlgebra, CairoMakie, MakiePublication

band = ["g", "r", "i", "z"]
for b in band
    lc1 = load_data("test/data/montano22_n1_g.txt")
    lc2 = load_data("test/data/montano22_n1_" * b * ".txt")
    lc1.time = lc1.time * 24 * 3600
    lc2.time = lc2.time * 24 * 3600

    select_time(lc1, lc2, "test/data/montano22_n1_" * b * "_binned.txt")
end

lc1 = load_data("test/data/montano22_n1_i_binned.txt", [1, 3, 4]; band = "i")
lc2 = load_data("test/data/montano22_n1_z_binned.txt", [1, 3, 4]; band = "z")

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

linec = Cycle([:color, :linestyle], covary=true)
scatterc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)



function myplot()
    # set_theme!(fonts = (; regular = "Times New Roman"))
    set_theme!(fonts = (; regular = "/Users/suzhenbo/opt/anaconda3/lib/python3.8/site-packages/smplotlib/ttf/AVHersheyComplexMedium.ttf"))

    fig = Figure()
    ax = Axis(fig[1, 1],aspect = 1, xscale=log10, yscale=log10, 
    
    xminorticks = IntervalsBetween(5),
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

    l1 = lines!(ax, t_fit_1, sf_fit_1)
    l2 = lines!(ax, t_fit_2, sf_fit_2)
    s1 = scatter!(ax, 10 .^ binsf1.x, binsf1.y)
    s2 = scatter!(ax, 10 .^ binsf2.x, binsf2.y)

    e1 = errorbars!(ax, 10 .^ binsf1.x, binsf1.y, binsf1.yerr, binsf1.yerr, whiskerwidth = 0.2)
    e2 = errorbars!(ax, 10 .^ binsf2.x, binsf2.y, binsf2.yerr, binsf2.yerr, whiskerwidth = 0.2)
    # vspan!([10 ^ t_min_1], [10 ^ t_break_1], color = (:blue))

    band!(10 ^ t_min_1:10^t_break_1, 0.006, 0.1, color= (:blue, 0.1) )
    # hspan!(-1.1, -0.9, color = (:blue, 0.5))

    xlims!(80, 3e4)
    ylims!(0.006, 0.1)

    resize_to_layout!(fig)
    # fig[1, 1] = ax
    # return fig
end

myplot()

with_theme(myplot,
theme_web(
    width=400, 
    colors=MakiePublication.tableau_10(),
    linestyles=[:dash, :dash],
    ishollowmarkers=[true, true],
    markers=[:circle, :diamond],
    linecycle=linec,
    scattercycle=scatterc,
    markerstrokewidth=1, heightwidthratio=0.9)
)


save("test/fig/plot_sf.pdf", fig, px_per_unit=4)
# t_break_1