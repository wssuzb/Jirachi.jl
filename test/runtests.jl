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

using UncertainData, Random, Peaks, LaTeXStrings, DelimitedFiles, LsqFit, Interpolations, StatsBase, HDF5, Printf, LinearAlgebra, CairoMakie, MakiePublication

include("../src/fitting.jl")

band = ["g", "r", "i", "z"]

for b in band
    lc1 = load_data("./test/data/montano22_n1_g.txt")
    lc1.time = lc1.time * 3600 * 24
    lc2 = load_data("./test/data/montano22_n1_" * b * ".txt")
    lc2.time = lc2.time * 3600 * 24
    select_time(lc1, lc2, "./test/data/montano22_n1_" * b * "_binned.txt")
end

# lc1 = load_data("/Users/suzhenbo/Mylibrary/Projects/lib_julia_external/VariabilityTools/test_data/lc/montano22_n1_i_binned.txt", [1, 3, 4]; band = "i")
# lc2 = load_data("/Users/suzhenbo/Mylibrary/Projects/lib_julia_external/VariabilityTools/test_data/lc/montano22_n1_z_binned.txt", [1, 3, 4]; band = "z")


lc1 = load_data("./test/data/montano22_n1_i_binned.txt", [1, 3, 4]; band = "i")
# lc1.time =  lc1.time* 3600 * 24

lc2 = load_data("./test/data/montano22_n1_z_binned.txt", [1, 3, 4]; band = "z")
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
  
# function log_tick_formatter(values)
#     return map(v -> Makie.UnicodeFun.to_superscript(round(Int64, v)), values)
# end

function mysf()
    # set_theme!(fonts = (; regular = "Times New Roman"))
    # set_theme!(fonts = (; regular = "/Users/suzhenbo/opt/anaconda3/lib/python3.8/site-packages/smplotlib/ttf/AVHersheyComplexMedium.ttf"))

    fig = Figure()  

    ax = Axis(fig[1, 1],
    xlabel=L"$\Delta t$ [sec]",
    ylabel="SF [mag]",
    # xtickformat = ,
    xscale=log10, 
    yscale=log10, 
    xticksize=8,
    yticksize=8,
    # xtickformat=log_tick_formatter,
    # ytickformat=log_tick_formatter,
    # xminorticks = IntervalsBetween(10),
    # yminorticks = IntervalsBetween(10),
    xticksmirrored = true, yticksmirrored = true, 
    # xminorticksvisible = true,
    # xminorgridvisible = false,      
    # yminorticksvisible = true,
    # yminorgridvisible = false,  
    # xminortickalign = 1,
    # yminortickalign = 1,
    # xtickalign = 1,
    # ytickalign = 1,
    # xgridvisible = false,
    # ygridvisible = false
    )

    l1 = lines!(ax, t_fit_1, sf_fit_1, )
    l2 = lines!(ax, t_fit_2, sf_fit_2)
    s1 = scatter!(ax, 10 .^ binsf1.x, binsf1.y)
    s2 = scatter!(ax, 10 .^ binsf2.x, binsf2.y)

    e1 = errorbars!(ax, 10 .^ binsf1.x, binsf1.y, binsf1.yerr, binsf1.yerr, linewidth=0.5, whiskerwidth = 0.2)
    e2 = errorbars!(ax, 10 .^ binsf2.x, binsf2.y, binsf2.yerr, binsf2.yerr,  linewidth=0.5,  whiskerwidth = 0.2; transparency=true)
    # vspan!([10 ^ t_min_1], [10 ^ t_break_1], color = (:blue))

    band!(ax, 10 ^ t_min_1:10^t_break_1, 0.003, 0.1, color= (:red, 0.1) )
    # hspan!(-1.1, -0.9, color = (:blue, 0.5))

    xlims!(80, 3e4)
    # ylims!(4e-3, 2e-1)
    ylims!(5e-3, 0.1)
    

    resize_to_layout!(fig)
    # fig[1, 1] = ax
    return fig
end

llc = Cycle([:color, :linestyle], covary=true)
ssc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
tex_web = merge(theme_web(width=350,
                            colors=MakiePublication.tableau_10(),
                            linestyles=[nothing, :dash, :dash],
                            ishollowmarkers=[true, true, false],
                            markers=[:circle, :diamond, :rtriangle],
                            linecycle=llc,
                            scattercycle=ssc,
                            markerstrokewidth=0.8, 
                            heightwidthratio=0.9,
                            ), 
                            theme_latexfonts())

fig_sf = with_theme(tex_web) do 
    mysf()
end

save("./test/fig/plot_sf.svg", fig_sf, px_per_unit=4)
# with_theme(mysf, theme_web( width=300,
# colors=MakiePublication.tableau_10(),
# linestyles=[nothing, :dash, :dash],
# ishollowmarkers=[true, true, false],
# markers=[:circle, :diamond, :rtriangle],
# linecycle=llc,
# scattercycle=ssc,
# markerstrokewidth=0.8, heightwidthratio=0.8   # font="Times New Roman"
#     )
# )

function myhist!(fig::Figure)
    # set_theme!(fonts = (; regular = "Times New Roman"))
    # set_theme!(fonts = (; regular = "/Users/suzhenbo/opt/anaconda3/lib/python3.8/site-packages/smplotlib/ttf/AVHersheyComplexMedium.ttf"))

  
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
    fig[1, 1] = ax
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
    # yminorticks = IntervalsBetween(10),
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


function myplot()
    fig = Figure(figure_padding=(2, 6, 1, 6))
    ax = Axis(fig, xlabel=L"x", ylabel=L"f(x)")
    x = [1, 2, 3, 4, 5]

    y_A = [1, 2, 1, 2, 1]
    err_A = [0.1, 0.2, 0.1, 0.2, 0.3]
    
    y_B = [2, 3, 2, 3, 2]
    err_B = [0.2, 0.3, 0.2, 0.1, 0.1]
    
    scatter!(ax, x, y_A, label='A')
    errorbars!(ax, x, y_A, err_A, label='A', linewidth=0.3)
    
    scatter!(ax, x, y_B, label='B')
    errorbars!(ax, x, y_B, err_B, label='B', linewidth=0.3)
    fig[1,1] = ax
    fig
end



f = Figure()
ax = Axis(f[1, 1], xscale=log10,yscale=log10, limits = (10, 2e4, 2, 4e3))
# ax.limits=(0.0011,3,0.0001,20)
stephist!(ax, num_pos,  bin=  10 .^ collect(2:0.1:5),linewidth=2, color=:dodgerblue4)
stephist!(ax, num_all, bin=10 .^ collect(2:0.1:5),linewidth=2, color=:black)

# stephist(f[1, 2], data, bins = 20, color = :red, strokewidth = 1, strokecolor = :black)
# stephist(f[2, 1], data, bins = [-5, -2, -1, 0, 1, 2, 5], color = :gray)
# stephist(f[2, 2], data, normalization = :pdf)
f

# stephist(num_all)
# stephist()
stephist(num_cut,xlims=(90, 10 ^ 4.6), xaxis=(:log10), yaxis=(:log10),bin= 10 .^ sf_bin_edges, ylims=(2, 4e3), linewidth=2, color=:dodgerblue4)



tex_web = merge(theme_web(width=350,
colors=MakiePublication.seaborn_dark(),
linestyles=[nothing, :dash, :dash],
ishollowmarkers=[true, true, false],
markers=[:circle, :diamond, :rtriangle],
linecycle=llc,
scattercycle=ssc,
markerstrokewidth=0.8, heightwidthratio=0.9,), theme_latexfonts())

function mysf(sf1::binned_result, sf2::binned_result; )

    with_theme(tex_web) do 
        fig = Figure()  

        ax = Axis(fig[1, 1],
        xlabel=L"$\Delta t$ [sec]",
        ylabel="SF [mag]",
        xscale=log10, 
        yscale=log10, 
        xticksize=8,
        yticksize=8,
        xticksmirrored = true, yticksmirrored = true, 
        )

        l1 = lines!(ax, t_fit_1, sf_fit_1)
        l2 = lines!(ax, t_fit_2, sf_fit_2)
        s1 = scatter!(ax, 10 .^ sf1.x, sf1.y)
        s2 = scatter!(ax, 10 .^ sf2.x, sf2.y)

        e1 = errorbars!(ax, 10 .^ sf1.x, sf1.y, sf1.yerr, sf1.yerr, linewidth=0.3, whiskerwidth = 0.2)
        e2 = errorbars!(ax, 10 .^ sf2.x, sf2.y, sf2.yerr, sf2.yerr,  linewidth=0.3,  whiskerwidth = 0.2; transparency=true)
        # vspan!([10 ^ t_min_1], [10 ^ t_break_1], color = (:blue))

        band!(ax, 10 ^ t_min_1:10^t_break_1, 0.003, 0.1, color= (:red, 0.1) )
        # hspan!(-1.1, -0.9, color = (:blue, 0.5))

        xlims!(80, 3e4)
        # ylims!(4e-3, 2e-1)
        ylims!(5e-3, 0.1)

        resize_to_layout!(fig)
        save("./test/fig/plot_sf.png", fig, px_per_unit=4)
    end    
    # return fig_sf
    
end

mysf()



lc_bin_err(dataset::Vector{T}) where T =  sqrt.(sum(dataset .^ 2)) / length(dataset)




function bin_light_curve(lc::lightcurve; lc_edges::AbstractArray)

    _bin_edges = collect(lc_edges)
    _bin_width = round.([(_bin_edges[i+1] - _bin_edges[i]) / 2 for i=1: lastindex(_bin_edges)-1], digits=2)
    bin_center = round.(_bin_edges[2:end] .- _bin_width, digits=2)

    flux_bin = bin(mean, lc_edges, lc.time, lc.flux)
    err_bin = bin(lc_bin_err, lc_edges, lc.time, lc.err)
    return lightcurve(bin_center, flux_bin, err_bin, lc.band)
end

function get_common(lc1::lightcurve, lc2::lightcurve)
    
    idx = all(!isnan, hcat(lc1.flux, lc2.flux); dims=2) |> vec

    lc1.time, lc1.flux, lc1.err = lc1.time[idx], lc1.flux[idx], lc1.err[idx]
    lc2.time, lc2.flux, lc2.err = lc2.time[idx], lc2.flux[idx], lc2.err[idx]

    # lc1.time[t1_idx] ∩ lc2.time[t2_idx]
    # filter(row -> all(x -> !(x isa Number && isnan(x)), row), t_bin)
    return lc1, lc2
end

lc = lc1
_bin_edges = collect(lc_edges)
_bin_width = round.([(_bin_edges[i+1] - _bin_edges[i]) / 2 for i=1: lastindex(_bin_edges)-1], digits=2)
bin_center = round.(_bin_edges[2:end] .- _bin_width, digits=2)

flux_bin = bin(mean, lc_edges, lc.time, lc.flux)
err_bin = bin(lc_bin_err, lc_edges, lc.time, lc.err)


# function check_quasi_time()

# lc1 = load_data("/Users/suzhenbo/Mylibrary/Projects/archive_data/ngc4395/McHardy22/LightCurve_NotNormalised_mjy_g_8000.dat", [1, 2, 3]; band = "g")
# # lc1.time = round.(lc1.time * 24 * 3600, digits=2)
# lc2 = load_data("/Users/suzhenbo/Mylibrary/Projects/archive_data/ngc4395/McHardy22/LightCurve_NotNormalised_mjy_z_8000.dat", [1, 2, 3]; band = "r")
# # lc2.time = round.(lc2.time * 24 * 3600, digits=2)


# t_cad = mean(diff(lc1.time))
t_cad = 3
t_start =lc1.time[1] - (t_cad / 2)
t_end = lc1.time[end] + (t_cad /2 )
lc_edges = t_start:t_cad:t_end

lc1_bin = bin_light_curve(lc1; lc_edges=lc_edges)
lc2_bin = bin_light_curve(lc2;lc_edges = lc_edges)

lc1_bin, lc2_bin = get_common_lc(lc1_bin, lc2_bin)


fig = Figure(size=(500, 200))
ax = Axis(fig[1, 1])
scatter!(ax, lc1.time, lc1.flux, markersize=1, color=:black)
scatter!(ax, lc2.time, lc2.flux, markersize=1, color=:black)
errorbars!(ax, lc1.time, lc1.flux, lc1.err, color=:black)
errorbars!(ax, lc2.time, lc2.flux, lc2.err, color=:black, linewidth=1)

scatter!(ax, lc1_bin.time, lc1_bin.flux, markersize=1, color=:red)
scatter!(ax, lc2_bin.time, lc2_bin.flux, markersize=1, color=:red)
errorbars!(ax, lc1_bin.time, lc1_bin.flux, lc1_bin.err, color=:red)
errorbars!(ax, lc2_bin.time, lc2_bin.flux, lc2_bin.err, color=:red, linewidth=1)
xlims!(ax, 0, 400)
fig


lc1 = lc1_bin
lc2 = lc2_bin
sf_bin_edges=0:0.1:5
cv_bin_edges=0:0.2:5
nsigma=3
erron=true
nsim=10
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

f = Figure()
ax = Axis(f[1, 1])
lines!(ax, bincv_flux.x, bincv_flux.y)
errorbars!(ax, bincv_flux.x, bincv_flux.y, bincv_flux.yerr)
ylims!(0.5, 0.9)
f



# cv in mag-mag
cv_mag_res = color_variation(lc1, lc2, nsigma, erron, "mag")
cv_mag = cv_mag_res

bincv_mag = binned_color_variation(cv_mag, cv_bin_edges)
  