export theme_sf, theme_hist, theme_lc, IntegerTicks

llc = Cycle([:color, :linestyle], covary=true)
ssc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)

struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

theme_lc = merge(theme_web(width=1000,
                            colors=MakiePublication.tableau_10(),
                            linestyles=[nothing, :dash, :dash],
                            ishollowmarkers=[true, true, false],
                            markers=[:circle, :diamond, :rtriangle],
                            linecycle=llc,
                            scattercycle=ssc,
                            markerstrokewidth=0.8, 
                            heightwidthratio=0.3,
                            ), theme_latexfonts())

theme_sf = merge(theme_web(width=350, colors=MakiePublication.tableau_10(),linestyles=[nothing, :dash, :dash], ishollowmarkers=[true, true, false], markers=[:circle, :diamond, :rtriangle], linecycle=llc, scattercycle=ssc, markerstrokewidth=0.8, heightwidthratio=0.9,), theme_latexfonts())

theme_hist = merge(theme_web(width=350, colors=MakiePublication.tableau_10(),linestyles=[nothing, :dash, :dash], ishollowmarkers=[true, true, false], markers=[:circle, :diamond, :rtriangle], linecycle=llc, scattercycle=ssc, markerstrokewidth=0.8, heightwidthratio=0.9,), theme_latexfonts())



function mysf()
    fig_sf = with_theme(tex_web) do 
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
    end    
    save("./test/fig/plot_sf.svg", fig_sf, px_per_unit=4)
end

