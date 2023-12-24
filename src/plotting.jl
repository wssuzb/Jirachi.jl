export theme_sf, theme_hist, theme_lc, IntegerTicks


struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)


theme_lc = merge(theme_web(width=1000,
                            colors=MakiePublication.tableau_10(),
                            linestyles=[nothing, :dash, :dash],
                            ishollowmarkers=[true, true, false],
                            markers=[:circle, :diamond, :rtriangle],
                            linecycle=Cycle([:color, :linestyle], covary=true),
                            scattercycle=Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true),
                            markerstrokewidth=0.8, 
                            heightwidthratio=0.3,
                            ), theme_latexfonts())


theme_hist = merge(theme_web(width=350, colors=MakiePublication.tableau_10(),linestyles=[nothing, :dash, :dash], ishollowmarkers=[true, true, false], markers=[:circle, :diamond, :rtriangle], linecycle=Cycle([:color, :linestyle], covary=true), scattercycle=Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true), markerstrokewidth=0.8, heightwidthratio=0.9,), theme_latexfonts())

# theme_llc = merge(theme_web(width=1000,
#                             colors=MakiePublication.COLORS[5],
#                             linestyles=[:dash, :dash],
#                             ishollowmarkers=[true, true],
#                             markers=[:circle, :diamond],
#                             linecycle=Cycle([:color, :linestyle], covary=true),
#                             scattercycle=Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true),
#                             markerstrokewidth=0.8, 
#                             heightwidthratio=0.3,
#                             ), theme_latexfonts())



function plotlc(lc::lightcurve...; label=[], xlim=[], ms=2, lw=0.3, width=1000, hwratio=0.5, ft=20, lc_edges::AbstractArray=[])

    theme_llc = merge(
        theme_web(
        width=width,
        colors=MakiePublication.tol_bright(),
        linestyles=[:dash, :dash],
        ishollowmarkers=[true, false],
        markers=[:circle, :diamond],
        linecycle=Cycle([:color, :linestyle], covary=true),
        scattercycle=Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true),
        markerstrokewidth=0.8, 
        heightwidthratio=hwratio
        ), 
        theme_latexfonts()
    )

    isempty(lc_edges) ? (lc_edges = IntervalsBetween(9)) : lc_edges
    
    isempty(label) ? (label = [string(i) for i=1: length(lc)]) : label
    length(label) == length(lc) ? label :  (label = [string(i) for i=1: length(lc)])

    kwargs = (; xminorticks = lc_edges)

    fig = with_theme(theme_llc) do
        f = Figure(fontsize=ft)

        ax = Axis(f[1, 1],
        xlabel="Time",
        ylabel=L"f",
        xminorgridvisible = true, 
        xticksmirrored = true,
        xminorticksvisible = false, # xgridvisible = true,
        yticksmirrored = true; kwargs...)

        for (i, x) in enumerate(lc)
            sca = scatter!(ax, x.time, x.flux, markersize=ms, label=label[i])
            ebar = errorbars!(ax, x.time, x.flux, x.err, linewidth=lw)
        end

        axislegend(ax, nbanks = 2)

        isempty(xlim) ? println(" ")  : xlims!(ax, xlim[1], xlim[2])

        f
    end

    # savefig("./figure/lc_test.pdf", fig)
    return fig
end



function mysf()
    theme_sf = merge(theme_web(width=350, colors=MakiePublication.tableau_10(),linestyles=[nothing, :dash, :dash], ishollowmarkers=[true, true, false], markers=[:circle, :diamond, :rtriangle], linecycle=Cycle([:color, :linestyle], covary=true), scattercycle=Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true), markerstrokewidth=0.8, heightwidthratio=0.9,), theme_latexfonts())

    fig_sf = with_theme(theme_sf) do 
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

