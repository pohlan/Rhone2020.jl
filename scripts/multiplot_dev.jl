#dat = "1308"
#figure()
#plot(ctd265[dat][:press] .- ctd309[dat][:press][1:length(ctd265[dat][:press])])
#figure()
#plot(ctd265[dat][:press][idx_plot[dat]] .- ctd309[dat][:press][idx_plot[dat]])
#figure()
#plot(ctd265[dat][:t], ctd265[dat][:press] .- ctd309[dat][:press][1:length(ctd265[dat][:press])])


### remove later ###
boxcar = R.boxcar
####################

idx_plot = Dict("0808" => 15000:24000, # which indices to plot the hydraulic gradient; exclude parts that only disturb the scale of the y-axis
                "0908" => 8000:28000,  # [8000:19100, 12300:27700]
                "1008" => 12400:20750,
                "1108" => 8300:20700,
                "1308" => 13600:15900,
                "2108" => 1:28640);
idx_gaps = Dict("0908" => 11100:11500,
                "1108" => 8800:10000)




function multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)
    dates = ["0808", "0908", "1008", "1108", "1308", "2108"]
    props = [:dphi_dz, :Q, :v, :S, :f, :n_manning, :Re]
    panel_labs = [L"\bf{a}", L"\bf{b}", L"\bf{c}", L"\bf{d}", L"\bf{e}", L"\bf{f}", L"\bf{g}"]
    ylabels = [#L"$p_w\,\mathrm{(mH_2O)}$",
               L"$\partial\phi/\partial z$ $\mathrm{(mH_2O\,m^{-1})}$",
               L"$Q$ $\mathrm{(m^3\,s^{-1})}$", # or in l/s ???
               L"$v$ $\mathrm{(m\,s^{-1})}$",
               L"$S$ $\mathrm{(m^2)}$",
               L"$f$",
               L"$n'$ $(\mathrm{s\,m^{-1/3}})$",
               L"$Re$"]
    nprops = length(props)

    # font properties
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    fs = 20 # font size for figure numbering

    # time axis format
    # for the Locator it is necessary to define this function,
    # otherwise the xticks from previous plots are removed
    # https://stackoverflow.com/questions/55010236/set-major-locator-removes-x-ticks-and-labels-from-previous-subplots

    function format_xaxis(ax, loc)
        fmt = matplotlib.dates.DateFormatter("%H:%M")
        loc = matplotlib.dates.HourLocator(;loc...)
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(loc)
    end

    # width ratios for left subfigure
    widths = []
    for date in dates[1:end-1]
        append!(widths, length(idx_plot[date]))
    end
    w_ratios = widths ./ mean(widths)
    w_space = 0.04
    grid_dict_left = Dict(:width_ratios => widths,
                          :wspace => w_space, # horizontal space between panels
                          :hspace => 0.)  # vertical space between panels

    # grid_dict for right subfigure
    grid_dict_right = Dict(:hspace => 0.)

    # style of the spine line at x-axis breakpoints
    spine_line = (0, (2, 4))

    # for drawing breakpoints
    break_point_length = 0.017
    dxs = break_point_length ./ w_ratios
    kwargs = Dict(:transform => nothing,
                  :color     => "k",
                  :linewidth => 1.0,
                  :clip_on   => false
                  )

    # draw subplots
    fig = plt.figure(figsize=(20,15))
    subfigs = fig.subfigures(1, 2, width_ratios=[sum(widths)+mean(widths)*4*w_space, length(idx_plot["2108"])], wspace = 0.0)
    axesleft  = subfigs[1].subplots(length(props), 5, sharey="row", sharex="col", gridspec_kw=grid_dict_left)
    axesright = subfigs[2].subplots(length(props), 1, sharex="col", gridspec_kw=grid_dict_right)

    for (nd, date) in enumerate(dates)
        i = idx_plot[date]

        # time axes injections
        t_inj = mid_309_265[date][:t_inj][pick[date]]

        # pressure
        # time windows over which to average pressures, in seconds
        if date == "2108"
            window = 900
        else
            window = 60
        end
        presstop = boxcar(ctd309[date][:press][i], window)
        pressbot = boxcar(ctd265[date][:press][i], window)

        # hydraulic gradient
        dz = mid_309_265[date][:dz]
        dpress = Particles.(n_partcl, Normal.(pressbot, e_p*300)) .-
                 Particles.(n_partcl, Normal.(presstop, e_p*100))
        dp_smooth = Particles.(n_partcl, Normal.(mean.(dpress), std.(dpress))) # smooth the pressure difference
        dphi_dz = dp_smooth ./ dz .- rhow*g # compute hydraulic gradient, Pa/m

         # conversion to mH2O
        presstop = presstop ./ (rhow * g)
        pressbot = pressbot ./ (rhow * g)
        dphi_dz = dphi_dz ./ (rhow*g)

        # remove data points where we moved CTDs up and down or blocked the water inflow
        if haskey(idx_gaps, date)
            presstop[idx_gaps[date]] .= NaN
            pressbot[idx_gaps[date]] .= NaN
            dphi_dz[idx_gaps[date]] .= NaN
        end

        for (row, (prop, ylab)) in enumerate(zip(props, ylabels))
            if date == "2108"
                ax = axesright[row]
                nd = 1
            else
                ax = axesleft[nprops*(nd-1) + row]
            end
            #if row == 1
            #    ax.plot(ctd309[date][:t][i], presstop, "k", linewidth=0.5)
            #    ax.plot(ctd309[date][:t][i], pressbot, linewidth=0.5)
            if row == 1 # first row is for hydraulic gradient, manually
                ax.plot(ctd309[date][:t][i], mean.(dphi_dz), "k", linewidth=0.5)
                ax.fill_between(ctd309[date][:t][i], mean.(dphi_dz) .+ std.(dphi_dz), mean.(dphi_dz) .- std.(dphi_dz), color="grey")
            else
                data = mid_309_265[date][prop][pick[date]]
                ax.errorbar(t_inj, mean.(data), yerr=std.(data), fmt="k+", markersize=11)
            end

            # logarithmic y-scale for some parameters
            if any(prop .== [:f, :n_manning, :Re])
                ax.set_yscale("log")
            end

            # title
            if row == 1
                ax.set_title(date[1:2] * "-Aug")
            end

            # y-label and panel label
            if date == "0808"
                ax.set_ylabel(ylab, labelpad=25., wrap=true) #, rotation="horizontal", va="center")
                text(0.1, 0.6, panel_labs[row], fontsize=fs, transform=ax.transAxes, ha="left", va="top")
            end

            # make intermediate spines dashed lines and remove axis ticks
            if date !== "2108"
                dx = dxs[nd]
                dy = break_point_length * 3 # needs the factor because plots longer than high, factor chosen randomly
                if nd !== length(dates)-1 # adjust right side of panel
                    ax.spines["right"].set_linestyle(spine_line)
                    ax.tick_params(right=false)
                end
                if nd !== 1 # adjust left side of panel
                    ax.spines["left"].set_linestyle(spine_line)
                    ax.tick_params(which="both", left=false) # both -> remove also minor ticks in logscale plots, default is major
                end
            end
            if row !== nprops
                ax.tick_params(bottom=false)
            end

            # draw breakpoints
            kwargs[:transform] = ax.transAxes
            #if row == 1
            #    if nd !== 1 # draw in upper left corner
            #        ax.plot((-dx,+dx), (1-dy,1+dy); kwargs...)
            #    end
            #    if nd !== length(dates)-1 # draw in upper right corner
            #        ax.plot((1-dx,1+dx),(1-dy,1+dy); kwargs...)
            #    end
            # elseif row == nprops
            if row == nprops && date !== "2108"
                if nd !== 1 # draw in lower left corner
                    ax.plot((-dx,+dx), (-dy,+dy); kwargs...)
                end
                if nd !== length(dates)-1 # draw in lower right corner
                    ax.plot((1-dx,1+dx),(-dy,+dy); kwargs...)
                end
            end
            # adjust format of time axis
            #if date == "1008"
            #    format_xaxis(ax, Dict(:interval=>1)) # interval of x-tick labels 1 hour
            #if date == "1108"
            #    format_xaxis(ax, Dict(:byhour=>(12, 13)))
            if date == "1308"
                format_xaxis(ax, Dict(:byhour=>14))
            else
                format_xaxis(ax, Dict(:interval=>2))
            end
            ax.tick_params("x", pad=10., labelrotation=30.)
            ax.margins(y=0.2) # so that ylabels don't overlap each other
        end
    end

    subfigs_pos = Dict(:left => 0.18, :right => 0.97, :bottom => 0.12, :top => 0.92)
    subfigs[1].subplots_adjust(;subfigs_pos...)
    subfigs[2].subplots_adjust(;subfigs_pos...)
    subfigs[1].suptitle(L"\bf{AM15}", y=0.98)
    subfigs[2].suptitle(L"\bf{AM13}", y=0.98)
    subfigs[1].supxlabel("Time", y=0.02)

#gcf()
savefig("figure3.png")

end

multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)

