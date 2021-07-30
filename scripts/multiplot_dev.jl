
### remove later ###
boxcar = R.boxcar
####################

idx_plot = Dict("0808" => 14900:26600, # which indices to plot the hydraulic gradient; exclude parts that only disturb the scale of the y-axis
                "0908" => 8000:28000,  # [8000:19100, 12300:27700]
                "1008" => 12400:20750,
                "1108" => 8000:25600,
                "1308" => 13330:15975,
                "2108" => 1:28640);
idx_gaps = Dict("0908" => 11100:11500)




function multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)
    dates = ["0808", "0908", "1008", "1108", "1308"]
    props = [:pw, :dphi_dz, :Q, :v, :S, :f, :n_manning, :Re]
    ylabels = [L"$p_w\,\mathrm{(mH_2O)}$",
               L"$\partial \phi /\partial z\,\mathrm{(mH_2O\,m^{-1})}$",
               L"$Q\,\mathrm{(m^3\,s^{-1})}$", # or in l/s ???
               L"$v\,\mathrm{(m\,s^{-1})}$",
               L"$S\,\mathrm{(m^2)}$",
               L"$f$",
               L"$n'\,(\mathrm{s\,m^{-1/3}})$",
               L"$Re$"]
    nprops = length(props)

    # time windows over which to average pressures, in seconds
    window_am15 = 60  # AM15
    window_am13 = 900 # AM13

    # time axis format
    # for the Locator it is necessary to define this function,
    # otherwise the xticks from previous plots are removed
    # https://stackoverflow.com/questions/55010236/set-major-locator-removes-x-ticks-and-labels-from-previous-subplots

    function format_xaxis(axis)
        fmt = matplotlib.dates.DateFormatter("%H:%M")
        loc = matplotlib.dates.HourLocator(interval=2)
        axis.xaxis.set_major_formatter(fmt)
        axis.xaxis.set_major_locator(loc)
    end

    # width ratios
    widths = []
    for date in dates
        append!(widths, length(idx_plot[date]))
    end
    w_ratios = widths ./ mean(widths)
    grid_dict = Dict(:width_ratios => widths,
                     :wspace => 0.02, # horizontal space between panels
                     :hspace => 0.1)  # vertical space between panels

    # style of the spine line at x-axis breakpoints
    spine_line = (0, (2, 4))

    # for drawing breakpoints
    break_point_length = 0.017
    dxs = break_point_length ./ w_ratios
    kwargs = Dict(:transform => nothing,
                  :color     => "k",
                  :clip_on   => false)

    # draw subplots
    f, axes = plt.subplots(8, 5, sharey="row", sharex="col", gridspec_kw=grid_dict)

    for (nd, date) in enumerate(dates)
        i = idx_plot[date]

        # time axes injections
        t_inj = mid_309_265[date][:t_inj][pick[date]]

        # pressure
        presstop = boxcar(ctd309[date][:press][i], window_am15)
        pressbot = boxcar(ctd265[date][:press][i], window_am15)

        # hydraulic gradient
        dz = mid_309_265[date][:dz]
        #dpress = Particles.(n_partcl, Normal.(pressbot, e_p*300)) .-
        #         Particles.(n_partcl, Normal.(presstop, e_p*100))
        #dp_smooth = Particles.(n_partcl, Normal.(mean.(dpress), std.(dpress))) # smooth the pressure difference
        #dphi_dz = dp_smooth ./ dz .- rhow*g # compute hydraulic gradient, Pa/m

         # conversion to mH2O
        presstop = presstop ./ (rhow * g)
        pressbot = pressbot ./ (rhow * g)
        #dphi_dz = dphi_dz ./ (rhow*g)

        # remove data points where we moved CTDs up and down or blocked the water inflow
        if haskey(idx_gaps, date)
            presstop[idx_gaps[date]] .= NaN
            pressbot[idx_gaps[date]] .= NaN
            #dphi_dz[idx_gaps[date]] .= NaN
        end

        # plot
        # pressure
        axes[nprops*(nd-1)+1].plot(ctd309[date][:t][i], presstop, "k", linewidth=0.5)
        axes[nprops*(nd-1)+1].plot(ctd309[date][:t][i], pressbot, linewidth=0.5)
        # hydraulic gradient
        #axes[nprops*(nd-1)+2].plot(ctd309[date][:t][i], mean.(dphi_dz), "k", linewidth=0.5)
        #axes[nprops*(nd-1)+2].fill_between(ctd309[date][:t][i], mean.(dphi_dz) .+ std.(dphi_dz), mean.(dphi_dz) .- std.(dphi_dz), color="grey")

        dx = dxs[nd]
        dy = break_point_length * 3 # needs the factor because plots longer than high, factor chosen randomly

        for (row, (prop, ylab)) in enumerate(zip(props, ylabels))
            ax = nprops*(nd-1) + row
            if row >= 3 # first two rows are for pressure and hydraulic gradient, manually
                data = mid_309_265[date][prop][pick[date]]
                axes[ax].errorbar(t_inj, mean.(data), yerr=std.(data), fmt="k_")
            end

            # logarithmic y-scale for some parameters
            if any(prop .== [:f, :n_manning, :Re])
                axes[ax].set_yscale("log")
            end

            # y-label
            if nd == 1
                axes[ax].set_ylabel(ylab)
            end

            # make intermediate spines dashed lines and remove axis ticks
            if nd !== length(dates) # adjust right side of panel
                axes[ax].spines["right"].set_linestyle(spine_line)
                axes[ax].tick_params(right=false)
            end
            if nd !== 1 # adjust left side of panel
                axes[ax].spines["left"].set_linestyle(spine_line)
                axes[ax].tick_params(left=false)
            end
            if row !== nprops
                axes[ax].tick_params(bottom=false)
            end

            # draw breakpoints
            kwargs[:transform] = axes[ax].transAxes
            #if row == 1
            #    if nd !== 1 # draw in upper left corner
            #        axes[ax].plot((-dx,+dx), (1-dy,1+dy); kwargs...)
            #    end
            #    if nd !== length(dates) # draw in upper right corner
            #        axes[ax].plot((1-dx,1+dx),(1-dy,1+dy); kwargs...)
            #    end
            # elseif row == nprops
            if row == nprops
                if nd !== 1 # draw in lower left corner
                    axes[ax].plot((-dx,+dx), (-dy,+dy); kwargs...)
                end
                if nd !== length(dates) # draw in lower right corner
                    axes[ax].plot((1-dx,1+dx),(-dy,+dy); kwargs...)
                end
            end
            format_xaxis(axes[ax]) # adjust format of time axis
        end
    end

gcf()

end

multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)

