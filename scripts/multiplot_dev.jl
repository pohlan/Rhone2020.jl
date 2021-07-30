
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

    function draw_breakpoints(axis, kwargs, dx, dy, spin) # draw breakpoints to interrupt x-axis
        if spin == "upper right" # draw breakpoints at right spines of panels
            axis.plot((1-dx,1+dx), (-dy,+dy); kwargs...)
            axis.plot((1-dx,1+dx),(1-dy,1+dy); kwargs...)
        elseif spin == "left" # draw breakpoints at left spines of panels
            axis.plot((-dx,+dx), (1-dy,1+dy); kwargs...)
            axis.plot((-dx,+dx), (-dy,+dy); kwargs...)
        end
    end



    # width ratios
    widths = []
    for date in dates
        append!(widths, length(idx_plot[date]))
    end
    w_ratios = widths ./ mean(widths)
    grid_dict = Dict(:width_ratios => widths,
                     :wspace => 0.05, # horizontal space between panels
                     :hspace => 0.0)  # vertical space between panels

    f, axes = plt.subplots(5, 5, sharey="row", sharex="col", gridspec_kw=grid_dict)

    break_point_length = 0.017 # needs adjustment because figures not equally high and wide
    ds = break_point_length ./ w_ratios
    kwargs = Dict(:transform => nothing,
                  :color     => "k",
                  :clip_on   => false)

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

        # other properties
        props = Dict(:pw => [],
                     :dphi_dz => [],
                     :Q => mid_309_265[date][:Q][pick[date]], # m^3 / s !!!!
                     :v => mid_309_265[date][:v][pick[date]],
                     :S => mid_309_265[date][:S][pick[date]])
        nprops = length(props)

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

        dx = ds[nd]
        dy = break_point_length * 2
        ax = 0
        for (row, key) in enumerate([:pw, :dphi_dz, :Q, :v, :S])
            ax = length(props)*(nd-1) + row

            if row >= 3 # first two rows are for pressure and hydraulic gradient, manually
                axes[ax].errorbar(t_inj, mean.(props[key]), yerr=std.(props[key]), fmt="k_")
            end
            if nd == 1 # left panel
                # make intermediate spines dashed lines
                axes[ax].spines["right"].set_linestyle((0, (3, 10)))
                # remove axis ticks
                axes[ax].tick_params(right=false)
                # draw break point lines
                kwargs[:transform] = axes[ax].transAxes
                #axes[ax].plot((1-dx,1+dx), (-dy,+dy); kwargs...)
                #axes[ax].plot((1-dx,1+dx),(1-dy,1+dy); kwargs...)
            elseif nd == length(dates) # right panel
                axes[ax].spines["left"].set_linestyle((0, (3, 10)))
                axes[ax].tick_params(left=false)
                #kwargs[:transform] = axes[ax].transAxes
                #axes[ax].plot((-dx,+dx), (1-dy,1+dy); kwargs...)
                #axes[ax].plot((-dx,+dx), (-dy,+dy); kwargs...)
            else # panels in between
                axes[ax].spines["left"].set_linestyle((0, (3, 10)))
                axes[ax].spines["right"].set_linestyle((0, (3, 10)))
                axes[ax].tick_params(left=false, right=false)
                #kwargs[:transform] = axes[ax].transAxes
                #axes[ax].plot((1-dx,1+dx), (-dy,+dy); kwargs...)
                #axes[ax].plot((1-dx,1+dx),(1-dy,1+dy); kwargs...)
                #axes[ax].plot((-dx,+dx), (1-dy,1+dy); kwargs...)
                #axes[ax].plot((-dx,+dx), (-dy,+dy); kwargs...)
            end

            # draw breakpoints
            kwargs[:transform] = axes[ax].transAxes
            if row == 1
                if nd !== 1 # draw in upper left corner
                    axes[ax].plot((-dx,+dx), (1-dy,1+dy); kwargs...)
                elseif nd !== length(dates) # draw in upper right corner
                    axes[ax].plot((1-dx,1+dx),(1-dy,1+dy); kwargs...)
                end
            elseif row == nprops
                if nd !== 1 # draw in lower left corner
                    axes[ax].plot((-dx,+dx), (-dy,+dy); kwargs...)
                elseif nd !== length(dates) # draw in lower right corner
                    axes[ax].plot((1-dx,1+dx),(-dy,+dy); kwargs...)
                end
            end
        end
        format_xaxis(axes[ax])
    end



# ax1.set_xlim(0, xend)

gcf()

end

multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)

