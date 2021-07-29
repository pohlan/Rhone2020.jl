
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

    # width ratios
    widths = []
    for date in dates
        append!(widths, length(idx_plot[date]))
    end
    w_ratios = widths ./ mean(widths)
    grid_dict = Dict(:width_ratios => widths,
                     :wspace => 0.05,
                     :hspace => 0.1)

    f, axes = plt.subplots(5, 5, sharey="row", gridspec_kw=grid_dict)

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
        presstop = ctd309[date][:press][i] ./ (rhow * g) # mH20
        pressbot = ctd265[date][:press][i] ./ (rhow * g) # mH20

        # hydraulic gradient
        dz = mid_309_265[date][:dz]
        dpress = Particles.(n_partcl, Normal.(ctd265[date][:press][i], e_p*300)) .-
                 Particles.(n_partcl, Normal.(ctd309[date][:press][i], e_p*100))
        dp_smooth = Particles.(n_partcl, Normal.(boxcar(mean.(dpress), window_am15), (boxcar(std.(dpress), window_am15)))) # smooth the pressure difference
        dphi_dz = dp_smooth ./ dz .- rhow*g # compute hydraulic gradient, Pa/m
        dphi_dz = dphi_dz ./ (rhow*g) # convert it to mH2O/m

        # other properties
        props = Dict(:pw => [],
                     :dphi_dz => [],
                     :Q => mid_309_265[date][:Q][pick[date]], # m^3 / s !!!!
                     :v => mid_309_265[date][:v][pick[date]],
                     :S => mid_309_265[date][:S][pick[date]])

        # remove data points where we moved CTDs up and down or blocked the water inflow
        if haskey(idx_gaps, date)
            presstop[idx_gaps[date]] .= NaN
            pressbot[idx_gaps[date]] .= NaN
            dphi_dz[idx_gaps[date]] .= NaN
        end

        # plot
        # pressure
        axes[length(props)*(nd-1)+1].plot(ctd309[date][:t][i], presstop, "k", linewidth=0.5)
        axes[length(props)*(nd-1)+1].plot(ctd309[date][:t][i], pressbot, linewidth=0.5)
        # hydraulic gradient
        axes[length(props)*(nd-1)+2].plot(ctd309[date][:t][i], mean.(dphi_dz), "k", linewidth=0.5)
        axes[length(props)*(nd-1)+2].fill_between(ctd309[date][:t][i], mean.(dphi_dz) .+ std.(dphi_dz), mean.(dphi_dz) .- std.(dphi_dz), color="grey")

        dx = ds[nd]
        d  = break_point_length
        for (row, key) in enumerate([:pw, :dphi_dz, :Q, :v, :S])
            ax = length(props)*(nd-1) + row

            if row >= 3 # first two rows are for pressure and hydraulic gradient, manually
                axes[ax].errorbar(t_inj, mean.(props[key]), yerr=std.(props[key]), fmt="k+")
            end
            if nd == 1 # left panel
                # make intermediate spines dashed lines
                axes[ax].spines["right"].set_linestyle((0, (3, 10)))
                # remove axis ticks
                axes[ax].tick_params(right=false)
                # draw break point lines
                kwargs[:transform] = axes[ax].transAxes
                axes[ax].plot((1-dx,1+dx), (-d,+d); kwargs...)
                axes[ax].plot((1-dx,1+dx),(1-d,1+d); kwargs...)
            elseif nd == length(dates) # right panel
                axes[ax].spines["left"].set_linestyle((0, (3, 10)))
                axes[ax].tick_params(left=false)
                kwargs[:transform] = axes[ax].transAxes
                axes[ax].plot((-dx,+dx), (1-d,1+d); kwargs...)
                axes[ax].plot((-dx,+dx), (-d,+d); kwargs...)
            else # panels in between
                axes[ax].spines["left"].set_linestyle((0, (3, 10)))
                axes[ax].spines["right"].set_linestyle((0, (3, 10)))
                axes[ax].tick_params(left=false, right=false)
                kwargs[:transform] = axes[ax].transAxes
                axes[ax].plot((1-dx,1+dx), (-d,+d); kwargs...)
                axes[ax].plot((1-dx,1+dx),(1-d,1+d); kwargs...)
                axes[ax].plot((-dx,+dx), (1-d,1+d); kwargs...)
                axes[ax].plot((-dx,+dx), (-d,+d); kwargs...)
            end
        end
    end



# ax1.set_xlim(0, xend)

gcf()

end

multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)

