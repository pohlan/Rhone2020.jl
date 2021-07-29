
function multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)
    f, axes = plt.subplots(1, 5, sharey=true)
    subplots_adjust(# left = 0.1, -> position
                    wspace=0.05, # horizontal space
                    hspace=0.1) # vertical space

    # plot
    dates = ["0808", "0908", "1008", "1108", "1308"]
    d = .015
    kwargs = Dict(:transform => nothing,
                  :color     => "k",
                  :clip_on   => false)
    for (nd, date) in enumerate(dates)
        axes[nd].plot(ctd309[date][:t], ctd309[date][:press])

        if nd == 1 # left panel
            # make intermediate spines dashed lines
            axes[nd].spines["right"].set_linestyle((0, (3, 10)))
            # remove axis ticks
            axes[nd].tick_params(right=false)
            # draw break point lines
            kwargs[:transform] = axes[nd].transAxes
            axes[nd].plot((1-d,1+d), (-d,+d); kwargs...)
            axes[nd].plot((1-d,1+d),(1-d,1+d); kwargs...)
        elseif nd == length(dates) # right panel
            axes[nd].spines["left"].set_linestyle((0, (3, 10)))
            axes[nd].tick_params(left=false)
            kwargs[:transform] = axes[nd].transAxes
            axes[nd].plot((-d,+d), (1-d,1+d); kwargs...)
            axes[nd].plot((-d,+d), (-d,+d); kwargs...)
        else # panels in between
            axes[nd].spines["left"].set_linestyle((0, (3, 10)))
            axes[nd].spines["right"].set_linestyle((0, (3, 10)))
            axes[nd].tick_params(left=false, right=false)
            kwargs[:transform] = axes[nd].transAxes
            axes[nd].plot((1-d,1+d), (-d,+d); kwargs...)
            axes[nd].plot((1-d,1+d),(1-d,1+d); kwargs...)
            axes[nd].plot((-d,+d), (1-d,1+d); kwargs...)
            axes[nd].plot((-d,+d), (-d,+d); kwargs...)
        end
    end



# ax1.set_xlim(0, xend)

gcf()

end

multi_plot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)

