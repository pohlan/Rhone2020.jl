# # Reproduce the figures from the paper and the supplements

using PyPlot
figuredir = "../../products/paper_figures/"

# #### Fig.1: Map
include("../scripts/overview-map.jl")
savefig(figuredir * "map.png")

# #### Fig. 3: Plot of derived quantities
idx_plot = Dict("0808" => 14900:28640, # which indices to plot the hydraulic gradient; exclude parts that only disturb the scale of the y-axis
                "0908" => 7200:28640,
                "1008" => 12500:20750,
                "1108" => 8000:22300,
                "1308" => 13330:15975,
                "2108" => 1:28640);
R.multiplot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot)
savefig(figuredir * "masterplot.png")

# #### Fig. 4: Evolution of the cross-sectional area
R.plot_model_outputs(mid_309_265, model_runs)
savefig(figuredir * "size_evolution.png")

# #### Fig. 5: heat transfer parameters
R.plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured)
savefig(figuredir * "heat-transfer-numbers.png")

# ### Figures from supllementary material
#
# #### Fig. S1: Closure rate
R.plot_opening_closure(mid_309_265, model_runs)
savefig(figuredir * "closure_rates.png")

# #### Fig. S2: heat transfer parameters extended, with estimated surface temperature
R.plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured; T_surf)
savefig(figuredir * "heat-transfer-numbers-extended.png")
