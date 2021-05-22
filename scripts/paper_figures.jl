# # Reproduce the figures from the paper and the supplements

using PyPlot
figuredir = "../../products/paper_figures/"
tabledir = "../../products/derived_data/"
mkpath(figuredir)
mkpath(tabledir)

# #### Fig.1: Map
include("../scripts/overview-map.jl")
savefig(figuredir * "figure1.png")

# Fig2 is a photo and sketch, i.e. not produced here.

# #### Fig. 3: Plot of derived quantities
idx_plot = Dict("0808" => 14900:28640, # which indices to plot the hydraulic gradient; exclude parts that only disturb the scale of the y-axis
                "0908" => 7200:28640,
                "1008" => 12500:20750,
                "1108" => 8000:22300,
                "1308" => 13330:15975,
                "2108" => 1:28640);
R.multiplot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot)
savefig(figuredir * "figure3.png")

# #### Fig. 4: Evolution of the cross-sectional area
R.plot_model_outputs(mid_309_265, model_runs)
savefig(figuredir * "figure4.png")

# #### Fig. 5: heat transfer parameters
fig = R.plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured)[1]
savefig(figuredir * "figure5.png")

# ### Figures from supllementary material
#
# #### Fig. S1: Closure rate
R.plot_opening_closure(mid_309_265, model_runs)
savefig(figuredir * "figure_s1.png")

# #### Fig. S2: heat transfer parameters extended, with estimated surface temperature
fig, table = R.plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured; T_surf)
savefig(figuredir * "figure_s2.png")
# save tale to text file
write(tabledir * "table_s1.tex", "% Contens of Table S1\n" * table)

# Also save table S2 (produced in size_evolution_models.jl):
write(tabledir * "table_s2.tex", "% Contens of Table S2\n" * table_s2)

# CSV files
using DelimitedFiles
for (fl, daystr) in [("SQ_09_08_2020.csv", "0908"), ("SQ_21_08_2020.csv", "2108")]
    writedlm(tabledir * fl, ["Time" "S [m^2]" "std_S" "Q [m^3/s]" "std_Q" "dphi/dz [Pa/m]" "std dphi/dz" "f []" "std f";
                             mid_309_265[daystr][:t_inj] mean.(mid_309_265[daystr][:S]) std.(mid_309_265[daystr][:S]) mean.(ctd265[daystr][:Q]) std.(ctd265[daystr][:Q]) mean.(mid_309_265[daystr][:dphi_dz]) std.(mid_309_265[daystr][:dphi_dz]) mean.(mid_309_265[daystr][:f]) std.(mid_309_265[daystr][:f])],
             ", ")
end