# # Reproduce the figures from the paper and the supplements

figuredir = "../../products/paper_figures/"
tabledir = "../../products/derived_data/"
mkpath(figuredir)
mkpath(tabledir)

# #### Fig.1: Map
include("overview-map.jl")
savefig(figuredir * "figure1.png")

# Fig2 is a photo and sketch, i.e. not produced here.

# #### Fig. 3: Plot of derived quantities
idx_plot = Dict("0808" => 17000:24000, # which indices to plot the hydraulic gradient; exclude parts that only disturb the scale of the y-axis
                "0908" => 8000:28000,
                "1008" => 12400:20750,
                "1108" => 8300:20700,
                "1308" => 13600:15900,
                "2108" => 1:28640);
idx_gaps = Dict("0908" => 11100:11500,
                "1108" => 8800:10000)

R.multiplot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)
savefig(figuredir * "figure3.png")

# #### Fig. 4: Evolution of the cross-sectional area
R.plot_model_outputs(mid_309_265, model_runs)
savefig(figuredir * "figure4.png")

# #### Fig. 5: heat transfer parameters
fig = R.plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured)[1]
savefig(figuredir * "figure5.png")

# ### Figures and Tables from supllementary material
#
# Table S1 was produced manually from elements of :depth in
# ctd309/ctd265 dictionaries, e.g. ctd309["0908"][:depth]

# ### Fig. S1: water pressure and temperature
R.plot_pw_Tw_supp(mid_309_265, pick, ctd309, ctd265, e_p, e_T, idx_plot, idx_gaps)
savefig(figuredir * "figure_s1.png")

# #### Fig. S2: Opening rates
R.plot_opening(mid_309_265, model_runs)
savefig(figuredir * "figure_s2.png")

# #### Fig. S3: Closure rate
R.plot_closure(mid_309_265, model_runs)
savefig(figuredir * "figure_s3.png")

# #### Fig. S4: heat transfer parameters extended, with estimated surface temperature
fig, table = R.plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured; T_surf)
savefig(figuredir * "figure_s4.png")
# save table to text file
write(tabledir * "table_s2.tex", "% Contens of Table S2\n" * table)

# #### Fig. S5: heat transfer parameters of each experiment
R.plot_heat_params_timeresolved([z_eq, tau_eq, tau_w],
                                [L"z_\mathrm{eq}\,\mathrm{(m)}", L"\tau_\mathrm{eq}\,(\mathrm{°C})", L"\tau_{w}\,(\mathrm{°C})"],
                                [nothing, nothing, tauw_measured])
savefig(figuredir * "figure_s5.png")

# Also save table S3 (produced in size_evolution_models.jl):
write(tabledir * "table_s3.tex", "% Contens of Table S3\n" * table_s3)

# CSV files
using DelimitedFiles
for (fl, daystr) in [("SQ_09_08_2020.csv", "0908"), ("SQ_21_08_2020.csv", "2108")]
    writedlm(tabledir * fl, ["Time" "S [m^2]" "std_S" "Q [m^3/s]" "std_Q" "dphi/dz [Pa/m]" "std dphi/dz" "f []" "std f";
                             mid_309_265[daystr][:t_inj] mean.(mid_309_265[daystr][:S]) std.(mid_309_265[daystr][:S]) mean.(ctd265[daystr][:Q]) std.(ctd265[daystr][:Q]) mean.(mid_309_265[daystr][:dphi_dz]) std.(mid_309_265[daystr][:dphi_dz]) mean.(mid_309_265[daystr][:f]) std.(mid_309_265[daystr][:f])],
             ", ")
end
