# # Reproduce the figures from the paper and the supplements

figuredir = "../../products/paper_figures/"
tabledir = "../../products/derived_data/"
mkpath(figuredir)
mkpath(tabledir)

# #### Fig.1: Map
include("overview-map.jl") # also saves Fig. S1

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

# Fig. S1: Aerial image, is procuded in "overview-map.jl" (see above)

# ### Fig. S2: water pressure and temperature
R.plot_pw_Tw_supp(mid_309_265, pick, ctd309, ctd265, e_p, e_T, idx_plot, idx_gaps)
savefig(figuredir * "figure_s2.png")

# #### Fig. S3: Opening rates
R.plot_opening(mid_309_265, model_runs)
savefig(figuredir * "figure_s3.png")

# #### Fig. S4: Closure rate
R.plot_closure(mid_309_265, model_runs)
savefig(figuredir * "figure_s4.png")

# #### Fig. S5: heat transfer parameters extended, with estimated surface temperature
fig, table = R.plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured; T_surf)
savefig(figuredir * "figure_s5.png")
# save table to text file
write(tabledir * "table_s2.tex", "% Contens of Table S2\n" * table)

# #### Fig. S6: heat transfer parameters of each experiment
dtau = Dict()
for k in keys(tau_w)
    dtau[k] = Dict()
    for corr in keys(tau_w[k])
        dtau[k][corr] = tau_w[k][corr] .- tau_eq[k][corr]
    end
end
R.plot_heat_params_timeresolved([Nus, z_eq, tau_eq, tau_w, dtau, T_surf],
                                [L"Nu", L"z_\mathrm{eq}\,\mathrm{(m)}", L"\tau_\mathrm{eq}\,(\mathrm{°C})", L"\tau_{w}\,(\mathrm{°C})", L"\tau_{w}-\tau_\mathrm{eq}\,(\mathrm{°C})", L"T_{surf}\,(\mathrm{°C})"],
                                [nothing, nothing, nothing, tauw_measured, nothing, nothing])
savefig(figuredir * "figure_s6.png")

# Also save table S3 (produced in size_evolution_models.jl):
write(tabledir * "table_s3.tex", "% Contens of Table S3\n" * table_s3)

# CSV files
using DelimitedFiles
for (fl, daystr) in [("SQ_09_08_2020.csv", "0908"), ("SQ_21_08_2020.csv", "2108")]
    writedlm(tabledir * fl, ["Time" "S [m^2]" "std_S" "Q [m^3/s]" "std_Q" "dphi/dz [Pa/m]" "std dphi/dz" "f []" "std f";
                             mid_309_265[daystr][:t_inj] pmean.(mid_309_265[daystr][:S]) pstd.(mid_309_265[daystr][:S]) pmean.(ctd265[daystr][:Q]) pstd.(ctd265[daystr][:Q]) pmean.(mid_309_265[daystr][:dphi_dz]) pstd.(mid_309_265[daystr][:dphi_dz]) pmean.(mid_309_265[daystr][:f]) pstd.(mid_309_265[daystr][:f])],
             ", ")
end
