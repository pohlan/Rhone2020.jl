# # Results of all experiments
#
## Run the calculations #hide
include("../read_measurements.jl"); #hide
include("../derived_quantities.jl"); #hide

# Filled symbols: experiments picked for results in the paper, criteria..
#
# - all sensors were below the water table in the borehole (if a sensor was not in the water column it was recording atmospheric pressure)
# - discharge values at the two sensors are consistent
#
# Relative position of the CTDs, depth increasing: CTD-145, CTD-309, CTD-265 (see `raw_data.pdf` for more details)
#
#
# Hydraulic gradients are computed from pressure smoothed with moving average (time window 2 minutes for AM15 and 30 minutes for AM13).
# Uncertainties of the hydraulic gradient are not shown in the plots, they would hardly be visible.
#
# ### AM13, 8 August to 19 August
#+
# In AM13 only CTD-309 and CTD-265 were operating.
#+ #hide
R.plot_data_unselected("0808", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("0908", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1008", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1108", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1208", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1308", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1408", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1708", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1808", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("1908", ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide

# ### AM15, 20 and 21 August
#+
# In theory the third CTD in AM15 allowed the characterisation of a second channel section between the upper two CTDs.
# However, the discharges between the upper two sensors did not agree to within their uncertainties.
# Either there was an inflow in between or the discharge at the uppermost sensor was erroneous, potentially because the salt tracer had not mixed across the cross-section or was not completely dissolved yet.
# In both cases our results would not be meaningful, hence we decided to focus on the measurements between the two lower sensors.
#
# On 20 August there was hardly any discharge, which was already visible in the field.
#+ #hide
R.plot_data_unselected("2008",ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
#+ #hide
R.plot_data_unselected("2108",ctd309, ctd265, ctd145, mid_surf_309, mid_309_265, mid_145_309, pick, e_p) #hide
