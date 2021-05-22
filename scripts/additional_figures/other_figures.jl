# # Figures neither in the paper nor in the supplement
#

using PyPlot #hide

include("../read_measurements.jl"); #hide
include("../derive_quantities.jl"); #hide
use_progress_meter = false #hide
include("../size_evolution_models.jl"); #hide
include("../heat_transfer.jl"); #hide

# ## Moody plot
#
# The Moody plot shows the relation between the Reynolds number and the friction factor.
# Between AM15 on 9 August and AM13 on 21 August the friction factor does not change even though the Reynolds numbers are considerably different.
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams") #hide
rcParams["font.size"] = 9 # general #hide
figure(figsize=(9,5)) #hide
for date in keys(pick) #hide
    errorbar(mean.(mid_309_265[date][:Re][pick[date]]),mean.(mid_309_265[date][:f][pick[date]]),xerr=std.(mid_309_265[date][:Re][pick[date]]),yerr=std.(mid_309_265[date][:f][pick[date]]),fmt="_",label=join(["Aug-",date[1:2]])) #hide
    xlabel("Reynolds number") #hide
    ylabel("Friction coefficient") #hide
end #hide
legend() #hide
xscale("log") #hide
yscale("log") #hide
gcf() #hide

# ## Thermodynamic variables
#
# These plots are based on the same data as Figure 5 in the paper but instead of daily averages the values for the individual tracer experiments are shown.

# ### Equilibrating length scale
R.plot_Nu_params(z_eq, L"z_\mathrm{eq}\,\mathrm{(m)}") #hide

# ### Equilibrium offset-temperature
R.plot_Nu_params(tau_eq, L"\tau_\mathrm{eq}\,(\mathrm{°C})") #hide

# ### Water offset-temperature
R.plot_Nu_params(tau_w, L"\tau_{w}\,(\mathrm{°C})", tauw_measured) #hide

# ### Water offset-temperature as a histogram
#
# This plot shows the same data as Figure 5d in the paper but plotted as a histogram instead of an errorbar plot.
R.plot_tauw_hist(tau_w, tauw_measured) #hide
