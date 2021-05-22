# first make sure that current working directory is scripts and not Rhone2020
cd(@__DIR__)

## settings for uncertainty propagation and MCMC:
ind = 2 # set to 1 for faster processing, to 2 for paper
n_iterations = [5*10^4, 10^6][ind] # at least 10^6 for nice distribution in plot
set_n_partcl([1000, 20000][ind]) # fig. 4 used 20000 to make a nice looking histogram.  But that takes time and memory to run!

# Plotting.  Set to true for interactive plots
pygui(false)

# run all scripts for computations
include("read_measurements.jl") # load all data
include("derive_quantities.jl") # compute hydraulic gradient, discharge, flow speed, cross-sectional area, friction factor
include("size_evolution_models.jl") # run ct-model and Bayesian-model
include("heat_transfer.jl") # comput thermodynamic variables (equilibrium offset-temperature etc.)
include("paper_figures.jl") # plot all figures of the paper and supplements and save them in the products folder

####### TODO #########
# write to files
#using DelimitedFiles
#writedlm("SQ_09_08_2020.txt", ["S [m^2]" "std_S" "Q [m^3/s]" "std_Q" "indices";mean.(mid_309_265["0908"][:S]) std.(mid_309_265["0908"][:S]) mean.(ctd265["0908"][:Q]) std.(ctd265["0908"][:Q]) t0["0908"]'], ", ")
#writedlm("SQ_21_08_2020.txt", ["S [m^2]" "std_S" "Q [m^3/s]" "std_Q" "indices";mean.(mid_309_265["2108"][:S]) std.(mid_309_265["2108"][:S]) mean.(ctd265["2108"][:Q]) std.(ctd265["2108"][:Q]) t0["2108"]'], ", ")
