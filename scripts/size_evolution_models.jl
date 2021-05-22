# ------------------------------------------ #
#                   Importing                #
# ------------------------------------------ #

cd(@__DIR__)
using KissMCMC, LaTeXStrings


# ------------------------------------------ #
#           Loop over the two days           #
# ------------------------------------------ #

model_runs = Dict()
if !@isdefined n_iterations
    n_iterations = 5*10^4 # at least 10^6 for nice distribution in plot
end

for date in ["0908", "2108"]

    out = Dict()


    # ------------------------------------------ #
    #            Bayesian-model, MCMC            #
    # ------------------------------------------ #

    # initial point of walker
    if date == "0908"
        gradT0 = -5e-4;
    elseif date == "2108"
        gradT0 = -2e-3;
    end
    Sinit_0 = mid_309_265[date][:S][1];
    Sinit_std = std.(mid_309_265[date][:S][1]);
    theta0 = [gradT0, mean.(Sinit_0)];

    # prior information
    dTmin, dTmax = -1e-2, -1e-6;
    # logarithmic posterior density function
    logposterior = R.make_logpdf(date, mid_309_265, [dTmin, dTmax]);

    # emcee MCMC sampler:
    thetas, accept_ratioe, logdensities, blobs = emcee(logposterior, make_theta0s(theta0, [mean.(Sinit_0), Sinit_std], logposterior, 6, hasblob=true), niter=n_iterations, hasblob=true);
    thetas, accept_ratioe, logdensities, blobs = squash_walkers(thetas, accept_ratioe, logdensities, blobs); # puts all walkers into one

    out[:dTdz_MCMC] = [thetas[k][1] for k in 1:length(thetas)];
    # out[:Sinit] = [thetas[k][2] for k in 1:length(thetas)];

    out[:S_MCMC] = [blobs[k][1] for k in 1:length(blobs)];
    out[:dSdt_MCMC] = [blobs[k][2] for k in 1:length(blobs)];
    out[:dSdt_sensible_MCMC] = [blobs[k][3] for k in 1:length(blobs)];


    # ------------------------------------------ #
    #                   ct-model                 #
    # ------------------------------------------ #

    dp = mean(mid_309_265[date][:dp], dims=1) # mean pressure difference over the tracer experiments
    out[:dTdz_ct] = dp * ct / mid_309_265[date][:dz]
    out[:S_ct],
        out[:dSdt_ct],
        out[:dSdt_sensible_ct] = R.model_S(mid_309_265, date, out[:dTdz_ct], Sinit_0)


    # ------------------------------------------ #
    #                closure rate                #
    # ------------------------------------------ #

    out[:closure] = R.closure_rate(ctd309, ctd265, mid_309_265, date)


    # ------------------------------------------ #
    #            dT/dz from measurements         #
    # ------------------------------------------ #

    i = indices[date][1][1]:indices[date][4][end] # indices from appearance of first tracer signal until disappearance of last one
    out[:dTdz_measured] = (Particles.(n_partcl, Normal.(
                            mean(ctd265[date][:temp_PT][i]), e_T))
                        - Particles.(n_partcl, Normal.(
                            mean(ctd309[date][:temp_PT][i]), e_T))
                        ) / mid_309_265[date][:dz]


    # ------------------------------------------ #
    #              Save in dictionary            #
    # ------------------------------------------ #

    model_runs[date] = out

end

# print latex table of total opening rate, opening rate due to sensible heat and closure rate
# for (MCMC, c_t, description) in zip([:dSdt_MCMC, :dSdt_sensible_MCMC],
#                                  [:dSdt_ct, :dSdt_sensible_ct],
#                                  [L"Total opening rate ($\mathrm{m^2\,s^{-1}}$)", L"Opening rate due to sensible heat ($\mathrm{m^2\,s^{-1}}$)"])
#     for date in ["0908", "2108"]
#         if date == "0908"
#             o = description
#         else
#             o = ""
#         end
#         values_ct = only(mean(model_runs[date][c_t], dims=1))
#         values_MCMC = only(mean(
#                     [ bootstrap(
#                         Particles([model_runs[date][MCMC][i][k] for i in 1:length(model_runs[date][MCMC])]),
#                         n_partcl)
#                       for k in 1:length(model_runs[date][MCMC][1]) ]
#                     , dims=1))
#         if date == "0908"
#             o = o * " & " * string(mean(values_ct) ± std(values_ct)) * " & " * string(mean(values_MCMC) ± std(values_MCMC)) * "\\\\ \n"
#         else
#             o = o * " & \\textbf{" * string(mean(values_ct) ± std(values_ct)) * "} & \\textbf{" * string(mean(values_MCMC) ± std(values_MCMC)) * "} \\\\ \n"
#         end
#         print(o)
#     end
# end
# print("\\hline \n")
# closure_rate = only(mean(model_runs["0908"][:closure], dims=1))
# o = L"Closure rate ($\mathrm{m^2\,s^{-1}}$)" * " & \\multicolumn{2}{" * string(mean(closure_rate) ± std(closure_rate)) * "} \\\\ \n"
# print(o)
# closure_rate = only(mean(model_runs["2108"][:closure], dims=1))
# o = " & \\multicolumn{2}{\\textbf{" * string(mean(closure_rate) ± std(closure_rate)) * "}} \\\\ \n"
# print(o)
