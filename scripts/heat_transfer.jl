# ------------------------------------------ #
#                  Constants                 #
# ------------------------------------------ #
Pr = 13.5
k = 0.57  # thermal conductivity of water, W / (K*m)
nu = mu/rhow # kin. viscosity


# ------------------------------------------ #
#       Nusselt number parameterisations     #
# ------------------------------------------ #

# Where A, a and b are defined: Dittus-Boelter correlation A * Pr^a * Re^b
coefs = (
    standard = (A=0.023, a=2/5, b=0.8),
    sommers = (A=0.0025, a=0.4, b=0.95),
    lunardini = (A=0.0078, a=1/3, b=0.927),
    vincent = (A=0.332, a=1/3, b=0.74),
    ogier = (A=1.78, a=1/3, b=0.58),
    gnielinski = ()
)

# ------------------------------------------ #
#                  Formulas                  #
# ------------------------------------------ #

r(S) = sqrt(S/pi)
Re(Q,S) = Q/S * 2 * r(S) / nu
Nu_Dittus_Boelter(Re, p) = p.A * Pr^p.a .* Re.^p.b
Nu_Gnielinski(f, Re) = f./8 .* (Re .- 1000) .* Pr ./
                       (1 .+ 12.7 .* (f./8).^0.5 .* (Pr^(2/3) .- 1))
# select between the two depending on parameter:
Nu_fn(f, Re, p::Tuple{}) = Nu_Gnielinski(f, Re)
Nu_fn(f, Re, p) = Nu_Dittus_Boelter(Re, p)

melt_energy_eq(Q, dphi_dz, dp_dz) = -Q * ( dphi_dz + rhow*cw*ct*dp_dz ) # Energy

# from Sommers and RAJARAM 2020, Eq. 29
z_eq_fn(Q, Nu) = rhow * cw * Q ./ (pi * k * Nu) # e-folding length, m
tau_eq_fn(Q, dphi_dz, dp_dz, Nu) = melt_energy_eq(Q, dphi_dz, dp_dz) / (k * pi * Nu ) # equilibrium offset-temperature: T_eq - T_melt, °C

# same paper, modification of Eq. 21
tau_w_diff(z_eq, dtauw_dz, tau_eq) = tau_eq - z_eq * dtauw_dz # water offset-temperature, T_w - T_melt, °C, calculated from differential equation
tau_w_sol(z, tau_eq, t0, z_eq) = tau_eq + (t0-tau_eq)*exp(-z/z_eq) # water offset-temperature, T_w - T_melt, °C, calculated from solution of differential equation as a function of z


# ------------------------------------------ #
#          Initialise dictionaries           #
# ------------------------------------------ #

z_eq = Dict()
tau_eq = Dict()
tau_w = Dict()
tau_diff = Dict() # Tw - Teq = tau_w - tau_eq
T_surf = Dict()
tauw_measured = Dict()
dtauw_dz = Dict()
Nus = Dict()


# ------------------------------------------ #
#    Calculate heat transfer parameters      #
# ------------------------------------------ #

days = ["0908", "2108"]

for date in days

    # Quantities derived from field experiments
    Qs = mid_309_265[date][:Q]
    Ss = mid_309_265[date][:S]
    f = mid_309_265[date][:f]
    dphi_dzs = mid_309_265[date][:dphi_dz]
    dp_dzs = mid_309_265[date][:dp] ./ mid_309_265[date][:dz]
    T_ct = ct * mid_309_265[date][:p_mean]

    # Data using output of free-gradient model: dtauw_dz = dTw_dz - dTmelt_dz, the first term comes from the MCMC runs
    local dtauw_dz = bootstrap(Particles(model_runs[date][:dTdz_MCMC]), n_partcl) .-  model_runs[date][:dTdz_ct]  # first term has to be converted from array to particles

    # Temperature measurements
    i = indexin(DateTime.(mid_309_265[date][:t_inj]), DateTime.(ctd309[date][:t]))  # extract the indices where there was a tracer experiment
    tauw_measured[date] = Particles(n_partcl, Normal(                             # compute the mean value over the date since it does not vary much
            mean(
            [mean(ctd309[date][:temp_PT][i]),
             mean(ctd265[date][:temp_PT][i])]), e_T)
            ) .- mean(T_ct, dims=1)

    # mean depth of test-section
    d = mean([ctd309[date][:depth], ctd265[date][:depth]], dims=1)

    # compute Reynolds number
    Res = Re.(Qs, Ss)

    # e-folding length
    z_eq[date] = Dict()
    for n = 1:length(coefs)
        c = coefs[n]
        z_eq[date][keys(coefs)[n]] = z_eq_fn.(Qs, Nu_fn(f, Res, c))
    end

    # equilibrium offset-temperature and Nusselt number
    tau_eq[date] = Dict()
    Nus[date] = Dict()
    for n = 1:length(coefs)
        c = coefs[n]
        tau_eq[date][keys(coefs)[n]] = tau_eq_fn.(Qs, dphi_dzs, dp_dzs, Nu_fn(f, Res, c))
        Nus[date][keys(coefs)[n]] = Nu_fn(f, Res, c)
    end

    # water offset-temperature, its difference to tau_eq and the estimated surface temperature
    tau_w[date] = Dict() # Tw - Tmelt
    tau_diff[date] = Dict()
    T_surf[date] = Dict()
    for k = keys(coefs)
        c = coefs[k]
        tau_w[date][k] = tau_w_diff.(z_eq[date][k], dtauw_dz, tau_eq[date][k])
        tau_diff[date][k] = tau_w[date][k] .- tau_eq[date][k]
        T_surf[date][k] = tau_w_sol.(-d, tau_eq[date][k], tau_w[date][k], z_eq[date][k])
    end
end
