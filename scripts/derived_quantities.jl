# ------------------------------------------ #
#                 Importing                  #
# ------------------------------------------ #

# packages
using Statistics


# ------------------------------------------ #
#           Initialise dictionaries          #
# ------------------------------------------ #

# dictionaries for quantities of probed section (between the sensors)
mid_309_265 = Dict("0808" => Dict(), # holds all quantities averaged between ctds 205309 and 207265
                   "0908" => Dict(),
                   "1008" => Dict(),
                   "1108" => Dict(),
                   "1208" => Dict(),
                   "1308" => Dict(),
                   "1408" => Dict(),
                   "1708" => Dict(),
                   "1808" => Dict(),
                   "1908" => Dict(),
                   "2008" => Dict(),
                   "2108" => Dict()
)
mid_145_309 = Dict("2008" => Dict(), # holds all quantities averaged between ctds 205145 and 205309
                   "2108" => Dict()
)
mid_surf_309 = Dict("0808" => Dict(), # holds all quantities averaged betweed the surface and ctd 205309
                    "0908" => Dict(),
                    "1008" => Dict(),
                    "1108" => Dict(),
                    "1208" => Dict(),
                    "1308" => Dict(),
                    "1408" => Dict(),
                    "1708" => Dict(),
                    "1808" => Dict(),
                    "1908" => Dict(),
                    "2008" => Dict(),
                    "2108" => Dict()
)
mid_surf_145 = Dict("2008" => Dict(), # holds all quantities averaged betweed the surface and ctd 205145
                    "2108" => Dict()
)


# ------------------------------------------ #
#            Derive all quantities           #
# ------------------------------------------ #

vtype = [:v, :v_HM][1] # choose whether to take v from peak (1) or harmonic mean (2) time difference

for date in keys(indices)
    # between CTD-309 and CTD-265
    # --------------------------- #

    # vertical distance between two ctds, m
    mid_309_265[date][:dz] = ctd265[date][:depth] - ctd309[date][:depth]

    # discharge, in m^3/s
    ctd309[date][:Q] = R.discharge(ctd309, date, indices[date][1:2], mass, e_m, e_cond)
    ctd265[date][:Q] = R.discharge(ctd265, date, indices[date][3:4], mass, e_m, e_cond)
    mid_309_265[date][:Q] = mean([ctd309[date][:Q],ctd265[date][:Q]],dims=1)[1] # mean discharge between sensors, used for further calculations

    # flow speed, in m/s, according to time difference of conductivity peaks (:v) and difference of harmonic mean times (:v_HM)
    # we use (:v); in AM15 (8-19 August) it doesn't make a significant difference, in AM13 (20 and 21 August) f and S are much more consistent for (:v)
    mid_309_265[date][:t_inj],
        mid_surf_309[date][:v],
        mid_309_265[date][:v],
        mid_surf_309[date][:v_HM],
        mid_309_265[date][:v_HM] = R.flow_speed(ctd309, ctd265, date, indices[date][1:4], t0, e_dt)
    mid_surf_309[date][:t_inj] = mid_309_265[date][:t_inj]

    # smoothed pressure difference, Pa, hydraulic gradient, Pa/m, and mean pressure, Pa
    mid_309_265[date][:dphi_dz],
        mid_309_265[date][:dp],
        mid_309_265[date][:p_mean] = R.hydraulic_gradient(ctd309, ctd265, mid_309_265, date, indices[date], e_p)

    # cross-sectional area S, m^2
    mid_309_265[date][:S] = mid_309_265[date][:Q] ./ mid_309_265[date][vtype]

    # Darcy-Weisbach friction factor f, unitless, and manning roughness n', s/m^{1/3}
    mid_309_265[date][:f],
        mid_309_265[date][:n_manning] = R.friction(mid_309_265, date, vtype)

    # Reynolds number, unitless
    mid_309_265[date][:Re] = R.Reynolds_number(mid_309_265, date, vtype)


    #   between CTD-145 and CTD-309, only on 20 and 21 August
    # ------------------------------------------------------- #

    if haskey(ctd145, date)
        # vertical distance between two ctds, m
        mid_145_309[date][:dz] = ctd309[date][:depth] - ctd145[date][:depth]

        # discharge, m^3/s
        ctd145[date][:Q] = R.discharge(ctd145,date,indices[date][5:6],mass,e_m,e_cond)
        mid_145_309[date][:Q] = mean([ctd145[date][:Q],ctd309[date][:Q]],dims=1)[1] # mean discharge between sensors, used for further calculations

        # flow speed, m/s
        mid_145_309[date][:t_inj],
            mid_surf_145[date][:v],
            mid_145_309[date][:v],
            mid_surf_145[date][:v_HM],
            mid_145_309[date][:v_HM] = R.flow_speed(ctd145, ctd309, date, indices[date][[5,6,1,2]], t0, e_dt)
        mid_surf_145[date][:t_inj] = mid_145_309[date][:t_inj]

        # smoothed pressure difference, Pa, hydraulic gradient, Pa/m, and mean pressure, Pa
        mid_145_309[date][:dphi_dz],
            mid_145_309[date][:dp],
            mid_145_309[date][:p_mean] = R.hydraulic_gradient(ctd145, ctd309, mid_145_309, date, indices[date], e_p)

        # cross-sectional area S, m^2
        mid_145_309[date][:S] = mid_145_309[date][:Q] ./ mid_145_309[date][vtype]

        # Darcy-Weisbach friction factor f, unitless, and manning roughness n', s/m^{1/3}
        mid_145_309[date][:f],
            mid_145_309[date][:n_manning] = R.friction(mid_145_309, date, vtype)

        # Reynolds number, unitless
        mid_145_309[date][:Re] = R.Reynolds_number(mid_145_309, date, vtype)
    end
end


# ------------------------------------------ #
#          Select trusted experiments        #
# ------------------------------------------ #

# pick out the data points where 1.) both sensors where in the water, 2.) discharges between the two sensors are consistent
pick = Dict("0808" => [3,4],
            "0908" => 1:12,
            "1008" => [2,3],
            "1108" => 4:9,
            "1308" => [3],
            "2108" => 1:17
);
