# ------------------------------------------ #
#  Import functions defined in Rhone2020.jl  #
# ------------------------------------------ #

cd(@__DIR__);
using Pkg;
Pkg.activate("..");
using Rhone2020;
const R = Rhone2020;

using Dates, MonteCarloMeasurements;


# ------------------------------------------ #
#    Load the data files into dictionaries   #
# ------------------------------------------ #

ctd145, ctd309, ctd265 = R.load_DCX22_CTD();
e_cond = 200*0.025 # error of conductivity in μS/cm, 2.5% of range FS
e_p = rhow*g*0.0005 # error of pressure, 0.05% FS; has to be multiplied with max. water column, e.g. *100 for ctd309-100mH20


# ------------------------------------------ #
#       Pick start/end points of traces      #
# ------------------------------------------ #

# plot conductivity to pick start and end points of individual injection signals
#using PyPlot
#pygui(true) # interactive window for figure
#plot(ctd309["1008"][:cond])
#plot(ctd265["1008"][:cond])

# manually determined start and end points of individual injection signals
# one column per injection
# 1st row: first index of ctd309
# 2nd row: last index of ctd309
# 3rd row: first index of ctd265
# 4th row: last index of ctd265
# 5th row: first index of ctd145
# 6th row: last index of ctd145

indices = Dict("0808" => ([16882 18408 19638 21934],  # traces 1 and 6-9: supersaturation
                          [17071 19079 19741 22000],
                          [17088 19355 19726 21992],
                          [17339 19474 19869 22080]),
               "0908" => ([8313 10072 11810 13638 15494 17204 18709 21049 22606 24239 26059 27914], # injection 13 too late for ctd265
                          [8386 10120 11900 13710 15560 17270 18784 21130 22682 24300 26140 27990],
                          [8339 10097 11833 13657 15514 17224 18732 21073 22630 24263 26086 27944],
                          [8440 10180 11950 13750 15590 17310 18815 21160 22716 24350 26170 28037]),
               "1008" => ([11829 15747 19984],
                          [11930 15795 20010],
                          [12070 15800 20040],
                          [12320 15880 20090]),
               "1108" => ([5307 7100 8213 9386 10902 11964 13122 18468 19575 25237],
                          [5337 7124 8246 9446 10984 12016 13152 18498 19604 25261],
                          [5329 7139 8273 9459 10960 12006 13164 18519 19631 25290],
                          [5367 7184 8330 9532 11050 12085 13213 18568 19693 25339]),
               "1208" => ([4134 5529 7557 8772 10018 25506],
                          [4181 5574 7592 8811 10050 25562],
                          [4148 5556 7592 8804 10052 25544],
                          [4218 5622 7662 8866 10104 25641]),
               "1308" => ([5437 7555 14943],
                          [5533 7606 14991],
                          [5441 7566 15017],
                          [5549 7641 15110]),
               "1408" => ([4642 5841 7764 9688 19782],
                          [4692 5893 7817 9729 19827],
                          [4649 5857 7805 9725 19822],
                          [4733 5955 7894 9830 19906]),
               "1708" => ([10072 12345 23422],
                          [10149 12437 23563],
                          [10101 12369 23444],
                          [10220 12514 23661]),
               "1808" => ([2416 21321],
                          [2530 21407],
                          [2421 21359],
                          [2540 21526]),
               "1908" => ([1732 2288 27037 28252],
                          [1925 2424 27187 28377],
                          [1740 2301 27091 28286],
                          [1939 2529 27318 28527]),
               "2008" => ([9960],[25900],[15160],[25900],[6213],[10571]), # not usable for calculations
               "2108" => ([722 2554 4313 6102 7894 9744 11763 13359 15143 16950 18758 20472 22224 23230 24193 26416 27983], # 309 and 265: end point = start point of next trace
                          [2554 4313 6102 7894 9744 11763 13359 15143 16950 18758 20472 22224 23230 24193 26416 27983 length(ctd309["2108"][:cond])],
                          [841 2670 4422 6195 7995 9879 11932 13499 15296 17101 18908 20634 22378 23389 24336 26564 28118],
                          [2670 4422 6195 7995 9879 11932 13499 15296 17101 18908 20634 22378 23389 24336 26564 28118 length(ctd265["2108"][:cond])],
                          [622 2446 4223 6021 7821 9627 11656 13233 15031 16830 18629 20363 22104 23122 24080 26307 27893],
                          [688 2512 4310 6100 7893 9760 11770 13355 15150 16943 18737 20461 22208 23228 24158 26400 27980]));

#plot to check if indices are ok
#using PyPlot
#pygui(true)
#date = "2108"
#figure()
#for tr in 1:length(indices[date][1])
#    plot(ctd309[date][:cond][indices[date][1][tr]:indices[date][2][tr]])
#end
#figure()
#for tr in 1:length(indices[date][1])
#    plot(ctd265[date][:cond][indices[date][3][tr]:indices[date][4][tr]])
#end


# ------------------------------------------ #
#            Indices of injections           #
# ------------------------------------------ #

# time of salt injection (translated to indices of ctd arrays)
t0 = Dict("0808" => ([16806 18186 19606 21906]),
          "0908" => ([8256 10026 11766 13596 15456 17166 18666 21006 22566 24196 26016 27876]),
          "1008" => ([11586 15441 19701]),
          "1108" => ([5291 7081 8191 9301 10816 11896 13066 18436 19521 25216]),
          "1208" => ([4126 5521 7546 8761 10006 25501]),
          "1308" => ([5431 7546 14926]),
          "1408" => ([4631 5831 7756 9676 19771]),
          "1708" => ([10066 12346 23411]),
          "1808" => ([2411 21316]),
          "1908" => ([1726 2281 27031 28246]),
          "2008" => ([4441]),
          "2108" => ([601 2416 4201 6001 7801 9601 11626 13201 15001 16801 18606 20341 22081 23101 24061 26281 27871]));
e_dt = 3 # estimated error for time difference in s between watch and ctd clock; dt error between two ctds is 0.5 s


# ------------------------------------------ #
#                Tracer masses               #
# ------------------------------------------ #

# masses of injected salt in kg
mass = Dict("0808" => [0.1;0.1;0.1;0.1],
            "0908" => [0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05],
            "1008" => [0.05;0.05;0.05],
            "1108" => [0.05;0.05;0.07;0.07;0.1;0.1;0.1;0.1;0.1;0.1],
            "1208" => [0.05;0.1;0.1;0.1;0.1;0.06],
            "1308" => [0.1;0.1;0.1],
            "1408" => [0.1;0.1;0.1;0.1;0.1],
            "1708" => [0.1;0.1;0.1],
            "1808" => [0.1;0.1],
            "1908" => [0.1;0.1;0.1;0.1],
            "2008" => [0.02],
            "2108" => [0.025;0.025;0.025;0.025;0.05;0.05;0.025;0.025;0.025;0.025;0.025;0.025;0.025;0.05;0.025;0.025;0.025]
            );

# repeatedly measured 100g with analogue suspended scale (as in the field) and compared it to kitchen scale:
m = [99.7,100.3,100.1,99.5,99.9,99.7,101.1,99.7,99.1,98.7,99.8,99.3,100.3,99.8,98.6,100.1,99.9,99.5,100.1,99.2,98.7,100.6];
e_m = 0.002 # error of mass in kg (small influence)


# ------------------------------------------ #
#              Depths of CTDs                #
# ------------------------------------------ #

e_dz = 0.5 # estimated error of vertical distance in m

# depths of ctds in m for each day [ctd309, ctd265, ctd145]
# upper 100m: add 0.04*(depth-100m),; 104m measured under tension
# lower 100m: add 0.01*depth; 101m measured under tension

ctd309["0808"][:depth], ctd265["0808"][:depth] = Particles.(n_partcl, Normal.((62.4, 164.0), e_dz))
ctd309["0908"][:depth], ctd265["0908"][:depth] = Particles.(n_partcl, Normal.((121.1, 171.7), e_dz))
ctd309["1008"][:depth], ctd265["1008"][:depth] = Particles.(n_partcl, Normal.((120.9, 171.7), e_dz))
ctd309["1108"][:depth], ctd265["1108"][:depth] = Particles.(n_partcl, Normal.((121.0, 171.7), e_dz))
ctd309["1208"][:depth], ctd265["1208"][:depth] = Particles.(n_partcl, Normal.((121.0, 171.7), e_dz))
ctd309["1308"][:depth], ctd265["1308"][:depth] = Particles.(n_partcl, Normal.((131.6, 171.9), e_dz))
ctd309["1408"][:depth], ctd265["1408"][:depth] = Particles.(n_partcl, Normal.((131.3, 171.7), e_dz))
ctd309["1708"][:depth], ctd265["1708"][:depth] = Particles.(n_partcl, Normal.((131.3, 171.7), e_dz))
ctd309["1808"][:depth], ctd265["1808"][:depth] = Particles.(n_partcl, Normal.((131.3, 171.7), e_dz))
ctd309["1908"][:depth], ctd265["1908"][:depth] = Particles.(n_partcl, Normal.((141.4, 171.7), e_dz))
ctd145["2008"][:depth], ctd309["2008"][:depth], ctd265["2008"][:depth] = Particles.(n_partcl, Normal.((10.4, 88.4, 189.9), e_dz))
ctd145["2108"][:depth], ctd309["2108"][:depth], ctd265["2108"][:depth] = Particles.(n_partcl, Normal.((10.4, 88.4, 189.9), e_dz))


# ------------------------------------------ #
#                Calibrations                #
# ------------------------------------------ #

# get functions converting conductivities to concentrations as well as shift of temperature sensors from 0°C
include("calibration.jl")
e_T = 0.05 # °C, error of temperature measurements introduced by calibration (otherwise it would be 0.1°C for PT sensor)


# ------------------------------------------ #
#           Subtract air pressure            #
# ------------------------------------------ #

stage = R.load_stage_pressure() # pressure sensor installed upstream, contains air pressure; sampling interval 2min
keller_p = R.load_DCX22_pressure() # pressure sensor installed in moulins over night; not used
for key in [collect(keys(indices));"1608"]
    if haskey(ctd309,key)
        local idx_stage309 = indexin(DateTime.(ctd309[key][:t]),DateTime.(stage[:t]))
        idx_stage309 = idx_stage309[idx_stage309.!==nothing]
        local idx_ctd309 = indexin(DateTime.(stage[:t]),DateTime.(ctd309[key][:t]))
        idx_ctd309 = idx_ctd309[idx_ctd309.!==nothing]
        local idx_stage265 = indexin(DateTime.(ctd265[key][:t]),DateTime.(stage[:t]))
        idx_stage265 = idx_stage265[idx_stage265.!==nothing]
        local idx_ctd265 = indexin(DateTime.(stage[:t]),DateTime.(ctd265[key][:t]))
        idx_ctd265 = idx_ctd265[idx_ctd265.!==nothing]
        ctd309[key][:dpress] = ctd309[key][:press][idx_ctd309]-stage[:airpress][idx_stage309]
        ctd309[key][:t_dpress] = ctd309[key][:t][idx_ctd309]
        ctd265[key][:dpress] = ctd265[key][:press][idx_ctd265]-stage[:airpress][idx_stage265]
        ctd265[key][:t_dpress] = ctd265[key][:t][idx_ctd265]
        ctd265[key][:idx_stage] = idx_stage265
    end
    if haskey(ctd145,key)
        local idx_stage145 = indexin(DateTime.(ctd145[key][:t]),DateTime.(stage[:t]))
        idx_stage145 = idx_stage145[idx_stage145.!==nothing]
        local idx_ctd145 = indexin(DateTime.(stage[:t]),DateTime.(ctd145[key][:t]))
        idx_ctd145 = idx_ctd145[idx_ctd145.!==nothing]
        ctd145[key][:dpress] = ctd145[key][:press][idx_ctd145]-stage[:airpress][idx_stage145]
        ctd145[key][:t_dpress] = ctd145[key][:t][idx_ctd145]
    end
    if haskey(keller_p,key) || key == "1608"
        local dt = keller_p[key][:t][1]-stage[:t]
        local i1_stage = findmin(abs.(dt))
        local iend_stage = min(i1_stage[2]+length(keller_p[key][:t])-1, length(stage[:t]))
        local iend_kellerp = iend_stage-i1_stage[2]+1
        keller_p[key][:dpress] = keller_p[key][:press][1:iend_kellerp]-stage[:airpress][i1_stage[2]:iend_stage]
        keller_p[key][:t_dpress] = keller_p[key][:t][1:iend_kellerp]
    end
end
stage[:dpress] = stage[:press]-stage[:airpress];
