# calibration readouts
# --------------------
bucketsize = 1.0 # calibration bucket size in liters
solution = 1.0 # calibration solution concentration (g/l)
# calis: total calibration ml solution, 1st column, vs sensor readout (μS/cm), 2nd column
calis = Dict( 309 => (# [0  1                    # 23.06.20 (only a test before going to the field)
                      # 1  2
                      # 6  10
                      # 16 25],
                      [ 0  1                    # 05.08.20
                        1 2
                        6 8
                        16 22
                        36 47
                        76 90
                        116 131
                        166 191
                        191 229],
                      [ 0 1                     # 10.08.20
                        1 2
                        6 8
                        16 27
                        36 48
                        76 95
                        116 144],
                      [ 0 1                     # 13.08.20
                        1 2
                        6 8
                        16 26
                        36 53
                        76 98
                        116 143
                        156 186
                        256 345],  # disturbs linear fit quite a lot
                      [ 0  1                    # 20.08.20
                        2  3
                        12 19
                        42 57
                        102 121
                        162 182
                        222 251
                        262 270]),
              265 => (# [ 0 1                     # 23.06.20  (only a test before going to the field)
                      #   1 3
                      #   6 11
                      #   16 26],
                      [ 0  1                    # 05.08.20
                        1 2
                        6 8
                        16 23
                        36 47
                        76 91
                        116 133
                        166 193
                        191 235],
                      [ 0 1                     # 10.08.20
                        1 2
                        6 8
                        16 27
                        36 54
                        76 100
                        116 147],
                      [ 0 1                     # 13.08.20
                        1 2
                        6 8
                        16 26
                        36 54
                        76 99
                        116 140
                        156 188
                        256 352],  # disturbs linear fit quite a lot
                      [ 0  1                    # 20.08.20
                        2  3
                        12 19
                        42 57
                        102 122
                        162 183
                        222 252
                        262 267]),
              145 => ([0  1                     # 23.06.20  (only a test before going to the field, but not enough data points for a proper calibration anyway)
                        1  2
                        6  10
                        16 25],
                       [0  0                     # 20.08.20
                        2  3
                        12 18
                        42 57
                        102 121
                        162 182
                        222 251])
            )


# fit a linear function through calibration readouts
ctd309["cond2conc"] = Rhone2020.fit_calibration(bucketsize, solution, e_cond, calis[309][1], calis[309][2], calis[309][3], calis[309][4])
ctd265["cond2conc"] = Rhone2020.fit_calibration(bucketsize, solution, e_cond, calis[265][1], calis[265][2], calis[265][3], calis[265][4])
ctd145["cond2conc"] = Rhone2020.fit_calibration(bucketsize, solution, e_cond, calis[145][1], calis[145][2])

# temperature calibration (see calibration_plots.md)
# only possible for lab calibrations of CTD-309 and CTD-265; take mean value between Christophe and Annegret
T309_shift = -0.035
T265_shift = -0.335
T145_shift = 0.033

for date in keys(indices)
  ctd309[date][:temp_PT] = ctd309[date][:temp_PT] .- T309_shift
  ctd265[date][:temp_PT] = ctd265[date][:temp_PT] .- T265_shift
  if haskey(ctd145, date)
    ctd145[date][:temp_PT] = ctd145[date][:temp_PT] .- T145_shift
  end
end
