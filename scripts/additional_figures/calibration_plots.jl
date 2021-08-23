# # Calibrations
#
# The three Keller DCX-22-CTD sensors 205309, 207265 and 205145 and are referred to as CTD-309, CTD-265 and CTD-145.

using Rhone2020  #hide
const R = Rhone2020  #hide

include("../read_measurements.jl"); #hide
using PyPlot, LaTeXStrings; #hide

# ## Conductivity and concentration
#
# The sensors were installed in a bucket where the salt concentration was increased successively. The corresponding conductivity readings were then used to fit a linear function in order to translate conductivity measurements to concentrations.
#

R.plot_conc_cali(calis,(ctd309, ctd265, ctd145), e_m, e_cond) #hide

# ## Zero degree calibration
#
# ### In the lab (23 June 2020)
#
# A box was filled with ice and the coldest possible tap water. After a while, the sensors were installed in the box such that they neither touched the box wall nor each other. The plots show the readings at equilibrium.
#

R.plot_temp_cali(   (ctd309, ctd265, ctd145), #hide
                    [280:600, 180:520, 284:600], #hide
                    "2306", #hide
                    dTOBs = ["-0.031 ± 0.1", # ctd-309 #hide
                             "0.152 ± 0.1", # ctd-265 #hide
                             "0.026 ± 0.1"], # ctd-145 #hide
                    dTs =   ["-0.049 ± 0.01", # ctd-309 #hide
                             "-0.351 ± 0.01", # ctd-265 #hide
                             "-0.010 ± 0.01"], # ctd-145 #hide
                    ) #hide

# ### In the field (20 August 2020)
#
# A bucket was filled with ice and meltwater from the stream. Sensors were installed in the bucket without touching each other, then the bucket was isolated by covering it with as much ice as possible. Unavoidably, some sensors touched the bucket wall. Even though several calibrations were conducted in the field, only the one from 20.08. is usable since only there the sensors stayed in the bucket for long enough to reach the equilibrium. It is advisable to leave the sensors in the bucket for at least 20 to 30 minutes.
#

R.plot_temp_cali(   (ctd309, ctd265, ctd145), #hide
                    [1:640, 1:630, 1:630], #hide
                    "2008", #hide
                    dTOBs = ["-0.06 ± 0.1", # ctd-309 #hide
                             "0.218 ± 0.1", # ctd-265 #hide
                             "0.273 ± 0.1"], # ctd-145 #hide
                    dTs =   ["-0.004 ± 0.01", # ctd-309 #hide
                             "-0.257 ± 0.01", # ctd-265 #hide
                             "?"], # ctd-145 #hide
                    ) #hide

# ### Summary of temperature calibrations
#
# Table shows the measurements in °C.
#
# CO = Christophe Ogier, AP = Annegret Pohle
#
#   | CTD number | Sensor | CO lab (08.10.19) | CO field (25.07., 20?) | AP lab (23.06.20)   | AP field (20.08.20) |
#   |------------|--------|-------------------|--------------------------|-------------------|---------------------|
#   | 205309     |  TOB1  |  0.018            |    0.054                 | -0.031            |    0.060            |
#   | 205309     |  T     | -0.021            |   -0.014                 | -0.049            |   -0.004            |
#   | 207265     |  TOB1  |  0.161            |    0.338                 |  0.152            |    0.218            |
#   | 207265     |  T     | -0.319            |   -0.097                 | -0.351            |   -0.257            |
#   | 205145     |  TOB1  |  0.002            |    0.109                 |  0.026            |    0.273            |
#   | 205145     |  T     |  0.076            |    0.125                 | -0.010            |       ?             |
#
# During field calibrations the temperatures seemed to be less stable. Therefore we take the mean value between CO's and AP's lab calibrations. For the paper we only use TOB1 sensors of CTD-309 and CTD-265.
