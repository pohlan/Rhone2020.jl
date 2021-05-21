# Rhone2020

This repository contains the processing of data from experiments in artificial moulins on Rhonegletscher in 2020. The folders are structured as follows.

### Scripts
- `discharge_flowspeed.jl`: computes discharges and flowspeeds; uses `fit_calibration.jl` for the conductivity-concentration conversion.
- `friction_cross_section.jl`: loads data from the other pressure sensors and calculates friction parameter, the cross-sectional area and the Reynolds number. `discharge_flowspeed.jl` needs to be run in advance.
- `master.jl`: computes all parameters by running `discharge_flowspeed.jl` and `fricion_cross_section.jl` (no plotting)
- `not_working.jl`: some failed tries to get a good stage-discharge relation for continuous discharge.

### Data
- `DCX22_CTD`: one file for each CTD on each day of measurements; the only sensors that are relevant for the paper
- `DCX22_pressure`: one file for each time the sensor was installed in the borehole (mostly at night)
- `Stage_pressure`: files for arbitrary time periods; sensor was installed upstream of the moulins
- `WTW`: handheld device for conductivity and temperature, installed a few times at 1 to 2 m depth

### Plots
This folder contains markdown (`.md`) files of all the plots and the corresponding julia (`.jl`) files to create them
- `raw_data.xx`: raw data of the CTDs (conductivity, temperature, pressure) and the other pressure sensors
- `calibration_plots.xx`: conductivity and temperature calibrations
- `results.xx`:
  - all the derived quantities for data points that were considered to be 'reasonable' (see `all_data.xx` for justification)
  - Moody plot
  - evolution of the cross-sectional area over time vs. calculations from simple model
- `all_data.xx`: shows results of all data points, also the ones which were not considered to be 'reasonable'
