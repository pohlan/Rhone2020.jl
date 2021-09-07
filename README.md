# Rhone2020

This repository contains the processing of data from experiments in artificial moulins on Rhonegletscher in 2020. The folders are structured as follows.

## Scripts
- `read_measurements.jl`: loads data from the sensors and creates dictionaries with manually measured quantities (such as mass of salt tracer); uses `calibration.jl` for the conductivity-concentration conversion and the temperature correction
- `derived_quantities.jl`: derives the hydraulic gradient, discharge, flow speed, cross-sectional area, friction factor and manning roughness
- `size_evolution_models.jl`: runs both the ct-model and the free-gradient model
- `heat_transfer.jl`: calculates all thermodynamic variables (e.g. equilibrium offset-temperature)
- `paper_figures.jl`: produces figures for paper and supplements as well as table and CSV files of derived quantities and saves them in the `products` folder); uses `overviewmap.jl` to produce Fig. 1 in the paper and Fig. S1 in the supplements
- `master.jl`: runs all the scripts

#### `additional_figures`
This folder contains julia files used to produce the `.pdf` files in the `products/additional_figures` folder. The figures included in these documents show raw data, results from intermediate steps and results plotted in different ways as in the paper.
- `raw_data.jl`: raw data of the CTDs (conductivity, temperature, pressure)
- `calibration_plots.jl`: conductivity and temperature calibrations
- `derived_quantities_extra.jl`: same quantities are plotted as Fig. 3 in the paper but for all tracer experiments not only the one we selected
- `other_figures.jl`: Moody plot and some plots that show heat transfer data (Figure 5 in paper) differently

- `make.jl`: produces the `.pdf` files and saves them in the `products/additional_figures` folder.

## Dependencies

All Julia dependencies are installed when installing this package with
(at terminal prompt)
```
$ git clone url_of_this_repo
$ cd Rhone2020.jl/scripts
$ julia --project
julia> ]
Rhone2020 pkg> instantiate
```

External dependencies which need to be installed by hand/via the Linux
package manager/homebrew/etc:
- pandoc
- unzip
- gdal (in particular gdal_merge.py)
