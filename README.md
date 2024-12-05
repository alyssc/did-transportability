# DID transportability

Code for "Transporting Difference-in-Differences Estimates to Assess Health Equity Impacts of Payment and Delivery Models"

Includes simulation scripts and plots used to generate paper results. 

* `HSR_simulation_fn.R`: Defines functions for generating simulations (eg. `make_regions`) and estimation (eg. `estimate_patt`) 
* `sim_scale.R`: Scales simulation; saves output as `.Rdata` files in data folder
* `HSR_plot_scaled_fn.R`: Defines functions to plot scaled simulation results
* `HSR_plot_scaled.R`: Plots scaled simulation results and saves to plots folder

Calibrating parameters to target values from literature (see Supplemental Materials, Table 1): 
* `HSR_calib_fn.R`: Helper functions
* `HSR_calib.R`: Script to calibrate parameters to target values. Done partly manually; see code for instructions. 
