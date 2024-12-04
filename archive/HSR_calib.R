source("HSR_simulation_fn.R", local = knitr::knit_global())
source("HSR_calib_fn.R", local = knitr::knit_global())

# Check current error
calib_avg(start_params,table=1)
calib_avg(start_params,table=2)
calib_avg(start_params,table=3)

# Calibration process
params1 <- calib_table1(start_params)
params2 <- calib_table2(params1)
params3 <- calib_table3(params2)
save(params3,file="calib_param.RData")

# Check error post-calibration
calib_avg(params3,table=1)
calib_avg(params3,table=2)
calib_avg(params3,table=3)


