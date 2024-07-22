source("HSR_simulation_fn.R", local = knitr::knit_global())
source("HSR_calib_fn.R", local = knitr::knit_global())

# Starting params; roughly from literature plus some manual error minimization (preliminary)
start_params <- data.frame(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= -.18, phi.2=-.06, 
  
  q= -1.38,
  om.1= -1,
  om.2= -1,
  om.3= -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1=-.03, 
  
  theta.P = -68.5, sigma.P = 69,
  gamma.3=-92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  
  beta.0 = -3.05,
  beta.3=0,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500 , 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200 # number of practices in a region
  
)

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


