source("HSR_simulation_fn.R", local = knitr::knit_global())
source("HSR_param_calib.R", local = knitr::knit_global())

# Rough starting params from literature; some are averages of a range
# Note that many of these are log-odds of values from literature

start_params <- data.frame(
  x1.r = -.6, 
  x2.r = -.7, 
  phi.1=.17, phi.2=.04, 
  
  q= -2.12, 
  om.1= 0, 
  om.2= 0, 
  om.3= -.3,
  
  H = 0, sigma.H = .02, 
  psi.1 = -.03, 
  
  theta.P = -76, sigma.P = 69, 
  gamma.3=-92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  
  beta.0 = -2.51, 
  beta.3=0, 
  beta.4 = .809,
  beta.5 = 1.17,
  beta.6 = .46, 
  
  alpha.0 = 10100,
  alpha.1 = 46500 , 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200
) 

# Check current error
calib_assist1(start_params)
calib_assist2(start_params)
calib_assist3(start_params)



calib_targets <- c(
  .146 (q)	(0.345, 0.353) (x.1r) 	(.33) (x.2r) 	0.1 (beta.0)
  .12 (om.3 (or om.1/om.2))	0.310 (phi.1)	0.316 (phi.2)	.18 (beta.6)
  .141 (beta.0,beta.6)	0.269 (beta.4)	0.121 (beta.0,beta.6)	0.309 (beta.5) 
  .11 
  0.069
  
  
)


test <- data.frame(
  a= .1,
  b= c(.1,2)
)

