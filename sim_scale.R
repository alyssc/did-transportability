## Scaling up the simulation
library(foreach)
library(doParallel)

# Parallelization
myCluster <- makeCluster(3, 
                         type = "PSOCK") # starts new R session here
registerDoParallel(myCluster)

source("HSR_simulation_fn.R")

nreps <- 1000
output <- foreach(r = c(1:nreps), .packages=c('tidyverse'), .combine = 'rbind') %dopar% {
  
  # Fixed and varied parameters
  global_params <- data.frame(expand.grid(
    x1.r = -.617,
    x2.r = -.715,
    phi.1= c(-.17,.17), # log odds ratio of Pr(SSP) in CPC+ vs non-CPC+
    phi.2= c(-.14,.14), # log odds ratio of Pr(system) in CPC+ vs non-CPC+

    q= -1.38,
    om.1= c(-.25,.25), #log odds ratio of Pr(Black) in SSP vs non-SSP
    om.2= c(-.23,.23), #log odds ratio of Pr(Black) in system vs indep
    om.3= -.235,

    H = 0, sigma.H = 0.02,
    psi.1=-.03,

    theta.P = -76, sigma.P =69,
    gamma.3=  -16,
    gamma.4 = 282,
    gamma.5 = -277,

    beta.0 = -3.05,
    beta.3=0,
    beta.4 = .839,
    beta.5 = 1.17,
    beta.6 = .778,

    alpha.0 = 10100,
    alpha.1 = 46500,
    alpha.2 = 600,
    alpha.3 = 12,

    P= 1200, # number of practices in a region
    replicate=r)) %>%
    mutate(scenario=rep(1:16))

  ## Simulate the data
  simdat <- global_params %>% group_by(scenario,replicate) %>%
    nest() %>% mutate(singlesim=map(data,make_regions)) %>%
    unnest(cols=c(data,singlesim)) %>%
    mutate(group=factor(S,levels=c('1','0'),labels=c('Sample','Target')),
           trt=factor(A,levels=c('0','1'),labels=c('Untreated','Treated'))) %>%
    group_by(scenario,replicate,S) %>%
    mutate(total.Black=sum(B)) %>% ungroup() %>%
    mutate(wt=B/total.Black) # proportion of all Black benes w/in {S}


  ## Truth and estimates of the SATT and PATT across simulation reps
  estimates <- simdat %>% group_by(scenario,replicate) %>%
    nest() %>% mutate(patt=map(data,estimate_patt),
                      true.patt=map(data,true_patt),
                      satt=map(data,lm.did)) %>%
    unnest(cols=c(patt,true.patt,satt)) %>% select(-data) %>%
    mutate(dr.patt.err=dr-true.patt,
           satt.err=est.satt-true.satt,
           true.diff=true.patt-true.satt,
           est.diff=dr - est.satt)
  return(estimates)
}

stopCluster(myCluster)




# Fixed params
global_params <- data.frame(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= -.18, phi.2=-.06, 
  
  q= -1.38,
  om.1= -.165,
  om.2= -3.2,
  om.3= -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1=-.03, 
  
  theta.P = -76, sigma.P =0, 
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

