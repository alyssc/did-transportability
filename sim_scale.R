## Scaling up the simulation
library(foreach)
library(doParallel)
source("HSR_simulation_fn.R")

nreps <- 1000

scenario_params <- data.frame(expand.grid(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= c(-.18,.18),  # log odds ratio of Pr(SSP) in CPC+ vs non-CPC+
  phi.2=c(-.06,.06),   # log odds ratio of Pr(system) in CPC+ vs non-CPC+
  
  q= -1.28,
  om.1 = c(-1,1),  # log odds ratio of Pr(Black) in SSP vs non-SSP
  om.2 = c(-1,1), # log odds ratio of Pr(Black) in system vs indep
  om.3 = -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1 = -.03, 
  
  theta.P = -68.5, sigma.P = 69, 
  gamma.3 = -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  
  beta.0 = -3.05,
  beta.3 = 0,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200)) %>%
  mutate(scenario=c(1:16))

# Helper function to rbind the lists
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

# Parallelization

myCluster <- makeCluster(100) 
registerDoParallel(myCluster)
start <- Sys.time()
output <- foreach(r = c(1:nreps), .packages=c('tidyverse'), .combine = 'comb',.multicombine = TRUE) %dopar% {
  ## Simulate the data
  simdat <- scenario_params %>% mutate(replicate=r) %>% group_by(scenario,replicate) %>%
    nest() %>% mutate(singlesim=map(data,make_regions)) %>%
    unnest(cols=c(data,singlesim)) %>%
    mutate(group=factor(S,levels=c('1','0'),labels=c('Sample','Target')),
           trt=factor(A,levels=c('0','1'),labels=c('Untreated','Treated'))) %>%
    group_by(scenario,replicate,S) %>%
    mutate(total.Black=sum(B)) %>% ungroup() %>%
    mutate(wt=B/total.Black) # proportion of all Black benes w/in {S}
  stats <- sumstats(simdat)

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
  return(list('sumstats'=stats,'ests'=estimates))
}
stop <- Sys.time()
stopCluster(myCluster)
stop-start

#quick and dirty repackaging of the results object
output$sumstats <- list('by.S'=do.call("rbind",output$sumstats[,'by.S']),
                        'by.X'=do.call("rbind",output$sumstats[,'by.X']),
                        'by.A'=do.call("rbind",output$sumstats[,'by.A']),
                        'by.AS'=do.call("rbind",output$sumstats[,'by.AS']),
                        'by.AS.long'=do.call("rbind",output$sumstats[,'by.AS.long']),
                        'diff.by.AS.long'=do.call("rbind",output$sumstats[,'diff.by.AS.long']))
save(output,file="2024-04-13_HSR_results.RData")
