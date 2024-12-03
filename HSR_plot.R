## Making figures from simulation output
library(tidyverse)
source("HSR_simulation_fn.R")
source("HSR_plot_fn.R")
NREP <- 25

#### Base case ####
## Final rated params
base_params <- data.frame(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1 = -.18, phi.2=-.06, 
  
  q = -1.28,
  om.1 = -1,
  om.2 = -1,
  om.3 = -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1 = -.03, 
  psi.2 =0,
  psi.3 = 0,
  
  theta.P = -68.5, sigma.P = 0, 
  gamma.3 = -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=0,
  gamma.7=0,
  
  beta.0 = -3.05,
  beta.3 = 0,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200,
  scenario=1,
  replicate= 1:NREP)

base_simdat <- simanalyze(base_params)
base_sumstats <- sumstats(base_simdat$data)
calibplots(base_sumstats,save.figs=T)
datplots(base_sumstats,save.figs=T,prefix="plots/base")
estplots(base_simdat$ests,base_sumstats,save.figs=T,prefix="plots/base")

#### 16 main scenarios ####
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
  psi.2 =0,
  psi.3 = 0,
  
  theta.P = -68.5, sigma.P =0, 
  gamma.3 = -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=0,
  gamma.7=0,
  
  beta.0 = -3.05,
  beta.3 = 0,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200,
  replicate = 1:NREP )) %>%
  mutate(scenario=rep(1:16,NREP))

scenario_simdat <- simanalyze(scenario_params)
scenario_sumstats <- sumstats(scenario_simdat$data)
datplots(scenario_sumstats,save.figs=T,prefix="plots/scenario")
estplots(scenario_simdat$ests,scenario_sumstats,save.figs=T,prefix="plots/scenario")

scenario_params %>% filter(replicate==1) %>% select(scenario,phi.1,phi.2,om.1,om.2)

# For reporting in the text:
# Scenario 1 and 16 are always labelled as the scenarios with the smallest and largest gaps in non-SSP,system
scenario_sumstats$by.AS.long %>% filter(name=="Non-SSP, system") %>% 
  group_by(scenario,group) %>% summarize(median=quantile(value,p=0.5)) %>%
  filter(scenario%in%c(1,16))

scenario_simdat$ests %>% group_by(scenario) %>% 
  summarize(true.diff=quantile(true.diff,0.5),
            est.diff=quantile(est.diff,0.5),
            true.patt=quantile(true.patt,0.5),
            true.satt=quantile(true.satt,0.5)) %>%
  filter(scenario%in%c(1,16))

#### Near and complete Violations of Assumption 10 ####
viol_10 <- data.frame(expand.grid(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= -.18, phi.2=-.06, 
  
  q= -1.28,
  om.1= -1,
  om.2= -1,
  om.3= -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1 = -.03, 
  psi.2 = 0,
  psi.3 = 0,
  
  theta.P = -68.5, sigma.P =0, 
  gamma.3= -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=c(0,200),
  gamma.7=c(0,200),
  
  beta.0 = -3.05,
  beta.3=0,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:4,NREP))

viol10_simdat <- simanalyze(viol_10)
viol10_sumstats <- sumstats(viol10_simdat$data, scenarios = 4)
datplots(viol10_sumstats,save.figs=T,prefix="plots/viol_10")
estplots(viol10_simdat$ests,viol10_sumstats,save.figs=T,prefix="plots/viol_10")

#### Near and complete Violations of Assumption 11 through U ####
viol_11 <- data.frame(expand.grid(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= -.18, phi.2=-.06, 
  
  q= -1.28,
  om.1= -1,
  om.2= -1,
  om.3= -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1 = c(0, .03,2), 
  psi.2 = 0,
  psi.3 = 0,
  
  theta.P = -68.5, sigma.P =0, 
  gamma.3= -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=0,
  gamma.7=0,
  
  beta.0 = -3.05,
  beta.3=.5,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.1b = 48000, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:3,NREP))

violation_simdat <- simanalyze(viol_11)
violation_sumstats <- sumstats(violation_simdat$data)
datplots(violation_sumstats,save.figs=T,prefix="plots/viol_11")
estplots(violation_simdat$ests,violation_sumstats,save.figs=T,prefix="plots/viol_11")

# Joint distribution of W,U by {S} across replicates and simulations
# ggplot(violation_simdat$data,aes(x=U,group=interaction(S,W))) + 
#   geom_density(aes(col=interaction(S,W))) + facet_wrap(S~scenario,ncol=3)
# baseline value =-0.03, .1 = little overlap, .2 = no overlap


#### Near and complete Violations of Assumption 10c ####
#### Near and complete Violations of Assumption 11 through U ####
viol_10c <- data.frame(expand.grid(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= -.18, phi.2=-.06, 
  
  q= -1.28,
  om.1= -1,
  om.2= -1,
  om.3= -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1 = 0, 
  psi.2 = 0,
  psi.3 = 0,
  
  theta.P = -68.5, sigma.P =0, 
  gamma.3= -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=0,
  gamma.7=0,
  gamma.8=c(0,50,100),
  
  beta.0 = -3.05,
  beta.3=.5,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.1b = 48000, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:3,NREP))

violation_simdat <- simanalyze(viol_10c)
violation_sumstats <- sumstats(violation_simdat$data)
datplots(violation_sumstats,save.figs=T,prefix="plots/viol_10c")
estplots(violation_simdat$ests,violation_sumstats,save.figs=T,prefix="plots/viol_10c")
