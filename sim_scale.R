## Making figures from simulation output
library(tidyverse)
source("HSR_simulation_fn.R")
theme_set(theme_minimal())
# Flag for saving figures to disk
save.figs <- TRUE

## Combine params that vary and params that are fixed
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
  
  theta.P = -68.5, sigma.P =69, 
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
  
  P= 100, # number of practices in a region
  replicate=1:50)) %>%
  mutate(scenario=rep(1:16,50))

# What are the scenarios?
if (save.figs) write_csv(global_params%>% filter(replicate==1) %>% select(scenario,phi.1,phi.2,om.1,om.2),
                         file="param_combos.csv")

## Simulate the data
simdat <- global_params %>% group_by(scenario,replicate) %>%
  nest() %>% mutate(simdat=map(data,make_regions)) %>%
  unnest(cols=c(data,simdat)) %>% 
  mutate(group=factor(S,levels=c('1','0'),labels=c('Sample','Target')),
         trt=factor(A,levels=c('0','1'),labels=c('Untreated','Treated'))) %>%
  group_by(scenario,replicate,S) %>%
  mutate(total.Black=sum(B)) %>% ungroup() %>%
  mutate(wt=B/total.Black) # proportion of all Black benes w/in {S}

sumstats.AS <- simdat %>%
  group_by(scenario,replicate,group,trt) %>%
  # Weighted by their importance to Black benes
  summarize(pct.Black=weighted.mean(b,w=wt),
            non.indep=weighted.mean(W==1,w=wt),
            non.sys=weighted.mean(W==2,w=wt),
            ssp.indep=weighted.mean(W==3,w=wt),
            ssp.sys=weighted.mean(W==4,w=wt))
sumstats.AS.long <- sumstats.AS %>% pivot_longer(5:9) %>% 
  mutate(name=factor(name,levels=c('pct.Black',
                                   'non.indep','non.sys','ssp.indep','ssp.sys'),
                     labels=c('% Black',
                              'Non-SSP, independent','Non-SSP, system',
                              'SSP, independent','SSP, system')))
diffstats.AS.long <- filter(sumstats.AS.long,trt=="Treated") %>% pivot_wider(names_from=group,values_from=value) %>% 
  mutate(diff=Target-Sample) %>% select(scenario,replicate,name,diff)

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