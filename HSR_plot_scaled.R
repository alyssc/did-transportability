## Making figures from simulation output
library(tidyverse)
source("HSR_plot_scaled_fn.R")

# Run all three together
load("2024-07-22_HSR_results.RData")
load("2024-07-22_HSR_results_viol10.RData")
scenario.output <- output

scenario.order <- (scenario.output$sumstats$diff.by.AS.long %>% filter(name=="Non-SSP, system") %>%
  group_by(scenario) %>% summarize(median=quantile(diff,p=0.5)) %>% arrange(median))$scenario
# Find the new number of the base scenario

viol_10.scenario.order <- c(1, 2, 4, 5, 3, 7, 6, 8, 9)

for(i in c(1:2)){
  if(i==1){
    output <- scenario.output
    prefix <- "plots/scenario"
    order <- scenario.order
  }else{
    output <- viol_10_output
    prefix <- "plots/viol_10"
    order <- viol_10.scenario.order
  }
  datplots(output$sumstats,save.figs=T,prefix=prefix,order=order)
  estplots(output$ests,output$sumstats,save.figs=T,prefix=prefix,order=order)
}

calibplots(scenario.output$sumstats,save.figs=T,which.base=which(order==1))


# For reporting in the text:
output$sumstats$by.AS.long %>% mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  filter(name=="Non-SSP, system",trt=="Treated") %>%
  group_by(scenario,group) %>% summarize(median=quantile(value,p=0.5)) %>%
  filter(scenario%in%c(1,16))

output$sumstats$diff.by.AS.long %>% 
  mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  filter(scenario==9) %>%
  group_by(name) %>% summarize(median=quantile(diff,p=0.5))

output$ests %>% mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  group_by(scenario) %>% 
  summarize(true.diff=quantile(true.diff,0.5),
            est.diff=quantile(est.diff,0.5),
            true.patt=quantile(true.patt,0.5),
            true.satt=quantile(true.satt,0.5)) %>%
  filter(scenario%in%c(1,9,16))

scenario_params <- data.frame(expand.grid(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= c(-0.188,0.188),  # log odds ratio of Pr(SSP) in CPC+ vs non-CPC+
  phi.2=c(-0.057,0.057),   # log odds ratio of Pr(system) in CPC+ vs non-CPC+
  
  q= -1.38,
  om.1 = c(-1.25,1.25),  # log odds ratio of Pr(Black) in SSP vs non-SSP
  om.2 = c(-1,1), # log odds ratio of Pr(Black) in system vs indep
  om.3 = -0.1389,
  
  H = 0, sigma.H = 0.02,
  psi.1 = -.03, 
  
  theta.P = -68.5, sigma.P = 69, 
  gamma.3 = -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  
  beta.0 = -3.063,
  beta.3 = 0,
  beta.4 = 1.02,
  beta.5 = 1.301,
  beta.6 = 0.6076, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200)) %>%
  mutate(scenario=c(1:16)) %>%
  mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order)))

write_csv(scenario_params %>% select(scenario,om.1,om.2,phi.1,phi.2),"scenario_params_varying.csv")



viol_10_params <- data.frame(expand.grid(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= -0.188,  # log odds ratio of Pr(SSP) in CPC+ vs non-CPC+
  phi.2= -0.057,   # log odds ratio of Pr(system) in CPC+ vs non-CPC+
  
  q= -1.38,
  om.1 = -1.25,  # log odds ratio of Pr(Black) in SSP vs non-SSP
  om.2 = -1, # log odds ratio of Pr(Black) in system vs indep
  om.3 = -0.1389,
  
  H = 0, sigma.H = 0.02,
  psi.1 = -.03, 
  
  theta.P = -68.5, sigma.P = 69, 
  gamma.3 = -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=c(0,50,200),
  gamma.7=c(0,50,200),
  
  beta.0 = -3.063,
  beta.3 = 0,
  beta.4 = 1.02,
  beta.5 = 1.301,
  beta.6 = 0.6076, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200)) %>%
  mutate(scenario=c(1:9)) %>%
  mutate(scenario=factor(scenario,levels=viol_10.scenario.order,labels=1:length(viol_10.scenario.order)))

write_csv(viol_10_params %>% select(scenario,gamma.6,gamma.7),"viol_10_params_varying.csv")
