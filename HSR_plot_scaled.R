## Making figures from simulation output
library(tidyverse)
source("HSR_plot_scaled_fn.R")
load("result_data/2024-07-23_HSR_results.RData")
load("result_data/2024-07-23_HSR_results_viol10.RData")

scenario.output <- output

scenario.order <- (scenario.output$sumstats$diff.by.AS.long %>% filter(name=="Non-SSP, system") %>%
  group_by(scenario) %>% summarize(median=quantile(diff,p=0.5)) %>% arrange(median))$scenario
# Find the new number of the base scenario

viol_10.scenario.order <- c(1, 2, 4, 3, 5, 6, 7, 8, 9)

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
scenario.output$sumstats$by.AS.long %>% mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  filter(name=="Non-SSP, system",trt=="Treated") %>%
  group_by(scenario,group) %>% summarize(median=quantile(value,p=0.5)) %>%
  filter(scenario%in%c(1,16))

scenario.output$sumstats$diff.by.AS.long %>% 
  mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  filter(scenario==8) %>%
  group_by(name) %>% summarize(median=quantile(diff,p=0.5))

scenario.output$ests %>% mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  group_by(scenario) %>% 
  summarize(true.diff=quantile(true.diff,0.5),
            est.diff=quantile(est.diff,0.5),
            true.patt=quantile(true.patt,0.5),
            true.satt=quantile(true.satt,0.5)) %>%
  filter(scenario%in%c(1,8,16))

# Put the parameters in the right order:
scenario_params <- read_csv("scenario_params_varying.csv") %>%
  mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  arrange(scenario)
write_csv(scenario_params,"scenario_params_sorted.csv")

viol_params <- read_csv("viol_10_params_varying.csv") %>%
  mutate(scenario=factor(scenario,levels=viol_10.scenario.order,labels=1:length(viol_10.scenario.order))) %>%
  arrange(scenario)
write_csv(viol_params,"viol_10_params_sorted.csv")
