## Making figures from simulation output
library(tidyverse)
source("HSR_plot_scaled_fn.R")
load("2024-04-13_HSR_results.RData")

scenario.order <- (output$sumstats$diff.by.AS.long %>% filter(name=="Non-SSP, system") %>%
  group_by(scenario) %>% summarize(median=quantile(diff,p=0.5)) %>% arrange(median))$scenario
# Find the new number of the base scenario:
which(scenario.order==1)

datplots(output$sumstats,save.figs=T,prefix="scenario",order=scenario.order)
estplots(output$ests,output$sumstats,save.figs=T,prefix="scenario",order=scenario.order)

# For reporting in the text:
output$sumstats$by.AS.long %>% mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  filter(name=="Non-SSP, system") %>%
  group_by(scenario,group) %>% summarize(median=quantile(value,p=0.5)) %>%
  filter(scenario%in%c(1,16))

output$ests %>% mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order))) %>%
  group_by(scenario) %>% 
  summarize(true.diff=quantile(true.diff,0.5),
            est.diff=quantile(est.diff,0.5),
            true.patt=quantile(true.patt,0.5),
            true.satt=quantile(true.satt,0.5)) %>%
  filter(scenario%in%c(1,16))
