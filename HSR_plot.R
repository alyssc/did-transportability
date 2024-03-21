## Making figures from simulation output
library(tidyverse)
source("HSR_simulation_fn.R")

## Combine params that vary and params that are fixed
global_params <- data.frame(expand.grid(
  theta.P = .2, sigma.P = .01, 
  H = 1, sigma.H = .1,
  x1.r = -.49, x2.r = -.5,
  phi.1=.14, phi.2=.16, 
  gamma.3=.1, gamma.4 = .1,
  beta.0 = -1.38, beta.3=0, beta.4 = .1, beta.5 = .1, 
  beta.6 = c(-1,0,1), # Vary the take-up in target versus sample regions
  alpha.0 = 10100, alpha.1 = 46500, alpha.2 = 600, alpha.3 = 12,
  om.1= .1, om.2= .1, om.3= -.13,
  q= -.92,
  P= 10)) %>%
mutate(scenario=1:9)

sumstats <- global_params %>% group_by(scenario) %>%
  nest() %>% mutate(simdat=map(data,make_regions)) %>%
  unnest(cols=c(data,simdat)) %>% 
  mutate(group=factor(S,levels=c('1','0'),labels=c('Sample','Target')),
  trt=factor(treated,levels=c('0','1'),labels=c('Untreated','Treated'))) %>%
  group_by(scenario,group,trt) %>%
  summarize(pct.Black=mean(b),participation=mean(treated),ssp.sys=mean(W==1),
            ssp.indep=mean(W==2),non.sys=mean(W==3),non.indep=mean(W==4))
sumstats.long <- sumstats %>% pivot_longer(4:9) %>% 
  mutate(name=factor(name,levels=c('pct.Black','participation',
                                   'ssp.sys','ssp.indep','non.sys','non.indep'),
                     labels=c('% Black','CPC+ participation',
                              'SSP, system','SSP, independent',
                              'Non-SSP, system','Non-SSP, independent')))
ggplot(sumstats.long,aes(x=value,y=name)) + 
  geom_point(aes(col=group,shape=trt),position=position_dodge(width=.5)) +
  xlab("Mean") + ylab("") + theme(legend.position="bottom")
