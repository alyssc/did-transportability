## Making figures from simulation output
library(tidyverse)
source("HSR_simulation_fn.R")

## Combine params that vary and params that are fixed
global_params <- data.frame(
  x1.r = -.6, x2.r = -.7, phi.1=.17, phi.2=.04, 
  
  q= -2.12, om.1= .1, om.2= .6, om.3= c(-.3,0,.3),
  
  H = 0, sigma.H = .02, psi.1=-.03, 
  
  theta.P = -76, sigma.P = 69,gamma.3=-92,gamma.4 = 206, gamma.5 = -87,
  
  beta.0 = -2.51, beta.3=0, beta.4 = -.809,
  beta.5 = -1.17, beta.6 = c(-.46,0,.46), 
  
  alpha.0 = 10100, alpha.1 = 46500 , 
  alpha.2 = 600, alpha.3 = 12,
  
  P= 50
) %>%
mutate(scenario=1:9)

sumstats <- global_params %>% group_by(scenario) %>%
  nest() %>% mutate(simdat=map(data,make_regions)) %>%
  unnest(cols=c(data,simdat)) %>% 
  mutate(group=factor(S,levels=c('1','0'),labels=c('Sample','Target')),
  trt=factor(A,levels=c('0','1'),labels=c('Untreated','Treated'))) %>%
  group_by(scenario,group,trt) %>%
  summarize(pct.Black=mean(b),participation=mean(A),ssp.sys=mean(W==1),
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
