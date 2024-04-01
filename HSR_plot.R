## Making figures from simulation output
library(tidyverse)
source("HSR_simulation_fn.R")
theme_set(theme_minimal())

## Combine params that vary and params that are fixed
global_params <- data.frame(expand.grid(
  x1.r = -.64, 
  x2.r = -.74,
  phi.1= c(.15,.19), # log odds ratio of Pr(SSP) in CPC+ vs non-CPC+
  phi.2= c(.03,.07), # log odds ratio of Pr(system) in CPC+ vs non-CPC+
  
  q= -2.12,
  om.1= c(-.165,.149), #log odds ratio of Pr(Black) in SSP vs non-SSP
  om.2= c(-.097,.138), #log odds ratio of Pr(Black) in system vs indep
  om.3= -.3,
  
  H = 0, sigma.H = .02,
  psi.1=-.03, 
  
  theta.P = -76, sigma.P = 69, 
  gamma.3=  -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  
  beta.0 = -2.51,
  #beta.1=0,
  #beta.2=0,
  beta.3 = 0,
  beta.4 = .809,
  beta.5 = 1.17,
  beta.6 = .46, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 50, # number of practices in a region
  replicate=1:50)) %>%
mutate(scenario=rep(1:16,50))

# What are the scenarios?
global_params%>% filter(replicate==1) %>% select(scenario,phi.1,phi.2,om.1,om.2)

## Simulate the data
simdat <- global_params %>% group_by(scenario,replicate) %>%
  nest() %>% mutate(simdat=map(data,make_regions)) %>%
  unnest(cols=c(data,simdat)) %>% 
  mutate(group=factor(S,levels=c('1','0'),labels=c('Sample','Target')),
         trt=factor(A,levels=c('0','1'),labels=c('Untreated','Treated'))) %>%
  group_by(scenario,replicate,S) %>%
  mutate(total.Black=sum(B)) %>% ungroup() %>%
  mutate(wt=B/total.Black) # proportion of all Black benes w/in {S}

## Figure 1: simulated practice summary statistics, comparing target CPC+ versus sample CPC+,
# Weighted by their importance to Black benes
## (weighted) avg values of W and B w/in {A x S} 
sumstats.AS <- simdat %>%
  group_by(scenario,replicate,group,trt) %>%
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
panelA <- ggplot(filter(sumstats.AS.long,trt=="Treated"),aes(x=value,y=name)) + 
  geom_boxplot(aes(col=group),position=position_dodge(width=1)) + facet_wrap(~scenario) +
  scale_shape_manual(values=c('Sample'=3,'Target'=19)) +
  scale_x_continuous("Mean") + ylab("") + theme(legend.position="bottom") +
  ggtitle('CPC+ practices, weighted')
panelA

# Re-express as the distribution of the difference between sample and target CPC+ practices:
# no visual subtraction needed
diffstats.AS.long <- filter(sumstats.AS.long,trt=="Treated") %>% pivot_wider(names_from=group,values_from=value) %>% 
  mutate(diff=Target-Sample) %>% select(scenario,replicate,name,diff)
panelB <- ggplot(filter(diffstats.AS.long),aes(x=diff,y=name)) + 
  geom_boxplot() + facet_wrap(~scenario) +
  scale_shape_manual(values=c('Sample'=3,'Target'=19)) +
  scale_x_continuous("Difference (target - sample)") + ylab("") + theme(legend.position="bottom") +
  ggtitle('CPC+ practices, weighted') + geom_vline(xintercept =0)
panelB


## Avg values of A and B across groups of {W X S}
sumstats.WS <- simdat %>%
  group_by(scenario,replicate,group,W) %>%
  summarize(pct.Black=mean(b),
            pct.CPC=mean(A)) %>%
  mutate(W=factor(W,levels=1:4,labels=c( 'Non-SSP, independent','Non-SSP, system',
                                         'SSP, independent','SSP, system')))
sumstats.WS.long <- sumstats.WS %>% pivot_longer(5:6) %>% 
  mutate(name=factor(name,levels=c('pct.Black','pct.CPC'),
                     labels=c('% Black','CPC+ participation')))
ggplot(sumstats.WS.long,aes(x=value,y=name)) + 
  geom_boxplot(aes(col=interaction(W,group)),position=position_dodge(width=1)) + facet_wrap(~scenario) +
  scale_x_continuous("Mean") + ylab("") + theme(legend.position="bottom") +
  ggtitle('All practices, unweighted')


### Figure 2: showing the PATT estimates across scenarios (and replications)
## Still TBD

