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

# Arrange scenarios by median difference in proportion of Non-SSP, system practices:
scenario.order <- (diffstats.AS.long %>% filter(name=="Non-SSP, system") %>% group_by(scenario) %>% summarize(median=quantile(diff,p=0.5)) %>% arrange(median))$scenario

# For reporting in the text:
View(sumstats.AS.long %>% filter(name=="Non-SSP, system") %>% 
  group_by(scenario,group) %>% summarize(median=quantile(value,p=0.5)) %>% arrange(factor(scenario,levels=scenario.order)))
View(estimates %>% group_by(scenario) %>% summarize(median=quantile(true.diff,0.5)) %>% arrange(factor(scenario,levels=scenario.order)))
View(estimates %>% group_by(scenario) %>% summarize(median=quantile(est.diff,0.5)) %>% arrange(factor(scenario,levels=scenario.order)))
View(estimates %>% group_by(scenario) %>% summarize(median=quantile(true.patt,0.5)) %>% arrange(factor(scenario,levels=scenario.order)))
View(estimates %>% group_by(scenario) %>% summarize(median=quantile(true.satt,0.5)) %>% arrange(factor(scenario,levels=scenario.order)))

## Practice characteristics in Sample and Target
ggplot(filter(sumstats.AS.long,trt=="Treated",name!="% Black"),aes(x=value,y=name)) + 
  geom_boxplot(aes(col=group),position=position_dodge(width=1)) + facet_wrap(~factor(scenario,levels=scenario.order)) +
  scale_x_continuous("Proportion") + ylab("") + theme(legend.position="bottom")
if (save.figs) ggsave("W_dist.png",width=8,height=6)

# Difference in practice characteristics between sample and target
ggplot(filter(diffstats.AS.long,name!="% Black"),
       aes(y=diff,x=factor(scenario,levels=scenario.order))) + 
  geom_boxplot() + facet_wrap(~name) +
  scale_y_continuous("Difference in Proportion (target - sample)") + xlab("Scenario") + geom_hline(yintercept =0)
if (save.figs) ggsave("W_diff_dist_by_W.png",width=8,height=6)

### Figure 2: showing the PATT estimates across scenarios (and replications)
# For plotting the error between the estimate and the truth (for PATT, SATT) 
err.plot.dat <- estimates %>% select(replicate,scenario,dr.patt.err,satt.err) %>% 
  pivot_longer(dr.patt.err:satt.err)
# For plotting the difference between the PATT and the SATT (both true and estimated)
diff.plot.dat <- estimates %>% select(replicate,scenario,true.diff,est.diff) %>%
  pivot_longer(true.diff:est.diff)
# For plotting the true PATT and SATT
true.plot.dat <- estimates %>% select(scenario,replicate,true.patt,true.satt) %>% 
  pivot_longer(true.patt:true.satt)

# SATT and PATT truths
ggplot(true.plot.dat,aes(x=factor(scenario,levels=scenario.order),y=value)) + geom_boxplot(aes(col=name)) + 
  geom_hline(yintercept=0) + ylab("True value") + xlab("Scenario") + 
  scale_color_discrete(labels=c('PATT','SATT')) + theme(legend.position="bottom")
if (save.figs) ggsave("true_PATT_SATT.png",width=6,height=4)

# Showing the estimation error
ggplot(err.plot.dat,aes(x=factor(scenario,levels=scenario.order),y=value)) + 
  geom_boxplot() + facet_wrap(~factor(name,levels=c('dr.patt.err','satt.err'),labels=c("DR PATT","OLS SATT"))) + 
  geom_hline(yintercept=0) + xlab("Scenario") + ylab("Error")
if (save.figs)  ggsave("est_vs_true.png",width=8,height=6)

# Showing the difference between target and sample (both estimated and true)
ggplot(diff.plot.dat,aes(x=factor(scenario,levels=scenario.order),y=value,group=interaction(name,scenario))) + 
  geom_boxplot(aes(col=name),position=position_dodge(width=.75)) + geom_hline(yintercept=0) +
  labs(x="Scenario",y="Target ATT - Sample ATT",color="") + scale_color_discrete(labels=c('Estimated','True'))
if (save.figs) ggsave("PATT_minus_SATT.png",width=8,height=6)
