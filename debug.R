test <- params %>% 
  group_by(scenario,replicate) %>%
  nest() %>% mutate(simdat=map(data,make_regions))


test$simdat

test1 <- sim_data %>% 
  group_by(scenario,replicate) %>%
  nest()
head(test1)
test2 <- test1 %>% mutate(patt=map(data,estimate_patt),
       true.patt=map(data,true_patt),
       satt=map(data,lm.did))

head(test2)

first_data <- test1$data[[1]]
estimate_patt(first_data)
true_patt(first_data)
lm.did(first_data)

head(first_data)
summary(first_data$Yb.post)



sumstats <- violation_sumstats

ests <-  ests %>%
  mutate(scenario=factor(scenario,levels=sumstats$scenario.order,
                         labels=1:length(sumstats$scenario.order)))
# SATT and PATT truths
true.plot.dat <- ests %>% select(scenario,replicate,true.patt,true.satt) %>% 
  pivot_longer(true.patt:true.satt)

panelA <- ggplot(true.plot.dat,aes(x=scenario,y=value)) + geom_boxplot(aes(col=name)) + 
  geom_hline(yintercept=0) + scale_color_discrete(labels=c('PATT','SATT')) + 
  labs(y="True ATT",x="Simulation Scenario",color="") + theme(legend.position="bottom")
if (save.figs) ggsave(paste0(prefix,"_","true_PATT_SATT.png"),width=6,height=4) else print(panelA)

# Error between the estimate and the truth (for PATT, SATT) 
err.plot.dat <- ests %>% select(replicate,scenario,dr.patt.err,satt.err) %>%
  pivot_longer(dr.patt.err:satt.err) %>%
  mutate(name=factor(name,levels=c('dr.patt.err','satt.err'),labels=c("DR PATT","OLS SATT")))

panelB <- ggplot(err.plot.dat,aes(x=scenario,y=value)) + geom_boxplot() + 
  facet_wrap(~name) + geom_hline(yintercept=0) + 
  labs(x="Simulation Scenario",y="Error (est - true)")
if (save.figs)  ggsave(paste0(prefix,"_","est_vs_true.png"),width=8,height=4) else print(panelB)

# Difference between target and sample (both estimated and true)
diff.plot.dat <- ests %>% select(replicate,scenario,true.diff,est.diff) %>%
  pivot_longer(true.diff:est.diff)

panelC <- ggplot(diff.plot.dat,aes(x=scenario,y=value,group=interaction(name,scenario))) + 
  geom_boxplot(aes(col=name),position=position_dodge(width=.75)) + geom_hline(yintercept=0) +
  labs(x="Simulation Scenario",y="PATT - SATT",color="") + scale_color_discrete(labels=c('Estimated','True')) +
  theme(legend.position="bottom")
if (save.figs) ggsave(paste0(prefix,"_","PATT_minus_SATT.png"),width=6,height=4) else print(panelC)
