datplots <- function(sumstats,save.figs,prefix,order){
  theme_set(theme_minimal())
  sumstats$by.AS.long <- sumstats$by.AS.long %>%
    mutate(scenario=factor(scenario,levels=order,labels=1:length(order)))
  
  sumstats$diff.by.AS.long <- sumstats$diff.by.AS.long %>%
    mutate(scenario=factor(scenario,levels=order,labels=1:length(order)))
  
  ## Practice characteristics in sample and target
  panelA <- ggplot(filter(sumstats$by.AS.long,trt=="Treated"),aes(x=value,y=name)) + 
    geom_boxplot(aes(col=group),position=position_dodge(width=1)) + facet_wrap(~scenario) +
    labs(x="Percent",y="",color="") + theme(legend.position="bottom")
  if (save.figs) ggsave(paste0(prefix,"_","W_dist.png"),panelA,width=8,height=6) else print(panelA)
  
  # Difference in practice characteristics between sample and target
  panelB <- ggplot(filter(sumstats$diff.by.AS.long),aes(y=diff,x=scenario,group=scenario)) + 
    geom_boxplot() + geom_hline(yintercept=0) + facet_wrap(~name) +
    labs(y="Difference in Percent (target - sample)",x="Simulation Scenario")
  if (save.figs) ggsave(paste0(prefix,"_","W_diff_dist_by_W.png"),width=8,height=6) else print(panelB)
}

estplots <- function(ests,sumstats,save.figs,prefix,order){
  theme_set(theme_minimal())
  ests <- ests %>% 
    mutate(scenario=factor(scenario,levels=order,labels=1:length(order)))
  
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
}

calibplots <- function(sumstats,save.figs,which.base){
  by.S.targets <- tibble(S=c('Target','Sample'),
                         'Black'=c(.146,.12),
                         'SSP'=c(.35,.31),
                         'sys'=c(.33,.316),
                         'A'=c(.1,.18)) %>% pivot_longer(2:5)
  panelA <- ggplot(sumstats$by.S %>% filter(scenario==which.base) %>% pivot_longer(4:7),
                   aes(x=name,y=value,group=interaction(S,name))) + 
    geom_boxplot(aes(col=S),position=position_dodge(width=1)) +
    geom_point(data=by.S.targets,aes(x=name,y=value,col=S,group=interaction(S,name)),
               shape=8,size=2,position=position_dodge(width=1)) + 
    #scale_y_continuous(limits=c(0,1)) + 
    labs(x="",y="Proportion",color="") + theme(legend.position="bottom")
  if(save.figs) ggsave("plots/calib_by_S.png",width=6,height=4) else print(panelA)
  
  by.X.targets <- tibble(X=c('non-SSP','SSP','independent','system'),
                         A=c(.141,.269,.121,.309))
  panelB <- ggplot(sumstats$by.X %>% filter(scenario==which.base) ,aes(x=X,y=A)) + 
    geom_boxplot(aes(col=X)) +
    geom_point(data=by.X.targets,aes(x=X,y=A,col=X,group=X),shape=8,size=2) +
    scale_y_continuous(limits=c(0,1)) + 
    scale_color_discrete(guide="none") + labs(x="",y="Proportion treated")
  if (save.figs) ggsave("plots/calib_by_X.png",width=6,height=4) else print(panelB)
  
  by.A.targets <- tibble(A=c('Untreated','Treated'),
                         Black=c(.131,.069))
  panelC <- ggplot(sumstats$by.A %>% filter(scenario==which.base),aes(x=A,y=Black)) + 
    geom_boxplot(aes(col=A)) +
    geom_point(data=by.A.targets,aes(x=A,y=Black,col=A),shape=8,size=2) +
    scale_y_continuous(limits=c(0,1)) + 
    scale_color_discrete(guide="none") + labs(x="",y="Proportion Black",col="")
  if (save.figs) ggsave("plots/calib_by_A.png",width=6,height=4) else print(panelC)
}
