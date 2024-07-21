simanalyze <- function(params){
  sim_data <- params %>% 
    group_by(scenario,replicate) %>%
    nest() %>% mutate(simdat=map(data,make_regions)) %>%
    unnest(cols=c(data,simdat)) %>% 
    mutate(group=factor(S,levels=c('1','0'),labels=c('Sample','Target')),
           trt=factor(A,levels=c('0','1'),labels=c('Untreated','Treated'))) %>%
    group_by(scenario,replicate,S) %>%
    mutate(total.Black=sum(B)) %>% ungroup() %>%
    mutate(wt=B/total.Black) # proportion of all Black benes w/in {S}
  
  ests <- sim_data %>% 
    group_by(scenario,replicate) %>%
    nest() %>% mutate(patt=map(data,estimate_patt),
                      true.patt=map(data,true_patt),
                      satt=map(data,lm.did)) %>%
    unnest(cols=c(patt,true.patt,satt)) %>% select(-data) %>%
    mutate(dr.patt.err=dr-true.patt,
           satt.err=est.satt-true.satt,
           true.diff=true.patt-true.satt,
           est.diff=dr - est.satt)
  
  return(list('data'=sim_data,'ests'=ests))
}

sumstats <- function(data, scenarios = 0){
  # Calibration targets
  # Pr(B|S), Pre(SSP|S), Pr(sys|S), Pr(A|S)
  by.S <- data %>% group_by(scenario,replicate,S) %>%
    # unweighted (practice-level) stats
    summarize(Black=mean(b),
              SSP=mean(X1),
              sys=mean(X2),
              A=mean(A)) %>%
    mutate(S=factor(S,levels=c(1,0),labels=c('Sample','Target')))
  # Pr(A|SSP=0,S=1), Pr(A|SSP=1,S=1), Pr(A|sys=0,S=1) , Pr(A=1,sys=1,S=1)
  by.X <- data %>% filter(S==1) %>% group_by(scenario,replicate,X1) %>%
    summarize(A=mean(A)) %>% rename(X=X1) %>% ungroup() %>%
    mutate(X=factor(X,levels=c(0,1),labels=c('non-SSP','SSP'))) %>%
    bind_rows(data %>% filter(S==1) %>% group_by(scenario,replicate,X2) %>%
                summarize(A=mean(A)) %>% rename(X=X2) %>% ungroup() %>%
                mutate(X=factor(X,levels=c(0,1),labels=c('independent','system')))) %>%
    mutate(X=factor(X,levels=c('non-SSP','SSP','independent','system')))
  # Pr(B|S=1,A)
  by.A <- data %>% filter(S==1) %>% group_by(scenario,replicate,A) %>%
    summarize(Black=mean(b)) %>%
    mutate(A=factor(A,levels=c(0,1),labels=c('Untreated','Treated')))
  
  # W by {AxS}
  by.AS <- data %>%
    group_by(scenario,replicate,group,trt)  %>%
    # Weighted by their importance to Black benes
    summarize(non.indep=weighted.mean(W==1,w=wt)*100,
              non.sys=weighted.mean(W==2,w=wt)*100,
              ssp.indep=weighted.mean(W==3,w=wt)*100,
              ssp.sys=weighted.mean(W==4,w=wt)*100)
  by.AS.long <- by.AS %>% pivot_longer(5:8) %>% 
    mutate(name=factor(name,levels=c('non.indep','non.sys','ssp.indep','ssp.sys'),
                       labels=c('Non-SSP, independent','Non-SSP, system',
                                'SSP, independent','SSP, system')))
  # Differences between Target and Sample
  diff.by.AS.long <- filter(by.AS.long,trt=="Treated") %>% 
    pivot_wider(names_from=group,values_from=value) %>% 
    mutate(diff=Target-Sample) %>% 
    select(scenario,replicate,name,diff)
  
  if(scenarios == 0){
    # Arrange scenarios by median difference in proportion of Non-SSP, system practices:
    scenario.order <- (diff.by.AS.long %>% filter(name=="Non-SSP, system") %>% 
                         group_by(scenario) %>% summarize(median=quantile(diff,p=0.5)) %>% 
                         arrange(median))$scenario
  }else{
    scenario.order <- c(1:scenarios)
  }
  
  by.S <- by.S %>%
    mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order)))
  by.X <- by.X %>%
    mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order)))
  by.A <- by.A %>% 
    mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order)))
  by.AS <- by.AS %>% 
    mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order)))
  by.AS.long <- by.AS.long %>%
    mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order)))
  diff.by.AS.long <- diff.by.AS.long %>%
    mutate(scenario=factor(scenario,levels=scenario.order,labels=1:length(scenario.order)))
  
  return(list('by.S'=by.S,'by.X'=by.X,'by.A'=by.A,
              'by.AS'=by.AS,'by.AS.long'=by.AS.long,'diff.by.AS.long'=diff.by.AS.long,
              'scenario.order'=scenario.order))
}

calibplots <- function(sumstats,save.figs){
  by.S.targets <- tibble(S=c('Target','Sample'),
                         'Black'=c(.146,.12),
                         'SSP'=c(.35,.31),
                         'sys'=c(.33,.316),
                         'A'=c(.1,.18)) %>% pivot_longer(2:5)
  panelA <- ggplot(sumstats$by.S %>% pivot_longer(4:7),aes(x=name,y=value,group=interaction(S,name))) + 
    geom_boxplot(aes(col=S),position=position_dodge(width=1)) +
    geom_point(data=by.S.targets,aes(x=name,y=value,col=S,group=interaction(S,name)),
               shape=8,size=2,position=position_dodge(width=1)) + 
    #scale_y_continuous(limits=c(0,1)) + 
    labs(x="",y="Proportion",color="") + theme(legend.position="bottom")
  if(save.figs) ggsave("plots/calib_by_S.png",width=6,height=4) else print(panelA)
  
  by.X.targets <- tibble(X=c('non-SSP','SSP','independent','system'),
                         A=c(.141,.269,.121,.309))
  panelB <- ggplot(sumstats$by.X,aes(x=X,y=A)) + geom_boxplot(aes(col=X)) +
    geom_point(data=by.X.targets,aes(x=X,y=A,col=X,group=X),shape=8,size=2) +
    #scale_y_continuous(limits=c(0,1)) + 
    scale_color_discrete(guide="none") + labs(x="",y="Proportion treated")
  if (save.figs) ggsave("plots/calib_by_X.png",width=6,height=4) else print(panelB)
  
  by.A.targets <- tibble(A=c('Untreated','Treated'),
                         Black=c(.131,.069))
  panelC <- ggplot(sumstats$by.A,aes(x=A,y=Black)) + geom_boxplot(aes(col=A)) +
    geom_point(data=by.A.targets,aes(x=A,y=Black,col=A),shape=8,size=2) +
    #scale_y_continuous(limits=c(0,1)) + 
    scale_color_discrete(guide="none") + labs(x="",y="Proportion Black",col="")
  if (save.figs) ggsave("plots/calib_by_A.png",width=6,height=4) else print(panelC)
}

datplots <- function(sumstats,save.figs,prefix){
  theme_set(theme_minimal())
  ## Practice characteristics in sample and target
  panelA <- ggplot(filter(sumstats$by.AS.long,trt=="Treated"),aes(x=value,y=name)) + 
    geom_boxplot(aes(col=group),position=position_dodge(width=1)) + facet_wrap(~scenario) +
    labs(x="Percent",y="",color="") + theme(legend.position="bottom")
  if (save.figs) ggsave(paste0(prefix,"_","W_dist.png"),panelA,width=8,height=6) else print(panelA)
  
  # Difference in practice characteristics between sample and target
  panelB <- ggplot(filter(sumstats$diff.by.AS.long),aes(y=diff,x=scenario)) + 
    geom_boxplot() + geom_hline(yintercept=0) + facet_wrap(~name) +
    labs(y="Difference in Percent (target - sample)",x="Simulation Scenario")
  if (save.figs) ggsave(paste0(prefix,"_","W_diff_dist_by_W.png"),width=8,height=6) else print(panelB)
}

estplots <- function(ests,sumstats,save.figs,prefix){
  theme_set(theme_minimal())
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
}
