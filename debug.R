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



## Question: Now including alpha.1b != alpha.1 and beta.3 = .5, why is SATT error ~600? 
# Step 1: checking breakdown of A,S
viol_11_db <- data.frame(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1= -.18, phi.2=-.06, 
  
  q= -1.28,
  om.1= -1,
  om.2= -1,
  om.3= -.235,
  
  H = 0, sigma.H = 0.02,
  psi.1 = -.03, 
  psi.2 = 1,
  psi.3 = 1,
  
  theta.P = -68.5, sigma.P =0, 
  gamma.3= -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=0,
  gamma.7=0,
  
  beta.0 = -3.05,
  beta.3=.5,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.1b = 48000, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200)

# Looks to be similar; there are A=0,1 in S=0,1
df1 <- make_regions(viol_11_db)
table(df1$A,df1$S)
df1 %>% filter(S==0) %>% with(plot(U,A))
df1 %>% filter(W==1) %>% with(plot(U,A))

df2 <- make_regions(global_params)
table(df2$A,df2$S)

lm <- lm.did(df1)
lm

## Paring down U to be minimum to satisfy assumptions

U_mod_params <- data.frame(expand.grid(
  x1.r = -.617, 
  x2.r = -.715,
  phi.1 = -.18, phi.2=-.06, 
  
  q = -1.28,
  om.1 = -1,
  om.2 = -1,
  om.3 = -.235,
  
  H = 0, sigma.H = 0.2,
  psi.1 = c(0,.3,2), 
  psi.2 =0,
  psi.3 = 0,
  
  theta.P = -68.5, sigma.P = 0, 
  gamma.3 = -92, 
  gamma.4 = 206, 
  gamma.5 = -87,
  gamma.6=0,
  gamma.7=0,
  
  beta.0 = -3.05,
  beta.3 = .5,
  beta.4 = .839,
  beta.5 = 1.17,
  beta.6 = .778, 
  
  alpha.0 = 10100,
  alpha.1 = 46500, 
  alpha.2 = 600,
  alpha.3 = 12,
  
  P= 1200,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:3,NREP))

U_mod_simdat <- simanalyze(U_mod_params)
U_mod_sumstats <- sumstats(U_mod_simdat$data)
calibplots(U_mod_sumstats,save.figs=T)
datplots(U_mod_sumstats,save.figs=T,prefix="plots/U_mod")
estplots(U_mod_simdat$ests,U_mod_sumstats,save.figs=T,prefix="plots/U_mod")


## Illustrating how psi.1 allows us to control distributional overlap of U between S=0,1

viol_11 <- U_mod_params %>% filter(scenario == 1, replicate == 1)
viol_11$psi.1 <- 3
df <- make_regions(viol_11)
df <- tibble::as_tibble(df)

# Looking at U between S=0, 1
p <- ggplot( df, aes(x=U, color=factor(S, levels = c("0","1")),fill=factor(S, levels = c("0","1")))) + 
  geom_density(alpha=.3)
p

