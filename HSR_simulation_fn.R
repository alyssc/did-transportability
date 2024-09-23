library(tidyverse)

# Default global params

# Defining functions for simulation study

## DATA GENERATING MECHANISMS ##

# Generate region parameters
# Including region_id, S
get_region_params <- function(id, global_params, in_sample = 0){
  region_params <- data.frame(
      region_id = id,
      S = in_sample
    )
  return(region_params)
} 

inv.logit <- function(p){
  return(exp(p)/(1+exp(p)))
}

rbern <- function(n, prob){
  rbinom(n,size=1,prob=prob)
}

#' Generate simulation of universe of regions and practices
make_regions <- function(global_params){
  
  # some parameters are new and may not be included in previously written params
  if(is.null(global_params$gamma.6)){
    global_params$gamma.6 = 0
  }
  if(is.null(global_params$gamma.7)){
    global_params$gamma.7 = 0
  }
  if(is.null(global_params$psi.2)){
    global_params$psi.2 = 0
  }
  if(is.null(global_params$psi.3)){
    global_params$psi.3 = 0
  }
  if(is.null(global_params$alpha.1b)){
    global_params$alpha.1b = global_params$alpha.1
  }

  n_regions <- 50
  n_sregions <- 18 # number of CPC+ regions
  
  df <- data.frame(region_id = integer(),
                   S = integer(),
                   # practice_id added later
                   b = numeric(),
                   B = numeric(),
                   W = integer(),
                   X1 = integer(),
                   X2 = integer(),
                   U = numeric(),
                   delta = numeric(),
                   A = integer(), 
                   Yb.pre = numeric(),
                   Yb.post = numeric(),
                   Yb.post0 = numeric(),
                   Yb.post1 = numeric())
                           
  # Add practices to dataframe, region by region
  for(region_id in 1:n_regions){
    S = 0
    
    # Adjust parameters if region is CPC+ region
    if(region_id <= n_sregions)
    {
      S = 1
    }
    
    params <- cbind(global_params, get_region_params(region_id, global_params, S))
    
    P <- params$P
    n <- rbinom(P, 840, .81304) #patients in practices
    
    X1 <- rbern(n=P, prob=inv.logit(params$x1.r+params$phi.1*S) )
    X2 <- rbern(n=P, prob=inv.logit(params$x2.r+params$phi.2*S) )
    W <- 1*(1-X1)*(1-X2) + 2*(1-X1)*X2 + 3*(1-X2)*X1 + 4*X1*X2
    
    B <- rbinom(P, n, inv.logit(params$q + params$om.1*X1 + params$om.2*X2 + params$om.3*S))
    b <- B/n
    
    U <- rnorm(P, mean = params$H + params$psi.1*S + params$psi.2*X1 + params$psi.3*X2 
               #+ params$psi.4*X1*S + params$psi.5*X2*S
               , sd = params$sigma.H) # unobserved H
    # U.post <- rnorm(P, mean = params$H + params$psi.2, sd = params$sigma.H) #post-period unobserved H
    
    # make sure at least within range of ATT from JAMA paper
    delta <- rnorm(P, params$theta.P + params$gamma.3*X1 + params$gamma.4*X2 + params$gamma.5*X1*X2 + 
                     params$gamma.6*X1*(S-1) + params$gamma.7*X2*(S-1), sd = params$sigma.P)
    
    betas <- c(params$beta.0, params$beta.3, params$beta.4, params$beta.5, params$beta.6)
    trt.prob <- inv.logit(betas %*% rbind(1, U, X1, X2, S))
    A <- rbern(n = P, prob = trt.prob)
    
    Yb.pre <-   params$alpha.0 + params$alpha.1*U + params$alpha.2 * X1 + params$alpha.3 * X2
    Yb.post <-  params$alpha.0 + params$alpha.1b*U + params$alpha.2 * X1 + params$alpha.3 * X2 + delta * A # observed outcome
    Yb.post0 <- params$alpha.0 + params$alpha.1b*U + params$alpha.2 * X1 + params$alpha.3 * X2 # untreated potential outcome
    Yb.post1 <- params$alpha.0 + params$alpha.1b*U + params$alpha.2 * X1 + params$alpha.3 * X2 + delta # treated potential outcome

    df <- rbind(df, data.frame(region_id, S, b, B, W,X1,X2, U, delta, A, Yb.pre, Yb.post, Yb.post0, Yb.post1))
  }

  # Practice-level weights (proportion of Black benes in {S} cared for by each practice)
  B0 <- sum(df$B[df$S==0])
  B1 <- sum(df$B[df$S==1])
  # For sample practices
  df$wt <- (df$B/B1)*as.numeric(df$S==1) + (df$B/B0)*as.numeric(df$S==0)
  
  df$W <- factor(df$W, levels = c(1:4))
  #colnames(df)[2] <- "S"
  return(df)
  
}


## ESTIMATION ##


## Returns vector of three estimates from G comp, IPW, DR estimators
estimate_patt <- function(df){
  m <- lm(Yb.post-Yb.pre ~ W*A*S, data = df) 
  m0 <- predict(m, newdata=df %>% mutate(A=0,S=1))
  m1 <- predict(m, newdata=df %>% mutate(A=1,S=1))
  
  I10 = 1*(df$A==1 & df$S==0)
  I11 = 1*(df$A==1 & df$S==1)
  I01 = 1*(df$A==0 & df$S==1)
  
  p10 = weighted.mean(df$A*(1-df$S),w=df$wt)
  
  gT <- lm(A ~ W*S, data = df)
  gS <- lm(S ~ W, data = df)
  
  g11 <- predict(gT, newdata=df %>% mutate(S=1)) * predict(gS, newdata=df)
  g01 <- (1-predict(gT, newdata=df %>% mutate(S=1))) * predict(gS, newdata=df)
  g10 <- predict(gT, newdata=df %>% mutate(S=0)) * (1-predict(gS, newdata=df))
  
  gcomp.val <- weighted.mean(I10/p10*(m1-m0),w=df$wt)
  # This isn't right yet: 
  ipw.val <- weighted.mean((I11*g10*(df$Yb.post-df$Yb.pre)/(p10*g11)) - 
                           (I01*g10*(df$Yb.post-df$Yb.pre)/(p10*g01)),w=df$wt)
  
  dr.val <- weighted.mean(I11*g10*((df$Yb.post-df$Yb.pre) - m1)/(p10*g11) -
               I01*g10*((df$Yb.post-df$Yb.pre) - m0)/(p10*g01) +
               I10*(m1 - m0)/p10,w=df$wt)
  
  return(tibble('gc'=gcomp.val, 'ipw'=ipw.val,'dr'= dr.val))
}

## True PATT
true_patt <- function(df){
  I10 = 1*(df$A==1 & df$S==0)
  p10 = weighted.mean(df$A*(1-df$S),w=df$wt)
  return(weighted.mean(I10*(df$Yb.post1 - df$Yb.post0)/p10,w=df$wt))
}


## In-sample DiD
lm.did <- function(df, plot = F){
  df_long <- gather(df, time, Y, Yb.pre:Yb.post)
  df_long$period <- ifelse(df_long$time=="Yb.pre",0,1)
  
  est.satt <- summary(lm(Y ~ A*period, data = df_long %>% filter(S==1) ,weights=wt))$coef['A:period','Estimate']
  true.satt <- (df %>% filter(S==1,A==1) %>% mutate(delta=Yb.post1-Yb.post0) %>%
    summarize(mean=weighted.mean(delta,w=wt)))$mean
  
  return(tibble('est.satt'=est.satt,'true.satt'=true.satt))
}

# Plot in-sample DiD
plot.did <- function(df){
  df_long <- gather(df, time, Y, Yb.pre:Yb.post)
  df_long$time[df_long$time=="Yb.pre"] <- 0
  df_long$time[df_long$time=="Yb.post"] <- 1
  
  plot_data <- df_long %>% 
    mutate(A = factor(A, labels = c("non-CPC+ participant", "CPC+ participant")),
           time = factor(time, labels = c("Pre-period", "Post-period"))) %>% 
    group_by(A,time) %>% 
    summarize(mean_Y = mean(Y),
              se_Y = sd(Y) / sqrt(n()),
              upper = mean_Y + (-1.96 * se_Y),
              lower = mean_Y + (1.96 * se_Y)) 
  
  return(ggplot(plot_data, aes(x = time, y = mean_Y, color = A)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 1) + 
    geom_line(aes(group = A)))
}


# Calculate estimated and true PATT across variation of a single parameter
sim_patt <- function(default_params,var_name,var_seq,nsim){
  x <- var_seq
  results <- matrix(NA,nrow = 0, ncol = 5)
  
  for (i in 1:length(x)){
    default_params[var_name] <- x[i]
    y_singlesims <- matrix(NA,nrow = 0, ncol = 5)
    
    for (j in 1:nsim){
      
      df <- make_regions(default_params)
      y_singlesims <- rbind(y_singlesims, c(x[i],true_patt(df),estimate_patt(df)))
    }
    
    results <- rbind(results, colMeans(y_singlesims))
    
  }
  results <- as.data.frame(results)
  colnames(results) <- c("psi.1","true", "gcomp", "ipw", "dr")
  return(results)
}

# Plots the results of sim_patt
plot_patt <- function(results, var_name){
  results_long <- melt(results, id = var_name)  #requires reshape2, not included right now
  plot_patt <- ggplot(results_long,
                      aes(x = get(var_name), 
                          y = value, 
                          color = variable)) +  geom_line() 
  return(plot_patt)
}

sumstats <- function(data){
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

  return(list('by.S'=by.S,'by.X'=by.X,'by.A'=by.A,
              'by.AS'=by.AS,'by.AS.long'=by.AS.long,'diff.by.AS.long'=diff.by.AS.long))
}


