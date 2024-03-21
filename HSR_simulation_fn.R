library(boot)
library(Rlab)
library(tidyverse)

# Defining parameters and functions for simulation study

## DATA GENERATING MECHANISMS ##

global_params <- data.frame(theta.P = .2, sigma.P = .01, 
                            H = .2, sigma.H = .1,
                            x1.r = -.49, 
                            x2.r = -.5,
                            phi.1=.14, phi.2=.16, 
                            gamma.2 = 0,
                            gamma.3=.1, 
                            gamma.4 = .1, 
                            psi.1=.1, 
                            beta.0 = 2,
                            beta.1=.1,
                            beta.2=0,
                            beta.3=0,
                            beta.4 = .1,
                            beta.5 = .1,
                            alpha.1 = 50, 
                            alpha.2 = 70,
                            alpha.3 = 60,
                            om.1= .1,
                            om.2= .1,
                            om.3= -.13,
                            q= -.92,
                            P= 10 # number of practices in a region
                            
                            )

# Generate region parameters
# Including region_id, S
get_region_params <- function(id, global_params, in_sample = 0){
  region_params <- data.frame(
      region_id = id,
      S = in_sample
    )
  return(region_params)
} 


#' Generate simulation of universe of regions and practices
make_regions <- function(global_params){
  n_regions <- 50
  n_sregions <- 18 # number of CPC+ regions
  
  df <- data.frame(region_id = integer(),
                   S = integer(),
                   # practice_id added later
                   b = numeric(),
                   W = integer(),
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
    n <- rbinom(P, 1000, .3) #patients in practices
    
    X1 <- rbern(n=P, prob=inv.logit(params$x1.r+params$phi.1*S) )
    X2 <- rbern(n=P, prob=inv.logit(params$x2.r+params$phi.2*S) )
    W <- 1*(1-X1)*(1-X2) + 2*(1-X1)*X2 + 3*(1-X2)*X1 + 4*X1*X2
    
    B <- rbinom(P, n, inv.logit(params$q + params$om.1*X1 + params$om.2*X2 + params$om.3*S))
    b <- B/n
    
    U <- rnorm(P, mean = params$H, sd = params$sigma.H + params$psi.1*S) # unobserved H
    # U.post <- rnorm(P, mean = params$H + params$psi.2, sd = params$sigma.H) #post-period unobserved H
    
    # make sure at least within range of ATT from JAMA paper
    delta <- rnorm(P, params$theta.P + params$gamma.2*U + params$gamma.3*X1 + params$gamma.4*X2, sd = params$sigma.P)
    
    betas <- c(params$beta.0, params$beta.1, params$beta.2, params$beta.3, params$beta.4, params$beta.5)
    trt.prob <- inv.logit(betas %*% rbind(1, n, b, U, X1, X2))
    A <- rbern(n = P, prob = trt.prob)
    
    Yb.pre <- params$alpha.1*U + params$alpha.2 * X1 + params$alpha.3 * X2
    Yb.post <- params$alpha.1*U + params$alpha.2 * X1 + params$alpha.3 * X2 + delta * A # observed outcome
    Yb.post0 <- params$alpha.1*U + params$alpha.2 * X1 + params$alpha.3 * X2 # untreated potential outcome
    Yb.post1 <- params$alpha.1*U + params$alpha.2 * X1 + params$alpha.3 * X2 + delta # treated potential outcome

    df <- rbind(df, data.frame(region_id, S, b, W, U, delta, A, Yb.pre, Yb.post, Yb.post0, Yb.post1))
    
  }
  
  df$W <- factor(df$W, levels = c(1:4))
  colnames(df)[2] <- "S"
  return(df)
  
}

## Displaying and testing parameters
# In the form of P(B=1 | S=1) etc. 



## ESTIMATION ##


## Returns vector of three estimates from G comp, IPW, DR estimators
estimate_patt <- function(df){
  m <- lm(Yb.post-Yb.pre ~ W*A*S, data = df)
  m0 <- predict(m, newdata=df %>% mutate(A=0,S=1))
  m1 <- predict(m, newdata=df %>% mutate(A=1,S=1))
  
  I10 = 1*(df$A==1 & df$S==0)
  I11 = 1*(df$A==1 & df$S==1)
  I01 = 1*(df$A==0 & df$S==1)
  
  p10 = mean(df$A*(1-df$S))
  
  gT <- lm(A ~ W*S, data = df)
  gS <- lm(S ~ W, data = df)
  
  g11 <- predict(gT, newdata=df %>% mutate(S=1)) * predict(gS, newdata=df)
  g01 <- (1-predict(gT, newdata=df %>% mutate(S=1))) * predict(gS, newdata=df)
  g10 <- predict(gT, newdata=df %>% mutate(S=0)) * (1-predict(gS, newdata=df))
  
  gcomp.val <- mean(I10*(m1-m0)/p10)
  ipw.val <- mean(I11*g10*(df$Yb.post-df$Yb.pre)/(p10*g11) -
                I01*g10*(df$Yb.post-df$Yb.pre)/(p10*g01))
  dr.val <- mean(I11*g10*((df$Yb.post-df$Yb.pre) - m1)/(p10*g11) -
               I01*g10*((df$Yb.post-df$Yb.pre) - m0)/(p10*g01) +
               I10*(m1 - m0)/p10)
  
  return(c(gcomp.val, ipw.val, dr.val))
}

## True PATT 
true_patt <- function(df){
  I10 = 1*(df$T==1 & df$S==0)
  return(mean(I10*(df$Yb.post1 - df$Yb.post0)))
}


## In-sample DiD







# BELOW ARE DRAFT FUNCTIONS WHICH ARE NOT CURRENTLY CONSISTENT WITH ABOVE NOTATION

#' Get Average Treatment Effect of Treated given a region's practice-level attributes
#' @param region Region's practice-level attributes from make_region()
#' @param type "t", "b", or "w" for total ATT, ATT for black beneficiaries, and ATT for white beneficiaries
get_ATT <- function(region, type = "t"){
  if(type == "t"){
    return(mean(region$delta[region$treated == 1]))
  } else if(type == "b"){
    return(mean(region$delta.b[region$treated == 1]))
  } else if(type == "w"){
    return(mean(region$delta.w[region$treated == 1]))
  }
}

get_ATE <- function(region, type = "t"){
  if(type == "t"){
    return(mean(region$delta))
  } else if(type == "b"){
    return(mean(region$delta.b))
  } else if(type == "w"){
    return(mean(region$delta.w))
  }
}

# Plot ATT.B - ATE.B across variation of a single parameter
vary_param_plot <- function(default_params,var_name,var_seq,nsim){
  x <- var_seq
  y <- rep(NA, length(x))
  
  for (i in 1:length(x)){
    default_params[var_name] <- x[i]
    y_singlesims <- rep(NA, nsim)
    
    for (j in 1:nsim){
      region <- make_region(default_params)
      y_singlesims[j] <- get_ATT(region,"b") - get_ATE(region,"b")
    }
    
    y[i] <- mean(y_singlesims)
    
  }
  
  plot(x, y, type="b", xlab=paste("Values of ",var_name), 
       ylab="ATT-ATE for Black patients")
}


