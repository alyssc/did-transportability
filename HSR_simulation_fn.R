library(boot)
library(Rlab)
library(tidyverse)
library(simstudy)
library(furrr)

global_params <- data.frame(theta.P = .2, sigma.P = .01, 
                            H = 1, sigma.H = .4, b.overall = .14, 
                            x.baseline = 0)

region_params <- data.frame(
                            gamma.1=.1, 
                            gamma.2=.1, 
                            psi.1=-.1, 
                            psi.2=.1, 
                            beta.1=.1,
                            beta.2=0,
                            beta.3=.01,
                            beta.4 = 1, # relationship btwn proportion Black and treatment variable
                            beta.5 = 0,
                            beta.6 = 0,
                            alpha.1 = 50, 
                            alpha.2 = 70,
                            alpha = .5, beta = .5, # distribution parameters for proportion of Black beneficiaries in practices in CPC+ regions
                            w = 0, c = .2, m = .1,
                            meanlog = 5, sdlog = 1, # mean and SD of FFS Medicare population in practices in CPC+ regions
                            P = 100, # number of practices in a CPC+ region
                            H.r  = .9 # baseline risk score of Black beneficiaries in a practice in CPC+ regions 
)

# Generate region parameters
# Including region_id, S

#' Generate practices within a region
#' @param params A dataframe of parameters, see global_params for example. 
make_region <- function(params){
  P <- params$P
  b <- rbeta(P, params$alpha, params$beta)
  n <- round(rlnorm(P, params$meanlog, params$sdlog))
  
  H.b0 <- rnorm(P, mean = params$H.r + params$psi.1*b, 
                sd = params$sigma.H) #pre-period H
  H.b1 <- rnorm(P, mean = params$H.r + params$psi.1*b + params$psi.2, 
                sd = params$sigma.H) #post-period H
  H.w0 <- rnorm(P, mean = H.b0 + params$c, sd = params$sigma.H)
  H.w1 <- rnorm(P, mean = H.b1 + params$c, sd = params$sigma.H)
  
  # make sure at least within range of ATT from JAMA paper
  delta.b <- rnorm(P, params$theta.P + params$gamma.1*b + params$gamma.2*H.b0, sd = params$sigma.P)
  delta.w <- rnorm(P, delta.b + params$m, sd = params$sigma.P)
  delta <- delta.b * b + delta.w * (1-b)
  
  betas <- c(params$beta.1, params$beta.2, params$beta.3, params$beta.4, params$beta.5, params$beta.6)
  trt.prob <- inv.logit(betas %*% rbind(delta.b, delta.w, n, b, H.b0, H.w0))
  treated <- rbern(n = P, prob = trt.prob)
  
  Yb.pre <- params$alpha.1*H.b0
  Yw.pre <- params$alpha.2*H.w0
  Yb.post <- params$alpha.1*H.b1 + delta.b * treated
  Yw.post <- params$alpha.2*H.w1 + delta.w * treated
  
  # Returns generated practice-level attributes, excluding intermediary values like trt.prob
  return(data.frame(b, n, H.b0, H.b1, H.w0, H.w1, delta.b, delta.w, delta, treated, 
                    Yb.pre, Yw.pre, Yb.post, Yw.post))
  
}

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

