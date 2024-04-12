## Helper functions for calibration

library(tidyverse)
setwd("/Users/alyssachen/Desktop/Projects/target-practice")
source("HSR_simulation_fn.R")


# Calibration helper
calib_assist <- function(new_params){
  df <- make_regions(new_params)
  df <- tibble::as_tibble(df)
  
  attributes_by_S <- df %>% 
    group_by(S) %>% 
    summarise("P(B=1 | S)" = mean(b),
              "P(SSP | S)" = mean(X1),
              "P(system | S)" = mean(X2),
              "P(A | S)" = mean(A)
    ) 
  
  # Now printing difference from target
  target <- data.frame(
    S=c(0,0),
    "P(B=1 | S)" = c(.146, .12),
    "P(SSP | S)" = c(0.35, 0.31),
    "P(system | S)" = c(.33, .316),
    "P(A | S)" = c(.1,.18)
  )
  return(attributes_by_S-target)
}

calib_assist1 <- function(new_params){
  df <- make_regions(new_params)
  df <- tibble::as_tibble(df)

  
  val <- c( df %>% 
              filter(S==1,X1==0) %>% 
              summarise(prop = mean(A)) %>%
              pull(prop), 
           df %>% 
             filter(S==1,X1==1) %>% 
             summarise(prop = mean(A)) %>%
             pull(prop),
           df %>% 
             filter(S==1,X2==0) %>% 
             summarise(prop = mean(A)) %>%
             pull(prop),
           df %>% 
             filter(S==1,X2==1) %>% 
             summarise(prop = mean(A)) %>%
             pull(prop)
           )
  
  
  target <- c(.141, .269, .121, .309)
  
  return(val-target)
}

# input is params
calib_assist2 <- function(new_params){
  df <- make_regions(new_params)
  df <- tibble::as_tibble(df)
  
  # proportion Black in CPC+ and non-CPC+ participants within CPC+ regions
  race_by_participation <- df %>% 
    filter(S==1) %>% 
    group_by(A) %>%
    summarise("P(B=1 | S=1,A)" = mean(b)
    ) 
  
  # Now returning difference from target
  target <- data.frame(
    A=c(0,0),
    "P(B=1 | S=1,A)" = c(.123, .069)
  )
  return(race_by_participation-target)
}

# Similar to above but input is df and output is error
calib_assist3 <- function(df){
  # proportion Black in CPC+ and non-CPC+ participants within CPC+ regions
  race_by_participation <- df %>% 
    filter(S==1) %>% 
    group_by(A) %>%
    summarise(race_prob = mean(b) # race_prob = P(B=1 | S=1,A)
    ) 
  
  # Now returning difference from target
  target <- data.frame(
    A=c(0,0),
    race_prob = c(.12, .069)
  )
  return(race_by_participation-target)
}


# nreps <- 1
# ## Combine params that vary and params that are fixed
# global_params <- data.frame(expand.grid(
#   x1.r = -.64, 
#   x2.r = -.74,
#   phi.1=.17, phi.2=.04,
#   
#   q= -2.12,
#   om.1= seq(-.1,.1, .1), #log odds ratio of Pr(Black) in SSP vs non-SSP
#   om.2= seq(-.1,.2,.1), #log odds ratio of Pr(Black) in system vs indep
#   om.3= -.3,
#   
#   H = 0, sigma.H = .02,
#   psi.1=-.03, 
#   
#   theta.P = -76, sigma.P = 69, 
#   gamma.3=  -92, 
#   gamma.4 = 206, 
#   gamma.5 = -87,
#   
#   beta.0 = -2.51,
#   #beta.1=0,
#   #beta.2=0,
#   beta.3 = 0,
#   beta.4 = seq(-1,1,1),
#   beta.5 = seq(-1,1,1),
#   beta.6 = seq(-1,1,1), 
#   
#   alpha.0 = 10100,
#   alpha.1 = 46500, 
#   alpha.2 = 600,
#   alpha.3 = 12,
#   
#   P= 100, # number of practices in a region
#   replicate=1:nreps)) %>%
#   mutate(scenario=rep(1:324,nreps))
# 
# # What are the scenarios?
# global_params%>% filter(replicate==1) %>% select(scenario,om.1,om.2,beta.4,beta.5,beta.6)
# 
# ## Simulate the data
# simdat <- global_params %>% group_by(scenario,replicate) %>%
#   nest() %>% mutate(singlesim=map(data,make_regions))
# 
# errs <- simdat %>% mutate(error=map(singlesim,calib_assist3))%>%
#   unnest(cols=c(error))%>% 
#   group_by(scenario, A) %>%
#   summarize(mean.err = mean(race_prob)) 
# 
# errs %>% filter(A==0, abs(mean.err) < .02)
# 
# global_params %>% filter(scenario == 85)
# 
