## Helper functions for calibration

library(tidyverse)
setwd("/Users/alyssachen/Desktop/Projects/target-practice")
source("HSR_simulation_fn.R")

# Run calib_assist multiple times to get average errors of calibration error
# Example call: calib_avg(global_params)
calib_avg <- function(new_params, nruns = 100, table = 1){
  
  if(table == 1){
    
  }
}

# Get standard errors of calibration error
calib_se <- function(){
  
}

# calib_assist1, calib_assist2, calib_assist3 correspond to the 3 calibration target tables
# Given parameters, simulates a universe and returns error of calibration targets
calib_assist1 <- function(new_params){
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

calib_assist2 <- function(new_params){
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


calib_assist3 <- function(new_params){
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
calib_assist_df <- function(df){
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



