## Helper functions for calibration

library(tidyverse)
source("HSR_simulation_fn.R")

# Run calib_assist multiple times to get average errors of calibration error
# Example call: calib_avg(global_params)
calib_avg <- function(new_params, nruns = 10, table = 1){
  total <- 0
  for(i in 1:nruns){
    if(table == 1){
      total <- total + calib_assist1(new_params)
    }
    else if(table == 2){
      total <- total + calib_assist2(new_params)
    }
    else{
      total <- total + calib_assist3(new_params)
    }
  }
  return(total/nruns)
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

# Calib table functions take a starting set of parameters
# and calibrates the parameters in 1 of 3 calibration target tables

calib_table1 <- function(start_params){
  new_params <- start_params
  param_names <- c("om.3", "phi.1", "phi.2", "beta.0", "beta.6")
  rows <- c(2,2,2,1,2)
  cols <- c("P(B=1 | S)", "P(SSP | S)", "P(system | S)", "P(A | S)","P(A | S)")
  done <- rep(0, length(param_names))
  table <- 1
  
  # Iterate through parameters; three rounds of calibration for each parameter
  for(j in 1:length(param_names)){
    param <- param_names[j]
    row <- rows[j]
    col <- cols[j]
    
    start <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[1]
    end <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[2]
    
    for(iter in 1:3){
      if(done[j]){
        break
      }
      
      new_param_vals <- seq(start, end, (end-start)/11)
      new_errs <- rep(NA, length(new_param_vals))
      
      for(i in 1:length(new_param_vals) ){
        new_params[param]= new_param_vals[i]
        new_errs[i] <- calib_avg(new_params, nruns = 10, table = table)[row,col]
      }
      
      m <- which.min(abs(new_errs))
      if(abs(new_errs[m]) <=.001){
        done[j] = 1
      }
      
      # If error is minimized on the edges of the tested param values, 
      # shift the param value range lower or higher correspondingly
      if(m == 1 | m == length(new_param_vals)){
        next_start <- new_param_vals[m] - ((end-start)/2)
        next_end <- new_param_vals[m] + ((end-start)/2)
        
      }else{
        # Otherwise, narrow the range of param values to test
        # around the minimized error value
        next_start <- new_param_vals[m-1]
        next_end <- new_param_vals[m+1]
      }
      start<-next_start
      end<-next_end
    }
    
    # update start_params with new param value; i.e. each subsequent parameter
    # is calibrated from an updated start_params
    start_params[[param]] <- new_param_vals[m]
  }
  return(start_params)
}


calib_table2 <- function(start_params){
  new_params <- start_params
  param_names <- c("beta.4", "beta.5")
  done <- rep(0, length(param_names))
  table <- 2
  
  # Iterate through parameters to calibrate; three rounds of calibration each
  for(j in 1:length(param_names)){
    param <- param_names[j]
    a <- 1
    if(param == "beta.5"){
      a <- 3
    }
    
    start <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[1]
    end <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[2]
    
    for(iter in 1:3){
      if(done[j]){
        break
      }
      
      new_param_vals <- seq(start, end, (end-start)/10)
      new_errs <- rep(NA, length(new_param_vals))
      
      for(i in 1:length(new_param_vals) ){
        new_params[param]= new_param_vals[i]
        avg_err <- calib_avg(new_params, nruns = 10, table = table)
        new_errs[i] <- abs(avg_err[a])+abs(avg_err[a+1])
      }
      
      m <- which.min(abs(new_errs))
      if(abs(new_errs[m]) <=.002){
        done[j] = 1
      }
      
      if(m == 1 | m == length(new_param_vals)){
        next_start <- new_param_vals[m] - ((end-start)/2)
        next_end <- new_param_vals[m] + ((end-start)/2) 
      }else{
        next_start <- new_param_vals[m-1]
        next_end <- new_param_vals[m+1]
      }
      start<-next_start
      end<-next_end
    }
    
    start_params[[param]] <- new_param_vals[m]
  }
  return(start_params)
}


calib_table3 <- function(start_params){
  new_params <- start_params
  param_names <- c("om.1", "om.2")
  done <- 0
  table <- 3
  
  # Calibrate both om.1 and om.2 simultaneously
  
  param <- "om.1"
  start1 <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[1]
  end1 <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[2]
  
  param <- "om.2"
  start2 <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[1]
  end2 <- sort(c(start_params[[param]]*1.5, start_params[[param]]*.5), decreasing = FALSE)[2]
  
  for(iter in 1:3){
    if(done){
      break
    }
    
    new_param_vals1 <- seq(start1, end1, (end1-start1)/4)
    new_param_vals2 <- seq(start2, end2, (end2-start2)/4)
    
    new_errs <- matrix(NA, nrow = length(new_param_vals1), ncol = length(new_param_vals2))
    
    for(i in 1:length(new_param_vals1) ){
      for(j in 1:length(new_param_vals2)){
        new_params$om.1= new_param_vals1[i]
        new_params$om.2 = new_param_vals2[i]
        
        avg_err <- calib_avg(new_params, nruns = 10, table = table)
        new_errs[i,j] <- abs(avg_err[1,2])+abs(avg_err[2,2])
        
      }
    }
    
    m <- which(abs(new_errs) == min(abs(new_errs)),arr.ind = TRUE) 
    
    if(abs(new_errs[m[1],m[2]]) <=.01){
      done = 1
    }
    
    if(m[1] == 1 | m[1] == length(new_param_vals1)){
      next_start1 <- new_param_vals1[m[1]] - ((end1-start1)/2)
      next_end1 <- new_param_vals1[m[1]] + ((end1-start1)/2)
    }else{
      next_start1 <- new_param_vals1[m[1]-1]
      next_end1 <- new_param_vals1[m[1]+1]
    }
    start1<-next_start1
    end1<-next_end1
    
    if(m[2] == 1 | m[2] == length(new_param_vals2)){
      next_start2 <- new_param_vals2[m[2]] - ((end2-start2)/2)
      next_end2 <- new_param_vals2[m[2]] + ((end2-start2)/2)
    }else{
      next_start2 <- new_param_vals2[m[2]-1]
      next_end2 <- new_param_vals2[m[2]+1]
    }
    start2<-next_start2
    end2<-next_end2
  }
  
  start_params$om.1 <- new_param_vals1[m[1]]
  start_params$om.2 <- new_param_vals2[m[2]]
  
  return(start_params)
}





