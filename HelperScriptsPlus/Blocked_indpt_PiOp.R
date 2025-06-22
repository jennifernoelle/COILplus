# Blocked sampler using
# Gibbs step for binary indicator
# MH step for probability
# Note that I'm supplying an initial value to make it easier to plug into the Gibbs sampler


Blockupdate_OpP  <- function(R = 2000, burn_in = 500,
                             mh_p_step = 0.1, mh_pprior_sd = 0.1, 
                             p_1to0 = 0.65, p_0to1 = 0.25, 
                             p_curr, occ_curr,
                             occur_prior_probs, probobs_curr, probobs_others, occur_others,
                             curr_inter, focus, A){
  
  # Progress bar 
  pb = txtProgressBar(min = 0, max = R + burn_in, initial = 0) 
  
  # Extract important structures
  p <- 2
  num_obs <- nrow(p_curr) # number of species of this type
  num_studies <- ncol(p_curr)
  pipj <- outer(probobs_curr, probobs_others)
  
  # To do later: think about good default switch probabilities
  
  # Don't replace occurrence indicators for species that were actually detected
  detected <- t(apply(A,2, function(x) colSums(x) > 0))*1 # species are rows as in curr_inter
  
  # Set up storage
  pi_samples <- array(data = NA, dim = c(R, num_obs, num_studies))
  occ_samples <- array(data = NA, dim = c(R, num_obs, num_studies))
  accepted = array(0, dim = c(R, num_obs, num_studies))
  
  
  for(r in 1:(burn_in + R)){
    
    setTxtProgressBar(pb,r)
    
    #### SAMPLE NEW PROBABILITY
    
    for(st in 1:num_studies){ # Loop over studies
      
      #### 1. Set up useful quantities for this study: all plants
      pi_prior_st <- occur_prior_probs[, st]
      p_curr_st  <- p_curr[,st]
      occur_curr_st <- occ_curr[,st]
      
      ### 2. Sample new p and new occurrence indicator for all nodes in this study, independently
      p_prop_st <- rtruncnorm(n = 1,  mean = p_curr_st, sd = mh_p_step, a = 0, b = 1) # Sample new p from trunc norm
      prop_probs_st <- ifelse(occur_curr_st == 1, 1-p_1to0, p_0to1)
      occur_prop_st <- rbinom(num_obs, 1, prop_probs_st)
      
      
      ### 3.  Accept or reject both p, occur for each node within study st
      
      # Set up useful values for both current and proposal values
      occ_others_st <- outer(rep(1,num_obs), occur_others[, st] )
      prod1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* occ_others_st), 1, prod)   # Take product over other species
      
      ## Proposal values
      log_prob1_prop <- loglik(p_prop_st, prod1, occur_prop_st) # p(occ | pi, ...)
      log_prob_pi_prop <- log(dtruncnorm(x = p_prop_st, a = 0, b = 1, mean = pi_prior_st, sd = mh_pprior_sd))   # p(pi | ...) just the trunc normal prior
      log_prop_pi_prop <- log(dtruncnorm(x = p_prop_st, a = 0, b = 1, mean = p_curr_st, sd = mh_p_step)) # Asymmetric proposal for pi (trunc norm centered at p_prev)
      log_prop_occ_prop <- log(prop_probs_st^occur_prop_st * (1-prop_probs_st)^(1-occur_prop_st)) # Asymmetric proposal for occ (indicator switch)
      
      ## Current values
      log_prob1_curr <- loglik(p_curr_st, prod1, occur_curr_st) # p(occ | pi, ...)
      log_prob_pi_curr <- log(dtruncnorm(x = p_curr_st, a = 0, b = 1, mean = pi_prior_st, sd = mh_pprior_sd))   # p(pi | ...) just the trunc normal prior
      log_prop_pi_curr <- log(dtruncnorm(x = p_curr_st, a = 0, b = 1, mean = p_prop_st, sd = mh_p_step)) # Asymmetric proposal for pi (trunc norm centered at p_prev)
      prop_probs_st_curr <- ifelse(occur_prop_st == 1, 1-p_1to0, p_0to1) # Prob of moving from proposal to current
      log_prop_occ_curr <- log(prop_probs_st_curr^occur_curr_st * (1-prop_probs_st_curr)^(1-occur_curr_st)) # Asymmetric proposal for occ (indicator switch)
      
      ## Acceptance probability
      AP <- log_prob1_prop + log_prob_pi_prop + log_prop_pi_prop + log_prop_occ_curr
      AP <- AP - (log_prob1_curr + log_prob_pi_curr + log_prop_pi_curr + log_prop_occ_prop) 
      update <- log(runif(num_obs)) < AP
      
      ## Update stored quantities for the selected species
      p_curr[update, st] <- p_prop_st[update]
      occ_curr[update, st] <- occur_prop_st[update]
      #accepted[update, st] <- 1 # change the dimension for single iteration
      
      if(r > burn_in){
        accepted[r-burn_in,update, st] <- 1 # change the dimension for single iteration 
      }
      
    } # End loop over studies st
    
    # Now, if a species was detected to interact in a study, set their occurrence
    # probability to 1.
    wh_detected <- which(detected == 1)
    p_curr[wh_detected] <- 1
    
    # Store the values after burn-in: for all studies at once
    if (r > burn_in) {
      pi_samples[r-burn_in,,] <- p_curr
      occ_samples[r - burn_in,,] <- occ_curr
    } 
    
    
  }# End loop over iterations r
  
  ret <- list(pi_samples = pi_samples, occ_samples = occ_samples, accepted = accepted)
  return(ret)
  
}

