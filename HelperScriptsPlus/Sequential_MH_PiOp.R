# Sequential sampler: unblocked
# Note that I'm supplying an initial value to make it easier to plug into the Gibbs sampler

MH_sequential <- function(R = 2000, burn_in = 500,
                          mh_occ_step = 0.1, mh_pprior_sd = 0.1, 
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
  
  # Don't replace occurrence indicators for species that were actually detected
  detected <- t(apply(A,2, function(x) colSums(x) > 0))*1 # species are rows as in curr_inter
  
  # AMH storage
  pi_samples <- array(data = NA, dim = c(R, num_obs, num_studies))
  occ_samples <- array(data = NA, dim = c(R, num_obs, num_studies))
  accepted_p = array(0, dim = c(R, num_obs, num_studies))
  
  
  for(r in 1:(burn_in + R)){
    
    setTxtProgressBar(pb,r)
    
    #### SAMPLE NEW PROBABILITY
    
    for(st in 1:num_studies){ # Loop over studies
      
      #### 1. Set up useful quantities for this study: all plants
      p_prior_st <- occur_prior_probs[, st]
      beta_st  <- beta[,st,]
      pi_curr_st  <- p_curr[,st]
      occur_curr_st <- occ_curr[,st]
      
      ### 2. Sample new p
      p_prop <- rtruncnorm(n = 1,  mean = pi_curr_st, sd = mh_occ_step, a = 0, b = 1)
      
      ### 3.  Accept or reject for each node within study st
      # Calculating the likelihood, prior, and proposal for proposed value for each proposed value: 
      log_prior_prop <- log(dtruncnorm(x = p_prop, a = 0, b = 1, mean = p_prior_st, sd = mh_pprior_sd))
      log_lik_prop <- log(p_prop^occur_curr_st * (1-p_prop)^(1-occur_curr_st))
      log_prop <- log(dtruncnorm(x = p_prop, a = 0, b= 1, mean = pi_curr_st, sd = mh_pprior_sd))  
      
      # Calculating the likelihood, prior, and proposal for proposed value for each current value: 
      log_prior_curr <- log(dtruncnorm(x = pi_curr_st, a = 0, b = 1, mean = p_prior_st, sd = mh_pprior_sd))
      log_lik_curr <- log(pi_curr_st^occur_curr_st * (1-pi_curr_st)^(1-occur_curr_st))
      log_curr <- log(dtruncnorm(x = pi_curr_st, a = 0, b = 1, mean = p_prop, sd = mh_pprior_sd))
      
      AP <- log_lik_prop + log_prior_prop + log_curr
      AP <- AP - (log_lik_curr + log_prior_curr + log_prop)
      
      update <- log(runif(num_obs)) < AP
      p_curr[update, st] <- p_prop[update]
      if(r>burn_in){
        accepted_p[r-burn_in, update, st] <- 1        
      }
      
    } 
    
    # Now, if a species was detected to interact in a study, set their occurrence
    # probability to 1.
    wh_detected <- which(detected == 1)
    p_curr[wh_detected] <- 1
    
    #### SAMPLE OCCURRENCE INDICATORS
    for (st in 1 : num_studies) {
      occ_others_st <- outer(rep(1,num_obs), occur_others[, st] )
      # Take product over other species
      prob1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* occ_others_st), 1, prod) 
      prob1 <- p_curr[, st] * prob1 # Multiply by occurrence probabilities to get prob occurrence == 1
      prob0 <- 1 - p_curr[, st] # Prob occurrence == 0
      prob1 <- prob1 / (prob1 + prob0) # Standardize
      occ_curr[, st] <- rbinom(num_obs, 1, prob1) # Sample from Bernoullis for each species in this study
    }
    
    
    
    # Now, if a species was detected to interact in a study, set their occurrence
    # indicator to 1.
    wh_detected <- which(detected == 1)
    occ_curr[wh_detected] <- 1
    
    # Store the values after burn-in: for all studies at once
    if (r > burn_in) {
      pi_samples[r-burn_in,,] <- p_curr
      occ_samples[r - burn_in,,] <- occ_curr
    } 
    
    
  }# End loop over iterations r
  
  ret <- list(pi_samples = pi_samples, occ_samples = occ_samples, accepted_p = accepted_p)
  return(ret)
  
}

