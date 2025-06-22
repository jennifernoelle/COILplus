
library(truncnorm)
library(mvtnorm)



# Define RWMH function
# Note that I'm supplying an initial value to make it easier to plug into the Gibbs sampler

MH_RW <- function(R = 2000, burn_in = 500, #s.extra = 5,
                  mh_pprior_sd = 0.1, mh_prop_sd =0.25, 
                  new_value_p, new_value_occ, new_value_theta, 
                  occur_prior_probs, probobs_curr, probobs_others, occur_others,
                  curr_inter, focus, A){
  
  # Progress bar 
  pb = txtProgressBar(min = 0, max = R + burn_in, initial = 0) 
  
  # Extract important structures
  p <- 2 # 2 parameters
  num_obs <- nrow(p_curr) # number of species of this type
  num_studies <- ncol(p_curr)
  pipj <- outer(probobs_curr, probobs_others)
  
  
  # Don't replace occurrence indicators for species that were actually detected
  detected <- t(apply(A,2, function(x) colSums(x) > 0))*1 # species are rows as in curr_inter
  
  # AMH storage
  pi_samples <- array(data = NA, dim = c(R, num_obs, num_studies))
  occ_samples <- array(data = NA, dim = c(R, num_obs, num_studies))
  accepted = array(0, dim = c(R, num_obs, num_studies))
  theta_prop_st <- matrix(data = NA, nrow = num_obs, ncol = 2)
  
  for(r in 1:(R + burn_in)){
    setTxtProgressBar(pb,r)
    
    for (st in 1 : num_studies) {
      
      p_prior_st <- occur_prior_probs[, st]
      theta_curr_st  <- new_value_theta[,st,]
      pi_curr_st  <- new_value_p[,st]
      occur_curr_st <- new_value_occ[,st]
      
      ### Generate proposals independently for each species and then transform all 
      for(i in 1:num_obs){
        cov_prop <- mh_prop_sd*diag(1,p) # Modify this after checking acceptance rate
        theta_prop_st[i,] <- rmvnorm(1, theta_curr_st[i,], cov_prop) 
      }
      
      pi_prop_st <- pnorm(theta_prop_st[,1])
      occur_prop_st <- ifelse(theta_prop_st[,2] > 0 , 1, 0)
      
      ### Compute acceptance probabilities
      
      ## Proposal 
      # p(occ | pi, ...)
      # Double check: I should be using new value of probs in here?
      occ_others_st <- outer(rep(1,num_obs), occur_others[, st] )
      prod1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* occ_others_st), 1, prod)   # Take product over other species
      prob1 <- pi_prop_st * prod1 # Multiply by occurrence probabilities to get prob occurrence == 1
      prob0 <- 1 - pi_prop_st # Prob occurrence == 0
      prob1 <- prob1 / (prob1 + prob0) # Standardize
      this_prob <- prob1*occur_prop_st + (1-prob1)*(1-occur_prop_st) # p = p1 if occ = 1, (1-p1) if occ = 0
      log_prob1_prop <- log(this_prob) # this gives a vector length num_obs
      
      # p(pi | ...) just the trunc normal prior
      log_prob_pi_prop <- log(dtruncnorm(x = pi_prop_st, a = 0, b = 1, mean = p_prior_st, sd = mh_pprior_sd))
      
      ## Current: can probably reuse some of this from above if i rename things
      # p(occ | pi, ...)
      #p_occ_others_st <- outer(rep(1,num_obs), p_occur_others[, st] )
      #prob1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* occ_others_st), 1, prod)   # Take product over other species
      prob1 <- pi_curr_st * prod1 # Multiply by current occurrence probabilities to get prob occurrence == 1, same product as above
      prob0 <- 1 - pi_curr_st # Prob occurrence == 0
      prob1 <- prob1 / (prob1 + prob0) # Standardize
      this_prob <- prob1*occur_curr_st + (1-prob1)*(1-occur_curr_st) # p = p1 if occ = 1, (1-p1) if occ = 0
      log_prob1_curr <- log(this_prob) # this gives a vector length num_obs
      
      
      # p(pi | ...) just the trunc normal prior
      log_prob_pi_curr <- log(dtruncnorm(x = pi_curr_st, a = 0, b = 1, mean = p_prior_st, sd = mh_pprior_sd))
      
      ## Accept or reject for each species (num_obs)
      AP <- log_prob_pi_prop + log_prob1_prop
      AP <- AP - (log_prob_pi_curr + log_prob1_curr)
      
      update <- log(runif(num_obs)) < AP
      
      ## Update stored quantities for the selected species
      new_value_p[update, st] <- pi_prop_st[update]
      new_value_occ[update, st] <- occur_prop_st[update]
      new_value_theta[update, st,] <- theta_prop_st[update,]
      #accepted[update, st] <- 1 # change the dimension for single iteration
      if(r>burn_in){
        accepted[r - burn_in,update, st] <- 1 # change the dimension for single iteration      
      }
      
    }
    
    # Now, if a species was detected to interact in a study, set their occurrence ########### THIS COULD HAPPEN AT THE END B/C BOTH P AND OCC ARE REPLACED, BUT POST MEAN WOULD BE OFF
    # indicator and probability to 1
    wh_detected <- which(detected == 1)
    new_value_occ[wh_detected] <- 1
    new_value_p[wh_detected] <- 1
    
    if(r>burn_in){
      pi_samples[r - burn_in,,] <- new_value_p
      occ_samples[r - burn_in ,,] <- new_value_occ    
    }
    
    
  }
  
  ret <- list(pi_samples = pi_samples, occ_samples = occ_samples, accepted = accepted)
  return(ret)
}
