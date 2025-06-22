# Define AMH function
# Note that I'm supplying an initial value to make it easier to plug into the Gibbs sampler

MH_Adaptive <- function(R = 2000, burn_in = 500, R0 = 10,
                        epsilon = 1e-6, s.extra = 5, mh_pprior_sd = 0.1, 
                        p_curr, occ_curr, beta, 
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
  beta_new <- beta # set up storage for proposals
  mu <- beta # Initialize running mean mu and set up storage
  
  # Adapative covariance setup
  S <- diag(epsilon, p) #Initial proposal matrix S
  Sigma <- diag(0, p) # Initial covariance matrix
  S_array <- Sigma_array <- S0_array <- array(0, dim = c(p,p,num_obs, num_studies))
  for(i in 1:num_obs){
    for(st in 1:num_studies){
      S_array[,,i,st] <- S # Will store the actual proposals
      Sigma_array[,,i,st] <- Sigma # Will store the recursive covariance computations
      S0_array[,,i,st] <- s.extra*2.38^2 * Sigma/p # C0: non-adaptive portion
    }
  }
  
  
  for(r in 1:(burn_in + R)){
    
    setTxtProgressBar(pb,r)
    
    for(st in 1:num_studies){ # Loop over studies
      
      #### 1. Set up useful quantities for this study: all plants
      p_prior_st <- occur_prior_probs[, st]
      beta_st  <- beta[,st,]
      pi_curr_st  <- p_curr[,st]
      occur_curr_st <- occ_curr[,st]
      
      #### 2. Update proposal distribution and draw sample for each node within this study
      for(i in 1:num_obs){ 
        
        # Updating the covariance matrix and mean based on the previous iteration
        if(r > 1){ 
          Sigma_array[,,i,st] <-  (r-2)/(r-1)*Sigma_array[,,i,st] + tcrossprod(beta[i,st,] - mu[i,st,])/r # recursive covar def when r = 2 we DO have 2 samples before drawing our new value due to initialization
          mu[i,st,] <- (r - 1) / r * mu[i,st,] + beta[i,st,] / r
          S_array[,,i,st] <- s.extra*2.38^2 * Sigma_array[,,i,st] / p + diag(epsilon, p) # slight difference from Haario here
        }  
        
        # Sample new value
        if(r < (R0+1)){
          beta_new[i,st,] <- rmvnorm(1, beta[i,st,], S0_array[,,i,st])       
        } else{
          beta_new[i,st,] <- rmvnorm(1, beta[i,st,], S_array[,,i,st])       
        }
        
      } # End loop over nodes i within study st
      
      ### 3.  Accept or reject for each node within study st
      
      # Transform proposal parameters
      pi_prop_st <- pnorm(beta_new[,st,1])
      occur_prop_st <- ifelse(beta_new[,st,2] > 0 , 1, 0)
      
      # Set up useful values for both current and proposal values
      occ_others_st <- outer(rep(1,num_obs), occur_others[, st] )
      prod1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* occ_others_st), 1, prod)   # Take product over other species
      
      ## Proposal values
      log_prob1_new <- loglik(pi_prop_st, prod1, occur_prop_st) # p(occ | pi, ...)
      log_prob_pi_new <- log(dtruncnorm(x = pi_prop_st, a = 0, b = 1, mean = p_prior_st, sd = mh_pprior_sd))   # p(pi | ...) just the trunc normal prior
      
      ## Current values
      log_prob1_curr <- loglik(pi_curr_st, prod1, occur_curr_st) # p(occ | pi, ...)
      log_prob_pi_curr <- log(dtruncnorm(x = pi_curr_st, a = 0, b = 1, mean = p_prior_st, 
                                         sd = mh_pprior_sd)) # p(pi | ...) just the trunc normal prior
      
      ## Accept or reject for each species (num_obs)
      AP <- log_prob_pi_new + log_prob1_new
      AP <- AP - (log_prob_pi_curr + log_prob1_curr)
      update <- log(runif(num_obs)) < AP
      
      ## Update stored quantities for the selected species
      p_curr[update, st] <- pi_prop_st[update]
      occ_curr[update, st] <- occur_prop_st[update]
      beta[update, st,] <- beta_new[update, st,]
      #accepted[update, st] <- 1 # change the dimension for single iteration
      
      if(r > burn_in){
        accepted[r-burn_in,update, st] <- 1 # change the dimension for single iteration 
      }
      
      
    } # End loop over studies st
    
    # Now, if a species was detected to interact in a study, set their occurrence ########### THIS COULD HAPPEN AT THE END B/C BOTH P AND OCC ARE REPLACED, BUT POST MEAN WOULD BE OFF
    # indicator and probability to 1
    wh_detected <- which(detected == 1) ###################### DEBUG HERE, SEEMS WRONG
    occ_curr[wh_detected] <- 1
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

