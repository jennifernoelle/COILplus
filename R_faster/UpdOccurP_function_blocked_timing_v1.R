UpdOccurP_blocked_timing <- function(mh_p_step = mh_p_step, mh_pprior_sd = mh_pprior_sd, 
                              p_1to0 = p1to0, p_0to1 = p_0to1, 
                              p_curr, occ_curr,
                              occur_prior_probs, probobs_curr, probobs_others, occur_others,
                              curr_inter, detected, 
                              focus_rows_idx_list, focus_cols_idx_list, focus_int_list) {
  t0 <- proc.time()
  
  # timing accumulators (seconds)
  t_setup   <- 0
  t_propose <- 0
  t_prod1   <- 0
  t_loglik  <- 0
  t_accept  <- 0
  t_write   <- 0
  t_enforce <- 0
  
  # sizes
  num_obs     <- nrow(p_curr)
  num_studies <- ncol(p_curr)
  
  # precompute once per call (pi/pj don't depend on study)
  pipj        <- outer(probobs_curr, probobs_others)
  log1m_pipj <- log1p(-outer(probobs_curr, probobs_others))  # n_i x n_j
  
  # # rows/cols with any focus per study
  # focus_rows_idx <- vector("list", num_studies)
  # focus_cols_idx <- vector("list", num_studies)
  # for (st in seq_len(num_studies)) {
  #   Fst <- focus[,,st]
  #   focus_rows_idx[[st]] <- which(rowSums(Fst != 0L) > 0L)
  #   focus_cols_idx[[st]] <- which(colSums(Fst != 0L) > 0L)
  # }
  # 
  # bookkeeping
  accepted <- matrix(0, nrow = num_obs, ncol = num_studies)
  
  t_setup <- t_setup + elapsed(t0)
  
  
  for (st in 1:num_studies) {
    
    ## --- propose (new p, new occ) ---
    t0 <- proc.time()
    pi_prior_st    <- occur_prior_probs[, st]
    p_curr_st      <- p_curr[, st]
    occur_curr_st  <- occ_curr[, st]
    
    # vectorized truncated normal proposal (still rtruncnorm; we’ll optimize if this dominates)
    p_prop_st      <- rtruncnorm(n = 1, mean = p_curr_st, sd = mh_p_step, a = 0, b = 1)
    prop_probs_st  <- ifelse(occur_curr_st == 1, 1 - p_1to0, p_0to1)
    occur_prop_st  <- rbinom(num_obs, 1, prop_probs_st)
    t_propose <- t_propose + elapsed(t0)
    
    ## --- prod1 build ---
    t0 <- proc.time()
    # columns that matter this study: occur==1 AND any focus in that column
    occ_idx <- intersect(which(occur_others[, st] == 1L), focus_cols_idx_list[[st]])
    if (length(occ_idx) == 0L) {
      log_prod1 <- rep(0, num_obs)
    } else {
      i_idx   <- focus_rows_idx_list[[st]]        # rows with any focus
      foc_st_i <- focus_int_list[[st]]            # pre-cast integer matrix
      log_prod1 <- row_logprod_mask_idx_slice(
        curr_inter, foc_st_i, as.integer(i_idx), as.integer(occ_idx), log1m_pipj
      )
    }
    t_prod1 <- t_prod1 + elapsed(t0)
    
    # ## --- prod1 build ---
    # 
    # t0 <- proc.time()
    # 
    # # columns that matter this study: occur==1 AND any focus in that column
    # occ_idx <- intersect(which(occur_others[, st] == 1L), focus_cols_idx[[st]])
    # if (length(occ_idx) == 0L) {
    #   log_prod1 <- rep(0, num_obs)
    # } else {
    #   i_idx <- focus_rows_idx[[st]]         # only rows with any focus this study
    #   # integer matrices for C++ (0/1)
    #   foc_st_i <- matrix(as.integer(focus[,,st]), nrow(curr_inter), ncol(curr_inter))
    #   # compute only on active rows/cols
    #   log_prod1 <- row_logprod_mask_idx_slice(
    #     curr_inter, foc_st_i, as.integer(i_idx), as.integer(occ_idx), log1m_pipj
    #   )
    # }
    # 
    # t_prod1 <- t_prod1 + elapsed(t0)
    
    
    # t0 <- proc.time()
    # # binary mask of exponents: L * focus_slice * occ_others_st
    # mask_st <- curr_inter * focus[,,st] * outer(rep(1, num_obs), occur_others[, st])
    # # log product = rowSums(mask * log(1 - p_i p_j)), then exp back
    # log_prod1 <- rowSums(mask_st * log1m_pipj)
    # prod1 <- exp(log_prod1)
    # t_prod1 <- t_prod1 + (proc.time() - t0)[3]
    
    
    # t0 <- proc.time()
    # occ_others_st <- outer(rep(1, num_obs), occur_others[, st])
    # mat  <- (1 - pipj) ^ (curr_inter * focus[, , st] * occ_others_st)
    # # if you have the Rcpp version, swap next line to: 
    # prod1 <- rowwise_prod(mat)
    # #prod1 <- exp(rowSums(log(pmax(mat, 1e-300))))
    # t_prod1 <- t_prod1 + (proc.time() - t0)[3]
    
    ## --- log-likelihood pieces ---
    t0 <- proc.time()
    log_prob1_prop <- fast_loglik_log(p_prop_st, log_prod1, occur_prop_st)
    log_prob1_curr <- fast_loglik_log(p_curr_st, log_prod1, occur_curr_st)
    # eta_prop <- qlogis(p_prop_st)
    # eta_curr <- qlogis(p_curr_st)
    # log_prob1_prop <- fast_loglik_eta_cpp(eta_prop, log_prod1, occur_prop_st)
    # log_prob1_curr <- fast_loglik_eta_cpp(eta_curr,  log_prod1, occur_curr_st)
    
    # priors and proposal densities (truncated normals)
    log_prob_pi_prop <- log(dtruncnorm(x = p_prop_st, a = 0, b = 1, mean = pi_prior_st, sd = mh_pprior_sd))
    log_prob_pi_curr <- log(dtruncnorm(x = p_curr_st,  a = 0, b = 1, mean = pi_prior_st, sd = mh_pprior_sd))
    log_prop_pi_prop <- log(dtruncnorm(x = p_prop_st, a = 0, b = 1, mean = p_curr_st,  sd = mh_p_step))
    log_prop_pi_curr <- log(dtruncnorm(x = p_curr_st,  a = 0, b = 1, mean = p_prop_st, sd = mh_p_step))
    
    prop_probs_st_curr <- ifelse(occur_prop_st == 1, 1 - p_1to0, p_0to1)
    log_prop_occ_prop  <- log(prop_probs_st^occur_prop_st * (1 - prop_probs_st)^(1 - occur_prop_st))
    log_prop_occ_curr  <- log(prop_probs_st_curr^occur_curr_st * (1 - prop_probs_st_curr)^(1 - occur_curr_st))
    t_loglik <- t_loglik + elapsed(t0)
    
    ## --- accept/reject + writeback ---
    t0 <- proc.time()
    AP <- (log_prob1_prop + log_prob_pi_prop + log_prop_pi_curr + log_prop_occ_curr) -
      (log_prob1_curr + log_prob_pi_curr + log_prop_pi_prop + log_prop_occ_prop)
    
    update <- (log(runif(num_obs)) < AP)
    t_accept <- t_accept + elapsed(t0)
    
    t0 <- proc.time()
    if (any(update)) {
      p_curr[update, st]   <- p_prop_st[update]
      occ_curr[update, st] <- occur_prop_st[update]
      accepted[update, st] <- 1
    }
    t_write <- t_write + elapsed(t0)
  }
  
  t0 <- proc.time()
  # enforce detected=1 ⇒ occurrence/probability = 1
  wh <- which(detected == 1)
  if (length(wh)) {
    p_curr[wh]   <- 1
    occ_curr[wh] <- 1
  }
  t_enforce <- t_enforce + elapsed(t0)
  
  list(
    p_curr   = p_curr,
    occ_curr = occ_curr,
    accepted = accepted,
    timings  = c(setup = t_setup, propose = t_propose, prod1 = t_prod1, loglik = t_loglik,
                 accept = t_accept, write = t_write, enforce = t_enforce)
  )
}
