UpdOccurP_blocked_timing_v3 <- function(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr, occ_curr,
    occur_prior_probs, probobs_curr, probobs_others, occur_others,
    curr_inter, 
    focus_rows_idx_list, focus_cols_idx_list, focus_int_list, 
    focus_rows_by_col_list, # NEW
    logZ_prior_mat = NULL, 
    detected   # matrix [num_obs x num_studies] (pass detected_P / detected_B)
    ){
  
  t_setup <- t_propose <- t_prod1 <- t_loglik <- t_accept <- t_write <- t_enforce <- 0
  t0 <- proc.time()
  
  # sizes
  num_obs     <- nrow(p_curr)
  num_studies <- ncol(p_curr)
  
  # normalize detected to logical matrix (per-study), preserving original semantics
  if (!(is.matrix(detected) && all(dim(detected) == c(num_obs, num_studies))))
    stop("`detected` must be a matrix with dim = nrow(p_curr) x ncol(p_curr)")
  detected_mat <- detected != 0L
  
  # per-call things (depend only on pi/pj)
  log1m_pipj <- log1p(-outer(probobs_curr, probobs_others))  # n_i x n_j
  
  accepted <- matrix(0L, nrow = num_obs, ncol = num_studies)
  
  t_setup <- t_setup + elapsed(t0)
  
  for (st in seq_len(num_studies)) {

    
    ## ACTIVE MASK (per-study): rows we may update
    active_idx <- which(!detected_mat[, st])
    
    ## --- propose (new p, new occ) ---
    t0 <- proc.time()
    pi_prior_st    <- occur_prior_probs[, st]
    p_curr_st      <- p_curr[, st]
    occur_curr_st  <- occ_curr[, st]
    
    # draw for ALL rows (same as old behavior)
    p_prop_full <- rtruncnorm(n = 1, mean = p_curr_st, sd = mh_p_step, a = 0, b = 1)
    occ_prop_full <- rbinom(
      n = num_obs, size = 1L,
      prob = ifelse(occur_curr_st == 1L, 1 - p_1to0, p_0to1)
    )
    
    # then mask: inactive rows keep current values
    p_prop_st     <- p_curr_st
    occur_prop_st <- occur_curr_st
    if (length(active_idx)) {
      p_prop_st[active_idx]     <- p_prop_full[active_idx]
      occur_prop_st[active_idx] <- occ_prop_full[active_idx]
    }
    t_propose <- t_propose + elapsed(t0)
    
    ## --- prod1 (restrict rows to active) ---
    t0 <- proc.time()
    occ_idx <- intersect(which(occur_others[, st] == 1L), focus_cols_idx_list[[st]])
    if (length(occ_idx) == 0L) {
      log_prod1 <- numeric(num_obs)  # zeros
    } else {
      # i_idx    <- intersect(focus_rows_idx_list[[st]], active_idx)
      i_idx <- integer(0)
      if (length(occ_idx)) {
        i_idx <- unique(unlist(focus_rows_by_col_list[[st]][occ_idx], use.names = FALSE))
        i_idx <- intersect(i_idx, active_idx)
      }
      if (length(i_idx)) {
        log_prod1 <- row_logprod_mask_idx_slice(
          curr_inter, focus_int_list[[st]],
          as.integer(i_idx), as.integer(occ_idx), log1m_pipj
        )
      } else {
        log_prod1 <- numeric(num_obs)
      }
    }
    t_prod1 <- t_prod1 + elapsed(t0)
    
    ## --- log-likelihood & priors (only active rows) ---
    t0 <- proc.time()
    # logZ_prior     <- log(pnorm((1 - pi_prior_st)/mh_pprior_sd) - pnorm((-pi_prior_st)/mh_pprior_sd))
    logZ_prior <- if (is.null(logZ_prior_mat)) {
      log(pnorm((1 - pi_prior_st)/mh_pprior_sd) - pnorm((-pi_prior_st)/mh_pprior_sd))
    } else {
      logZ_prior_mat[, st]
    }
    logZ_step_curr <- log(pnorm((1 - p_curr_st)/mh_p_step)      - pnorm((-p_curr_st)/mh_p_step))
    logZ_step_prop <- log(pnorm((1 - p_prop_st)/mh_p_step)      - pnorm((-p_prop_st)/mh_p_step))
    
    APa <- NULL  # acceptance stat for active rows only
    if (length(active_idx)) {
      # (keep your existing lp_prop/lp_curr + truncated-normal + Bernoulli pieces)
      lp_prop <- fast_loglik_log(p_prop_st[active_idx], log_prod1[active_idx], occur_prop_st[active_idx])
      lp_curr <- fast_loglik_log(p_curr_st[active_idx],  log_prod1[active_idx], occur_curr_st[active_idx])
      
      log_prob_pi_prop <- log_dtrunc01(p_prop_st[active_idx],  pi_prior_st[active_idx], mh_pprior_sd, logZ_prior[active_idx])
      log_prob_pi_curr <- log_dtrunc01(p_curr_st[active_idx],  pi_prior_st[active_idx], mh_pprior_sd, logZ_prior[active_idx])
      log_prop_pi_prop <- log_dtrunc01(p_prop_st[active_idx],  p_curr_st[active_idx],   mh_p_step,    logZ_step_curr[active_idx])
      log_prop_pi_curr <- log_dtrunc01(p_curr_st[active_idx],  p_prop_st[active_idx],   mh_p_step,    logZ_step_prop[active_idx])
      
      eps <- 1e-12
      pp  <- p_0to1 + occur_curr_st[active_idx] * ((1 - p_1to0) - p_0to1)
      ppc <- p_0to1 + occur_prop_st[active_idx] * ((1 - p_1to0) - p_0to1)
      pp  <- pmin(pmax(pp,  eps), 1 - eps)
      ppc <- pmin(pmax(ppc, eps), 1 - eps)
      log_prop_occ_prop <- occur_prop_st[active_idx]*log(pp)  + (1-occur_prop_st[active_idx])*log1p(-pp)
      log_prop_occ_curr <- occur_curr_st[active_idx]*log(ppc) + (1-occur_curr_st[active_idx])*log1p(-ppc)
      
      # Diagnostics printing
      if (anyNA(p_prop_st[active_idx]) || anyNA(p_curr_st[active_idx])) {
        warning(sprintf("[UpdOccurP] NA in p_* (st=%d): prop=%d curr=%d",
                        st, sum(is.na(p_prop_st[active_idx])), sum(is.na(p_curr_st[active_idx]))))
      }
      if (any(!is.finite(log_prod1[active_idx]))) {
        warning(sprintf("[UpdOccurP] non-finite log_prod1 (st=%d): %d",
                        st, sum(!is.finite(log_prod1[active_idx]))))
      }
      if (!(mh_p_step > 0 && mh_pprior_sd > 0)) stop("mh_p_step and mh_pprior_sd must be > 0")
      
      
      APa <- (lp_prop + log_prob_pi_prop + log_prop_pi_curr + log_prop_occ_curr) -
        (lp_curr + log_prob_pi_curr + log_prop_pi_prop + log_prop_occ_prop)
      # Treat any non-finite AP as automatic reject
      APa[!is.finite(APa)] <- -Inf
      
    }
    t_loglik <- t_loglik + elapsed(t0)
    
    ## --- accept/reject + writeback ---
    # NOTE: `accepted` marks MH acceptances only on non-detected rows (detected==1 are skipped).
    # This is intentional for speed; detected rows would have been enforced to 1's anyway.
    # Optional: make that explicit in diagnostics:    accepted[detected_mat] <- NA_integer_
    
    t0 <- proc.time()
    u_full <- runif(num_obs)                 # draw for ALL rows (preserves RNG stream)
    update <- rep(FALSE, num_obs)
    if (length(active_idx)) {
      update[active_idx] <- is.finite(APa) & (log(u_full[active_idx]) < APa)
    }
    t_accept <- t_accept + elapsed(t0)
    
    t0 <- proc.time()
    if (any(update)) {
      p_curr[update, st]   <- p_prop_st[update]
      occ_curr[update, st] <- occur_prop_st[update]
      accepted[update, st] <- 1L
    }
    t_write <- t_write + elapsed(t0)
  }
  
  # no row-wise “enforce all studies” block here (preserves original sparsity)
  # t_enforce remains 0, included in timings for consistency
  
  list(
    p_curr   = p_curr,
    occ_curr = occ_curr,
    accepted = accepted,
    timings  = c(
      setup   = t_setup,
      propose = t_propose,
      prod1   = t_prod1,
      loglik  = t_loglik,
      accept  = t_accept,
      write   = t_write,
      enforce = t_enforce
    )
  )
}


# UpdOccurP_blocked_timing_v2 <- function(mh_p_step = mh_p_step, mh_pprior_sd = mh_pprior_sd, 
#                                      p_1to0 = p1to0, p_0to1 = p_0to1, 
#                                      p_curr, occ_curr,
#                                      occur_prior_probs, probobs_curr, probobs_others, occur_others,
#                                      curr_inter, detected, 
#                                      focus_rows_idx_list, focus_cols_idx_list, focus_int_list) {
#   
#   elapsed <- function(t0) (proc.time() - t0)[[3]]
#   
#   t0 <- proc.time()
#   # timing accumulators (seconds)
#   t_setup   <- 0
#   t_propose <- 0
#   t_prod1   <- 0
#   t_loglik  <- 0
#   t_accept  <- 0
#   t_write   <- 0
#   t_enforce <- 0
#   
#   # sizes
#   num_obs     <- nrow(p_curr)
#   num_studies <- ncol(p_curr)
#   
#   # precompute once per call (pi/pj don't depend on study)
#   pipj        <- outer(probobs_curr, probobs_others)
#   log1m_pipj  <- log1p(-outer(probobs_curr, probobs_others))  # n_i x n_j
#   
#   # rows we actually update (skip detected==1)
#   active_idx <- which(detected == 0L)
#   
#   # normalize `detected` -> logical matrix [num_obs x num_studies]
#   if (is.matrix(detected)) {
#     if (!identical(dim(detected), c(num_obs, num_studies)))
#       stop("`detected` dim mismatch")
#     detected_mat <- detected != 0L
#   } else if (length(detected) == num_obs) {
#     detected_mat <- matrix(detected != 0L, nrow = num_obs, ncol = num_studies)
#   } else {
#     stop("`detected` must be length nrow(p_curr) or a nrow x ncol matrix.")
#   }
#   
#   # bookkeeping
#   accepted <- matrix(0L, nrow = num_obs, ncol = num_studies)
#   
#   t_setup <- t_setup + elapsed(t0)
#   
#   for (st in seq_len(num_studies)) {
#     
#     # Masking for this study
#     active_idx <- which(detected[, st] == 0L)
#     
#     ## --- propose (new p, new occ) ---
#     t0 <- proc.time()
#     pi_prior_st    <- occur_prior_probs[, st]
#     p_curr_st      <- p_curr[, st]
#     occur_curr_st  <- occ_curr[, st]
#     
#     # rows we actually update (skip detected + any NA rows)
#     i_idx <- intersect(focus_rows_idx_list[[st]], active_idx)
#     
#     # start from current, then overwrite only active rows
#     p_prop_st     <- p_curr_st
#     occur_prop_st <- occur_curr_st
#     
#     if (length(active_idx)) {
#       p_prop_st[active_idx] <- rtruncnorm(
#         n    = length(active_idx),
#         mean = p_curr_st[active_idx],
#         sd   = rep_len(mh_p_step, length(active_idx)),
#         a = 0, b = 1
#       )
#       
#       # proposal probs for occ indicator (clamped & NA-safe)
#       prop_probs_st <- ifelse(occur_curr_st[active_idx] == 1L, 1 - p_1to0, p_0to1)
#       prop_probs_st <- pmin(pmax(prop_probs_st, 0), 1)       # keep in [0,1]
#       prop_probs_st[!is.finite(prop_probs_st)] <- 0          # guard
#       occur_prop_st[active_idx] <- rbinom(length(active_idx), 1L, prop_probs_st)
#     }
#     t_propose <- t_propose + elapsed(t0)
#     
#     ## --- log-lik / priors only on active rows; build AP full length ---
#     t0 <- proc.time()
#     AP <- rep(-Inf, num_obs)
#     if (length(active_idx)) {
#       eps <- 1e-12
#       eta_prop <- qlogis(pmin(pmax(p_prop_st[active_idx], eps), 1 - eps))
#       eta_curr <- qlogis(pmin(pmax(p_curr_st[active_idx],  eps), 1 - eps))
#       
#       lp_prop <- fast_loglik_eta_cpp(eta_prop, log_prod1[active_idx], as.integer(occur_prop_st[active_idx]))
#       lp_curr <- fast_loglik_eta_cpp(eta_curr,  log_prod1[active_idx], as.integer(occur_curr_st[active_idx]))
#       
#       dtn <- function(x, m, s) log(dtruncnorm(x, 0, 1, mean = m, sd = s))
#       log_prob_pi_prop <- dtn(p_prop_st[active_idx],  pi_prior_st[active_idx], mh_pprior_sd)
#       log_prob_pi_curr <- dtn(p_curr_st[active_idx],  pi_prior_st[active_idx], mh_pprior_sd)
#       log_prop_pi_prop <- dtn(p_prop_st[active_idx],  p_curr_st[active_idx],  mh_p_step)
#       log_prop_pi_curr <- dtn(p_curr_st[active_idx],  p_prop_st[active_idx],  mh_p_step)
#       
#       pp   <- ifelse(occur_curr_st[active_idx]==1L, 1 - p_1to0, p_0to1)
#       ppc  <- ifelse(occur_prop_st[active_idx]==1L, 1 - p_1to0, p_0to1)
#       pp   <- pmin(pmax(pp,  0), 1)
#       ppc  <- pmin(pmax(ppc, 0), 1)
#       log_prop_occ_prop <- occur_prop_st[active_idx]*log(pp)  + (1-occur_prop_st[active_idx])*log1p(-pp)
#       log_prop_occ_curr <- occur_curr_st[active_idx]*log(ppc) + (1-occur_curr_st[active_idx])*log1p(-ppc)
#       
#       AP[active_idx] <- (lp_prop + log_prob_pi_prop + log_prop_pi_curr + log_prop_occ_curr) -
#         (lp_curr + log_prob_pi_curr + log_prop_pi_prop + log_prop_occ_prop)
#     }
#     t_loglik <- t_loglik + elapsed(t0)
#     
#     ## --- accept/reject + writeback (ensure full-length logical) ---
#     t0 <- proc.time()
#     update <- rep(FALSE, num_obs)
#     if (length(active_idx)) {
#       upd_active <- is.finite(AP[active_idx]) & (log(runif(length(active_idx))) < AP[active_idx])
#       update[active_idx] <- upd_active
#     }
#     t_accept <- t_accept + elapsed(t0)
#     
#     t0 <- proc.time()
#     if (any(update)) {
#       p_curr[update, st]   <- p_prop_st[update]
#       occ_curr[update, st] <- occur_prop_st[update]
#       accepted[update, st] <- 1L
#     }
#     t_write <- t_write + elapsed(t0)
#     
#     t0 <- proc.time()
#     # enforce detected=1 ⇒ occurrence/probability = 1
#     wh <- which(detected[, st] == 1L)
#     if (length(wh)) {
#       p_curr[wh, st]   <- 1
#       occ_curr[wh, st] <- 1L
#     }
#     # sanity check
#     stopifnot(all(occ_curr[wh, st] == 1L))
#     t_enforce <- t_enforce + elapsed(t0)
#   }
# 
#   
#   list(
#     p_curr   = p_curr,
#     occ_curr = occ_curr,
#     accepted = accepted,
#     timings  = c(
#       setup   = t_setup,
#       propose = t_propose,
#       prod1   = t_prod1,
#       loglik  = t_loglik,
#       accept  = t_accept,
#       write   = t_write,
#       enforce = t_enforce
#     )
#   )
# }