p_curr = this_occur_P
occ_curr = this_O_P
occur_prior_probs = p_occur_P
probobs_curr = this_pj
probobs_others = this_pi
occur_others = this_O_V
curr_inter = t(this_L)
detected = detected_P
focus_rows_idx_list = focus_rows_idx_P
focus_cols_idx_list = focus_cols_idx_P
focus_int_list      = focus_int_P
logZ_prior_mat = logZ_prior_P

bench_parallel <- function(cores_vec = c(1,2,4,6,8)) {
  out <- data.frame(cores = integer(), sec = numeric())
  for (k in cores_vec) {
    t0 <- proc.time()
    # call your occur update with use_parallel= (k>1), ncores=k
    invisible( UpdOccurP_blocked_parallel(
      mh_p_step = mh_p_step, mh_pprior_sd = mh_pprior_sd,
      p_1to0 = p_1to0, p_0to1 = p_0to1,
      p_curr = this_occur_P, occ_curr = this_O_P,
      occur_prior_probs = p_occur_P,
      probobs_curr = this_pj, probobs_others = this_pi,
      occur_others = this_O_V,
      curr_inter = t(this_L),
      detected = detected_P,
      focus_rows_idx_list = focus_rows_idx_P,
      focus_cols_idx_list = focus_cols_idx_P,
      focus_int_list      = focus_int_P,
      logZ_prior_mat = logZ_prior_P, 
      use_parallel = (k>1), ncores = k) )
    out <- rbind(out, data.frame(cores = k, sec = (proc.time()-t0)[3]))
  }
  out
}


# 3) run it
RNGkind("L'Ecuyer-CMRG"); set.seed(1)   # reproducible across cores
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")

bench_parallel(cores_vec = c(1,2,4,6,8))



#### Compare to old version

## ---- one-time: keep BLAS single-threaded for fair scaling
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")

## helper to time a single call
.time_one <- function(expr) {
  gc(FALSE); t0 <- proc.time(); force(expr); (proc.time()-t0)[[3]]
}

## run one P-block update with fresh copies (so timing is comparable)
.run_occurP_once <- function(parallel=FALSE, ncores=6L) {
  UpdOccurP_blocked_parallel(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = this_occur_P, occ_curr = this_O_P,
    occur_prior_probs = p_occur_P,
    probobs_curr = this_pj, probobs_others = this_pi,
    occur_others = this_O_V,
    curr_inter = t(this_L),
    focus_rows_idx_list = focus_rows_idx_P,
    focus_cols_idx_list = focus_cols_idx_P,
    focus_int_list      = focus_int_P,
    detected = detected_P,
    use_parallel = parallel, ncores = ncores
  )
  invisible(NULL)
}


## run one P-block update with fresh copies (so timing is comparable)
.run_occurP_old_once <- function(parallel=FALSE, ncores=6L) {
  UpdOccurP_blocked_timing_v2(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = this_occur_P, occ_curr = this_O_P,
    occur_prior_probs = p_occur_P,
    probobs_curr = this_pj, probobs_others = this_pi,
    occur_others = this_O_V,
    curr_inter = t(this_L),
    focus_rows_idx_list = focus_rows_idx_P,
    focus_cols_idx_list = focus_cols_idx_P,
    focus_int_list      = focus_int_P,
    detected = detected_P
  )
  invisible(NULL)
}


## run one P-block update with fresh copies (so timing is comparable)
.run_occurP_new_once <- function(parallel=FALSE, ncores=6L) {
  UpdOccurP_blocked_timing_v3(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = this_occur_P, occ_curr = this_O_P,
    occur_prior_probs = p_occur_P,
    probobs_curr = this_pj, probobs_others = this_pi,
    occur_others = this_O_V,
    curr_inter = t(this_L),
    focus_rows_idx_list = focus_rows_idx_P,
    focus_cols_idx_list = focus_cols_idx_P,
    focus_int_list      = focus_int_P,
    focus_rows_by_col_list = focus_rows_by_col_P,   # NEW
    detected = detected_P
  )
  invisible(NULL)
}


## benchmark both modes
bench_occurP_compare <- function(nreps=5L, ncores_parallel=6L, seed=1L) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  
  serial_times_old <- replicate(nreps, .time_one(.run_occurP_old_once(FALSE, 1L)))
  serial_times_new <- replicate(nreps, .time_one(.run_occurP_new_once(FALSE, 1L)))
  serial_times <- replicate(nreps, .time_one(.run_occurP_once(FALSE, 1L)))
  parallel_times <- replicate(nreps, .time_one(.run_occurP_once(TRUE, ncores_parallel)))
  
  data.frame(
    mode   = rep(c("serial_old", "serial_new",  "serial_fast","parallel"), each = nreps),
    cores  = rep(c(1L, 1L, 1L, ncores_parallel), each = nreps),
    rep    = rep(seq_len(nreps),4),
    sec    = c(serial_times_old, serial_times_new, serial_times, parallel_times)
  )
}

## ---- run it
res <- bench_occurP_compare(nreps = 10, ncores_parallel = 6)
aggregate(sec ~ mode + cores, data = res,
          FUN = function(x) c(median=median(x), min=min(x), max=max(x)))
res


## 

mean_occ_cols <- mean(sapply(seq_len(ncol(this_O_V)), function(st)
  sum(this_O_V[, st] == 1L & (seq_len(nrow(this_O_V)) %in% focus_cols_idx_P[[st]]))
))
mean_i_rows <- mean(sapply(seq_len(ncol(this_O_V)), function(st)
  length(intersect(focus_rows_idx_P[[st]], which(!detected_P[, st])))
))
cat("avg occ_cols per study =", mean_occ_cols,
    "  avg active_rows per study =", mean_i_rows, "\n")
