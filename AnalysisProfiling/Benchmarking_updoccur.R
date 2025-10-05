# This code benchmarks the Upd_occur_p function variants

## --- freeze a snapshot of the current state ---
snap <- list(
  p_curr   = this_occur_P,
  occ_curr = this_O_P,
  p_prior  = p_occur_P,
  pj       = this_pj,
  pi       = this_pi,
  occ_Other= this_O_V,
  L        = this_L,
  det      = detected_P,
  # focus helpers
  rows_idx = focus_rows_idx_P,
  cols_idx = focus_cols_idx_P,
  focus_i  = focus_int_P
  #rows_by_col = focus_rows_by_col_P  # only used by new version; set to NULL if not
)

.time <- function(expr) { t0 <- proc.time(); force(expr); (proc.time()-t0)[[3]] }

run_old_once <- function() {
  # fresh copies so runs are comparable
  UpdOccurP_blocked(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = snap$p_curr, occ_curr = snap$occ_curr,
    occur_prior_probs = snap$p_prior,
    probobs_curr = snap$pj, probobs_others = snap$pi,
    occur_others = snap$occ_Other,
    curr_inter = t(snap$L),
    focus = aperm(focus, c(2,1,3)),     # your old call’s focus
    detected = snap$det
  )
}

run_new_once <- function() {
  UpdOccurP_blocked_timing_v2(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = snap$p_curr, occ_curr = snap$occ_curr,
    occur_prior_probs = snap$p_prior,
    probobs_curr = snap$pj, probobs_others = snap$pi,
    occur_others = snap$occ_Other,
    curr_inter = t(snap$L),
    focus_rows_idx_list = snap$rows_idx,
    focus_cols_idx_list = snap$cols_idx,
    focus_int_list      = snap$focus_i,
   # focus_rows_by_col_list = snap$rows_by_col,  # if you added this optimization
    detected = snap$det,
    #use_parallel = FALSE, ncores = 1L
  )
}

bench_occur_versions <- function(n=5L, seed=1L) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  invisible(run_old_once()); invisible(run_new_once())  # warmup
  
  old <- replicate(n, .time(run_old_once()))
  new <- replicate(n, .time(run_new_once()))
  data.frame(
    version = c("old","new"),
    median  = c(median(old), median(new)),
    min     = c(min(old), min(new)),
    max     = c(max(old), max(new))
  )
}

## run it
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
bench_occur_versions()



### Now also benchmark v3

## Build rows-by-col once (if you haven't already)
## using your existing per-study integer focus slices
build_rows_by_col <- function(focus_int_list) {
  lapply(seq_along(focus_int_list), function(st) {
    Fst <- focus_int_list[[st]] != 0L
    apply(Fst, 2, function(col) which(col))  # list-of-int per column
  })
}

## One-shot A/B/C on identical state (old vs v2 vs v3)
bench_occur_versions3 <- function(n = 5L, seed = 1L) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
             MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
  
  # snapshot current inputs (P-block)
  snap <- list(
    p_curr   = this_occur_P,
    occ_curr = this_O_P,
    p_prior  = p_occur_P,
    pj       = this_pj,
    pi       = this_pi,
    occ_Other= this_O_V,
    L        = this_L,
    det      = detected_P,
    rows_idx = focus_rows_idx_P,
    cols_idx = focus_cols_idx_P,
    focus_i  = focus_int_P
  )
  rows_by_col_P <- build_rows_by_col(snap$focus_i)  # for v3
  
  .time <- function(f) { t0 <- proc.time(); r <- f(); list(sec=(proc.time()-t0)[[3]], res=r) }
  
  run_old <- function() UpdOccurP_blocked(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = snap$p_curr, occ_curr = snap$occ_curr,
    occur_prior_probs = snap$p_prior,
    probobs_curr = snap$pj, probobs_others = snap$pi,
    occur_others = snap$occ_Other,
    curr_inter = t(snap$L),
    focus = aperm(focus, c(2,1,3)),
    detected = snap$det
  )
  
  run_v2 <- function() UpdOccurP_blocked_timing_v2(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = snap$p_curr, occ_curr = snap$occ_curr,
    occur_prior_probs = snap$p_prior,
    probobs_curr = snap$pj, probobs_others = snap$pi,
    occur_others = snap$occ_Other,
    curr_inter = t(snap$L),
    focus_rows_idx_list = snap$rows_idx,
    focus_cols_idx_list = snap$cols_idx,
    focus_int_list      = snap$focus_i,
    detected = snap$det,
    #use_parallel = FALSE, ncores = 1L
  )
  
  run_v3 <- function() UpdOccurP_blocked_timing_v3(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = snap$p_curr, occ_curr = snap$occ_curr,
    occur_prior_probs = snap$p_prior,
    probobs_curr = snap$pj, probobs_others = snap$pi,
    occur_others = snap$occ_Other,
    curr_inter = t(snap$L),
    focus_rows_idx_list   = snap$rows_idx,
    focus_cols_idx_list   = snap$cols_idx,
    focus_int_list        = snap$focus_i,
    focus_rows_by_col_list= rows_by_col_P,      # <-- key difference
    detected = snap$det
  )
  
  invisible(.time(run_old)); invisible(.time(run_v2)); invisible(.time(run_v3))  # warmup
  
  old <- replicate(n, .time(run_old), simplify = FALSE)
  v2  <- replicate(n, .time(run_v2),  simplify = FALSE)
  v3  <- replicate(n, .time(run_v3),  simplify = FALSE)
  
  vec <- function(lst, key) vapply(lst, function(x) x[[key]], numeric(1))
  tab <- data.frame(
    version = c("old","v2","v3"),
    median  = c(median(vec(old,"sec")), median(vec(v2,"sec")), median(vec(v3,"sec"))),
    min     = c(min(vec(old,"sec")),    min(vec(v2,"sec")),    min(vec(v3,"sec"))),
    max     = c(max(vec(old,"sec")),    max(vec(v2,"sec")),    max(vec(v3,"sec")))
  )
  
  # micro breakouts for v2/v3 (old has no timings)
  get_timing <- function(lst, name) median(sapply(lst, function(x) x$res$timings[[name]]))
  micros <- data.frame(
    version   = c("v2","v3"),
    prod1_med = c(get_timing(v2,"prod1"),  get_timing(v3,"prod1")),
    loglik_med= c(get_timing(v2,"loglik"), get_timing(v3,"loglik"))
  )
  
  list(summary = tab, micros = micros)
}

## Run
bench <- bench_occur_versions3(n = 5, seed = 1)
bench$summary
bench$micros

############### EQUIVALENCE ###################################################
## One-shot equivalence check on a frozen snapshot (P-block)
check_occur_equivalence <- function(seed = 1L) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  
  # snapshot current inputs
  snap <- list(
    p_curr   = this_occur_P,
    occ_curr = this_O_P,
    p_prior  = p_occur_P,
    pj       = this_pj,
    pi       = this_pi,
    occ_Other= this_O_V,
    L        = this_L,
    det      = detected_P,
    rows_idx = focus_rows_idx_P,
    cols_idx = focus_cols_idx_P,
    focus_i  = focus_int_P
  )
  
  run_old <- function() UpdOccurP_blocked(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = snap$p_curr, occ_curr = snap$occ_curr,
    occur_prior_probs = snap$p_prior,
    probobs_curr = snap$pj, probobs_others = snap$pi,
    occur_others = snap$occ_Other,
    curr_inter = t(snap$L),
    focus = aperm(focus, c(2,1,3)),
    detected = snap$det
  )
  
  run_v2 <- function() UpdOccurP_blocked_timing_v2(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr = snap$p_curr, occ_curr = snap$occ_curr,
    occur_prior_probs = snap$p_prior,
    probobs_curr = snap$pj, probobs_others = snap$pi,
    occur_others = snap$occ_Other,
    curr_inter = t(snap$L),
    focus_rows_idx_list = snap$rows_idx,
    focus_cols_idx_list = snap$cols_idx,
    focus_int_list      = snap$focus_i,
    detected = snap$det
  )
  
  old <- run_old()
  set.seed(seed)               # same RNG stream for v2
  v2  <- run_v2()
  
  # compare ONLY where not detected in that study
  M <- (snap$det == 0L)
  
  max_abs_dp <- max(abs(old$p_curr[M] - v2$p_curr[M]))
  occ_mm     <- sum(old$occ_curr[M] != v2$occ_curr[M])
  acc_mm     <- sum(old$accepted[M] != v2$accepted[M])
  
  cat(sprintf(
    "max|Δp| (non-detected) = %.3e\nocc mismatches (non-detected) = %d\naccepted mismatches (non-detected) = %d\n",
    max_abs_dp, occ_mm, acc_mm
  ))
  
  # simple pass/fail (tweak tol if you like)
  tol <- 1e-12
  ok <- (max_abs_dp <= tol) && occ_mm == 0 && acc_mm == 0
  if (ok) message("✅ Equivalent on non-detected rows.")
  else    message("⚠️ Differences found on non-detected rows.")
  invisible(list(max_abs_dp=max_abs_dp, occ_mm=occ_mm, acc_mm=acc_mm))
}

## run once
check_occur_equivalence(seed = 1)

