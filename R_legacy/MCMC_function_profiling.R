# This versions tests out the *alt version of the blocked sampler which corrects possible issue in asymmetric proposal
# And also saves unthinned quantities of interest for traceplots including logL


#' MCMC for bipartite network model with trait information and unrecorded
#' interactions
#' 
#' Perform MCMC to acquire samples from the posterior distribution of model
#' parameters for a model for the true interaction matrix among different sets
#' of species with available trait information.
#' 
#' @param obs_A Observed interaction array. Contains values of 0 and 1. The
#' array is of three dimensions. The first dimension represents the first set
#' of species, the second dimension is the second set of species, and the third
#' dimension is the different studies.
#' @param focus Array of three dimensions corresponding to the two sets of
#' species, and the different studies. Values are 0 or 1. This array represents
#' whether the current study would be willing to record the interactions of a
#' given species. The value should be 1 except for studies that are animal or
#' plant-oriented, where the value for species that are not of focus will be 0.
#' @param p_occur_V Matrix with rows corresponding to the first set of species
#' and columns to studies. Values represent the prior probability of occurence, and
#' they should be between 0 and 1.
#' @param p_occur_P Same as occur_B but for the second set of species.
#' @param obs_X Matrix of observed covariates for the first set of species.
#' Rows correspond to species, columns to covariates. Continuous covariates
#' should be first, and binary covariates should follow.
#' @param obs_W Same structure as obs_X but for the second set of species.
#' @param Cu Phylogenetic correlation matrix for the first set of species.
#' @param Cv Phylogenetic correlation matrix for the second set of species.
#' @param Nsims Number of posterior samples to be kept.
#' @param burn Number of samples to be burnt in the beginning of the MCMC.
#' @param thin MCMC sampling thinning.
#' @param use_H Number of latent factors to be used. Defaults to 10.
#' @param use_shrinkage Logical. Whether the shrinkage prior on the variance 
#' terms should be used. Defaults to TRUE.
#' @param bias_cor Logical. Whether the model should aim to correct for
#' geographical and taxonomical bias. If set to FALSE, the MCMC is simplied
#' since the observed matrix of interactions is considered the same as the
#' true matrix of interactions. Defaults to TRUE.
#' @param theta_inf Variance value theta infinity in the spike part of the
#' prior distribution for increasing shrinkage. Defaults to 0.01.
#' @param mh_n_pis Parameter n in the Beta proposal for updating pi. Defaults
#' to 100.
#' @param mh_n_pjs Same as mh_n_pis but for the update of pj.
#' @param mh_n_rho Same as mh_n_pis but for the update of rho.
#' @param mh_p_step Step size in updating occurrence probs under blocked and unblocked sampling
#' @param mh_pprior_sd Prior sd in trunc normal prior for probs under blocked and unblocked sampling
#' @param p_1to0 Stepsize, e.g. probability of flipping 1 to 0 for occurrence indicator
#' (blocked sampling) 
#'  @param p_0to1 Stepsize, e.g. probability of flipping 0 to 1 for occurrence indicator
#' (blocked sampling) 
#' @param stick_alpha Alpha value in the increasing shrinkage prior. Defaults
#' to 5.
#' @param prior_theta Hyperparameters of the inverse gamma distribution in the
#' slab part of the increasing shrinkage prior. Defaults to (1, 1).
#' @param pruir_tau Hyperarameters of the inverse gamma prior for the
#' additional variance parameter in the coefficients. Defaults to (5, 5).
#' @param prior_rho Hyperparameters of the beta prior for the weight rho in 
#' the correlation matrix for latent factors. Defaults to (5, 5).
#' @param prior_mu0 Mean of the normal prior distribution for all intercepts.
#' Defaults to 0.
#' @param prior_sigmasq0 Variance of the normal prior distribution for all
#' intercepts. Defaults to 10.
#' @param prior_sigmasq Hyperparameters of the inverse gamma prior on the
#' variance terms of continuous traits. Defaults to (1, 1).
#' @param start_values List that can include starting values for any subset of
#' parameters. Starting values for parameters that are not specified in
#' start_values will be sampled. Defaults to NULL.
#' @param sampling List specifying which parameters should be sampled by
#' setting the value to TRUE, and which should not and should be kept at
#' their starting value by setting the corresponding list element to FALSE.
#' If sampling is set to FALSE for a parameter, it is recommended that the
#' corresponding start_value is specified to be set to the parameter's true
#' value. Defaults to NULL, and when set to NULL all parameters are sampled.
#' @param cut_feed Logical. Whether we should cut the feedback from the
#' interaction and detectability submodels into the latent factors. Defaults to
#' FALSE, in which case the full posterior is considered.
#' 
#' @export
#' 
#' 
#' 

MCMC <- function(obs_A, focus, p_occur_V, p_occur_P, obs_X, obs_W, Cu, Cv,
                 Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
                 bias_cor = TRUE, theta_inf = 0.01,
                 mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
                 mh_p_step = 0.1, mh_pprior_sd = 0.1, 
                 p_1to0 = 0.65, p_0to1 = 0.25, 
                 stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                 prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                 prior_sigmasq = c(1, 1), start_values = NULL,
                 sampling = NULL, 
                 cut_feed = FALSE, block_sampleOccP = FALSE, 
                 save_logL = FALSE, save_rhos = FALSE) {

  
  # -------------------- PART 0 ------------------- #
  # ------- Specifying what will be sampled. ------- #
  
  if (is.null(sampling)) {  # All parameters will be updated.
    sampling <- list(L = TRUE, lambda = TRUE, tau = TRUE, beta = TRUE,
                     gamma = TRUE, sigmasq = TRUE, sigmasq_p = TRUE,
                     delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                     z = TRUE, theta = TRUE, pis = TRUE, pjs = TRUE, rU = TRUE,
                     rV = TRUE, miss_X = TRUE, miss_W = TRUE, O_V = TRUE,
                     O_P = TRUE, p_OV = TRUE, p_OP = TRUE)
  }
  
  # Return error if block updating is specified, but occurrence probs are not sampled
  if(block_sampleOccP==TRUE & (sampling$p_OV == FALSE & sampling$p_OP == FALSE)){
    stop("Error: one or both of sampling$p_OV and sampling$p_OP must be TRUE in
            order to use block updating of occurrence indicators and probabilities")
  }
  
  # Without bias adjustment, some parameters are not updated:
  if (!bias_cor) {
    cat('Without bias correction a number of parameters will not be sampled.', fill = TRUE)
    sampling$L <- FALSE
    sampling$sigmasq_p <- FALSE
    sampling$delta <- FALSE
    sampling$zeta <- FALSE
    sampling$pis <- FALSE
    sampling$pjs <- FALSE
    sampling$O_V <- FALSE
    sampling$O_P <- FALSE
  }
  
  # If the shrinkage prior is NOT used, some parameters will not be updated.
  if (!use_shrinkage) {
    cat('Without the shrinkage prior, some parameters are not updated.', fill = TRUE)
    sampling$v <- FALSE
    sampling$z <- FALSE
  }
  
  # If the probabilities of occurence for the first set of species is either 0
  # or 1, we never have to update their occurrence indicator or occurrence probability.
  
  if(sum(p_occur_V==0) + sum(p_occur_V==1) == length(p_occur_V)){
    sampling$O_V <- sampling$p_OV <- FALSE
  }
  
  # Similarly for the second set of species:
  if(sum(p_occur_P==0) + sum(p_occur_P==1) == length(p_occur_P)){
    sampling$O_P <- sampling$p_OP <- FALSE
  }


  # ---------------- PART 1 ------------- #
  # Getting the parameters that we use throughout:
  
  # Sample size
  nB <- dim(obs_A)[1]
  nP <- dim(obs_A)[2]
  nstudies <- dim(obs_A)[3]
  
  cat('MCMC on', nB, 'x', nP, 'number of species.', fill = TRUE)
  
  # Getting the combined network:
  comb_A <- apply(obs_A, c(1, 2), sum)
  comb_A <- (comb_A > 0) * 1
  
  # Which species in the combined network have unrecorded interactions.
  index_A0 <- which(comb_A == 0)
  quant_A0 <- length(index_A0)
  
  # Getting whether the species were detected in each study.
  # Note that in some cases, species are detected in the area but not in interactions,
  # Particularly if the data is subset for analysis
  # In this case the probability of occurrence is 1
  detected_B <- (apply(obs_A, c(1, 3), sum) > 0) * 1
  detected_P <- (apply(obs_A, c(2, 3), sum) > 0) * 1
  
  detected_B_occ <- ifelse(p_occur_V == 1, 1, 0)
  detected_P_occ <- ifelse(p_occur_P == 1, 1, 0)
  
  if(sum(detected_B > detected_B_occ)){warning("Error in p_occur_V: All detected species should have occurrence probability = 1.")}
  if(sum(detected_P > detected_P_occ)){warning("Error in p_occur_P: All detected species should have occurrence probability = 1.")}
  
  if(sum(detected_B < detected_B_occ)){
    message("Note on p_occur_V: Some species with occurrence probability = 1 do not have observed interactions (likely due to subsetting)")
    }
  if(sum(detected_P < detected_P_occ)){
    message("Note on p_occur_P: Some species with occurrence probability = 1 do not have observed interactions (likely due to subsetting)")
    }
  
  detected_B[detected_B_occ == 1] <- 1
  detected_P[detected_P_occ == 1] <- 1
  
  # Continuous and binary covariates.
  entries_X <- apply(obs_X, 2, function(x) length(unique(x[!is.na(x)]))) # unique non-NA values 
  pB <- c(sum(entries_X > 2), sum(entries_X == 2)) # number continuous, number binary
  entries_W <- apply(obs_W, 2, function(x) length(unique(x[!is.na(x)])))
  pP <- c(sum(entries_W > 2), sum(entries_W == 2))
  # Making sure covariates are ordered as continuous first.
  if (any(c(entries_X[1 : pB[1]] == 2, entries_W[1 : pP[1]] == 2))) {
    stop('Reorder covariates')
  }
  
  # If the covariates have missing values, we check for it here.
  any_X_miss <- any(apply(obs_X, 2, function(x) sum(is.na(x))) > 0)
  any_W_miss <- any(apply(obs_W, 2, function(x) sum(is.na(x))) > 0)
  
  # Indices with missing values:
  if (any_X_miss) {
    miss_X_ind <- apply(obs_X, 2, function(x) which(is.na(x)))
  }
  if (any_W_miss) {
    miss_W_ind <- apply(obs_W, 2, function(x) which(is.na(x)))
  }
  
  
  # ---------------- PART 2 ------------- #
  # Creating the elements where the MCMC samples will be saved.
  # Reduce the dimension of stored elements
  
  Ls <- array(NA, dim = c(Nsims, nB, nP))  # True interaction matrix.
  mod_pL1s <- array(NA, dim = c(Nsims, nB, nP))  # Probability of interaction without correction.
  #pL1s <- array(NA, dim = c(Nsims, nB, nP))  # Probability of interaction. #TRIM
  
  # Don't store latent factors or marginal probs for now 
  #Us <- array(NA, dim = c(Nsims, nB, use_H))  # Bird latent factors. #TRIM
  #Vs <- array(NA, dim = c(Nsims, nP, use_H))  # Plant latent factors. #TRIM
  U_cumsum <- array(0, dim = c(nB, use_H)) #TRIM
  V_cumsum <- array(0, dim = c(nP, use_H)) #TRIM
  
  lambdas <- array(NA, dim = c(Nsims, use_H + 1))  # Coefficients, network model.
  taus_lambda <- array(NA, dim = c(Nsims, use_H))  # Extra variance, network model.
  betas <- array(NA, dim = c(Nsims, sum(pB), use_H + 1))  # Coefficients, bird traits.
  gammas <- array(NA, dim = c(Nsims, sum(pP), use_H + 1))  # Coefficients, plant traits.
  taus_beta <- array(NA, dim = c(Nsims, sum(pB), use_H))  # Extra variance, bird traits.
  taus_gamma <- array(NA, dim = c(Nsims, sum(pP), use_H))  # Extra variance, plant traits.
  sigmasq_m <- array(NA, dim = c(Nsims, pB[1]))  # Residual variance, bird traits.
  sigmasq_l <- array(NA, dim = c(Nsims, pP[1]))  # Residual variance, plant traits.
  deltas <- array(NA, dim = c(Nsims, use_H + 1))  # Coefficients, observing bird.
  taus_delta <- array(NA, dim = c(Nsims, use_H))  # Extra variance, observing bird.
  zetas <- array(NA, dim = c(Nsims, use_H + 1))  # Coefficients, observing plant.
  taus_zeta <- array(NA, dim = c(Nsims, use_H))  # Extra variance, observing plant.
  thetas <- array(NA, dim = c(Nsims, use_H))  # Shrinking variance for loadings.
  vs <- array(NA, dim = c(Nsims, use_H))  # Stick breaking prior for shrinking variance.
  omegas <- array(NA, dim = c(Nsims, use_H))  # Length of broken sticks.
  zs <- array(NA, dim = c(Nsims, use_H))  # Which stick break the H components belong to.
  #pis <- array(NA, dim = c(Nsims, nB))  # Probabilities of observing bird.#TRIM
  pi_cumsum <- rep(0, nB) #TRIM
  pi_accepted <- rep(0, nB)  # Number of MCMC iterations for which the bird proposal was accepted.
  #pjs <- array(NA, dim = c(Nsims, nP))  # Probabilities of observing plant.#TRIM
  pj_accepted <- rep(0, nP)  # Number of MCMC iterations for which the plant proposal was accepted.
  pj_cumsum <- rep(0, nP) #TRIM
  sigmasq_pB <- rep(NA, Nsims)  # Residual variance, observing bird.
  sigmasq_pP <- rep(NA, Nsims)  # Residual variance, observing plant.
 
  p_OV_accepted <- array(0, dim = c(nB,nstudies)) # Number of MCMC iteractions for which for bird occurrence prob proposal was accepted
  p_OV_cumsum <- array(0, dim = c(nB, nstudies))
  occP_B_accepted <- array(0, dim = c(nB, nstudies))
  #p_OVs <- array(NA, dim = c(Nsims, nB, nstudies)) # Trim this after debugging
  p_OP_accepted <- array(0, dim = c(nP,nstudies)) # Number of MCMC iteractions for which for plant occurrence prob proposal was accepted
  p_OP_cumsum <- array(0, dim = c(nP, nstudies))
  occP_P_accepted <- array(0, dim = c(nP, nstudies))
  #p_OPs <- array(NA, dim = c(Nsims, nP, nstudies)) # Trim this after debugging
  OV_cumsum <- array(0, dim = c(nB, nstudies))
  OP_cumsum <- array(0, dim = c(nP, nstudies))
  
  # Setup storage for imputed NA covariates at each iteration
  if (any_X_miss) {
    Xs <- lapply(miss_X_ind, function(x) matrix(NA, nrow = Nsims, ncol = length(x)))
  }
  if (any_W_miss) {
    Ws <- lapply(miss_W_ind, function(x) matrix(NA, nrow = Nsims, ncol = length(x)))
  }
  
  # Setup storage for trace quantities
  rU <- rep(NA, Nsims*thin + burn)  # Correlation, latent factors for frugivores.
  ru_accepted <- 0
  rV <- rep(NA, Nsims*thin + burn)  # Correlation, latent factors for plants.
  rv_accepted <- 0
  
  logL <- rep(NA, Nsims*thin + burn)
  
  # ---------------- PART 3 ------------- #
  # Setting starting values for the parameters.
  
  this_L <- matrix(rbinom(nB * nP, 1, 1 / 2), nB, nP)
  this_L[comb_A == 1] <- 1
  this_U <- matrix(rnorm(nB * use_H), nB, use_H)
  this_V <- matrix(rnorm(nP * use_H), nP, use_H)
  this_lambda <- rnorm(use_H + 1, 0, 1)
  this_tau_lambda <- 1 / rgamma(use_H, 10, 10)
  this_beta <- matrix(rnorm(sum(pB) * (use_H + 1)), sum(pB), use_H + 1)
  this_gamma <- matrix(rnorm(sum(pP) * (use_H + 1)), sum(pP), use_H + 1)
  this_tau_beta <- matrix(1 / rgamma(sum(pB) * use_H, 10, 10), sum(pB), use_H)
  this_tau_gamma <- matrix(1 / rgamma(sum(pP) * use_H, 10, 10), sum(pP), use_H)
  this_sigmasq_m <- 1 / rgamma(pB[1], 10, 10)
  this_sigmasq_l <- 1 / rgamma(pP[1], 10, 10)
  this_delta <- rnorm(use_H + 1, 0, 1)
  this_tau_delta <- 1 / rgamma(use_H, 10, 10)
  this_zeta <- rnorm(use_H + 1, 0, 1)
  this_tau_zeta <- 1 / rgamma(use_H, 10, 10)
  this_theta <- sort(1 / rgamma(use_H, 5, 5), decreasing = TRUE)
  this_v <- c(rbeta(use_H - 1, 1, 3), 1)
  this_sigmasq_pB <- 1 / rgamma(1, 10, 10)
  this_sigmasq_pP <- 1 / rgamma(1, 10, 10)
  this_pi <- runif(nB, 0, 1)
  this_pj <- runif(nP, 0, 1)
  this_ru <- rbeta(1, 5, 5)
  this_rv <- rbeta(1, 5, 5)
  this_X <- obs_X
  this_W <- obs_W
  this_mod_pL1 <- matrix(NA, nrow = nB, ncol = nP)
  this_pL1 <- matrix(NA, nrow = nB, ncol = nP)

  if(sampling$p_OV == TRUE){
    this_occur_B <- matrix(rtruncnorm(n = nB*nstudies, mean = p_occur_V, sd = mh_pprior_sd, a = 0, b = 1), 
                           nrow = nB, ncol = nstudies) # Note naming convention is off, can fix later if wanted (this_occur_B)
    this_occur_B[p_occur_V == 1] <- 1
  }else {this_occur_B <- p_occur_V} # New so that we don't randomly initialize the probs when we aren't sampling them

  if(sampling$p_OP == TRUE){
    this_occur_P <- matrix(rtruncnorm(n = nP*nstudies, mean = p_occur_P, sd = mh_pprior_sd, a = 0, b = 1), 
                           nrow = nP, ncol = nstudies) 
    this_occur_P[p_occur_P == 1] <- 1    
  } else {this_occur_P <- p_occur_P}  # New so that we don't randomly initialize the probs when we aren't sampling them

  
  # ------- Dealing with individual study focus and co-occurence.
  # Indicator of occurrence for each set of species:
  this_O_V <- matrix(rbinom(nB * nstudies, 1, prob = p_occur_V), nB, nstudies) # New: initialize with prior probs
  this_O_P <- matrix(rbinom(nP * nstudies, 1, prob = p_occur_P), nP, nstudies) # New: initialize with prior probs 
  # Indicator of co-occurence:
  this_O <- as.numeric(do.call(abind::abind, c(lapply(seq(nstudies), function(ss)
    outer(this_O_V[,ss], this_O_P[,ss])), along = 3)))
  
  # Number of studies with focus and occurence with recorded interaction
  this_AFO <- apply(obs_A * focus * this_O, c(1, 2), sum)
  # Number of studies with focus and occurence without recorded interaction
  this_1_AFO <- apply((1 - obs_A) * focus * this_O, c(1, 2), sum)
  # Total number of studies with focus and occurence
  this_FO <- this_AFO + this_1_AFO
  
  # Draw starting values for missing covariates from observed distribution.
  if (any_X_miss) {
    # Continuous X covariates
    if(pB[1] != 0){
    for (jj in 1 : pB[1]) {
      use_stats <- c(mean(obs_X[, jj], na.rm = TRUE), sd(obs_X[, jj], na.rm = TRUE))
      this_X[miss_X_ind[[jj]], jj] <- rnorm(length(miss_X_ind[[jj]]),
                                            use_stats[1], use_stats[2])
    }
    }
    # Binary X covariates
    if(pB[2] != 0){
      for (jj in (pB[1] + 1) : sum(pB)) {
      use_stats <- mean(obs_X[, jj], na.rm = TRUE)
      this_X[miss_X_ind[[jj]], jj] <- rbinom(length(miss_X_ind[[jj]]), 1, use_stats)
        }
    }
  }
  
  if (any_W_miss) {
    # Continuous W covariates
    if(pP[1] != 0){
    for (jj in 1 : pP[1]) {
      use_stats <- c(mean(obs_W[, jj], na.rm = TRUE), sd(obs_W[, jj], na.rm = TRUE))
      this_W[miss_W_ind[[jj]], jj] <- rnorm(length(miss_W_ind[[jj]]),
                                            use_stats[1], use_stats[2])
    }
    }
    # Binary W covariates
    if(pP[2] != 0){
    for (jj in (pP[1] + 1) : sum(pP)) { 
      use_stats <- mean(obs_W[, jj], na.rm = TRUE)
      this_W[miss_W_ind[[jj]], jj] <- rbinom(length(miss_W_ind[[jj]]), 1, use_stats)
    }
    }
  }
  
  # If the starting values have been specified in start_values, use those. 
  # ERROR: START VALUES ARE NAMED WRONG BY DEFAULT AND WON'T BE USED
  # If you don't sample, they will just be fixed at the random initialization above
  # I didn't correct, I just renamed my list of start values
  if (!is.null(start_values)) {
    for (pp in 1 : length(start_values)) {
      assign(x = names(start_values)[pp], start_values[[pp]])
    }
  }
  
  this_Su <- this_ru * Cu + (1 - this_ru) * diag(nB) # DEBUG: THIS ISN'T ACTUALLY USING START VALUES IF THEY'RE SET AS IN DOCUMENTATION
  this_Sv <- this_rv * Cv + (1 - this_rv) * diag(nP)
  
  # If we do not perform bias correction, change the starting values:
  if (!bias_cor) {
    this_L <- comb_A
    this_delta <- rep(NA, use_H + 1)
    this_tau_delta <- rep(NA, use_H)
    this_zeta <- rep(NA, use_H + 1)
    this_tau_zeta <- rep(NA, use_H)
    this_sigmasq_pB <- NA
    this_sigmasq_pP <- NA
    this_pi <- rep(1, nB)
    this_pj <- rep(1, nP)
    this_AFO <- matrix(NA, nB, nP)  # Does not inform anything.
    this_1_AFO <- matrix(NA, nB, nP)  # Does not inform anything.
    this_FO <- matrix(NA, nB, nP)  # Does not inform anything.
  }
  
  # If we do not use the shrinkage prior, change the starting values:
  if (!use_shrinkage) {
    this_v <- rep(0, use_H)
  }
  
  this_omega <- OmegaFromV(v_val = this_v)
  if (sum(this_omega) > 0) {  # Shrinkage prior is used.
    this_z <- sample(c(1 : use_H), use_H, replace = TRUE, prob = this_omega)
  } else {  # When shrinkage prior NOT used: thetas from slab, z not updated.
    this_z <- rep(use_H + 1, use_H)
  }
  
  # ---------------- PART 4 ------------- #
  # Performing the MCMC.
  
  keep_index <- 1
  timings <- list()  # container for timing
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  track_time <- function(name, expr) {
    t0 <- Sys.time()
    out <- eval(expr)
    timings[[name]] <<- c(timings[[name]] %||% numeric(),
                          as.numeric(difftime(Sys.time(), t0, units = "secs")))
    return(out)
  }
  
  cat('Total number of iterations:', Nsims * thin + burn, fill = TRUE)
  
  for (ss in 1:(Nsims * thin + burn)) {
    if (ss %% 100 == 0) print(ss)
    
    # ------ L: Update true interactions -------- #
    if (sampling$L | sampling$lambda | sampling$U | sampling$V) {
      track_time("pL1_calc", quote({
        logit_pijL <- matrix(this_lambda[1], nB, nP)
        for (hh in 1:use_H) {
          lat_prod <- matrix(this_U[, hh], ncol = 1) %*% matrix(this_V[, hh], nrow = 1)
          logit_pijL <- logit_pijL + this_lambda[hh + 1] * lat_prod
        }
        this_mod_pL1 <- expit(logit_pijL)
      }))
    }
    
    if (sampling$L) {
      track_time("L_update", quote({
        pipj <- outer(this_pi, this_pj)
        pL1 <- this_mod_pL1 * (1 - pipj) ^ this_FO
        pL0 <- 1 - this_mod_pL1
        this_pL1 <- pL1 / (pL0 + pL1)
        this_L <- matrix(1, nB, nP)
        this_L[index_A0] <- rbinom(n = quant_A0, 1, prob = this_pL1[index_A0])
      }))
    }
    
    # ------- lambda update -------
    if (sampling$lambda | sampling$U) {
      omega_L <- track_time("omega_L", quote(
        BayesLogit::rpg(nB * nP, h = rep(1, nB * nP), z = as.numeric(logit_pijL))
      ))
    }
    
    if (sampling$lambda) {
      track_time("lambda_update", quote({
        des_mat <- matrix(1, nB * nP)
        for (hh in 1:use_H) {
          lat_prod <- matrix(this_U[, hh], ncol = 1) %*% matrix(this_V[, hh], nrow = 1)
          des_mat <- cbind(des_mat, as.numeric(lat_prod))
        }
        prior_m <- c(prior_mu0, rep(0, use_H))
        prior_S_inv <- diag(1 / c(prior_sigmasq0, this_tau_lambda * this_theta))
        new_S <- sweep(t(des_mat), 2, FUN = '*', omega_L) %*% des_mat
        new_S <- chol2inv(chol(new_S + prior_S_inv))
        new_m <- t(des_mat) %*% matrix(as.numeric(this_L - 1 / 2), ncol = 1)
        new_m <- new_S %*% (new_m + prior_S_inv %*% prior_m)
        this_lambda <- mvnfast::rmvn(1, new_m, new_S)
      }))
    }
    
    # ------- U latent factors -------
    if (sampling$U) {
      track_time("U_update", quote({
        this_Su_inv <- chol2inv(chol(this_Su))
        this_U <- UpdLatFac(
          latfac = this_U, latfac_others = this_V, probobs = this_pi,
          coefs_probobs = this_delta, var_probobs = this_sigmasq_pB,
          obs_covs = this_X, omega_obs_covs = omega_obsX,
          num_covs = pB, coefs_covs = this_beta,
          var_covs = this_sigmasq_m,
          curr_inter = this_L, coefs_inter = this_lambda,
          omega_inter = matrix(omega_L, nrow = nB, ncol = nP),
          prior_S_inv = this_Su_inv,
          cut_feed = cut_feed
        )
      }))
    }
    
    # ------- V latent factors -------
    if (sampling$V) {
      track_time("V_update", quote({
        this_Sv_inv <- chol2inv(chol(this_Sv))
        this_V <- UpdLatFac(
          latfac = this_V, latfac_others = this_U, probobs = this_pj,
          coefs_probobs = this_zeta, var_probobs = this_sigmasq_pP,
          obs_covs = this_W, omega_obs_covs = omega_obsW,
          num_covs = pP, coefs_covs = this_gamma,
          var_covs = this_sigmasq_l,
          curr_inter = t(this_L), coefs_inter = this_lambda,
          omega_inter = t(matrix(omega_L, nrow = nB, ncol = nP)),
          prior_S_inv = this_Sv_inv,
          cut_feed = cut_feed
        )
      }))
    }
    
    # ------- Occurrence indicators & probs -------
    if (sampling$O_V | sampling$O_P | sampling$p_OV | sampling$p_OP) {
      track_time("Occur_update", quote({
        # leave your full existing occurrence update block here unchanged
      }))
    }
    
    # ------- Log likelihood -------
    if (save_logL) {
      track_time("logL_calc", quote({
        this_pipj <- outer(this_pi, this_pj)
        log_this_pipj <- log(this_pipj)
        log_1_minus_this_pipj <- log(1 - this_pipj)
        this_FO_bin <- focus * this_O
        thislogL <- 0
        for (st in 1:nstudies) {
          thislogL <- thislogL +
            sum(obs_A[, , st] * log_this_pipj) +
            sum((1 - obs_A[, , st]) * log_1_minus_this_pipj)
        }
        logL[ss] <- thislogL
      }))
    }
    
    # ------- Save results -------
    if ((ss - burn) %% thin == 0 & ss > burn) {
      track_time("save_results", quote({
        Ls[keep_index, , ] <- this_L
        lambdas[keep_index, ] <- this_lambda
        keep_index <- keep_index + 1
      }))
    }
  }
  
  # ---------- PART 5 ---------- #
  # Returning the results:
  
  r <- list(timings = timings,
            Ls = Ls, 
            logL = logL,
            #pL1s = pL1s, #TRIM
            mod_pL1s = mod_pL1s, 
            lambdas = lambdas, taus_beta = taus_beta,
            taus_gamma = taus_gamma, taus_lambda = taus_lambda,
            taus_delta = taus_delta, taus_zeta = taus_zeta,
            betas = betas, gammas = gammas, sigmasq_m = sigmasq_m,
            sigmasq_l = sigmasq_l, sigmasq_pB = sigmasq_pB,
            sigmasq_pP = sigmasq_pP, deltas = deltas, zetas = zetas,
            # Us = Us, Vs = Vs, #TRIM
            U_mean = U_mean, V_mean = V_mean, #TRIM
            vs = vs, 
            omegas = omegas, zs = zs,
            thetas = thetas, 
            p_OV_mean = p_OV_mean, # running mean of occurrence probabilities
            p_OP_mean = p_OP_mean,
            OP_mean = OP_mean, # running mean of occurrence indicators
            OV_mean = OV_mean,
            pi_mean = pi_mean, pj_mean = pj_mean, #TRIM
            rU = rU, rV = rV,
            pi_accepted = pi_accepted / (Nsims * thin + burn),
            pj_accepted = pj_accepted / (Nsims * thin + burn), 
            ru_accepted = ru_accepted / (Nsims * thin + burn),
            rv_accepted = rv_accepted / (Nsims * thin + burn), 
            p_OV_accepted = p_OV_accepted/(Nsims * thin + burn), 
            p_OP_accepted = p_OP_accepted/(Nsims * thin + burn)            
            )
  
  if (any_X_miss) {
    r$Xs <- Xs
  }
  if (any_W_miss) {
    r$Ws <- Ws
  }
  
  return(r)
  
}