# This file runs COIL+ across four chains using: 
  # Expert-defined occurrence probabilities
  # Blocked sampling of occurrence probabilities and indicators

# ---------------------------------- TO DO ------------------------------------#


  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, along = 4)
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
  
  # Running mean of detection probabilities
  p_detect <- list(pis = mcmc$pi_mean, pjs = mcmc$pj_mean)
  
  # Running mean of latent factors
  factors <- list(U = mcmc$U_mean, V = mcmc$V_mean)
  
  # Imputed values of missing covariates
  Xs <- mcmc$Xs
  Ws <- mcmc$Ws
  
  # Occurrence indicators, probabilities, acceptance rates (last two not relevant if we don't sample P)
  occ_plants <- list(OP_mean = mcmc$OP_mean, p_OPs = mcmc$p_OP_mean, p_accept = mcmc$p_OP_accepted) 
  occ_verts <- list(OV_mean = mcmc$OV_mean, p_VPs = mcmc$p_OV_mean, p_accept = mcmc$p_OV_accepted) 
  
  # Log likelihood
  logL <- mcmc$logL
  
  # Combining the results we are interested in to a list and saving:
  res <- list(all_pred = all_pred, logL = logL, correlations = correlations, 
              p_detect = p_detect, factors = factors, occ_plants = occ_plants, 
              occ_verts = occ_verts,  Xs = Xs, Ws = Ws)
  
 # save(res, file = each_filepath)
  
  rm(res)



