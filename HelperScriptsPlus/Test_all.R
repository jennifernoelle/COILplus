# This script tests AMH RWMH on simulated data
# With and without prior misspecification

# Debugging code
# X Check that recursive computation matches iterative
# X Write functions for logpost, etc
# X Fill in AMH skeleton, can delete iterative covar
# X Write function to simulate data separately
# X Debug: why does my new accept/reject code fail - compare to old way
# Modify code so we aren't sampling known occurrences
# X Look at which plants-sites get stuck at zero
# X Experiment with MH parameters to get better acceptance rate
# X Add non-adaptive phase
# X Turn AMH into general function
# X Turn RW into general function
# X Fix p_occur_others vs occur_others confusion
# Add code so we aren't sampling occurrence for detected species to all three functions <- DEUBG this, we have suspicious posterior 1's
# X Create unblocked version as a general function
# Try Jondrow modification 
# X Create blocked version with independent sampling of pi, occ
# Create simulations with prior misspecification 
# Run with more realistic simulations (e.g. larger networks and sample size, sparsity)
# Modify for full sampler:
# Only running one iteration at each gibbs step
# We don't sample known 0/1 occurrences

# Main difference between the blocked and unblocked samplers
# Unblocked sampler has a Gibbs step for occurrence indicators and trunc normal for pi
# Blocked samplers are on the probit scale for both (pnorm for prob, indicator for occ),
# so it seems easier to get stuck at 0/1
# Unblocked/sequential sampler is much faster

# Tuning questions
# Acceptance rate is high, but if I increase the proposal variance I'm afraid of 
# getting stuck at 0/1


library(pROC)
#--------------------------- NO PRIOR MISSPECIFICATION ------------------------#

#  We simulate from the prior, but note that in the simulation function, 
#  Prior probs are reset to 1 for detected species for simplification
#  This is as it would be in the input data

### 0. Load functions

# This script tests AMH RWMH on simulated data
# With and without prior misspecification

# Debugging code
# X Check that recursive computation matches iterative
# X Write functions for logpost, etc
# X Fill in AMH skeleton, can delete iterative covar
# X Write function to simulate data separately
# X Debug: why does my new accept/reject code fail - compare to old way
# Modify code so we aren't sampling known occurrences
# X Look at which plants-sites get stuck at zero
# X Experiment with MH parameters to get better acceptance rate
# X Add non-adaptive phase
# X Turn AMH into general function
# X Turn RW into general function
# X Fix p_occur_others vs occur_others confusion
# Add code so we aren't sampling occurrence for detected species to all three functions <- DEUBG this, we have suspicious posterior 1's
# X Create unblocked version as a general function
# Try Jondrow modification 
# X Create blocked version with independent sampling of pi, occ
# Create simulations with prior misspecification 
# Run with more realistic simulations (e.g. larger networks and sample size, sparsity)
# Modify for full sampler:
  # Only running one iteration at each gibbs step
  # We don't sample known 0/1 occurrences

# Think about if comparative results will differ when embedded in a larger gibbs sampler
# vs fixing L, etc at truth

# Main difference between the blocked and unblocked samplers
  # Unblocked sampler has a Gibbs step for occurrence indicators and trunc normal for pi
  # Blocked samplers are on the probit scale for both (pnorm for prob, indicator for occ),
  # so it seems easier to get stuck at 0/1
  # Unblocked/sequential sampler is much faster

# Tuning questions
  # Acceptance rate is high, but if I increase the proposal variance I'm afraid of 
  # getting stuck at 0/1


library(pROC)

source_path <- 'HelperScriptsPlus/' 
source(paste0(source_path, "Utils.R"))
source(paste0(source_path, "Sequential_MH_PiOp.R"))
source(paste0(source_path, "AMH_PiOp.R"))
source(paste0(source_path, "RW_blocked_PiOp.R"))
source(paste0(source_path, "Blocked_indpt_PiOp.R"))

#--------------------------- NO PRIOR MISSPECIFICATION ------------------------#

#  We simulate from the prior, but note that in the simulation function, 
#  Prior probs are reset to 1 for detected species for simplification
#  This is as it would be in the input data

### 1. Simulate data without prior misspecification

# Set up key external function inputs: these come from inputs and current state
set.seed(1234)
nS <- 100
nM <- 10
nP <- 100
pi_L <- 0.25 # Sets sparsity of the true interactions matrix
mh_pprior_sd <- 0.1

sim.data <- simulate_plants(nS = nS, nM = nM, nP = nP, mh_pprior_sd = mh_pprior_sd, pi_L = pi_L)
probobs_curr <- sim.data$q # Detection probs for this species type
probobs_others <- sim.data$p # Detection probs for the other species type
curr_inter = t(sim.data$L)
focus = aperm(sim.data$Focus, c(2, 1, 3))
occur_prior_probs <- sim.data$prior_probs
p_occur_others <- sim.data$p_occur_others
occur_others <- sim.data$occur_others # actually the same as p above in my sims, but good to keep separate
p_curr <- sim.data$p_curr
occ_curr <- sim.data$occ_curr
beta <- sim.data$theta_curr

# For validation in simulations
Pi_p <- sim.data$Pi_p
O_p <- sim.data$O_p
A <- sim.data$A
A.presence <- t(apply(A,2, function(x) colSums(x) > 0)) # plants are rows, studies are columns

### 2. Fit with each sampler
R <- 2000
burn_in <- 500

# 1. Fit with AMH
t0 <- Sys.time()
AMH <- MH_Adaptive(R = R, burn_in = burn_in, R0 = 10,
                         epsilon = 1e-6, s.extra = 5, mh_pprior_sd = 0.1, 
                         p_curr, occ_curr, beta, 
                         occur_prior_probs, probobs_curr, probobs_others, occur_others, 
                         curr_inter, focus, A)
time.amh <- Sys.time() - t0

# 2. Fit with RWMH
t0 <- Sys.time()
RMH <- MH_RW(R = R, burn_in = burn_in,
                         mh_pprior_sd = 0.1, mh_prop_sd = 2, 
                         new_value_p = p_curr, new_value_occ = occ_curr, new_value_theta = beta, 
                         occur_prior_probs, probobs_curr, probobs_others, occur_others,
                         curr_inter, focus, A)
time.rmh <- Sys.time() - t0

# 3. Fit with sequential sampler
t0 <- Sys.time()
MHS <- MH_sequential(R = R, burn_in = burn_in,
                          mh_occ_step = 0.1, mh_pprior_sd = 0.1, 
                          p_curr = p_curr, occ_curr = occ_curr, 
                          occur_prior_probs, probobs_curr, probobs_others, occur_others,
                          curr_inter, focus, A)
time.mhs <- Sys.time() - t0

# 4. Fit with blocked sampler
t0 <- Sys.time()
BMH <- Blockupdate_OpP(R = R, burn_in = burn_in,
                          mh_p_step = 0.1, mh_pprior_sd = 0.1, 
                          p_1to0 = 0.5, p_0to1 = 0.5, 
                          p_curr, occ_curr,
                          occur_prior_probs, probobs_curr, probobs_others, occur_others,
                          curr_inter, focus, A)
time.bmh <- Sys.time() - t0

time.amh
time.rmh
time.mhs
time.bmh

### 4. Compare performance

### Accuracy
# Accuracy for pi is better under sequential sampling, then RW blocked, then adaptive
png(file = paste0(save_path, "AccuracyP_all_noPM.png"), width = 1000, height = 500)
par(mfcol = c(2,4))

post.pi.mhs <- apply(MHS$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.mhs) , main = "Posterior vs Prior: Sequential", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.mhs), main = "Posterior vs True: Seq", xlab = "True Pi",
     ylab = "Post Mean Pi") 
mse.mhs <- mean((Pi_p - post.pi.mhs)^2)

post.pi.bmh <- apply(BMH$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.bmh) , main = "Posterior vs Prior: Blocked No DA", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.bmh), main = "Posterior vs True: BMH", xlab = "True Pi",
     ylab = "Post Mean Pi") 
mse.bmh <- mean((Pi_p - post.pi.bmh)^2)

post.pi.rmh <- apply(RMH$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.rmh) , main = "Posterior vs Prior: Blocked DA", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.rmh), main = "Posterior vs True: RW", xlab = "True Pi",
     ylab = "Post Mean Pi")
mse.rmh <- mean((Pi_p - post.pi.rmh)^2)

post.pi.amh <- apply(AMH$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.amh) , main = "Posterior vs Prior: AMH", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.amh), main = "Posterior vs True: AMH", xlab = "True Pi",
     ylab = "Post Mean Pi") 
mse.amh <- mean((Pi_p - post.pi.amh)^2)

mtext("Pi Accuracy: Simulations from Prior", side = 3, line = -1.25, outer = TRUE)
dev.off()

mse.amh
mse.rmh
mse.mhs
mse.bmh

# Accuracy for O_p is best under sequential sampling, then adaptive, then RW blocked
post.occ.amh <- apply(AMH$occ_samples, c(2,3), mean)
post.occ.rmh <- apply(RMH$occ_samples, c(2,3), mean)
post.occ.s <- apply(MHS$occ_samples, c(2,3), mean)
post.occ.bmh <- apply(BMH$occ_samples, c(2,3), mean)
postbin.occ.amh <- rbinom(n = nS*nP, size = 1, prob = post.occ.amh)
postbin.occ.rmh <- rbinom(n = nS*nP, size = 1, prob = post.occ.rmh)
postbin.occ.s <- rbinom(n = nS*nP, size = 1, prob = post.occ.s)
postbin.occ.bmh <- rbinom(n = nS*nP, size = 1, prob = post.occ.bmh)
table(c(postbin.occ.amh), c(O_p))
table(c(postbin.occ.rmh), c(O_p))
table(c(postbin.occ.s), c(O_p))
table(c(postbin.occ.bmh), c(O_p))

png(file = paste0(save_path, "AccuracyO_all_noPM.png"), width = 1000, height = 500)
par(mfrow = c(1,4), oma=c(0,0,3,0))

roc.s <- roc(c(O_p), c(post.occ.s))
plot(roc.s, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.s$auc, 3)))
title(main = "Sequential", line = 2.5, cex = 1.5)

roc.bmh <- roc(c(O_p), c(post.occ.bmh))
plot(roc.bmh, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.bmh$auc, 3)))
title(main = "Blocked No DA", line = 2.5, cex = 1.5)

roc.rmh <- roc(c(O_p), c(post.occ.rmh))
plot(roc.rmh, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.rmh$auc, 3)))
title(main = "Blocked DA", line = 2.5, cex = 1.5)
  
roc.amh <- roc(c(O_p), c(post.occ.amh))
plot(roc.amh, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.amh$auc, 3)))
title(main = "AMH", line = 2.5, cex = 1.5)

mtext("Occ Indicator Accuracy: Simulations from Prior", side = 3, line = 1, outer = TRUE)
dev.off()
  
# Acceptance rate
summary(apply(AMH$accepted, 1, mean))
summary(apply(RMH$accepted, 1, mean))
summary(apply(MHS$accepted_p, 1, mean))
summary(apply(BMH$accepted, 1, mean))

# Mixing
png(file = paste0(save_path, "Mixing_all_noPM.png"), width = 1000, height = 500)
par(mfrow = c(3,4))
plot.indices <- arrayInd(sample(which(A.presence==0), 3), dim(A.presence))  
for(i in 1:nrow(plot.indices)){
  plot(MHS$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l",        
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title("Sequential")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue")
  
  plot(BMH$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l",        
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title("Blocked: No DA")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue") 

  plot(RMH$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l",        
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title("Blocked: DA")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue")
  
  plot(AMH$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l", 
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title(main = "AMH")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue")
  
}

mtext("Mixing: Simulations from Prior", side = 3, line = -1.25, outer = TRUE)
dev.off()


#---------------------------  PRIOR MISSPECIFICATION ------------------------#

## We use the same data above, but now we shuffle some of the non-1 prior probs

bad.prop <- 0.25 # proportion of non-1 prior probs we replace
occur_prior_probs_bad <- sim.data$prior_probs
bad.indices <- arrayInd(sample(which(A.presence==0), nS*nP*bad.prop), dim(A.presence))  
occur_prior_probs_bad[bad.indices] <- runif(nS*nP*bad.prop)

### 2. Fit with each sampler

# 1. Fit with AMH
AMH.pm <- MH_Adaptive(R = R, burn_in = burn_in, R0 = 10,
                   epsilon = 1e-6, s.extra = 5, mh_pprior_sd = 0.1, 
                   p_curr, occ_curr, beta, 
                   occur_prior_probs_bad, probobs_curr, probobs_others, occur_others, 
                   curr_inter, focus, A)

# 2. Fit with RWMH
RMH.pm <- MH_RW(R = R, burn_in = burn_in, #s.extra = 5,
             mh_pprior_sd = 0.1, mh_prop_sd = 2,
             new_value_p = p_curr, new_value_occ = occ_curr, new_value_theta = beta, 
             occur_prior_probs_bad, probobs_curr, probobs_others, occur_others,
             curr_inter, focus, A)

# 3. Fit with sequential sampler
MHS.pm <- MH_sequential(R = R, burn_in = burn_in,
                     mh_occ_step = 0.1, mh_pprior_sd = 0.1, 
                     p_curr = p_curr, occ_curr = occ_curr, 
                     occur_prior_probs_bad, probobs_curr, probobs_others, occur_others,
                     curr_inter, focus, A)

# 4. Fit with blocked indpt
BMH.pm <- Blockupdate_OpP(R = R, burn_in = burn_in,
                          mh_p_step = 0.1, mh_pprior_sd = 0.1, 
                          p_1to0 = 0.65, p_0to1 = 0.25, 
                          p_curr, occ_curr,
                          occur_prior_probs_bad, probobs_curr, probobs_others, occur_others,
                          curr_inter, focus, A)

### 4. Compare performance

### Accuracy
# Accuracy for pi is better under sequential sampling, then RW blocked, then adaptive
png(file = paste0(save_path, "AccuracyP_all_PM.png"), width = 1000, height = 500)
par(mfcol = c(2,4))

post.pi.mhs <- apply(MHS.pm$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.mhs) , main = "Posterior vs Prior: Sequential", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.mhs), main = "Posterior vs True: Seq", xlab = "True Pi",
     ylab = "Post Mean Pi") 
mse.mhs <- mean((Pi_p - post.pi.mhs)^2)

post.pi.bmh <- apply(BMH.pm$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.bmh) , main = "Posterior vs Prior: Blocked No DA", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.bmh), main = "Posterior vs True: BMH", xlab = "True Pi",
     ylab = "Post Mean Pi") 
mse.bmh <- mean((Pi_p - post.pi.bmh)^2)

post.pi.rmh <- apply(RMH.pm$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.rmh) , main = "Posterior vs Prior: Blocked DA", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.rmh), main = "Posterior vs True: RW", xlab = "True Pi",
     ylab = "Post Mean Pi")
mse.rmh <- mean((Pi_p - post.pi.rmh)^2)

post.pi.amh <- apply(AMH.pm$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.amh) , main = "Posterior vs Prior: AMH", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
plot(c(Pi_p), c(post.pi.amh), main = "Posterior vs True: AMH", xlab = "True Pi",
     ylab = "Post Mean Pi") 
mse.amh <- mean((Pi_p - post.pi.amh)^2)

mtext("Pi Accuracy: Simulations with Prior Miss-specification", side = 3, line = -1.25, outer = TRUE)
dev.off()


mse.amh
mse.rmh
mse.mhs
mse.bmh

# Accuracy for O_p is best under sequential sampling or blocked 
post.occ.amh <- apply(AMH.pm$occ_samples, c(2,3), mean)
post.occ.rmh <- apply(RMH.pm$occ_samples, c(2,3), mean)
post.occ.s <- apply(MHS.pm$occ_samples, c(2,3), mean)
post.occ.bmh <- apply(BMH.pm$occ_samples, c(2,3), mean)
postbin.occ.amh <- rbinom(n = nS*nP, size = 1, prob = post.occ.amh)
postbin.occ.rmh <- rbinom(n = nS*nP, size = 1, prob = post.occ.rmh)
postbin.occ.s <- rbinom(n = nS*nP, size = 1, prob = post.occ.s)
postbin.occ.bmh <- rbinom(n = nS*nP, size = 1, prob = post.occ.bmh)
table(c(postbin.occ.amh), c(O_p))
table(c(postbin.occ.rmh), c(O_p))
table(c(postbin.occ.s), c(O_p))
table(c(postbin.occ.bmh), c(O_p))

png(file = paste0(save_path, "AccuracyO_all_PM.png"), width = 1000, height = 500)
par(mfrow = c(1,4), oma=c(0,0,3,0))

roc.s <- roc(c(O_p), c(post.occ.s))
plot(roc.s, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.s$auc, 3)))
title(main = "Sequential", line = 2.5, cex = 1.5)

roc.bmh <- roc(c(O_p), c(post.occ.bmh))
plot(roc.bmh, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.bmh$auc, 3)))
title(main = "Blocked No DA", line = 2.5, cex = 1.5)

roc.rmh <- roc(c(O_p), c(post.occ.rmh))
plot(roc.rmh, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.rmh$auc, 3)))
title(main = "Blocked DA", line = 2.5, cex = 1.5)

roc.amh <- roc(c(O_p), c(post.occ.amh))
plot(roc.amh, print.auc = FALSE, auc.polygon = TRUE, auc.polygon.col = "skyblue", 
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) 
text(0.4, 0.4, paste("AUC:", round(roc.amh$auc, 3)))
title(main = "AMH", line = 2.5, cex = 1.5)

mtext("Occ Indicator Accuracy: Simulations with Prior Miss-Specification", side = 3, line = 1, outer = TRUE)
dev.off()


# Acceptance rate
summary(apply(AMH.pm$accepted, 1, mean))
summary(apply(RMH.pm$accepted, 1, mean))
summary(apply(MHS.pm$accepted_p, 1, mean))
summary(apply(BMH.pm$accepted, 1, mean))


# Mixing
png(file = paste0(save_path, "Mixing_all_PM.png"), width = 1000, height = 500)
par(mfrow = c(3,4))
plot.indices <- arrayInd(sample(which(A.presence==0), 3), dim(A.presence))  
for(i in 1:nrow(plot.indices)){
  plot(MHS.pm$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l",        
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title("Sequential")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue")
  
  plot(BMH.pm$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l",        
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title("Blocked: No DA")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue") 
  
  plot(RMH.pm$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l",        
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title("Blocked: DA")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue")
  
  plot(AMH.pm$pi_samples[,plot.indices[i, 1], plot.indices[i,2]], type = "l", 
       ylab = paste0("Pi[,", plot.indices[i,1], ",", plot.indices[i,2], "]"))
  title(main = "AMH")
  abline(h = Pi_p[plot.indices[i,1], plot.indices[i,2]], col = "blue")
  
}

mtext("Mixing: Simulations from Prior", side = 3, line = -1.25, outer = TRUE)
dev.off()



