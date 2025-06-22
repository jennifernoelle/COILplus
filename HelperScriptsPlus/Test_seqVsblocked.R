# This script tests Blocked and sequential samplers on simulated data
# Doesn't include adaptive MH because too slow
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
# Experiment with 0/1 swich prob in blocked version
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

library(pROC)

source_path <- 'HelperScriptsPlus/' 
save_path <- 'HelperScriptsPlus/'
source(paste0(source_path, "Utils.R"))
source(paste0(source_path, "Sequential_MH_PiOp.R"))
source(paste0(source_path, "Blocked_indpt_PiOp.R"))


#--------------------------- NO PRIOR MISSPECIFICATION ------------------------#

#  We simulate from the prior, but note that in the simulation function, 
#  Prior probs are reset to 1 for detected species for simplification
#  This is as it would be in the input data

### 1. Simulate data without prior misspecification

# Set up key external function inputs: these come from inputs and current state
set.seed(1234)
nS <- 200
nM <- 30
nP <- 100
pi_L <- 0.05 # Sets sparsity of the true interactions matrix
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

# A. Fit with sequential sampler
t0 <- Sys.time()
MHS <- MH_sequential(R = R, burn_in = burn_in,
                     mh_occ_step = 0.1, mh_pprior_sd = 0.1, 
                     p_curr = p_curr, occ_curr = occ_curr, 
                     occur_prior_probs, probobs_curr, probobs_others, occur_others,
                     curr_inter, focus, A)
time.mhs <- Sys.time() - t0

# B. Fit with blocked sampler
t0 <- Sys.time()
BMH <- Blockupdate_OpP(R = R, burn_in = burn_in,
                       mh_p_step = 0.1, mh_pprior_sd = 0.1, 
                       p_1to0 = 0.25, p_0to1 = 0.25, 
                       p_curr, occ_curr,
                       occur_prior_probs, probobs_curr, probobs_others, occur_others,
                       curr_inter, focus, A)
time.bmh <- Sys.time() - t0

time.mhs
time.bmh

### 4. Compare performance


### Accuracy
png(file = paste0(save_path, "AccuracyP_SeqVBlock_noPM.png"), width = 1000, height = 500)
par(mfcol = c(2,2))

post.pi.mhs <- apply(MHS$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.mhs) , main = "Posterior vs Prior: Seq", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
abline(a = 0, b = 1)
plot(c(Pi_p), c(post.pi.mhs), main = "Posterior vs True: Seq", xlab = "True Pi",
     ylab = "Post Mean Pi") 
abline(a = 0, b = 1)
mse.mhs <- mean((Pi_p - post.pi.mhs)^2)


post.pi.bmh <- apply(BMH$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.bmh) , main = "Posterior vs Prior: BMH", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
abline(a = 0, b = 1)
plot(c(Pi_p), c(post.pi.bmh), main = "Posterior vs True: Blocked No DA", xlab = "True Pi",
     ylab = "Post Mean Pi") 
abline(a = 0, b = 1)
mse.bmh <- mean((Pi_p - post.pi.bmh)^2)

mtext("Pi Accuracy: Simulations from Prior", side = 3, line = -1.25, outer = TRUE)
dev.off()

mse.mhs
mse.bmh

# Accuracy is similar

post.occ.s <- apply(MHS$occ_samples, c(2,3), mean)
post.occ.bmh <- apply(BMH$occ_samples, c(2,3), mean)

postbin.occ.s <- rbinom(n = nS*nP, size = 1, prob = post.occ.s)
postbin.occ.bmh <- rbinom(n = nS*nP, size = 1, prob = post.occ.bmh)

table(c(postbin.occ.s), c(O_p))
table(c(postbin.occ.bmh), c(O_p))


png(file = paste0(save_path, "AccuracyO_SeqV_noPM.png"), width = 1000, height = 500)
par(mfrow = c(1,2), oma=c(0,0,3,0))

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

mtext("Occ Indicator Accuracy: Simulations from Prior", side = 3, line = 1, outer = TRUE)
dev.off()

# Acceptance rate
summary(apply(MHS$accepted_p, 1, mean))
summary(apply(BMH$accepted, 1, mean))


# Mixing
png(file = paste0(save_path, "Mixing_SeqVBlock_noPM.png"), width = 1000, height = 500)
par(mfrow = c(3,2))
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

}
mtext("Mixing: Simulations from Prior", side = 3, line = -1.25, outer = TRUE)
dev.off()





#---------------------------  PRIOR MISSPECIFICATION ------------------------#

## We use the same data above, but now we shuffle some of the non-1 prior probs

bad.prop <- 0.25 # proportion of non-1 prior probs we replace
occur_prior_probs_bad <- sim.data$prior_probs
bad.indices <- arrayInd(sample(which(A.presence==0), nS*nP*bad.prop), dim(A.presence))  
occur_prior_probs_bad[bad.indices] <- runif(nS*nP*bad.prop)

### 1. Fit with each sampler

# A. Fit with sequential sampler
MHS.pm <- MH_sequential(R = 1000, burn_in = 500,
                        mh_occ_step = 0.1, mh_pprior_sd = 0.1, 
                        p_curr = p_curr, occ_curr = occ_curr, 
                        occur_prior_probs_bad, probobs_curr, probobs_others, occur_others,
                        curr_inter, focus, A)

# B. Fit with sequential sampler
BMH.pm <- Blockupdate_OpP(R = 1000, burn_in = 500,
                          mh_p_step = 0.1, mh_pprior_sd = 0.1, 
                          p_1to0 = 0.5, p_0to1 = 0.5, 
                          p_curr, occ_curr,
                          occur_prior_probs_bad, probobs_curr, probobs_others, occur_others,
                          curr_inter, focus, A)

### 2. Compare performance

### Accuracy
png(file = paste0(save_path, "AccuracyP_SeqVBlock_PM.png"), width = 1000, height = 500)
par(mfcol = c(2,2))

post.pi.mhs <- apply(MHS.pm$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.mhs) , main = "Posterior vs Prior: Seq", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
abline(a = 0, b = 1)
plot(c(Pi_p), c(post.pi.mhs), main = "Posterior vs True: Seq", xlab = "True Pi",
     ylab = "Post Mean Pi") 
abline(a = 0, b = 1)
mse.mhs <- mean((Pi_p - post.pi.mhs)^2)


post.pi.bmh <- apply(BMH.pm$pi_samples, c(2,3), mean) 
plot(c(occur_prior_probs), c(post.pi.bmh) , main = "Posterior vs Prior: BMH", 
     xlab = "Prior Pi", ylab = "Post Mean Pi")
abline(a = 0, b = 1)
plot(c(Pi_p), c(post.pi.bmh), main = "Posterior vs True: Blocked No DA", xlab = "True Pi",
     ylab = "Post Mean Pi") 
abline(a = 0, b = 1)
mse.bmh <- mean((Pi_p - post.pi.bmh)^2)

mtext("Pi Accuracy: Prior Miss-specification", side = 3, line = -1.25, outer = TRUE)
dev.off()

mse.mhs
mse.bmh

# Accuracy for O_p is best under sequential sampling or blocked 
post.occ.s <- apply(MHS.pm$occ_samples, c(2,3), mean)
post.occ.bmh <- apply(BMH.pm$occ_samples, c(2,3), mean)

postbin.occ.s <- rbinom(n = nS*nP, size = 1, prob = post.occ.s)
postbin.occ.bmh <- rbinom(n = nS*nP, size = 1, prob = post.occ.bmh)

table(c(postbin.occ.s), c(O_p))
table(c(postbin.occ.bmh), c(O_p))



png(file = paste0(save_path, "AccuracyO_SeqVBlock_PM.png"), width = 1000, height = 500)
par(mfrow = c(1,2), oma=c(0,0,3,0))

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

mtext("Occ Indicator Accuracy: Simulations from Prior", side = 3, line = 1, outer = TRUE)
dev.off()

# Acceptance rate
summary(apply(MHS.pm$accepted_p, 1, mean))
summary(apply(BMH.pm$accepted, 1, mean))

# Mixing
png(file = paste0(save_path, "Mixing_SeqVBlock_PM.png"), width = 1000, height = 500)
par(mfrow = c(3,2))
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
  
}
mtext("Mixing: Prior Miss-Specification", side = 3, line = -1.25, outer = TRUE)
dev.off()
