# Function for block updating of occurrence indicator and occurrence probability
# Random walk version
# This file simulates and tests the function


#' @param detected Matrix. Rows are species, columns are studies. Indicator of
#' whether the species were detected to be present in the given study. Derived
#' from the provided occurrence matrix, and checked against obs_A. Note that
#' while most studies only report interacting species, some studies do report
#' present but not interacting species. This often occurs when subsetting. 
#' @param p_occur Matrix. Rows are species, columns are studies. It includes the
#' prior probabilities of occurrence.
#' @param p_occur_others. Matrix Rows are species, columns are studies. It includes
#' the current values of occurrence for the other set of species, the ones that
#' are not being updated. 
#' @param probobs_curr Vector of current values for the probability of
#' observation for the set of species we are updating.
#' @param probobs_others Vector of current values for the probability of
#' observation for the other set of species, the ones that are not being
#' updated.
#' @param curr_inter Matrix with 0,1 entries for the current posterior samples
#' of possible interactions. The rows correspond to the species we are updating
#' and the columns correspond to the other set of species.
#' @param focus Array of three dimensions corresponding to the two sets of
#' species, and the different studies. Values are 0 or 1. This array represents
#' whether the current study would be willing to record the interactions of a
#' given species. The value should be 1 except for studies that are animal or
#' plant-oriented, where the value for species that are not of focus will be 0.


library(truncnorm)
library(mvtnorm)

### TO DO: 

## Run RW version for many iterations to test parameter recovery
## Code Adaptive version - do we need MH, I forget
## Test it for recovery as above
## Incorporate adaptive version into the MCMC function, maybe with an option?

#------------------------ 1. SIMULATE FROM THE BASIC MODEL -----------------------#
nS <- 70
nM <- 55
nP <- 500

set.seed(1234)

# True interaction matrix
L <- matrix(rbinom(nM*nP, 1, 0.5), nM, nP)

# Detection probabilities
p <- runif(nM, 0, 1)
q <- runif(nP, 0, 1)
pq <- outer(p, q)

# Study focus: all studies are animal-focused
Focus <- array(0, dim = c(nM, nP, nS))
p_occur_others <- matrix(0, nM, nS) # probability of occurrence for other species type (mammals, so in this case it will be 0/1)
for(s in 1:nS){
  # First focus animal
  these.mammals <- sample(1:nM, 2)
  Focus[these.mammals[1],,s] <- 1
  p_occur_others[these.mammals[1], s] <- 1
  # Second focus animal (sometimes)
  u <- runif(1)
  if(u > 0.85){
    Focus[these.mammals[2],,s] <- 1
    p_occur_others[these.mammals[2],s] <- 1
  }
}

# Plant occurrences and occurrence probabilities
probs = seq(0,1, 0.1) # just some random probs to sample from
prior_probs <- matrix(sample(probs, size = nP*nS, replace = TRUE), nP, nS) # Prior occurrence probabilities to use as input and to simulate true probs
mh_pprior_sd <- 0.1
Pi_p <- matrix(rtruncnorm(n = nP*nS,  mean = prior_probs, sd = mh_pprior_sd, a = 0, b = 1), nP, nS) # True occurrence probabilities for this species type
O_p <- matrix(rbinom(nP*nS, 1, Pi_p), nP, nS) # True occurrence indicators

# Observed interactions
A <- array(data = NA, dim = c(nM, nP, nS))
for(s in 1:nS){
  this.prob <- L * Focus[,,s] * pq *  matrix(1, nM, 1) %*% O_p[,s]
  A[,,s] <- matrix(rbinom(nM*nP, 1, this.prob), nM, nP)
}


# Randomly initialize current states to be stored and provided as input
theta_curr <- array(rnorm(nP*nS*2,0,1), dim = c(nP,nS,2)) # Store in mcmc function: must retrieve from previous iters, Tricky because of probit link
p_occur_curr <- pnorm(theta_curr[,,1]) # the current state of occurence prob to start from, these should converge to the true values in Pi_p
occur_curr <- ifelse(theta_curr[,,2] >0 , 1, 0) # the current state of occurrence ind to start from, these should converge to the true values in O_p

new_value_theta <- theta_curr # maybe we should just use these vars instead, if so, delete above but copy comments
new_value_p <-p_occur_curr
new_value_occ <- occur_curr

############# now add these into the sampler instead of _curr or assign _curr_st based on?

#----------------------------- Step 2: Sampler --------------------------------#

# Idea: we want to use Adaptive Metropolis to sample the occurrence indicators and 
# probabilities in one block. Use probit data augmentation to do it, so we'll have
# pi = phi(y), y ~ N(0,e)
# O_p = 1(z>0), z ~ N(0,e)
# start with e = 1


#################### Things the function actually needs: these are inputs not computed quantities

#p_occur <- Pi_p # Current occurrence probs for this species type; nope - Pi_p is true values

probobs_curr <- q # Detection probs for this species type
probobs_others <- p # Detection probs for the other species type
#occur_ind <- O_p # nope, these are the true values
curr_inter = t(L)
focus = aperm(Focus, c(2, 1, 3))
occur_prior_probs <- prior_probs


## Compute these
num_obs <- nrow(p_occur_curr) # number of species of this type
num_studies <- ncol(p_occur_curr)

# Setup storage etc
d <- 2
cov_prop <- diag(1,d)
new_O <- matrix(NA, nrow = num_obs, ncol = num_studies)

# Return something like this to keep track of acceptance?
#accepted = matrix(0, nrow = num_obs, ncol = num_studies)
n_iter <- 100
accepted = array(0, dim = c(n_iter, num_obs, num_studies))

# If the species were observed to interact in the study, the posterior for
# occurrence is always equal to 1. First, I will draw the occurence indicator
# as if the species were not observed to interact. Then, I will assign the
# ones observed interacting to 1.

# Updating the interaction indicator as if the species was not detected in
# the study.

# Questions and to do
  # Rename occur to be more informative, it's actually a probability
  # separate cov for each species?
  # we want prob, occ to be blocked for each species, but not for all species
  # so cov updates don't involve different species
  # ?? Is it faster to use a big mvtnorm with kronecker product covariance or to do a loop?
  # storing theta? we don't have to store all samples, just update with recursion

### Proposal: start with RW Metropolis

### Acceptance probability
## Compute prob1 = p(occ = 1 | pi, ....)
# For the probability that indicator is equal to 1:

pipj <- outer(probobs_curr, probobs_others)

# Not sure if we need these, store in new_value_theta instead?
# theta_star <- array(NA, dim = c(num_obs, num_studies, d)) # store proposal thetas here, INPUTS not computed quantities
# pi_occ_star <- array(NA, dim = c(num_obs, num_studies, d)) # store proposal transformed thetas here

# Debugging only: in the real MCMC we'll just do one iteration per step 
# X check mixing
# X check that we're updating theta etc
# X seems odd that occ ind isn't in logL - oops, need more complex logL for 0/1 cases
# ??? Pattern: pi shoots to 0/1 then we reject indefinitely
# This only started happening after I "fixed" the logL
# Maybe this happens when we have studies with 0 observed interactions?



# Set up some storage used only inside of the loop

#pi_samples <- matrix(data = NA, nrow = length())
pi_samples <- array(data = NA, dim = c(n_iter, num_obs, num_studies))
occ_samples <- array(data = NA, dim = c(n_iter, num_obs, num_studies))
theta_prop_st <- matrix(data = NA, nrow = num_obs, ncol = 2)

pb = txtProgressBar(min = 0, max = n_iter, initial = 0) 

for(r in 1:n_iter){
  setTxtProgressBar(pb,r)
  start <- Sys.time()
  
for (st in 1 : num_studies) {
  
  #p_curr <- p_occur_curr[, st]
  p_prior_st <- occur_prior_probs[, st]
  theta_curr_st  <- new_value_theta[,st,]
  pi_curr_st  <- new_value_p[,st]
  occur_curr_st <- new_value_occ[,st]
  #theta_star_st <- theta_star[,st,]
  
  ### Generate proposals independently for each species and then transform all 
  for(i in 1:num_obs){
    cov_prop <- 0.1*diag(1,d) # Since this version isn't adaptive, this is the same for all species, but that will change
    theta_prop_st[i,] <- rmvnorm(1, theta_curr_st[i,], cov_prop) 
  }
  

  #Sys.time() - start
  pi_prop_st <- pnorm(theta_prop_st[,1])
  occur_prop_st <- ifelse(theta_prop_st[,2] > 0 , 1, 0)
  
  ### Compute acceptance probabilities
  
  ## Proposal 
  # p(occ | pi, ...)
  # Double check: I should be using new value of probs in here?
  p_occ_others_st <- outer(rep(1,num_obs), p_occur_others[, st] )
  prod1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* p_occ_others_st), 1, prod)   # Take product over other species
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
  
  # DEBUGGING when p_curr = 1 but occur_cur_st = 0, we get log(0) = -inf for that component
  
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
  accepted[r,update, st] <- 1 # change the dimension for single iteration
}
  
  pi_samples[r,,] <- new_value_p
  occ_samples[r,,] <- new_value_occ
  #cat("\n Iteration ", r, " complete: that took ", Sys.time() - start)

}

plot(c(new_value_p), c(Pi_p))
table(c(new_value_occ), c(O_p))
apply(accepted, 1, mean)

plot(pi_samples[,1,1], type = "l")
plot(pi_samples[,2,1], type = "l")

# Look at plant 6, study 3, jumps around iteration 27 - 28
# Debugging: we're sampling 0 for occurrence each time, but pi_occ goes to 1
# This plant has L = 1 for most things, but only one interaction partner present
# But A = 0 for everything, in fact no interacitons observed in this study
# When P is high, this results in a likelihood approx 0, and a lower value of prior
# So we should be rejecting higher P. Why is P climbing?
# So why is logL prop > logL current?
# Both prior prob and l seem off? Prior should be greater for curr because curr is closer to prior than prop
pi_curr_st <- pi_samples[3,,3]
pi_prop_st <- pi_samples[4,,3] 
p_prior_st <-  occur_prior_probs[, 3]
p_prior_st[6]
pi_curr_st[6]
pi_prop_st[6] # Proposal is higher and further from prior than current
log(dtruncnorm(x = pi_curr_st[6], a = 0, b = 1, mean = p_prior_st[6], sd = mh_pprior_sd)) # much higher as expected
log(dtruncnorm(x = pi_prop_st[6], a = 0, b = 1, mean = p_prior_st[6], sd = mh_pprior_sd))

occur_curr_st <- occ_samples[3,,3]
occur_prop_st <- occ_samples[4,,3] # this is getting reassigned somewhere in my compute logL code
occur_curr_st[6]
occur_prop_st[6]

log_prob1_curr[6] # Significantly higher as expected
log_prob1_star[6]

############## left off here: occur is constantly 0 in my samples but not in occur_curr_st? am i subsetting wrong? saving samples wrong?
# this could be why logL is wrong and higher probs for occ = 0 are shooting up?

plot(pi_samples[,6,3], type = "l")
plot(occ_samples[,6,3], type = "l")
head(pi_samples[,6,3])

# Can skip this for starters
# Now, if a species was detected to interact in a study, set their occurrence
# indicator to 1.
wh_detected <- which(detected == 1)
new_O[wh_detected] <- 1












