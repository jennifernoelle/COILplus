#' Function that updates the indicator of occurence. Corrects errors in original version
#' 
#' Each species might occur in the study area or not. Even though probabilities
#' of occurence can be provided, the detected interactions also inform us of
#' whether the species occur.
#' 
#' @param detected Matrix. Rows are species, columns are studies. Indicator of
#' whether the species were detected to be present in the given study. Derived
#' from the provided occurrence matrix, and checked against obs_A. Note that
#' while most studies only report interacting species, some studies do report
#' present but not interacting species. This often occurs when subsetting. 
#' @param occur Matrix. Rows are species, columns are studies. It includes the
#' prior probabilities of occurence.
#' @param occur_others. Matrix Rows are species, columns are studies. It includes
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
#' 
UpdOccur <- function(detected, occur, occur_others, probobs_curr, probobs_others, curr_inter,
                     focus) {
  
  num_obs <- nrow(occur)
  num_studies <- ncol(occur)
  
  new_O <- matrix(NA, nrow = num_obs, ncol = num_studies)
  
  # If the species were observed to interact in the study, the posterior for
  # occurence is always equal to 1. First, I will draw the occurence indicator
  # as if the species were not observed to interact. Then, I will assign the
  # ones observed interacting to 1.
  
  # Updating the interaction indicator as if the species was not detected in
  # the study.
  
  # For the probability that indicator is equal to 1:
  pipj <- outer(probobs_curr, probobs_others)
  
  for (st in 1 : num_studies) {
    occ_others_st <- outer(rep(1,num_obs), occur_others[, st] )
    # Take product over other species
    prob1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* occ_others_st), 1, prod) ###### ERROR st was previously ss, should also include occur indicators for other species
    #prob1old <- apply((1 - pipj) ^ (curr_inter * focus[, , st]), 1, prod)
    #cat("\nstudy = ", st, sum(prob1old == prob1)) # if animal focused, this will be unchanged for plants and we don't care about occurrences for animals b/c focus
    prob1 <- occur[, st] * prob1 # Multiply by occurrence probabilities to get prob occurrence == 1
    prob0 <- 1 - occur[, st] # Prob occurrence == 0
    prob1 <- prob1 / (prob1 + prob0) # Standardize
    new_O[, st] <- rbinom(num_obs, 1, prob1) # Sample from Bernoullis for each species in this study
  }
  
  # Now, if a species was detected to interact in a study, set their occurence
  # indicator to 1.
  wh_detected <- which(detected == 1)
  new_O[wh_detected] <- 1
  
  return(new_O)
}
