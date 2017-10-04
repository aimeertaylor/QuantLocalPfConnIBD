################################################################################################
# Script to estimate Fst using WGS data and estimator outlined in supplementary of Reich et al. 2009
# See subheading for timings assuming 1000 repeats
# Only calculates for VeryLate (since insufficient data for earlier years)
# Refer to permute as perturb accidentally throughout
################################################################################################

# Load packages, data and function scripts
rm(list = ls())
require('foreach') # For running in parallel 
require('doMC') # For running in parallel 
require('rngtools') # For running in parallel 
require('abind') # For abind
require('plyr') # For aplyr

# Load sequencing data 
source('./calculate_pairwise_Fst.R')
load('./Data_store_WGS.RData')
load('./geo_dist_info.RData')

attach(Data_store)
attach(geo_dist_info, warn.conflicts = FALSE)
numComparisons <- length(geo_order)
Pair_wise_site_comparisons <- pairwise_site_distance[,c(1,2)]

# Magic numbers
registerDoMC(cores = 8) 
sample_size_threshold <- 3
n_repeat <- 1000

# ==================================================================================
# Create stores for Fst estimates
# ==================================================================================
stages <- as.character(unique(MetaData$stage))

Pair_wise_site_comparisons_Fst <- matrix(ncol = nrow(Pair_wise_site_comparisons), nrow = 1,
                                         dimnames = list(NULL, geo_order))
Pair_wise_site_comparisons_Fst_perSNP <- matrix(ncol = nrow(Pair_wise_site_comparisons),
                                                nrow = numSNPs,
                                                dimnames = list(NULL, geo_order))
Pair_wise_site_comparisons_Fst_perstage <- matrix(ncol = nrow(Pair_wise_site_comparisons), 
                                                  nrow = length(stages), 
                                                  dimnames = list(stages, geo_order))
Pair_wise_site_comparisons_Fst_perSNP_perstage <- array(dim = c(numSNPs,nrow(Pair_wise_site_comparisons), length(stages)),
                                                        dimnames = list(NULL, geo_order, stages))


# ==================================================================================
# Calculate pairwise FST estimates (~ 5 seconds)
# ==================================================================================
system.time(for(i in 1:nrow(Pair_wise_site_comparisons)){
  
  # Fst estimates averaging over years: indices, calculate, save 
  P1 <- (MetaData$collection_location == as.character(Pair_wise_site_comparisons[i,1]))
  P2 <- (MetaData$collection_location == as.character(Pair_wise_site_comparisons[i,2]))
  X <- calculate_pairwise_reich(P1, P2, SNPDataBinary = SNPData)
  Pair_wise_site_comparisons_Fst[i] <- X$F_st
  Pair_wise_site_comparisons_Fst_perSNP[,i] <- X$F_st_snp
  
  # Fst estimates per stage
  for(stage in 'VeryLate'){
    P1 <- (MetaData$collection_location == as.character(Pair_wise_site_comparisons[i,1]) & as.character(MetaData$stage) == stage)
    P2 <- (MetaData$collection_location == as.character(Pair_wise_site_comparisons[i,2]) & as.character(MetaData$stage) == stage)
    if(sum(P1) > sample_size_threshold & sum(P2) > sample_size_threshold){
      X <- calculate_pairwise_reich(P1, P2, SNPDataBinary = SNPData)  
      Pair_wise_site_comparisons_Fst_perstage[stage, i] <- X$F_st
      Pair_wise_site_comparisons_Fst_perSNP_perstage[,i,stage] <- X$F_st_snp
    }
  }
})

# Store altogether in a list 
Fst_WGS <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst, 
                Pair_wise_site_comparisons_Fst_perSNP = Pair_wise_site_comparisons_Fst_perSNP, 
                Pair_wise_site_comparisons_Fst_perstage = Pair_wise_site_comparisons_Fst_perstage, 
                Pair_wise_site_comparisons_Fst_perSNP_perstage = Pair_wise_site_comparisons_Fst_perSNP_perstage)


# ==================================================================================
# Calculate bootstrap over loci Fst estimates  (~ 4151 sec)
# ==================================================================================
system.time(Fst_WGS_bootstrap <- foreach(n = 1:n_repeat) %dopar%{
  
  # Resample 
  bootstrap_snps <- sample(numSNPs, replace = TRUE)
  Y <- SNPData[,bootstrap_snps]
  Z <- MetaData
  
  # Memory store
  store <- matrix(ncol = nrow(Pair_wise_site_comparisons), nrow = 1, dimnames = list(NULL, geo_order))
  store_perstage <- matrix(ncol = nrow(Pair_wise_site_comparisons), nrow = length(stages), dimnames = list(stages, geo_order))
  
  # Calculate pairwise thetas 
  for(i in 1:nrow(Pair_wise_site_comparisons)){
    P1 <- (Z$collection_location == as.character(Pair_wise_site_comparisons[i,1]))
    P2 <- (Z$collection_location == as.character(Pair_wise_site_comparisons[i,2]))
    X <- calculate_pairwise_reich(P1, P2, SNPDataBinary = Y)
    store[i] <- X$F_st
    
    for(stage in 'VeryLate'){
      P1 <- (Z$collection_location == as.character(Pair_wise_site_comparisons[i,1]) & as.character(Z$stage) == stage)
      P2 <- (Z$collection_location == as.character(Pair_wise_site_comparisons[i,2]) & as.character(Z$stage) == stage)
      if(sum(P1) > sample_size_threshold & sum(P2) > sample_size_threshold){
        X <- calculate_pairwise_reich(P1, P2, SNPDataBinary = Y)
        store_perstage[stage, i] <- X$F_st
      }
    }
  }
  stores <- list(store = store, store_perstage = store_perstage)
  return(stores)
})

# Pull out boostrap CIs result in an array, rather than in a list and save
A <- do.call(rbind, sapply(Fst_WGS_bootstrap, FUN = function(x){x['store']})) # Pull out values
A_deltas <- apply(A, 1, FUN = function(x){x[geo_order] - Fst_WGS$Pair_wise_site_comparisons_Fst[1,geo_order]}) # Calculate differences
A_percentiles <- apply(A_deltas, 1, quantile, probs = c(0.025, 0.975)) # Pull out quantiles
Fst_WGS$Pair_wise_site_comparisons_Fst_CIs <- apply(A_percentiles, 1, FUN = function(x){Fst_WGS$Pair_wise_site_comparisons_Fst[1,geo_order] - x}) # Calculate CI values

B <- sapply(Fst_WGS_bootstrap, FUN = function(x){x['store_perstage']}) # Pull out values
B_deltas <- lapply(B, FUN = function(x){x[stages,geo_order] - Fst_WGS$Pair_wise_site_comparisons_Fst_perstage[stages,geo_order]}) # Calculate differences
B_deltas2 <- aperm(do.call(abind, args = list(B_deltas, along = 3)), c(3,1,2)) # Reorder dimensions 
B_percentiles <- alply(apply(B_deltas2, c(2,3), quantile, probs = c(0.025, 0.975), na.rm = TRUE), 
                       1, .dims = TRUE) 
Fst_WGS$Pair_wise_site_comparisons_Fst_perstage_CIs <- lapply(B_percentiles, FUN = function(x){Fst_WGS$Pair_wise_site_comparisons_Fst_perstage[stages,geo_order] - x[stages,geo_order]})


# ==================================================================================
# Calculate perturbed Fst estimates (~ 6408 sec) # 47084
# ==================================================================================
system.time(Fst_WGS_peturbed <- foreach(n = 1:n_repeat) %dopar%{
  
  # Setup
  Location_true <- MetaData$collection_location # make copy of unperturbed values
  perturb <- function(x){x[sample(length(x))]} # function to perturb
  Pair_wise_site_comparisons_Fst <- matrix(ncol = nrow(Pair_wise_site_comparisons), # results store
                                           nrow = 1,
                                           dimnames = list(NULL, geo_order))
  
  for(i in 1:nrow(Pair_wise_site_comparisons)){

    # Perturb location labels 
    # We only want to perturb the labels of clinics A and B, not clincs A, B, C and D
    Location_temp <- Location_true
    ClinicA <- as.character(Pair_wise_site_comparisons[i,1])   
    ClinicB <- as.character(Pair_wise_site_comparisons[i,2])
    ClinicABInd <- Location_true == ClinicA | Location_true == ClinicB 
    Location_temp[ClinicABInd] <- perturb(Location_true[ClinicABInd]) 
    
    # Indices for populations to compare, calculate, save
    P1 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,1]))
    P2 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,2]))
  
    X <- calculate_pairwise_reich(P1, P2, SNPDataBinary = SNPData)
    Pair_wise_site_comparisons_Fst[i] <- X$F_st
    rm(Location_temp)
  }
  rm(Location_true)
  
  for(stage in 'VeryLate'){
    
    ind <- as.character(MetaData$stage) == stage
    SNPDataVeryLate <- SNPData[ind,]
    Location_true <- MetaData$collection_location[ind] # make new copy of unperturbed values
    
    for(i in 1:nrow(Pair_wise_site_comparisons)){
      # Perturb location labels of A and B only 
      Location_temp <- Location_true
      ClinicA <- as.character(Pair_wise_site_comparisons[i,1])   
      ClinicB <- as.character(Pair_wise_site_comparisons[i,2])
      ClinicABInd <- Location_true == ClinicA | Location_true == ClinicB 
      Location_temp[ClinicABInd] <- perturb(Location_true[ClinicABInd]) 
      
      # Calculate fst and then remove location_temp
      P1 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,1]))
      P2 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,2]))
      X <- calculate_pairwise_reich(P1, P2, SNPDataBinary = SNPDataVeryLate)  
      Pair_wise_site_comparisons_Fst_perstage[stage, i] <- X$F_st
      rm(Location_temp)
    }
  }  
  
  Fst_perturbed <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst, 
                        Pair_wise_site_comparisons_Fst_perstage = Pair_wise_site_comparisons_Fst_perstage)
  return(Fst_perturbed)
})

# Re-structure result in an array, rather than in a list and save
X <- do.call(rbind, sapply(Fst_WGS_peturbed, FUN = function(x){x['Pair_wise_site_comparisons_Fst']}))
Y <- do.call(abind, args = list(sapply(Fst_WGS_peturbed, FUN = function(x){x['Pair_wise_site_comparisons_Fst_perstage']}), along = 3))
Fst_WGS$Pair_wise_site_comparisons_Fst_perturbed = X
Fst_WGS$Pair_wise_site_comparisons_Fst_perstage_perturbed = Y


# ==================================================================================
# Save all 
#save(Fst_WGS, file = './Fst_WGS.RData')
# ==================================================================================

