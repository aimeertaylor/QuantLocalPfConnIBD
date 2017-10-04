##############################################################################################
# Script to generate fst estimates
# Perturbed should read permute 
# Commented permuted and bootstap W&C and inferred since no-longer using
# Also commented Years aside from "08" "09" "10" since insufficient data
##############################################################################################

rm(list = ls())
require('foreach') # For running in parallel 
require('doMC') # For running in parallel 
require('rngtools') # For running in parallel 
require('abind') 
require('plyr') # For aplyr

source('./calculate_pairwise_Fst.R')
load('./Data_store_barcode.RData')
load('./geo_dist_info.RData')

attach(geo_dist_info, warn.conflicts = FALSE)
Pair_wise_site_comparisons <- pairwise_site_distance[,c(1,2)]
numComparisons <- nrow(Pair_wise_site_comparisons)

# Magic no.s/variables
registerDoMC(cores = 8) # For foreach
n_repeat <- 1000
sample_size_threshold <- 3 # For per year
sample_size_threshold_freqs <- 10 # For per year
estimators <- c('WandC','Reich','Inferred') # Estimators (no-longer using )

# ============================================================================
# Function to calculate pairwise Fst estimates 
# =============================================================================
fst_calculations <- function(SNPDataBinary, Data){
  
  numSNPs <- ncol(SNPDataBinary)
  numSamples <- nrow(SNPDataBinary)
  
  #----------------------------------------------------------------------------
  # Between sites averaging over years 
  #----------------------------------------------------------------------------
  Pair_wise_site_comparisons_Fst <- array(dim = c(numComparisons, length(estimators)),
                                          dimnames = list(geo_order,estimators))
  Pair_wise_site_comparisons_Fst_perSNP <- array(dim = c(numSNPs, numComparisons,
                                                         length(estimators)),
                                                 dimnames = list(NULL, geo_order,estimators))
  Sample_size_min <- rep(NA, length = (numComparisons))
  
  for(i in 1:numComparisons){
    
    # Indices for populations to compare
    P1 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,1]))
    P2 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,2]))
    P1_name <- as.character(Pair_wise_site_comparisons[i,1])
    P2_name <- as.character(Pair_wise_site_comparisons[i,2])
    
    # Store sample size
    Sample_size_min[i] <- min(sum(P1), sum(P2))
    
    # Calculate Fst
    X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBinary)
    # X_WandC <- calculate_pairwise_WandC(P1, P2, SNPDataBinary)
    # X_inferred <- calculate_pairwise_inferred(P1_name, P2_name, 1:numSNPs) 
    
    # Save to matrix
    Pair_wise_site_comparisons_Fst[i, 'Reich'] <- X_Reich$F_st
    Pair_wise_site_comparisons_Fst_perSNP[,i, 'Reich'] <- X_Reich$F_st_snp
    # Pair_wise_site_comparisons_Fst[i, 'WandC'] <- X_WandC$theta_loci_hat_combined
    # Pair_wise_site_comparisons_Fst_perSNP[,i, 'WandC'] <- X_WandC$theta_hat
    # Pair_wise_site_comparisons_Fst[i, 'Inferred'] <- X_inferred$Fst   
    # Pair_wise_site_comparisons_Fst_perSNP[,i,'Inferred'] <- X_inferred$Fst_snp
  }
  
  #----------------------------------------------------------------------------
  # Between sites for each year 
  #----------------------------------------------------------------------------
  Years <- c("08", "09", "10") #as.character(unique(Data[,1]))
  
  # Allocate memory
  Pair_wise_site_comparisons_Fst_peryr <- array(dim = c(length(Years),
                                                        numComparisons, length(estimators)),
                                                dimnames = list(Years, geo_order, estimators))
  Pair_wise_site_comparisons_Fst_perSNP_peryr <- array(dim=c(numSNPs,length(Years), 
                                                             numComparisons, length(estimators)),
                                                       dimnames = list(NULL, Years, 
                                                                       geo_order, estimators))
  Sample_size_min_peryr <- array(dim = c(length(Years), numComparisons, length(estimators)), 
                                 dimnames = list(Years, geo_order, estimators))
  
  # Calculate pairwise thetas
  for(i in 1:numComparisons){
    for(year in Years){
      
      # Indices for populations to compare
      P1 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,1]) & Data$Year == year)
      P2 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,2]) & Data$Year == year)
      P1_name <- paste(as.character(Pair_wise_site_comparisons[i,1]), year, sep = ":", collapse = '')
      P2_name <- paste(as.character(Pair_wise_site_comparisons[i,2]), year, sep = ":", collapse = '')
      
      if(sum(P1) > sample_size_threshold & sum(P2) > sample_size_threshold){
        
        # Store sample size minimum
        Sample_size_min_peryr[year, geo_order[i], c('WandC', 'Reich')] <- min(sum(P1), sum(P2))
        
        # Calculate Fst
        X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBinary)
        # X_WandC <- calculate_pairwise_WandC(P1, P2, SNPDataBinary)
        
        Pair_wise_site_comparisons_Fst_peryr[year,i, 'Reich'] <- X_Reich$F_st
        Pair_wise_site_comparisons_Fst_perSNP_peryr[,year,i,'Reich'] <- X_Reich$F_st_snp
        # Pair_wise_site_comparisons_Fst_peryr[year, i,'WandC'] <- X_WandC$theta_loci_hat_combined
        # Pair_wise_site_comparisons_Fst_perSNP_peryr[,year,i,'WandC'] <- X_WandC$theta_hat
        
        # if(sum(P1) > sample_size_threshold_freqs & sum(P2) > sample_size_threshold_freqs){
        #   Sample_size_min_peryr[year, geo_order[i], 'Inferred'] <- min(sum(P1), sum(P2))  
        #   X_inferred <- calculate_pairwise_inferred(P1_name, P2_name, 1:numSNPs) 
        #   Pair_wise_site_comparisons_Fst_peryr[year,i,'Inferred'] <- X_inferred$Fst   # Save to matrix
        #   Pair_wise_site_comparisons_Fst_perSNP_peryr[,year,i,'Inferred'] <- X_inferred$Fst_snp
        # } 
      }
    }
  }
  
  Fst_barcode <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst, 
                      Pair_wise_site_comparisons_Fst_perSNP = Pair_wise_site_comparisons_Fst_perSNP,
                      Sample_size_min = Sample_size_min,
                      Pair_wise_site_comparisons_Fst_peryr = Pair_wise_site_comparisons_Fst_peryr,
                      Pair_wise_site_comparisons_Fst_perSNP_peryr = Pair_wise_site_comparisons_Fst_perSNP_peryr, 
                      Sample_size_min_peryr = Sample_size_min_peryr) 
  
  return(Fst_barcode)
}


# =============================================================================
# Function to calculate peturbed estimates for significance of estimates
# =============================================================================
fst_perturbed <- function(SNPDataBinary, Data){
  
  numSNPs <- ncol(SNPDataBinary)
  numSamples <- nrow(SNPDataBinary)
  Years <- c("08", "09", "10") #as.character(unique(Data[,1])) # For comparisons per year
  perturb <- function(x){x[sample(length(x))]}
  
  #----------------------------------------------------------------------------
  # Between sites averaging over years 
  #----------------------------------------------------------------------------
  foreach_return <- foreach(n = 1:n_repeat) %dopar% {
    
    # Allocate memory
    Pair_wise_site_comparisons_Fst <- array(dim = c(numComparisons, length(estimators)),
                                            dimnames = list(geo_order,estimators))
    Pair_wise_site_comparisons_Fst_peryr <- array(dim = c(length(Years),numComparisons, length(estimators)),
                                                  dimnames = list(Years, geo_order, estimators))
    
    # Preserve copy of unperturbed
    Location_true <- Data$Location.code
    
    for(i in 1:numComparisons){
      
      # We only want to perturb the labels of clinics A and B, not clincs A, B, C and D
      Location_temp <- Location_true
      ClinicA <- as.character(Pair_wise_site_comparisons[i,1])   
      ClinicB <- as.character(Pair_wise_site_comparisons[i,2])
      ClinicABInd <- Location_true == ClinicA | Location_true == ClinicB 
      Location_temp[ClinicABInd] <- perturb(Location_true[ClinicABInd]) 
      
      # Indices for populations to compare
      P1 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,1]))
      P2 <- (Location_temp == as.character(Pair_wise_site_comparisons[i,2]))
      P1_name <- as.character(Pair_wise_site_comparisons[i,1])
      P2_name <- as.character(Pair_wise_site_comparisons[i,2])
      
      # Calculate Fst
      X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBinary)
      # X_WandC <- calculate_pairwise_WandC(P1, P2, SNPDataBinary)
      # X_inferred <- calculate_pairwise_inferred(P1_name, P2_name, 1:numSNPs) 
      
      # Save to matrix
      Pair_wise_site_comparisons_Fst[i, 'Reich'] <- X_Reich$F_st
      # Pair_wise_site_comparisons_Fst[i, 'WandC'] <- X_WandC$theta_loci_hat_combined
      # Pair_wise_site_comparisons_Fst[i, 'Inferred'] <- X_inferred$Fst
      rm(Location_temp)
    }
    rm(Location_true)
    
    # For comparisons over years
    for(year in Years){
      
      ind <- Data$Year == year
      SNPDataYear <- SNPDataBinary[ind,]
      Location_true <- Data$Location.code[ind] # make new copy of unperturbed values
    
      for(i in 1:numComparisons){
        
        # We only want to perturb the labels of clinics A and B, not clincs A, B, C and D
        Location_temp <- Location_true
        ClinicA <- as.character(Pair_wise_site_comparisons[i,1])   
        ClinicB <- as.character(Pair_wise_site_comparisons[i,2])
        ClinicABInd <- Location_true == ClinicA | Location_true == ClinicB 
        Location_temp[ClinicABInd] <- perturb(Location_true[ClinicABInd]) 
        
        # Indices for populations to compare
        P1 <- Location_temp == as.character(Pair_wise_site_comparisons[i,1])
        P2 <- Location_temp == as.character(Pair_wise_site_comparisons[i,2]) 
        P1_name <- paste(as.character(Pair_wise_site_comparisons[i,1]), year, sep = ":", collapse = '')
        P2_name <- paste(as.character(Pair_wise_site_comparisons[i,2]), year, sep = ":", collapse = '')
        
        # Calculate estimates WandC and Reich
        # Perturbation of the inferred fst values would require re-estimation of all the frequencies upon every perturbation, which was not done 
        X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataYear)
        Pair_wise_site_comparisons_Fst_peryr[year,i, 'Reich'] <- X_Reich$F_st
        # X_WandC <- calculate_pairwise_WandC(P1, P2, SNPDataBinary)
        # Pair_wise_site_comparisons_Fst_peryr[year, i,'WandC'] <- X_WandC$theta_loci_hat_combined
        
      }
      
    }
    Fst_perturbed <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst, 
                          Pair_wise_site_comparisons_Fst_peryr = Pair_wise_site_comparisons_Fst_peryr)
    return(Fst_perturbed)
  }
  
  # Re-structure result in an array, rather than in a list
  X <- do.call(abind, args = list(sapply(foreach_return, FUN = function(x){x['Pair_wise_site_comparisons_Fst']}), along = 3))
  Y <- do.call(abind, args = list(sapply(foreach_return, FUN = function(x){x['Pair_wise_site_comparisons_Fst_peryr']}), along = 4))
  Fst_perturbed <- list(Pair_wise_site_comparisons_Fst = X, Pair_wise_site_comparisons_Fst_peryr = Y)
  return(Fst_perturbed)
}


# =============================================================================
# Function to calculate bootstrap estimates for confidence intervals
# =============================================================================
fst_bootstrapped <- function(SNPDataBinary, Data){
  
  numSNPs <- ncol(SNPDataBinary)
  numSamples <- nrow(SNPDataBinary)
  Years <- c("08", "09", "10") #as.character(unique(Data[,1])) # For comparisons per year
  
  # Between sites averaging over years -----------------------------------------
  foreach_return <- foreach(n = 1:n_repeat) %dopar% {
    
    # Allocate memory
    Pair_wise_site_comparisons_Fst <- array(dim = c(numComparisons, length(estimators)),
                                            dimnames = list(geo_order,estimators))
    Pair_wise_site_comparisons_Fst_peryr <- array(dim = c(length(Years),numComparisons, length(estimators)),
                                                  dimnames = list(Years, geo_order, estimators))
    
    # Bootstrap SNPs
    bootstrap_snps <- sample(numSNPs, replace = TRUE)
    SNPDataBootStrap <- SNPDataBinary[, bootstrap_snps]
    
    for(i in 1:numComparisons){
      
      # Indices for populations to compare
      P1 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,1]))
      P2 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,2]))
      P1_name <- as.character(Pair_wise_site_comparisons[i,1])
      P2_name <- as.character(Pair_wise_site_comparisons[i,2])
      
      # Calculate pairwise thetas Weir and Cockerham
      X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBootStrap)
      # X_WandC <- calculate_pairwise_WandC(P1, P2, SNPDataBootStrap)
      # X_inferred <- calculate_pairwise_inferred(P1_name, P2_name, bootstrap_snps)  
      
      # Save to matrix
      Pair_wise_site_comparisons_Fst[i, 'Reich'] <- X_Reich$F_st
      # Pair_wise_site_comparisons_Fst[i, 'WandC'] <- X_WandC$theta_loci_hat_combined
      # Pair_wise_site_comparisons_Fst[i, 'Inferred'] <- X_inferred$Fst
      
      # For comparisons over years
      for(year in Years){
        
        # Indices for populations to compare
        P1 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,1]) & Data$Year == year)
        P2 <- (Data$Location.code == as.character(Pair_wise_site_comparisons[i,2]) & Data$Year == year)
        P1_name <- paste(as.character(Pair_wise_site_comparisons[i,1]), year, sep = ":", collapse = '')
        P2_name <- paste(as.character(Pair_wise_site_comparisons[i,2]), year, sep = ":", collapse = '')
        
        if(sum(P1) > sample_size_threshold & sum(P2) > sample_size_threshold){
          
          # Calculate pairwise thetas Weir and Cockerham and save to matrix
          X_Reich <- calculate_pairwise_reich(P1, P2, SNPDataBootStrap)
          Pair_wise_site_comparisons_Fst_peryr[year,i, 'Reich'] <- X_Reich$F_st
          # X_WandC <- calculate_pairwise_WandC(P1, P2, SNPDataBootStrap)
          # Pair_wise_site_comparisons_Fst_peryr[year, i,'WandC'] <- X_WandC$theta_loci_hat_combined
          
          # # If there are frequencies available for said combination, calculate inferred and save
          # if(sum((Data[,4] == as.character(Pair_wise_site_comparisons[i,1]) & Data$Year == year)) > sample_size_threshold_freqs &
          #    sum((Data[,4] == as.character(Pair_wise_site_comparisons[i,2]) & Data$Year == year)) > sample_size_threshold_freqs){
          #   X_inferred <- calculate_pairwise_inferred(P1_name, P2_name, bootstrap_snps)
          #   Pair_wise_site_comparisons_Fst_peryr[year, i, 'Inferred'] <- X_inferred$Fst 
          # }
        }
      }
    }
    Fst_bootstrapped <- list(Pair_wise_site_comparisons_Fst = Pair_wise_site_comparisons_Fst, 
                             Pair_wise_site_comparisons_Fst_peryr = Pair_wise_site_comparisons_Fst_peryr)
    
    return(Fst_bootstrapped)
  }
  
  # Pull out CIs and re-structure result in an array
  X <- do.call(abind, args = list(sapply(foreach_return, FUN = function(x){x['Pair_wise_site_comparisons_Fst']}), along = 3))
  Y <- do.call(abind, args = list(sapply(foreach_return, FUN = function(x){x['Pair_wise_site_comparisons_Fst_peryr']}), along = 4))
  Fst_bootstrapped <- list(Pair_wise_site_comparisons_Fst = X, Pair_wise_site_comparisons_Fst_peryr = Y)
  return(Fst_bootstrapped)
}


# =============================================================================
# Generate Fst estimates using functions above
# =============================================================================
# < 1 sec
system.time(Fst_barcode <- fst_calculations(SNPDataBinary = Data_store$SNPData_no_multiclonal,
                                            Data = Data_store$Data_no_multiclonal))

# Add perturbed (120.412 sec)
system.time(Fst_barcode_perturbed <- fst_perturbed(SNPDataBinary = Data_store$SNPData_no_multiclonal,
                                                   Data = Data_store$Data_no_multiclonal))
Fst_barcode$Pair_wise_site_comparisons_Fst_perturbed <- Fst_barcode_perturbed$Pair_wise_site_comparisons_Fst
Fst_barcode$Pair_wise_site_comparisons_Fst_peryr_perturbed <- Fst_barcode_perturbed$Pair_wise_site_comparisons_Fst_peryr

# Calculate CIs (116.760 sec)
system.time(Fst_barcode_bootstrapped <- fst_bootstrapped(SNPDataBinary = Data_store$SNPData_no_multiclonal, 
                                                         Data = Data_store$Data_no_multiclonal))
A <- alply(Fst_barcode_bootstrapped$Pair_wise_site_comparisons_Fst, 3) 
A_deltas <- lapply(A, FUN = function(x){x[geo_order, estimators] - Fst_barcode$Pair_wise_site_comparisons_Fst[geo_order, estimators]}) # Calculate differences
A_deltas2 <- aperm(do.call(abind, args = list(A_deltas, along = 3)), c(3,1,2)) # Reorder dimensions 
A_percentiles <- alply(apply(A_deltas2, c(2,3), quantile, probs = c(0.025, 0.975), na.rm = TRUE), 1, .dims = TRUE) 
Fst_barcode$Pair_wise_site_comparisons_Fst_CIs <- lapply(A_percentiles, FUN = function(x){Fst_barcode$Pair_wise_site_comparisons_Fst[geo_order, estimators] - x[geo_order, estimators]})

Years <- c('08', '09', '10') #as.character(unique(Data_store$Data_no_multiclonal[,1]))
B <- alply(Fst_barcode_bootstrapped$Pair_wise_site_comparisons_Fst_peryr, 4)
B_deltas <- lapply(B, FUN = function(x){x[Years, geo_order, estimators] - Fst_barcode$Pair_wise_site_comparisons_Fst_peryr[Years, geo_order, estimators]}) # 3d array minus 3d array is valid!
B_deltas2 <- do.call(abind, args = list(B_deltas, along = 0)) # Reorder dimensions 
B_percentiles <- alply(apply(B_deltas2, c(2,3,4), quantile, probs = c(0.025, 0.975), na.rm = TRUE), 1, .dims = TRUE) 
Fst_barcode$Pair_wise_site_comparisons_Fst_peryr_CIs <- lapply(B_percentiles, FUN = function(x){Fst_barcode$Pair_wise_site_comparisons_Fst_peryr[Years, geo_order, estimators] - x[Years, geo_order, estimators]})

# Save in one big list
save(Fst_barcode, file ='./Fst_barcode.RData')



