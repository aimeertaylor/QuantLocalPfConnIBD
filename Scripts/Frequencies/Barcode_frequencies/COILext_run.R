#####################################################################
# MCMC to estimate allele frequencies and MOIs jointly accross barcode 
#
# Some very low acceptance rates on COIs but perhaps indicative of highly concentrated posterior
# Liklihood = 1 is a manual inside coil_extension (repeated each time code changed)
# Note when liklihood = 1, f accpetance much higher -> use my_sd=0.46
# and run for many iterations
#
#
# Activity log: 
# July 6th: generated frequencies per year and per site. 
# Wed 20th Aug: Used the average over all posterior estimates to run HMM on barcode in series
# Oct 18th 2016: 
# 1) TM_frequency_results_site.RData etc. were generated using data not inc. missing (deleted Aug 2017)
# 2) I ran all three partitions inc. missing data with 3 chains each in parallel
# Aug 2017:
# Moved directory, deleted not inc. missing RData, moved others and changed paths according (in case ever re-run), 
# but did not re-run any results, as they were not included in final manuscript draft (deemed unnecessary)
#####################################################################

rm(list = ls())
require(MCMCpack)
source('./coil_extension.R')
load('../../RData/Data_store_Barcode.RData')

# Filter missing data
MetaData <- Data_store_list$Data
SNPData <- Data_store_list$SNPData

# Hyperparameters
psi <- 1
phi <- 1
alpha <- 1
beta <- 1
no_iterations <- 10000
no_chains <- 3 
cores_max <- 3

# Per site/year (for Fst) -----------------------------
years <- unique(MetaData$Year)
sites <- unique(MetaData$Location.code)

results_store <- vector('list', 0)
for(site in sites){
  for(year in years){

    Data <- SNPData[MetaData$Year == year & MetaData$Location.code == site, ]

    # Number multiclonal
    no_multiclonal <- sum(apply(Data,1,FUN = function(x){any(x == 0.5)}), na.rm = TRUE)
    print(sprintf('%s, year %s: %s samples, %s of which are multiclonal', site, year, nrow(Data), no_multiclonal))

    if(nrow(Data) > 9 & no_multiclonal > 1){
      set.seed(1)
      results <- coil_extension(Data, psi, phi, alpha, beta,
                                no_iterations, no_chains, cores_max)
      results_store[[paste(site, year, sep = ':')]] <- results
    }
  }
}
#save(results_store, file = '../../RData/COILext_results_site_year.RData')


# Per site (for Fst) --------------------------------------------
results_store <- vector('list', 0)

for(site in sites){
  
  Data <- SNPData[MetaData$Location.code == site, ]
  no_multiclonal <- sum(apply(Data,1,FUN = function(x){any(x == 0.5)}), na.rm = TRUE)  # Number multiclonal
  print(sprintf('%s, %s samples, %s of which are multiclonal', site, nrow(Data), no_multiclonal))

  if(nrow(Data) > 9 & no_multiclonal > 1){
    set.seed(1)
    results <- coil_extension(Data, psi, phi, alpha, beta,
                              no_iterations, no_chains, cores_max)
    results_store[[site]] <- results
  }
}
save(results_store, file = '../../RData/COILext_results_site.RData')


# All sites and years combined -----------------------------------------
Data <- SNPData
no_multiclonal <- sum(apply(Data,1,FUN = function(x){any(x == 0.5)}),na.rm = TRUE); print(no_multiclonal) # Number multiclonal
set.seed(1)
results_store <- coil_extension(Data, psi, phi, alpha, beta, no_iterations, no_chains, cores_max)
#save(results_store, file = '../../RData/COILext_results_all.RData')

