##############################################################################################
# Script to estimate frequencies by proportions using barcode data 
# Missing data are ommitted if present
# Estimates are generated with minimum of one obs. 
##############################################################################################

rm(list = ls())
load('../../../RData/Data_store_Barcode.RData')
RawData <- Data_store$Data_no_multiclonal
SNPData <- Data_store$SNPData_no_multiclonal
years <- unique(RawData$Year)
sites <- unique(RawData$Location.code)
nSNPs <- ncol(SNPData)

# Site and year
site_years_store <- array(NA, dim = c(length(years), length(sites), nSNPs), dimnames = list(years, sites, NULL))
for(SNP in 1:nSNPs){
  for(year in years){
    for(site in sites){
      X <- SNPData[year == RawData$Year & site == RawData$Location.code,SNP]
      if(length(X) > 0){ 
        site_years_store[year, site, SNP] <- mean(X, na.rm = TRUE)
      }
    }
  }
}

# Sites only 
site_store <- array(NA, dim = c(length(sites), nSNPs), dimnames = list(sites, NULL))
for(SNP in 1:nSNPs){
  for(year in years){
    for(site in sites){
      X <- SNPData[site == RawData$Location.code,SNP]
      if(length(X) > 0){
        site_store[site, SNP] <- mean(X, na.rm = TRUE)
      }
    }
  }
}

# Unstratified
all_store <- colMeans(SNPData, na.rm = TRUE)

# Save results
FreqResultsStore_counting <- list(site_years = site_years_store, 
                                  site = site_store, 
                                  all = all_store)

save(FreqResultsStore_counting, file = '../../../RData/Barcode_frequencies_no_multiclonal.RData')




