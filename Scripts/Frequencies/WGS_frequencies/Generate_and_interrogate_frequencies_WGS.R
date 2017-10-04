##################################################################################################
# Script inc. plot of allele frequencies; 
# tests of differences in frequencies between clinics and years at pooled and SNP level
# using binary data and KW test as well as frequency point estimates with regression. 
# See section on "Significance by anova method" for result that features in the manuscript. 
# Script also sums of frequencies that fall with site specific CIS, where available, 
# counts of private alleles with and without comparisons to pooled null population. 
##################################################################################################

rm(list = ls())
load(file = '../../../RData/Data_store_WGS.RData')
attach(Data_store)
sites <- unique(MetaData$collection_location)
load(file = '../../../RData/Barcode_frequencies_no_multiclonal.RData')


#==================================================================================
# Plot of frequency estimates (as hist, since too many SNPs as barplot)
# and compare to barcode
#==================================================================================
WGS_freq <- colMeans(SNPData, na.rm = TRUE)
hist(sort(pmin(1-WGS_freq, WGS_freq)), col = 'grey', 
     xlab = 'Allele frequency', ylab = 'Density', main = '', freq = FALSE) 
hist(1-FreqResultsStore_counting$all, add = TRUE, freq = FALSE, 
     col = adjustcolor('blue', alpha.f = 0.5)) # Plot as minor
legend('topright', fill = c('grey', adjustcolor('blue', alpha.f = 0.5)), 
       bty = 'n', legend = c('WGS', 'Barcode'))


#==================================================================================
# Kruskall Wallis tests using SNP data (0s and 1s)
#==================================================================================

# Site indices
MLA <- (MetaData$collection_location == 'MLA')
MKT <- (MetaData$collection_location == 'MKT')
MKK <- (MetaData$collection_location == 'MKK')
WPA <- (MetaData$collection_location == 'WPA')

# ---------------------------------------------------------------------------------
# 1) Do pooled frequencies differ between sites? Yes but...
# ---------------------------------------------------------------------------------
DataList <- list(x = SNPData[MLA,],  
                 y = SNPData[MKT,],
                 z = SNPData[MKK,],
                 h = SNPData[WPA,])
kruskal.test(lapply(DataList, colMeans,TRUE))

# ---------------------------------------------------------------------------------
# 2) How many individual snps differ between sites? 257/34911 (<1%)
# ---------------------------------------------------------------------------------
pvalue_store <- rep(NA, numSNPs)
for(i in 1:numSNPs){
  DataList <- list(x = SNPData[MLA,i], 
                   y = SNPData[MKT,i],
                   z = SNPData[MKK,i],
                   h = SNPData[WPA,i])  
  X <- try(kruskal.test(DataList), TRUE)
  if(class(X) ==  "htest"){
    pvalue_store[i] <- X$p.value
  }
}
sum(pvalue_store < (0.05/numSNPs), na.rm = TRUE)  

# ---------------------------------------------------------------------------------
# 3) Do pooled frequencies differ over year? Yes, but ...
# ---------------------------------------------------------------------------------
# unique years: 14  8 12 11  1  2  4  3
Yr14 <- (MetaData$year == 14)
Yr12 <- (MetaData$year == 12)
Yr11 <- (MetaData$year == 11)
Yr8 <- (MetaData$year == 8)
Yr1 <- (MetaData$year == 1)
Yr2 <- (MetaData$year == 2)
Yr3 <- (MetaData$year == 3)
Yr4 <- (MetaData$year == 4)
DataList <- list(SNPData[Yr14,], 
                 SNPData[Yr12,], 
                 SNPData[Yr11,],
                 SNPData[Yr8,],
                 SNPData[Yr1,],
                 SNPData[Yr2,],
                 SNPData[Yr3,],
                 SNPData[Yr4,])
kruskal.test(lapply(DataList, colMeans, na.rm = TRUE))

# ---------------------------------------------------------------------------------
# 4) How many SNPs differ between years? 504/34911 (1.4%)
# ---------------------------------------------------------------------------------
pvalue_store <- rep(NA, numSNPs)
for(i in 1:numSNPs){
  DataList <- list(SNPData[Yr14,i], 
                    SNPData[Yr12,i], 
                    SNPData[Yr11,i],
                    SNPData[Yr8,i],
                    SNPData[Yr1,i],
                    SNPData[Yr2,i],
                    SNPData[Yr3,i],
                    SNPData[Yr4,i])
  X <- try(kruskal.test(DataList), TRUE)
  if(class(X) ==  "htest"){
    pvalue_store[i] <- X$p.value
  }
}
sum(pvalue_store < (0.05/numSNPs), na.rm = TRUE) 



#==================================================================================
# Regression tests over frequency point estimates 
# (since frequencies are used in the IBD model) 
#==================================================================================

# ---------------------------------------------------------------------------------
# Regression prep: generate frequencies by site and times
# ---------------------------------------------------------------------------------
Stage <- TRUE # choose between stage or year
if(Stage){
  MetaData$time <- MetaData$stage} else {
    MetaData$time <- MetaData$year
  }
times <- unique(MetaData$time)
sites <- unique(MetaData$collection_location)
site_times_store <- array(NA, dim = c(length(times), length(sites), numSNPs), 
                          dimnames = list(as.character(times), sites, NULL))

for(SNP in 1:numSNPs){ # Toggle between MetaData$stage and MetaData$year
  for(time in times){
    for(site in sites){
      X <- SNPData[time == MetaData$time & site == MetaData$collection_location, SNP]
      if(length(X) > 0){
        site_times_store[as.character(time), site, SNP] <- mean(X, na.rm = TRUE)
      }
    }
  }
}

# ---------------------------------------------------------------------------------
# Significance by anova 
# ---------------------------------------------------------------------------------
anova_fit_indv <- apply(as.matrix(Data_wide_format[,-(1:2)]), 2, function(x){
  anova(lm(x ~ Data_wide_format$site + Data_wide_format$times))
})
p_value_store <- t(sapply(anova_fit_indv, function(x){
  x[c('Data_wide_format$site', 'Data_wide_format$times'),'Pr(>F)']
}))
100*(apply(p_value_store, 2, function(x){sum(x < 0.05, na.rm = TRUE)})/numSNPs) # Before bonferroni correction
100*(apply(p_value_store, 2, function(x){sum(x < 0.05/numSNPs, na.rm = TRUE)})/numSNPs) # After bonferroni correction



