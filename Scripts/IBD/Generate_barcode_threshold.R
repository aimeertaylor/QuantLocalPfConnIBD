##############################################################
# Script to generate Barcode matrix with distance measures
##############################################################
rm(list = ls())
remove_clones <- TRUE # Remove clones from within clinics
require(gtools) # For combinations
require(proxy)
source('../../FunctionFiles/proportions.R')
load('../../RData/geo_dist_info.RData')
load('../../RData/Data_store_Barcode.RData')
SNPData <- Data_store$SNPData_no_multiclonal # single-genotype only
MetaData <- Data_store$MetaData[rownames(SNPData), ]

attach(geo_dist_info, warn.conflicts = FALSE)

# Remove clones from within sites: (removes 250 samples)
if(remove_clones){
  SNPData <- unique(cbind(MetaData$Location.code, SNPData))[,-1]
  MetaData <- MetaData[rownames(SNPData),]
}

# How many repeats per clinic? (not as many as 250!)
for(site in c("MLA", "MKK", "MKT", "WPA")){
  site_ind <- MetaData$Location.code == site
  allnames <- rownames(SNPData[site_ind, ])
  uniquenames <- rownames(unique(SNPData[site_ind, ]))
  print(sum(!allnames %in% uniquenames)) # No. excluded per clinic
  
  # Another way to count the same thing
  site_ind <- MetaData$Location.code == site
  Tab <- table(apply(SNPData[site_ind, ], 1, paste, collapse = ''))
  print(sum(Tab[Tab > 1]-1)) # minus one because we keep on barcode for each repeat motif
}

# How manys samples after removing all-but-one identical barcode per clinic? 
c('All' = nrow(MetaData),'2008' = sum(MetaData$Year == '08'),
  '2009' = sum(MetaData$Year == '09'),'2010' = sum(MetaData$Year == '10')) 

# How many single-genotype samples after removing all-but-one identical barcode per clinic? 
ind <- MetaData$Clonality == 'S'
c('All' = nrow(MetaData[ind,]),'2008' = sum(MetaData[ind,]$Year == '08'),
  '2009' = sum(MetaData[ind,]$Year == '09'),'2010' = sum(MetaData[ind,]$Year == '10')) 

# Returns all r choose 2 combinations for nSamples where combinations() fails (because nSamples too large) 
nSamples <- nrow(SNPData)
ind_j <- rep(1:nSamples, (nSamples:1)-1)
ind_i <- c()
for(i in 2:nSamples){
  ind_i <- c(ind_i, i:nSamples)
}

# Create Barcode matrix with metadata
Barcode <- data.frame(Sitei = MetaData$Location.code[ind_i]) 
Barcode$Sitej <- MetaData$Location.code[ind_j]
Barcode$Site_comparison <- apply(cbind(as.character(Barcode$Sitei), as.character(Barcode$Sitej)),
                                 1, FUN = function(x){paste(sort(x), collapse = '_')})
Barcode$geo_dist <- pairwise_site_distance_all[Barcode$Site_comparison]
Barcode$SampleIDi <- MetaData$Sample.ID[ind_i]
Barcode$SampleIDj <- MetaData$Sample.ID[ind_j]
Barcode$Sample_comparison <- paste(Barcode$SampleIDj, Barcode$SampleIDi, sep = '_')
rownames(Barcode) <- as.character(Barcode$Sample_comparison) # Name rows
Barcode$Yeari <- as.numeric(as.character(MetaData$Year[ind_i]))
Barcode$Yearj <- as.numeric(as.character(MetaData$Year[ind_j]))

# Collection date
Collection_date <- strptime(x = as.character(MetaData$Collection.date), format = '%d-%b-%y') 
Barcode$Collection_datei <- Collection_date[ind_i]
Barcode$Collection_datej <- Collection_date[ind_j]
Barcode$time_dist <- abs(difftime(Barcode$Collection_datei, Barcode$Collection_datej))

# Seasonal indicator: Jan-Mar, April-Jun (high), July-Sep, Oct-Dec (low)
Barcode$Monthi <- format(Barcode$Collection_datei, format = '%b') 
Barcode$Monthj <- format(Barcode$Collection_datej, format = '%b') 
Barcode$Seasoni <- Barcode$Seasonj <- NA
Barcode$Seasoni[Barcode$Monthi %in% c("Jan", "Feb", "Mar")] <- 'winter'
Barcode$Seasonj[Barcode$Monthj %in% c("Jan", "Feb", "Mar")] <- 'winter'
Barcode$Seasoni[Barcode$Monthi %in% c("Apr", "May", "Jun")] <- 'spring'
Barcode$Seasonj[Barcode$Monthj %in% c("Apr", "May", "Jun")] <- 'spring'
Barcode$Seasoni[Barcode$Monthi %in% c("Jul", "Aug", "Sep")] <- 'summer'
Barcode$Seasonj[Barcode$Monthj %in% c("Jul", "Aug", "Sep")] <- 'summer'
Barcode$Seasoni[Barcode$Monthi %in% c("Oct", "Nov", "Dec")] <- 'autumn'
Barcode$Seasonj[Barcode$Monthj %in% c("Oct", "Nov", "Dec")] <- 'autumn'
Barcode$Season_comparison <- apply(Barcode[,c("Seasoni", "Seasonj")],
                                   1, FUN = function(x){paste(sort(x), collapse = '_')})
Barcode$Spring_Summer <- FALSE # Add season predictor
Barcode$Spring_Summer <- Barcode$Season_comparison %in% c('spring_spring','spring_summer','summer_summer')
Barcode$diff_weeks <- round(as.numeric(Barcode$time_dist, units = "weeks"),0) # Added 5th March 2017

# Add site variable in {1 (within), 0 (across)} ---------------------------------
Barcode$MKT <- ((Barcode$Sitei == 'MKT') & (Barcode$Sitej == 'MKT'))
Barcode$MKK <- ((Barcode$Sitei == 'MKK') & (Barcode$Sitej == 'MKK'))
Barcode$MLA <- ((Barcode$Sitei == 'MLA') & (Barcode$Sitej == 'MLA'))
Barcode$WPA <- ((Barcode$Sitei == 'WPA') & (Barcode$Sitej == 'WPA'))

# 93-SNP IBD  --------------------------------------------------------------------------
Barcode$ProbIBD_93 <- NA
barcode93 <- read.table('../../TxtData/barcode93.hmm_fract.txt', sep = '\t', header = TRUE, as.is = TRUE)
barcode93$pair <- apply(cbind(barcode93$sample1, barcode93$sample2), 1, paste, collapse = '_', sep = '') 
Barcode[as.character(barcode93$pair), 'ProbIBD_93'] <- barcode93$fract_sites_IBD

# 24-SNP IBD  --------------------------------------------------------------------------
Barcode$ProbIBD_24 <- NA
barcode24 <- read.table('../../TxtData/barcode24.hmm_fract.txt', sep = '\t', header = TRUE, as.is = TRUE)
barcode24$pair <- apply(cbind(barcode24$sample1, barcode24$sample2), 1, paste, collapse = '_', sep = '') 
Barcode[as.character(barcode24$pair), 'ProbIBD_24'] <- barcode24$fract_sites_IBD 

# Remove any excess barcode due to addition of IBD93$sites_shared_fb 
# Note indexing barcode93$fact_sites_IBD above causes R to crash, hence posthoc filter
Barcode <- Barcode[!is.na(Barcode$Sample_comparison),]

# Thresholding 
threshold_IBD <- 0.5
Barcode$ProbIBD_93_tail <- Barcode$ProbIBD_93 >= threshold_IBD 
Barcode$ProbIBD_24_tail <- Barcode$ProbIBD_24 >= threshold_IBD 


# ======================= Save barcode itself =============================
if(remove_clones){
  Barcode_nowithinclinicclones <- Barcode
  save(Barcode_nowithinclinicclones, file = '../../RData/Barcode_nowithinclinicclones_threshold.RData')
} else {
  save(Barcode, file = '../../RData/Barcode_threshold.RData')
}

# ==================== Calculate and save proportions  ======================
# 890.121 sec for 1000 # 49 sec for 50
load('../../RData/Barcode_threshold.RData')
system.time(proportion_results <- proportions(X = Barcode, nrep = 10, distances = "ProbIBD_93_tail"))
save(proportion_results, file = c('../../RData/Barcode_proportions.RData'))
