##############################################################
# Script to generate WGS matrix with all distance measures: 
# IBD and IBS whole genome
# IBD 93 and 24
# 22nd May 2014: military coup, Thailand 
##############################################################
rm(list = ls())
require(gtools) # For combinations
source('../../FunctionFiles/proportions.R')

load('../../RData/Data_store_WGS.RData')
load('../../RData/geo_dist_info.RData')
attach(Data_store, warn.conflicts = FALSE)
attach(geo_dist_info, warn.conflicts = FALSE)
Additional_metadata <- read.csv('../../RawData/clearance_rate.csv',nrows = 198)
row.names(Additional_metadata) <- as.character(Additional_metadata$Sample_ID)

# Returns all r choose 2 combinations for nSamples where combinations() fails (because nSamples too large) 
nSamples <- nrow(SNPData)
ind_j <- rep(1:nSamples, (nSamples:1)-1)
ind_i <- c()
for(i in 2:nSamples){
  ind_i <- c(ind_i, i:nSamples)
}

# Create WGS matrix with metadata --------------------------------------------------------
WGS <- data.frame(Sitei = as.character(MetaData$collection_location[ind_i]))
WGS$Sitej <- as.character(MetaData$collection_location[ind_j])
WGS$Site_comparison <- apply(cbind(as.character(WGS$Sitei), as.character(WGS$Sitej)),
                                 1, FUN = function(x){paste(sort(x), collapse = '_')})
WGS$geo_dist <- pairwise_site_distance_all[WGS$Site_comparison] # MaeRaMat introduces NAs here
WGS$SampleIDi <- as.character(MetaData$strain[ind_i])
WGS$SampleIDj <- as.character(MetaData$strain[ind_j])
WGS$Sample_comparison <- paste(WGS$SampleIDj, WGS$SampleIDi, sep = '_')
rownames(WGS) <- as.character(WGS$Sample_comparison) # Name rows

# Collection date (not available for sequencing)
Collection_date <- strptime(x = Additional_metadata[rownames(SNPData),'Sample_Date'], 
                            format = '%d-%b-%y') # Make sure Dates in the same order as SNPData

# ======================================================================
# Aside: Are 2014 samples before after military coup? 
ind_14 <- format(Collection_date, format = '%y') == '14'
mean(difftime(Collection_date[ind_14], "2014-05-22") > 0) # 80% after the coup
# ======================================================================

WGS$Collection_datei <- Collection_date[ind_i]
WGS$Collection_datej <- Collection_date[ind_j]
WGS$time_dist <- abs(difftime(WGS$Collection_datei, WGS$Collection_datej))

# Seasonal indicator: Jan-Mar, April-Jun (high), July-Sep, Oct-Dec (low)
WGS$Monthi <- format(WGS$Collection_datei, format = '%b') 
WGS$Monthj <- format(WGS$Collection_datej, format = '%b') 
WGS$Seasoni <- WGS$Seasonj <- NA
WGS$Seasoni[WGS$Monthi %in% c("Jan", "Feb", "Mar")] <- 'winter'
WGS$Seasonj[WGS$Monthj %in% c("Jan", "Feb", "Mar")] <- 'winter'
WGS$Seasoni[WGS$Monthi %in% c("Apr", "May", "Jun")] <- 'spring'
WGS$Seasonj[WGS$Monthj %in% c("Apr", "May", "Jun")] <- 'spring'
WGS$Seasoni[WGS$Monthi %in% c("Jul", "Aug", "Sep")] <- 'summer'
WGS$Seasonj[WGS$Monthj %in% c("Jul", "Aug", "Sep")] <- 'summer'
WGS$Seasoni[WGS$Monthi %in% c("Oct", "Nov", "Dec")] <- 'autumn'
WGS$Seasonj[WGS$Monthj %in% c("Oct", "Nov", "Dec")] <- 'autumn'
WGS$Season_comparison <- apply(WGS[,c("Seasoni", "Seasonj")],
                                   1, FUN = function(x){paste(sort(x), collapse = '_')})
WGS$Spring_Summer <- FALSE 
WGS$Spring_Summer <- WGS$Season_comparison %in% c('spring_spring','spring_summer','summer_summer')


# Add year akin to Barcode matrix
WGS$Yeari <- as.numeric(format(WGS$Collection_datei, format = '%y')) 
WGS$Yearj <- as.numeric(format(WGS$Collection_datej, format = '%y')) 

# WG-SNP IBD  --------------------------------------------------------------------------
WGS$ProbIBD <- NA
WGSf <- read.table('../../TxtData/WGS.hmm_fract.txt', sep = '\t', header = TRUE, as.is = TRUE)
WGSf$pair <- apply(cbind(WGSf$sample1, WGSf$sample2), 1, paste, collapse = '_', sep = '') 
WGS[as.character(WGSf$pair), 'ProbIBD'] <- WGSf$fract_sites_IBD

# 93-SNP IBD  --------------------------------------------------------------------------
WGS$ProbIBD_93 <- NA
WGS93 <- read.table('../../TxtData/WGS93.hmm_fract.txt', sep = '\t', header = TRUE, as.is = TRUE)
WGS93$pair <- apply(cbind(WGS93$sample1, WGS93$sample2), 1, paste, collapse = '_', sep = '') 
WGS[as.character(WGS93$pair), 'ProbIBD_93'] <- WGS93$fract_sites_IBD

# 24-SNP IBD  --------------------------------------------------------------------------
WGS$ProbIBD_24 <- NA
WGS24 <- read.table('../../TxtData/WGS24.hmm_fract.txt', sep = '\t', header = TRUE, as.is = TRUE)
WGS24$pair <- apply(cbind(WGS24$sample1, WGS24$sample2), 1, paste, collapse = '_', sep = '') 
WGS[as.character(WGS24$pair), 'ProbIBD_24'] <- WGS24$fract_sites_IBD

# Add clinics
WGS$MKT <- ((WGS$Sitei == 'MKT') & (WGS$Sitej == 'MKT'))
WGS$MKK <- ((WGS$Sitei == 'MKK') & (WGS$Sitej == 'MKK'))
WGS$MLA <- ((WGS$Sitei == 'MLA') & (WGS$Sitej == 'MLA'))
WGS$WPA <- ((WGS$Sitei == 'WPA') & (WGS$Sitej == 'WPA'))

# Add Year 2014 indicator
WGS$Year <- WGS$Yeari == 14 & WGS$Yearj == 14 
WGS$diff_weeks <- round(as.numeric(WGS$time_dist, units = "weeks"),0) # Added 5th March 2017

# Threshold
threshold_IBD <- 0.5
WGS$ProbIBD_tail <- WGS$ProbIBD > threshold_IBD
WGS$ProbIBD_93_tail <- WGS$ProbIBD_93 > threshold_IBD
WGS$ProbIBD_24_tail <- WGS$ProbIBD_24 > threshold_IBD 


# =======================  Save barcode itself =============================
save(WGS, file = '../../RData/WGS_threshold.RData')

# ==================== Calculate and save proportions  ======================
load('../../RData/WGS_threshold.RData')
system.time(proportion_results <- proportions(X = WGS, distances = "ProbIBD_tail")) 
save(proportion_results, file = c('../../RData/WGS_proportions.RData')) 


