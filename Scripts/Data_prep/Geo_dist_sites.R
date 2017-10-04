################################################################################################
# Script to calculate distances between clinics by radian lat/long using the Haversine formula 
################################################################################################

rm(list = ls())
source('../FunctionFiles/gcd.hf.R')
require(gtools)

# Import data
Lonlat <- read.csv('../RawData/SMRU_clinics_lonlat.csv')
sites <- Lonlat$clinic[-4]
site_combinations <- combinations(n = 4, r = 2, as.character(sites))
pairwise_site_distance <- data.frame(site_combinations)
pairwise_site_distance$distance <- NA

for(i in 1:nrow(site_combinations)){
  site1 <- site_combinations[i,1]
  site2 <- site_combinations[i,2]
  
  pairwise_site_distance[i, 'distance'] <- gcd.hf(A = as.matrix(Lonlat[Lonlat$clinic == site1, 3:2]), 
                                                  B = as.matrix(Lonlat[Lonlat$clinic == site2, 3:2]))
}

# Sort pairwise distance
pairwise_site_distance <- pairwise_site_distance[sort(pairwise_site_distance$distance, 
                                                      index.return = TRUE)$ix,]
# Extract geo_order
geo_order <- apply(pairwise_site_distance[,c(1,2)], 1, paste, sep = '_', collapse = '_')

# Make matrix with all distances 
pairwise_site_distance_all <- rbind(pairwise_site_distance, 
                                    data.frame(X1 = pairwise_site_distance$X2,
                                               X2 = pairwise_site_distance$X1,
                                               distance = pairwise_site_distance$distance), data.frame(X1 = sites, X2 = sites, distance = 0))

names_pairwise_site_distance_all <- apply(pairwise_site_distance_all[,c(1,2)], 1, paste, sep = '_', collapse = '_')
pairwise_site_distance_all <- pairwise_site_distance_all$distance
names(pairwise_site_distance_all) <- names_pairwise_site_distance_all

geo_dist_info <- list(geo_order = geo_order, 
                      pairwise_site_distance = pairwise_site_distance,
                      pairwise_site_distance_all = pairwise_site_distance_all)

save(geo_dist_info, file = '../RData/geo_dist_info.RData')


