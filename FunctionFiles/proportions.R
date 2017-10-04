# Calculate proportions of highly related parasite sample pairs
proportions <- function(X, nrep = 1000, distances = colnames(X)[grepl('tail', colnames(X))]){
  
  load('../../RData/geo_dist_info.RData')
  attach(geo_dist_info, warn.conflicts = FALSE)
  site_comparisons <- names(pairwise_site_distance_all[-(7:12)])
  proportion <- array(dim = c(length(site_comparisons),length(distances)), 
                      dimnames = list(site_comparisons, distances))
  noisolatepairs <- array(dim = c(length(site_comparisons),length(distances)), 
                          dimnames = list(site_comparisons, distances))
  delta_store <- array(dim = c(length(site_comparisons),length(distances),nrep), 
                       dimnames = list(site_comparisons, distances, NULL))
  no_pairs_perclinic <- table(X$Site_comparison) 
  
  
  for(k in site_comparisons){
    print(k)
    ind <- which(X$Site_comparison == k)
    for(distance in distances){
      proportion[k,distance] <- mean(X[ind, distance], na.rm = TRUE) 
      noisolatepairs[k,distance] <- sum(!is.na(X[ind, distance])) 
    }
    no_pairs_k <- no_pairs_perclinic[k]
    pb <- txtProgressBar(min = 0, max = nrep, style = 3)
    for(i in 1:nrep){
      setTxtProgressBar(pb, i)
      boostrap <- X[sample(ind, size = no_pairs_k, replace = TRUE), ]
      for(distance in distances){
        delta_store[k,distance,i] <- mean(boostrap[, distance], na.rm = TRUE) - proportion[k,distance] 
      }
    }
  }
  
  deltaCIs <- apply(delta_store, c(1,2), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  results <- list(proportion = proportion, noisolatepairs = noisolatepairs, deltaCIs = deltaCIs)
  return(results)
}
