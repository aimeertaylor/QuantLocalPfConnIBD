# This script extracts the median and 95% CI of pop-level MOIs and frequencies
# Aimee Taylor, July 2016 

rm(list = ls())
require(MCMCpack)

# load data to get snp names
load("../../RData/Data_store_Barcode.RData")
MetaData <- Data_store$Data_no_multiclonal_missing
snp_names <- colnames(MetaData[-(1:6)])


# Results based on all data combined -----------------------------------------------
load('../../RData/COILext_results_all.RData')
burnin <- 1:(length(results_store[[1]]$lambda_store)/2)

# Pop-level COI results
lambda_results <- summary(mcmc.list(as.mcmc(results_store[[1]]$lambda_store[-burnin]),
                       as.mcmc(results_store[[2]]$lambda_store[-burnin]),
                       as.mcmc(results_store[[3]]$lambda_store[-burnin])))
MOI_quants <- array(1/lambda_results$quantiles[c('50%', '2.5%', '97.5%')], dim = 3, 
                   dimnames = list(c('50%', '2.5%', '97.5%')))

# Frequency results
f_results <- summary(mcmc.list(as.mcmc(results_store[[1]]$f_vector_store[-burnin,]),
                       as.mcmc(results_store[[2]]$f_vector_store[-burnin,]),
                       as.mcmc(results_store[[3]]$f_vector_store[-burnin,])))

frequency_quants <- array(f_results$quantiles[,c('50%', '2.5%', '97.5%')], dim = c(93, 3), 
                         dimnames = list(NULL, c('50%', '2.5%', '97.5%')))

minor_inferred_freq_store <- pmin(1 - frequency_quants[,'50%'], frequency_quants[,'50%'])

# Save 
save(frequency_quants, file = '../../RData/COIext_f_extracted_all.RData')
save(MOI_quants, file = '../../RData/COIext_MOIpop_extracted_all.RData')

# The below file was used to run the HMM to see if it gave different results when using frequencies based on proportions 
# compared with frequencies inferred under COI extension model. However, it has since been archived for deletion
# write.table(t(minor_inferred_freq_store), file = '../../TxtData/COIext_minor_inferred_freq_store_withnames.txt',
#             sep = ',', col.names = c(snp_names), row.names = FALSE)


# Results based on sites/site_years -----------------------------------------------
for(partition in c('site_year','site')){ 
  
  load(sprintf('../../RData/COILext_results_%s.RData',partition))
  
  # Allocate space
  MOI_quants <- array(NA, dim = c(3, length(names(results_store))), 
                     dimnames = list(c('50%', '2.5%', '97.5%'), 
                                     names(results_store)))
  frequency_quants <- array(NA, dim = c(93, 3, length(names(results_store))), 
                           dimnames = list(NULL, c('50%', '2.5%', '97.5%'), 
                                           names(results_store)))
  
  for(dataset in names(results_store)){
    
    # Lambda results
    lambda_results <- summary(mcmc.list(as.mcmc(results_store[dataset][[1]][[1]]$lambda_store[-burnin]),
                                as.mcmc(results_store[dataset][[1]][[2]]$lambda_store[-burnin]),
                                as.mcmc(results_store[dataset][[1]][[3]]$lambda_store[-burnin]))) 
    MOI_quants[c('50%', '2.5%', '97.5%'), dataset] <- 1/lambda_results$quantiles[c('50%', '2.5%', '97.5%')]
    
    # Frequency results
    f_results <- summary(mcmc.list(as.mcmc(results_store[dataset][[1]][[1]]$f_vector_store[-burnin,]),
                           as.mcmc(results_store[dataset][[1]][[2]]$f_vector_store[-burnin,]),
                           as.mcmc(results_store[dataset][[1]][[3]]$f_vector_store[-burnin,])))
    frequency_quants[,c('50%','2.5%','97.5%'), dataset] <- f_results$quantiles[,c('50%','2.5%','97.5%')]
  }
  
  save(MOI_quants, file = sprintf('../../RData/COIext_MOIpop_extracted_%s.RData',partition))
  save(frequency_quants, file = sprintf('../../RData/COIext_f_extracted_%s.RData',partition))
}












