##################################################################
# Script to Sensitivity of trend to threshold 
# Takes 7602.061 sec with nrep = 100 on pro to do both barcode and WGS
##################################################################

rm(list = ls())
nrep <- 100
air <- TRUE
if(air){source('../../FunctionFiles/simtests.R')}else{source('./simtests.R')}
if(air){source('../../FunctionFiles/glm_trends.R')}else{source('./glm_trends.R')}
if(air){load('../../RData/formulas.RData')}else{load('./formulas.RData')}
if(air){load('../../RData/WGS_threshold.RData')}else{load('./WGS_threshold.RData')}
if(air){load('../../RData/Barcode_threshold.RData')}else{load('./Barcode_threshold.RData')}
if(air){path <- '../../RData/'}else{path <- './'}
attach(formulas, warn.conflicts = FALSE)
require(tictoc)

sensitivity_threshold <- function(X, glmformula, distances, nrep){
  
  threshold_IBD <- 0.5 # Fixed threshold for IBD
  translation <- 0.3
  no_points <- 8
  percentile_results <- vector('list', length = 2) # Store results 
  names(percentile_results) <- names(distances)
  
  for(i in distances){
    
    if(grepl('IBD', i)){ # If IBD, define upper and lower by distance around fixed threshold 
      probs <- c(threshold_IBD, seq(threshold_IBD - translation, min(threshold_IBD + translation, 0.999), length.out = no_points))
    } else { # If IBS, define upper and lower by symmetry
      dens <- density(X[,i])
      mode_ibs <- dens$x[dens$y == max(dens$y)]
      min_ibs <- min(dens$x)
      threshold_IBS <- mode_ibs + (mode_ibs - min_ibs)  
      probs <- c(threshold_IBS, seq(max(threshold_IBS - translation, min_ibs + 0.001), min(threshold_IBS + translation, 0.999), length.out = no_points))
    }
    
    percentile_results_store <- array(dim = c(length(probs), 2), dimnames = list(as.character(probs), c('beta', 'p-value')))
    
    for(prob in probs){
      
      # Re-set tail, fit glm and extract results
      X$Response <- X[,i] > prob 
      glm_results <- glm_trends(X, distances = sprintf('%s_tail',i), glmformula, nrep, years = NULL, response_specified = TRUE, speed = TRUE)
      percentile_results_store[as.character(prob), 'beta'] <- glm_results[[sprintf('All %s_tail',i)]]$obs['geo_dist']
      percentile_results_store[as.character(prob), 'p-value'] <- glm_results[[sprintf('All %s_tail',i)]]$pvalue['geo_dist'] 
      rm(list = c('glm_results'))
    }
    percentile_results[[i]] <- percentile_results_store
  }
  return(percentile_results)
}

tic()

#====================== Barcode ======================#
unadj <- sensitivity_threshold(Barcode, unadjusted, distances = 'ProbIBD_93', nrep = nrep)
adj <- sensitivity_threshold(X = Barcode, glmformula = adjusted, distances = 'ProbIBD_93', nrep = nrep)
sensitivity_results <- list(sensitivity_unadjusted = unadj, sensitivity_adjusted = adj)
save(sensitivity_results, file = sprintf('%sBarcode_sensitivity.RData', path))
     
#====================== WGS ======================#
unadj <- sensitivity_threshold(WGS, unadjusted, distances = 'ProbIBD', nrep)
adj <- sensitivity_threshold(WGS, adjusted, distances = 'ProbIBD', nrep)
adjyear <- sensitivity_threshold(WGS, adjustedyear, distances = 'ProbIBD', nrep)
sensitivity_results <- list(sensitivity_unadjusted = unadj, sensitivity_adjusted = adj, sensitivity_adjustedyear = adjyear)
save(sensitivity_results, file = sprintf('%sWGS_sensitivity.RData', path))

toc()