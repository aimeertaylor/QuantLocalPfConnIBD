glm_trends <- function(X, distances, glmformula, nrep = 1000, years, response_specified = TRUE, speed = FALSE){
  
  # Use either speedglm or glm, where former is ~ 5 sec faster for barcode data
  if(speed){
    require(speedglm)
    glmfunction <- speedglm
    glmclass <- 'speedglm'
  } else {
    glmfunction <- glm
    glmclass <- 'glm'}
  
  results_store <- list() 
  
  for(distance in distances){
    
    if(!response_specified){X$Response <- X[, distance]} # Response specified in sensitivity script
    X <- X[!is.na(X$Response), ] # Remove lines with NA results 
    
    for(year in c('All', years)){
      
      # Stratify data 
      if(year != 'All'){
        Z <- X[X$Yeari == as.numeric(year) & X$Yearj == as.numeric(year),]
      } else {
        Z <- X}
      
      if(length(unique(c(as.character(X$Sitei), as.character(X$Sitej)))) < 4){next()}
      
      # Observation
      fit_glm <- glmfunction(formula = as.formula(glmformula), family = binomial(logit), data = Z, na.action = na.exclude)
      obs <- fit_glm$coefficients  
      
      # Permutation
      perm <- matrix(0, ncol = nrep, nrow = length(obs), dimnames = list(names(obs), NULL))
      for(n in 1:nrep){
        fit_glm <- NA
        while(!any(class(fit_glm) == glmclass)){ # Skip perturbations for which there are singularity issues
          Z$Response <- Z[sample.int(nrow(Z)), distance] # Perturb the response column
          fit_glm <- try(glmfunction(formula = as.formula(glmformula), family = binomial(logit),
                                     data = Z, na.action = na.exclude), silent = TRUE)
        }
        betas <- fit_glm$coefficients  
        perm[names(betas),n] <- betas
      }
      
      # P-value
      results_store[[paste(year, distance)]] <- mcarlo_rtest(obs = obs, sim = perm)
    }
  }
  return(results_store)
}




