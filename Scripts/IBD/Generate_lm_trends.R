#=============================================================================
# Generate results of direct regression of proportion IBD
# Takes ~ 30 mins to run with nrep = 1000 
#=============================================================================
rm(list = ls())
nrep <- 1000
set.seed(1)
source('../../FunctionFiles/simtests.R')
load('../../RData//formulas.RData')
load('../../RData/Barcode_threshold.RData')
load('../../RData/WGS_threshold.RData')
filename <- '../../RData/Direct_IBDprop_reg_coeff_store.RData'
require(tictoc) # for timings

# Function below copied and adapted from calculate_glm_trends_allpvalues.R on Sept 13th 207 (removed years option and multiple distances)
lm_trends <- function(X, distance, myformula, nrep){
  
  # Modify data matrix
  X$Response <- X[, distance] # Set outcome of interest
  X <- X[!is.na(X$Response), ] # Remove lines with NA results to run faster
  if(length(unique(c(as.character(X$Sitei), as.character(X$Sitej)))) < 4){next()}
  
  # Fit model
  fit <- lm(formula = as.formula(myformula), data = X, na.action = na.exclude)
  obs <- fit$coefficients
  
  # Permute
  perm <- matrix(0, ncol = nrep, nrow = length(obs), dimnames = list(names(obs), NULL))
  for(n in 1:nrep){
    fit <- NA
    while(!any(class(fit) == 'lm')){ # Skip permutations that don't return lm class
      X$Response <- X[sample.int(nrow(X)), distance] # permute the response column
      fit <- try(lm(formula = as.formula(myformula), data = X, na.action = na.exclude), silent = TRUE)
    }
    betas <- fit$coefficients
    perm[names(betas),n] <- betas
  }
  
  # P-value
  results <- mcarlo_rtest(obs = obs, sim = perm)
  return(results)
}

tic()
#----------------------------------------------------------------------------
# Barcode (estimate significance here since not clashing with any results elsewhere)
#----------------------------------------------------------------------------
fit2B <- lm_trends(X=Barcode, distance = 'ProbIBD_93', myformula = formulas$unadjusted, nrep)
fit2_B <- cbind(fit2B$obs, fit2B$pvalue)
fit8B <- lm_trends(X=Barcode, distance = 'ProbIBD_93', formulas$adjusted, nrep)
fit8_B <- cbind(fit8B$obs, fit8B$pvalue)

#----------------------------------------------------------------------------
# WGS (estimate significance here since not clashing with any results elsewhere)
#----------------------------------------------------------------------------
fit2W <- lm_trends(WGS, distance = 'ProbIBD', formulas$unadjusted, nrep)
fit2_W <- cbind(fit2W$obs, fit2W$pvalue)
fit8W <- lm_trends(WGS, distance = 'ProbIBD', formulas$adjusted, nrep)
fit8_W <- cbind(fit8W$obs, fit8W$pvalue)
fit11W <- lm_trends(WGS, distance = 'ProbIBD', formulas$adjustedyear, nrep)
fit11_W <- cbind(fit11W$obs, fit11W$pvalue)

#----------------------------------------------------------------------------
# Create tables of coefficients
#----------------------------------------------------------------------------
coeff_store <- array(dim = c(11, 5), dimnames = list(rownames(fit11_W), c('Barcode unadjusted', 'Barcode adjusted',
                                                                          'WGS unadjusted', 'WGS adjusted', 'WGS adjusted year')))

# Function to round and format trailing zeros
fmt5 <- function(x){formatC(round(x, 5), format='f', digits=5)}
fmt3 <- function(x){formatC(round(x, 3), format='f', digits=3)}

coeff_store[rownames(fit2_B),1] <- apply(fit2_B, 1, function(x){paste(fmt5(x[1]), ' (', fmt3(x[2]), ')' , sep = '')})
coeff_store[rownames(fit8_B),2] <- apply(fit8_B, 1, function(x){paste(fmt5(x[1]), ' (', fmt3(x[2]), ')' , sep = '')})
coeff_store[rownames(fit2_W),3] <- apply(fit2_W, 1, function(x){paste(fmt5(x[1]), ' (', fmt3(x[2]), ')' , sep = '')})
coeff_store[rownames(fit8_W),4] <- apply(fit8_W, 1, function(x){paste(fmt5(x[1]), ' (', fmt3(x[2]), ')' , sep = '')})
coeff_store[rownames(fit11_W),5] <- apply(fit11_W, 1, function(x){paste(fmt5(x[1]), ' (', fmt3(x[2]), ')' , sep = '')})

#save(coeff_store, file = filename)

toc()
