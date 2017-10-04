# Diagnostic check of COIL extension run in parallel with multiple imputation 
# and automation of sd parameter 

# Summary (Oct 19th 2016): 
# Likelihood = 1: COILext_likelihood_one
# not optimal mixing, 
# f posteriors approaching prior 
# lambda posterior not converged
# COIs converged and mixing
# freq acceptance ~ 0.5 for all
# COI acceptance ~ 0.7 for all

# All data combined:
# frequencies converged 
# lambda converged
# COIs not mixing, most likely due to concentrated posterior
# freq acceptance >0 <0.2
# COI acceptance: 417 have zero acceptance.

# site_year:
# frequencies converged: albeit with more uncertainty due to less data 
# lambda converged: albeit with more uncertainty due to less data 
# COIs mostly not mixing, most likely due to concentrated posterior
# Didn't look at all acceptance

# site:
# frequencies converged
# lambda converged
# COIs mostly not mixing, most likely due to concentrated posterior
# Didn't look at all acceptance

rm(list = ls())
filename <- 'COILext_results_site' # 'COILext_results_site_year' 
load(sprintf('../../RData/%s.RData', filename))
require(MCMCpack)

for(dataset in names(results_store)){

  # Frequency results
  f_results <- mcmc.list(as.mcmc(results_store[dataset][[1]][[1]]$f_vector_store),
                         as.mcmc(results_store[dataset][[1]][[2]]$f_vector_store),
                         as.mcmc(results_store[dataset][[1]][[3]]$f_vector_store))
  plot(f_results[,sample(1:93, 8)])

  # Lambda results
  lambda_results <- mcmc.list(as.mcmc(results_store[dataset][[1]][[1]]$lambda_store),
                              as.mcmc(results_store[dataset][[1]][[2]]$lambda_store),
                              as.mcmc(results_store[dataset][[1]][[3]]$lambda_store))
  plot(lambda_results)

  # COI results (plot first 4)
  m_results <- mcmc.list(as.mcmc(results_store[dataset][[1]][[1]]$m_vector_store),
                         as.mcmc(results_store[dataset][[1]][[2]]$m_vector_store),
                         as.mcmc(results_store[dataset][[1]][[3]]$m_vector_store))
  plot(m_results[,1:4])
}



# All data combined: COILext_results_all ------------------------- 
rm(list = ls())
filename <- 'COILext_results_all'
load(sprintf('../../RData/%s.RData', filename))
require(MCMCpack)

# Frequency results
f_results <- mcmc.list(as.mcmc(results_store[[1]]$f_vector_store),
                       as.mcmc(results_store[[2]]$f_vector_store),
                       as.mcmc(results_store[[3]]$f_vector_store))
plot(f_results[,sample(1:93, 8)])

# Lambda results
lambda_results <- mcmc.list(as.mcmc(results_store[[1]]$lambda_store),
                            as.mcmc(results_store[[2]]$lambda_store),
                            as.mcmc(results_store[[3]]$lambda_store))
plot(lambda_results)

# COI results (plot random subset of 4/1731)
m_results <- mcmc.list(as.mcmc(results_store[[1]]$m_vector_store),
                       as.mcmc(results_store[[2]]$m_vector_store),
                       as.mcmc(results_store[[3]]$m_vector_store))
plot(m_results[,sample(1:1731, 4)])

# Frequency acceptance
X <- apply(results_store[[1]]$Acceptance_store_f[-1,], 2, cumsum)
no_it <- nrow(results_store[[1]]$Acceptance_store_f[-1,])
Z <- X/matrix(1:no_it, nrow = no_it, ncol = 93)
par(mfrow = c(1,1))
plot(Z[,1], type = 'l', ylim = c(0,1))
for(i in 2:93){
  lines(Z[,i], col = i)
}

# COI acceptance (random sample of 100/1731)
X <- apply(results_store[[1]]$Acceptance_store_m[-1,], 2, cumsum)
no_it <- nrow(results_store[[1]]$Acceptance_store_m[-1,])
Z <- X/matrix(1:no_it, nrow = no_it, ncol = 1731)
par(mfrow = c(1,1))
plot(Z[,1], type = 'l', ylim = c(0,0.01))
for(i in sample(2:1731, 100)){
  lines(Z[,i], col = i)
}

# Final acceptance rates
acceptance_at_end <- tail(Z, 1)
barplot(table(acceptance_at_end), las = 2, 
        cex.names = 0.5)
