# Needs optimising (super slow)
# Missing data are multiply imputed under assumption of ignorability

coil_extension <- function(y_matrix, psi, phi, alpha, beta, 
                           no_iterations, no_chains, cores_max){
  
  # Parallel implementation set up ---------------------------------------
  # Packages required to run code in parallel
  installed_packages <- rownames(installed.packages())
  desired_packages <- c('foreach','doMC','rngtools')
  uninstalled_desired_packages <- desired_packages[!desired_packages %in% installed_packages]
  if(length(uninstalled_desired_packages) > 0){
    install.packages(uninstalled_desired_packages)
  }
  # Load packages if not already
  loaded_packages <- search()
  unloaded_desired_packages <- desired_packages[!paste('package:', desired_packages, sep = '') %in% loaded_packages]
  if(length(unloaded_desired_packages) > 0){
    for(i in 1:length(unloaded_desired_packages)){
      require(unloaded_desired_packages[i], character.only = TRUE)
    }
  }
  
  
  # Allocate functions --------------------------------------------------
  # Function for mapping between real line and (0,1)
  logit <- function(x){
    z <- log(x/(1 - x))
  }
  invlogit <- function(x){
    z <-1/(1 + exp(-x))
    zmin <- pmax(z,.Machine$double.eps)
    zmax <- pmin(zmin,1-.Machine$double.eps)
    return(zmax)
  }
  
  # Function for drawing missing data
  missing_ind <- is.na(y_matrix)
  if(any(missing_ind)){
    rdraw_likelihood <- function(missing_ind, m, f){
      # Each i and j has its own prob
      prob1 <- exp(m %*% log(f))
      prob0 <- exp(m %*% log(1-f))
      prob0.5 <- 1-prob1-prob0
      probs <- cbind(prob1[missing_ind], prob0[missing_ind], prob0.5[missing_ind]) 
      y_imputed <- apply(probs, 1, FUN = function(x){
        sample(c(0,1,0.5), prob = pmax(0,x), size = 1) # pmax prevents underflow problems
      })
      return(y_imputed)
    }} else {
      rdraw_likelihood <- function(missing_ind, m, f){
      }
    }
  
  # Function to calculate log likelihood
  log_likelihood <- function(y, m, f){
    x <- array(NA, dim = dim(y))
    A <- m %*% log(f)
    B <- m %*% log(1-f)
    maxAB <- pmax(A,B)
    C <- log(pmax(0, 1 - exp(A) - exp(B)))
    x[y == 1] <- A[y == 1]
    x[y == 0] <- B[y == 0]
    x[y == 0.5] <- C[y == 0.5]

    #x[] <- 0 # Likelihood check
    return(x)
  }
  
  # MCMC function 
  # Inner function defined in order to automate mysd proposal parameter specification 
  coil_extension_inner <- function(y_matrix, psi, phi, alpha, beta, mysd, no_iterations){
    
    # Allocate memory for stores
    T <- no_iterations # No. of iterations
    n <- nrow(y_matrix)
    p <- ncol(y_matrix)
    m_vector_store <- matrix(NA, nrow = T, ncol = n)
    f_vector_store <- matrix(NA, nrow = T, ncol = p)
    lambda_store <- matrix(NA, nrow = T, ncol = 1)
    Acceptance_store_f <- matrix(NA, nrow = T, ncol = p)
    Acceptance_store_m <- matrix(NA, nrow = T, ncol = n)
    
    # Initial parameter values
    lambda_store[1,] <- rbeta(1, shape1 = psi, shape2 = phi)
    f_vector_store[1,] <- rbeta(p, shape1 = alpha, shape2 = beta)
    m_vector_temp <- rgeom(n, lambda_store[1,])
    while(any(m_vector_temp == 0)){m_vector_temp[m_vector_temp == 0] <- rgeom(sum(m_vector_temp == 0), lambda_store[1,])}
    m_vector_store[1,] <- m_vector_temp
    y_matrix[missing_ind] <- rdraw_likelihood(missing_ind, m = t(m_vector_store[1,,drop = FALSE]), 
                                              f = f_vector_store[1,,drop = FALSE])
    
    # Allocate space ahead of time 
    m_star_vector <- matrix(NA, nrow = n, ncol = 1)
    log_likelihood_m_star <- rep(NA, n) 
    log_likelihood_f_star <- rep(NA, p)
    
    # Allocate prior prob ahead of time (and overwrite if accepted) to 
    log_prior_m_current <- dgeom(m_vector_store[1,], lambda_store[1,], log = TRUE)
    log_likelihood_m_current <- rowSums(log_likelihood(y = y_matrix, 
                                                       m = as.matrix(m_vector_store[1,], ncol = 1), 
                                                       f = f_vector_store[1,,drop = FALSE]))
    log_prior_lambda_current <- dbeta(lambda_store[1,], shape1 = psi, shape2 = phi, log = TRUE)
    log_likelihood_f_current <- colSums(log_likelihood(y = y_matrix, 
                                                       m = as.matrix(m_vector_store[1,], ncol = 1), 
                                                       f = f_vector_store[1,,drop = FALSE]))
    log_prior_f_current <- dbeta(f_vector_store[1,], shape1 = alpha, shape2 = beta, log = TRUE)
    
    #pb <- txtProgressBar(min = 0, max = T, style = 3)
    for(t in 2:T){
      
      #setTxtProgressBar(pb, t)
      #------------------------------------------------------------------
      # 1) Update the MOIs (random walk on [1,N+) - different from notes!)
      m_vector_store[t,] <- m_vector_store[t-1,] # Copy then overwrite after MH step
      
      # Propose
      ind_one <- m_vector_store[t-1,] == 1
      ind_plus_one <- !ind_one
      m_star_vector[ind_one,] <- m_vector_store[t-1,ind_one] + 1
      m_star_vector[ind_plus_one,] <- m_vector_store[t-1,ind_plus_one] + sample(c(-1,1), 
                                                                                size = sum(ind_plus_one),
                                                                                replace = TRUE)
      
      # Calculate log densities
      # log likelihood ratio
      log_likelihood_m_star <- rowSums(log_likelihood(y = y_matrix, 
                                                      m = m_star_vector, 
                                                      f = f_vector_store[t-1,,drop = FALSE]))
      log_likelihood_m_ratio <- log_likelihood_m_star - log_likelihood_m_current
      
      # log prior ratio (normalising constant cancels since lambda^(t-1) constant in mi step)
      log_prior_m_star <- dgeom(m_star_vector, lambda_store[t-1,], log = TRUE)
      log_prior_m_ratio <-  log_prior_m_star - log_prior_m_current
      
      # log proposal ratio
      log_proposal_m_star_given_current <- rep(log(0.5), n)
      log_proposal_m_star_given_current[ind_one] <- log(1)
      log_proposal_m_current_given_star <- rep(log(0.5), n)
      log_proposal_m_current_given_star[m_star_vector == 1] <- log(1)
      log_proposal_m_ratio <- log_proposal_m_current_given_star - log_proposal_m_star_given_current
      
      # MH step
      log_MH_acceptance_ratio <- log_likelihood_m_ratio + log_prior_m_ratio + log_proposal_m_ratio
      Accept_or_not <- log(runif(n,0,1)) <= log_MH_acceptance_ratio
      m_vector_store[t,Accept_or_not] <- m_star_vector[Accept_or_not]
      Acceptance_store_m[t, ] <- Accept_or_not
      
      # Remember to overwrite online stores (update log_prior_m_current in next step)
      log_likelihood_m_current[Accept_or_not] <- log_likelihood_m_star[Accept_or_not]
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      # 2) Update lambda (exactly gibbs sample)
      lambda_store[t,] <- rbeta(1, shape1 = psi + n, shape2 = phi + sum(m_vector_store[t,]) - n) 
      log_prior_m_current <- dgeom(m_vector_store[t,], lambda_store[t,], log = TRUE)  # Replace all
      #------------------------------------------------------------------
      
      #------------------------------------------------------------------
      # 3) Update f 
      # Note that the full conditional for fj can be expressed as a sum of betas for latent variables
      # aj, but there is no standard distribution for a sum of beta's (note that the support of a beta is (0,1) and the
      # support for a sum of betas not (0, 1)
      f_vector_store[t,] <- f_vector_store[t-1,] # Copy then overwrite after MH step
      
      # Propose 
      f_vector_current <- f_vector_store[t-1,]
      theta_vector_current <- logit(f_vector_current)
      theta_vector_star <- rnorm(p, theta_vector_current, mysd) # iid given theta_vector_current
      f_vector_star <- invlogit(theta_vector_star)
      
      # Calculate densitites
      # log likelihood
      log_likelihood_f_star <- colSums(log_likelihood(y = y_matrix, 
                                                      m = t(m_vector_store[t,,drop = FALSE]), 
                                                      f = matrix(f_vector_star, nrow = 1)))
      log_likelihood_f_ratio <- log_likelihood_f_star - log_likelihood_f_current
      
      # Prior 
      log_prior_f_star <- dbeta(f_vector_star, shape1 = alpha, shape2 = beta, log = TRUE)
      log_prior_f_ratio <- log_prior_f_star - log_prior_f_current
      
      # Proposal (normals cancel)
      log_proposal_f_ratio <- 2 * log(exp(theta_vector_current) + 1) - theta_vector_current + theta_vector_star - 2 * log(exp(theta_vector_star) + 1)
      
      # MH step
      log_MH_acceptance_ratio <- log_likelihood_f_ratio + log_prior_f_ratio + log_proposal_f_ratio
      Accept_or_not <- log(runif(p,0,1)) <= log_MH_acceptance_ratio
      f_vector_store[t, Accept_or_not] <- f_vector_star[Accept_or_not]
      Acceptance_store_f[t, ] <- Accept_or_not
      
      # Overwrite online stores 
      log_likelihood_f_current[Accept_or_not] <-log_likelihood_f_star[Accept_or_not]
      log_prior_f_current[Accept_or_not] <- log_prior_f_star[Accept_or_not]
      
      
      #------------------------------------------------------------------
      # 4) Update missing 
      y_matrix[missing_ind] <- rdraw_likelihood(missing_ind, m = t(m_vector_store[t,,drop = FALSE]), 
                                                f = f_vector_store[t,,drop = FALSE])
    }
    
    results <- list(m_vector_store = m_vector_store, 
                    f_vector_store = f_vector_store,
                    lambda_store = lambda_store, 
                    Acceptance_store_f = Acceptance_store_f, 
                    Acceptance_store_m = Acceptance_store_m)
    
    return(results)
  }

  # Choose sd proposal parameter ------------------------------
  # Start with sd; decrease if acceptance too low; increase if too high  
  mysd <- 0.36
  no_iterations_sd <- 1000
  burnin <- 1:(0.5*no_iterations_sd)
  results <- coil_extension_inner(y_matrix, psi, phi, alpha, beta, mysd, no_iterations_sd)
  X <- apply(results$Acceptance_store_f[-burnin,], 2, cumsum) 
  mean_f_accep <- mean(tail(X,1)/(0.5*no_iterations_sd))
  sd_f_accep <- sd(tail(X,1)/(0.5*no_iterations_sd))
  print(c(mean_f_accep, mysd, sd_f_accep))
  
  f_accep_range <- c(0.1, 0.3) 
  while(mean_f_accep < min(f_accep_range) | mean_f_accep > max(f_accep_range)){
    if(mean_f_accep < min(f_accep_range)){mysd <- mysd * 0.75}
    if(mean_f_accep > max(f_accep_range)){mysd <- mysd * 1.25}
    results <- coil_extension_inner(y_matrix, psi, phi, alpha, beta, mysd, no_iterations_sd)
    X <- apply(results$Acceptance_store_f[-burnin,], 2, cumsum) 
    mean_f_accep <- mean(tail(X,1)/(0.5*no_iterations_sd))
    sd_f_accep <- sd(tail(X,1)/(0.5*no_iterations_sd))
    print(c(mean_f_accep, mysd, sd_f_accep))
  }
  
  # Run model in parallel ------------------------------
  registerDoMC(cores = min(no_chains, cores_max)) # Register number of cores
  rng <- RNGseq(no_chains, 1234352) # Pre-specify seed per chain
  
  # Initialise chains
  parallel_return <- foreach(chain = 1:no_chains) %dopar% {
    rngtools::setRNG(rng[[chain]]) # Launch seed per chain
    results <- coil_extension_inner(y_matrix, psi, phi, alpha, beta, mysd, no_iterations)
  }
  
  return(parallel_return)
}






