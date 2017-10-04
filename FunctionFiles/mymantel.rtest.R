##############################################################################################
# No-longer used as generates non-uniform p-values (see S1 Appendix)
# Functions to calculate significance of trend estimate under simple model y ~ b0 + b1x + eps
# and simple correlation. Based on mantel.rtest from ade4 package, removing extra checks and classes. 
# Takes a vector instead of a matrix
##############################################################################################
mymantel_lm <- function (m1, m2, nrepet = 99) 
{ 
  require(ade4) 
  
  # Make m2 into a dist object
  y <- 0.5*(1 + sqrt(8*length(m2) + 1)) # Work out size of matrix
  m <- matrix(0, nrow = y, ncol = y) # Create matrix
  m[lower.tri(m)] <- m2 # Populate matrix
  m2 <- as.dist(m)
  
  # Function to swap comlumns and rows
  permutedist <- function(m) {
    w0 <- sample.int(attr(m, "Size"))
    m <- as.matrix(m)
    return(as.dist(m[w0, w0]))
  }
  
  obs <- summary(lm(m1 ~ m2))$coefficients['m2', 'Estimate'] 
  FUN <- function(x){
    X <- summary(lm(m1 ~ unclass(permutedist(m2))))
    coeff <- X$coefficients['unclass(permutedist(m2))', 'Estimate']
    return(coeff)
  }
  
  if (nrepet == 0){return(obs)}
  perm <- matrix(0, nrow = nrepet, ncol = 1)
  perm <- apply(perm, 1, FUN)
  w <- as.rtest(obs = obs, sim = perm, call = match.call())
  return(w)
}

mymantel_cor <- function (m1, m2, nrepet = 99) 
{ 
  require(ade4) 
  
  # Make m2 into a dist object
  y <- 0.5*(1 + sqrt(8*length(m2) + 1)) # Work out size of matrix
  m <- matrix(0, nrow = y, ncol = y) # Create matrix
  m[lower.tri(m)] <- m2 # Populate matrix
  m2 <- as.dist(m)
  
  # Function to swap comlumns and rows
  permutedist <- function(m) {
    w0 <- sample.int(attr(m, "Size"))
    m <- as.matrix(m)
    return(as.dist(m[w0, w0]))
  }
  
  obs <- cor(m1, m2)
  FUN <- function(x){
    X <- cor(m1, unclass(permutedist(m2)))
    return(X)
  }
  
  if (nrepet == 0){return(obs)}
  perm <- matrix(0, nrow = nrepet, ncol = 1)
  perm <- apply(perm, 1, FUN)
  w <- as.rtest(obs = obs, sim = perm, call = match.call())
  return(w)
}



