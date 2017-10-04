extract_major_freq <- function(x){
  y <- table(x)
  num_not_missing <- sum(x != 'N') # missing are discarded SNP-wise
  majorSNP <- names(which(y == max(y)))
  majorSNPfreq <- max(y)/num_not_missing # Assumes biallelic 
  return(list(majorSNP, majorSNPfreq))
}
