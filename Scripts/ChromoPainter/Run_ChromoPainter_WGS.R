##################################################################################################
# Script to run ChromoPainter on WGS data
# Takes ~ 1142 sec to run linked and ~ 285 sec unlinked
##################################################################################################
rm(list = ls())
rho <- 7.4e-7 # Same as hmmIBD v2.0.0

# Load data
SequencingData <- read.delim('../../TxtData/WGS.txt')
MetaData <- read.table('../../TxtData/WGS_metadata.txt', header = TRUE)

# Remove 'chrom' and 'pos' columns from sequences and transpose s.t. strain per row
SNPData <- t(SequencingData[,-(1:2)]) 
SNPData[SNPData == -1] <- NA

# Name metadata by strain so can order by SNPdata 
rownames(MetaData) <- as.character(MetaData$strain) 
MetaData <- MetaData[rownames(SNPData), ] 

# Summarise data dimensions
numSamples <- nrow(SNPData)
numSNPs <- ncol(SNPData) # biallelic snps only 


# ====================================================================================
# Create files for fs
# For input help see https://people.maths.bris.ac.uk/~madjl/finestructure/manualse4.html#x7-90004.1
# ====================================================================================
# Restrict data for initial set up (checking program runs etc.)
numSNPs_capped <- numSNPs
numSamples_capped <- numSamples

# <idfile>
# Space separated id file: N lines, one per individual: <NAME> <POPULATION> <INCLUSION>, which are string, string, 0 or 1
idfile <- cbind(MetaData[, c('strain', 'collection_location')], 1)
write.table(idfile[1:numSamples_capped,], file = './WGS.ids', sep = " ", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# <phasefiles> 
# First: imput missing SNPs, since not supported by fs (assume indpendence)
frequencies <- colMeans(SNPData, na.rm = TRUE)
set.seed(1) # set seed for reproducability
SNPDataImputed <- SNPData
for(i in 1:numSamples){ # for each missing per row, draw from a Bernoulli with prob = allele frequency
  missing_ind <- is.na(SNPDataImputed[i,])
  n_missing <- sum(missing_ind)
  SNPDataImputed[i, missing_ind] <- rbinom(n = n_missing, size = 1, prob = frequencies[missing_ind])  
}

# Second check order of haplotypes against idfile
rownames(SNPDataImputed) == as.character(idfile[,1])

# Third collapse each samples haplotype into a character with no spaces
haps <- matrix(apply(SNPDataImputed, 1, function(x){paste(x[1:numSNPs_capped], collapse = '')}), ncol = 1)

# Create phase file
phasefile <- rbind(numSamples_capped, # First row is the number of haplotypes
                   numSNPs_capped, # Second row is number of SNPs
                   paste('P', paste(SequencingData$pos[1:numSNPs_capped], collapse = ' ')), # Third row contains "P" followd by bp pos of each SNP
                   haps[1:numSamples_capped,, drop = FALSE]) # Each additional line contains a haplotype. NO SPACES. NO MISSING VALUES. 
write.table(phasefile, file = './WGS.phase', sep = " ", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)


# <recombfiles>
# Space separated recombfile with header file, pos, and distance in M/bp
default <- getOption("scipen") # Get default scipen so can temporarily surpress scientific notation
options(scipen = 999)
recomrateperbp <- rep(rho, numSNPs)

# Put '-9' at last bp position of preceding chrom  
chrom <- SequencingData$chrom[-numSNPs]
chrom_plus1 <- SequencingData$chrom[-1]
recomrateperbp[(chrom - chrom_plus1) == -1] <- -9

# Create and write file
recombfile <- cbind(start.pos = SequencingData$pos, 
                    recom.rate.perbp = recomrateperbp) # M/bp [Miles et al, Genome Res 26:1288-1299 (2016)]
write.table(recombfile[1:numSNPs_capped, ], file = './WGS.recombfile', sep = " ", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)

# Restore scipen
options(scipen = default)


# =======================================================================================
# Run chromopainter only within fs (by specifying -dos2) 
# =======================================================================================
# Linked 1142 sec
system.time(
  system('../../../../fs-2.1.1/fs ./WGS_linked.cp -n -phasefiles ./WGS.phase -recombfiles ./WGS.recombfile -idfile ./WGS.ids -ploidy 1 -dos2')
)

# Unlinked 285.070 sec
system.time(
  system('../../../../fs-2.1.1/fs ./WGS_unlinked.cp -n -phasefiles ./WGS.phase -idfile ./WGS.ids -ploidy 1 -dos2')
)





