##################################################################################################
# Script to run ChromoPainter on Barcode data
# Takes 31 secs
##################################################################################################
rm(list = ls())

# load data in IBD file format
BarcodeData <- read.delim('../../TxtData/Barcode93.txt', sep = '\t')
load('../../RData/Data_store_Barcode.RData')
MetaData <- Data_store$MetaData

# Remove 'chrom' and 'pos' columns from sequences and transpose s.t. Sample.ID per row
SNPData <- t(BarcodeData[,-(1:2)]) 
SNPData[SNPData == -1] <- NA

# Name metadata by Sample.ID so can order by SNPdata 
rownames(MetaData) <- as.character(MetaData$Sample.ID) 
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
idfile <- cbind(MetaData[, c('Sample.ID', 'Location.code')], 1)
write.table(idfile[1:numSamples_capped,], file = './Barcode.ids', sep = " ", quote = FALSE, 
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
                   paste('P', paste(BarcodeData$pos[1:numSNPs_capped], collapse = ' ')), # Third row contains "P" followd by bp pos of each SNP
                   haps[1:numSamples_capped,, drop = FALSE]) # Each additional line contains a haplotype. NO SPACES. NO MISSING VALUES. 
write.table(phasefile, file = './Barcode.phase', sep = " ", quote = FALSE, 
            row.names = FALSE, col.names  = FALSE)

# <recombfiles>
# Space separated recombfile with header file, pos, and distance in M/bp
default <- getOption("scipen") # Get default scipen so can temporarily surpress scientific notation
options(scipen = 999)
recomrateperbp <- rep(7.4e-7, numSNPs)

# Put '-9' at last bp position of preceding chrom  
chrom <- BarcodeData$chrom[-numSNPs]
chrom_plus1 <- BarcodeData$chrom[-1]
recomrateperbp[(chrom - chrom_plus1) == -1] <- -9

# Create and write file
recombfile <- cbind(start.pos = BarcodeData$pos, 
                    recom.rate.perbp = recomrateperbp) # M/bp [Miles et al, Genome Res 26:1288-1299 (2016)]
write.table(recombfile[1:numSNPs_capped, ], file = './Barcode.recombfile', sep = " ", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)

# Restore scipen
options(scipen = default)


# ====================================================================================
# Run fs with a recomb file
# WARNING: No regions found, insufficient data for calculating c.  
# Try running chromopainter either with a smaller "-k" option, or run chromocombine with the "-C" option. 
# See http://www.paintmychromosomes.com (faq page) for a discussion of this issue.
#
# Run chromopainter alone by specifying -dos2. Takes ~ 31 secs. 
# ====================================================================================
system.time(
  system('../../../../fs-2.1.1/fs ./Barcode_unlinked.cp -n -phasefiles ./Barcode.phase -idfile ./Barcode.ids -ploidy 1 -dos2')
)

# Can also run this way (need to specify the -o flag)
# system('../../../../fs-2.1.1/fs cp -t ./Barcode.ids -g ./Barcode.phase -u -j -a 0 0 -o ./Test')


