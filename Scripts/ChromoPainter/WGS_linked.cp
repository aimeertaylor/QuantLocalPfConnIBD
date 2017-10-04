section:fsproject
# This file contains the information about the project. It has the following sections:
############## stage0 (data processing to create the required input files):
############## stage1 (chromopainter parameter estimation via EM)
############## stage2 (chromopainter full painting)
############## stage3 (finestructure mcmc)
############## stage4 (finestructure tree)
#
# You may modify this file by hand, but it is possible to break it if you are careless.  It is advisable instead to set parameters via the "fs" interface; either from the command line or by placing the set of parameters to set in a settings file and using -import.
# Each section is separated by "section:<sectionname>".  Within each section the parameters are defined in the format:
# paramname:value
# A '#' symbol indicates the end of the readable line, everything afterwards is treated as a comment.
# Everything from the colon to the end of the line (or first #) is treated as content. Parameters that are empty are assigned their default values as described in "fs -h".

############## Command history:
section:history
CMDLINE: -phasefiles ./WGS.phase -recombfiles ./WGS.recombfile -idfile ./WGS.ids -ploidy 1 -dos2 
CMD:-phasefiles ./WGS.phase
CMD:-recombfiles ./WGS.recombfile
CMD:-idfile ./WGS.ids
CMD:-ploidy 1
CMD:-countdata
CMD:-makes1
CMD:-dos1
CMD:-combines1
CMD:-makes2
CMD:-dos2

section:parameters

###################
stage:0
## Data preparation and conversion. Not yet implemented

validatedoutput:1,1,0,0,0  # Derived. Whether we have validated output from each stage of the analysis (0-4)

###################
stage:1
## Universal Stage1-4 properties of the analysis

exec:fs  # Finestructure command line. Set this to be able to use a specific version of this software. (default: fs)
hpc:0  # THIS IS IMPORTANT FOR BIG DATASETS! Set hpc mode. 0: Commands are run 'inline' (see 'numthreads' to control how many CPU's to use). 1: Stop computation for an external batch process, creating a file containing commands to generate the results of each stage. 2: Run commands inline, but create the commands for reference. (default: 0.)
numthreads:0  # Maximum parallel threads in 'hpc=0' mode. Default: 0, meaning all available CPUs.
ploidy:1  # Haplotypes per individual. =1 if haploid, 2 if diploid. (default: 2)
linkagemode:linked  # unlinked/linked. Whether we use the linked model. default: unlinked / linked if recombination files provided.
indsperproc:0  # Desired number of individuals per process (default: 0, meaning autocalculate: use 1 In HPC mode, ceiling(N/numthreads) otherwise. Try to choose it such that you get a sensible number of commands compared to the number of cores you have available.
outputlogfiles:1  # 1=Commands are written to file with redirection to log files. 0: no redirection. (default:1)
allowdep:1  # Whether dependency resolution is allowed. 0=no, 1=yes. Main use is for pipelining. (default:1).
## ChromoPainter Stage1-2 generic properties

s12inputtype:phase  # What type of data input (currently only "phase" supported)
idfile:./WGS.ids  # IDfile location, containing the labels of each individual. REQUIRED, no default (unless -createids is used).
s12args:  # arguments to be passed to Chromopainter (default: empty)
## Quantities observed from data

ninds:178  # Derived. number of individuals observed in the idfile
nindsUsed:178  # Derived. number of individuals retained for processing from the idfile
nsnps:34911  # Derived. number of SNPs in total, over all files
## ChromoPainter Stage1 (EM) properties

s1args:-in -iM --emfilesonly  # Arguments passed to stage1 (default:-in -iM --emfilesonly)
s1emits:10  # Number of EM iterations (chromopainter -i <n>, default: 10)
s1minsnps:10000  # Minimum number of SNPs for EM estimation (for chromopainter -e, default: 10000)
s1snpfrac:0.1  # fraction of genome to use for EM estimation. (default: 0.1)
s1indfrac:1  # fraction of individuals to use for EM estimation. (default: 1.0)
s1outputroot:./WGS_linked_tmp_EM_linked_haploid  # output file for stage 1 (default is autoconstructed from filename)
## ChromoPainter Stage2 properties inferred from Stage1

Neinf:1904.43  # Derived. Inferred `Effective population size Ne' (chromopainter -n).
muinf:0.00393868  # Derived. Inferred Mutation rate mu (chromopainter -M)

###################
stage:2
## ChromoPainter Stage2 (main run) properties

s2chunksperregion:-1  # number of chunks in a "region" (-ve: use default of 100 for linked, nsnps/100 for unlinked)
s2samples:0  # number of samples of the painting to obtain per recipient haplotype, for examining the details of the painting. (Populates <root>.samples.out; default 0. Warning: these file can get large)
s2args:  # Additional arguments for stage 2 (default: none, "")
s2outputroot:./WGS_linked_tmp_mainrun.linked_haploid  # Output file name for stage 2 (default: autoconstructed).
s2combineargs:  # Additional arguments for stage 2 combine (fs combine; default: none, "")
## FineSTRUCTURE Stage3 properties inferred from Stage2

cval:-1  # Derived. 'c' as inferred using chromopainter. This is only used for sanity checking. See s34 args for setting it manually.
cproot:  # The name of the final chromopainter output. (Default: <filename>, the project file name)
cpchunkcounts:  # the finestructure input file, derived name of the chunkcounts file from cproot.

###################
stage:3
## FineSTRUCTURE Stage3-4 generic properties

fsroot:  # The name of the finestructure output (Default: <filename>, the project file name).
s34args:  # Additional arguments to both finestructure mcmc and tree steps. Add "-c <val>" to manually override 'c'.
## FineSTRUCTURE Stage3 MCMC inference

s3iters:100000  # Number of TOTAL iterations to use for MCMC. By default we assign half to burnin and half to sampling. (default: 100000)
s3iterssample:-1  # Number of iterations to use for MCMC (default: -ve, meaning derive from s3iters)
s3itersburnin:-1  # Number of iterations to use for MCMC burnin (default: -ve, meaning derive from s3iters)
numskip:-1  # Number of mcmc iterations per retained sample; (default: -ve, meaning derive from maxretained)
maxretained:500  # Maximum number of samples to retain when numskip -ve. (default: 500)
nummcmcruns:2  # Number of *independent* mcmc runs. (default: 2)
fsmcmcoutput:  # Filename to use for mcmc output (default: autogenerated)
mcmcGR:  # Derived. Gelman-Rubin diagnostics obtained from combining MCMC runs, for log-posterior, K,log-beta,delta,f respectively
threshGR:1.3  # Threshold for the Gelman-Rubin statistic to allow moving on to the tree building stage. We always move on if thresGR<0. (Default: 1.3)

###################
stage:4
## FineSTRUCTURE Stage4 tree inference

s4args:  # Extra arguments to the tree building step. (default: none, "")
s4iters:100000  # Number of maximization steps when finding the best state from which the tree is built. (default: 100000)
fstreeoutput:  # Filename to use for finestructure tree output. (default: autogenerated)

###################
stage:1
## Vector quantities placed at the end for readability

phasefiles:./WGS.phase  # Comma or space separated list of all 'phase' files containing the (phased) SNP details for each haplotype. Required. Must be sorted alphanumerically to ensure chromosomes are correctly ordered. So don't use *.phase, use file{1..22}.phase. Override this with upper case -PHASEFILES.
recombfiles:./WGS.recombfile  # Comma or space separated list of all recombination map files containing the recombination distance between SNPs. If provided, a linked analysis is performed. Otherwise an 'unlinked' analysis is performed. Note that linkage is very important for dense markers!
nsnpsvec:34911  # Derived. Comma separated list of the number of SNPs in each phase file.
s1outputrootvec:./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind1,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind2,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind3,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind4,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind5,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind6,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind7,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind8,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind9,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind10,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind11,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind12,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind13,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind14,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind15,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind16,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind17,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind18,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind19,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind20,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind21,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind22,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind23,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind24,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind25,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind26,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind27,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind28,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind29,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind30,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind31,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind32,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind33,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind34,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind35,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind36,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind37,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind38,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind39,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind40,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind41,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind42,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind43,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind44,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind45,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind46,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind47,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind48,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind49,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind50,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind51,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind52,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind53,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind54,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind55,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind56,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind57,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind58,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind59,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind60,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind61,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind62,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind63,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind64,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind65,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind66,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind67,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind68,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind69,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind70,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind71,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind72,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind73,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind74,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind75,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind76,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind77,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind78,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind79,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind80,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind81,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind82,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind83,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind84,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind85,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind86,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind87,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind88,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind89,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind90,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind91,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind92,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind93,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind94,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind95,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind96,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind97,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind98,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind99,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind100,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind101,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind102,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind103,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind104,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind105,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind106,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind107,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind108,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind109,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind110,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind111,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind112,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind113,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind114,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind115,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind116,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind117,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind118,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind119,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind120,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind121,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind122,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind123,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind124,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind125,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind126,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind127,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind128,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind129,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind130,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind131,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind132,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind133,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind134,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind135,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind136,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind137,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind138,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind139,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind140,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind141,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind142,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind143,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind144,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind145,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind146,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind147,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind148,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind149,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind150,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind151,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind152,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind153,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind154,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind155,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind156,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind157,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind158,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind159,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind160,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind161,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind162,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind163,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind164,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind165,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind166,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind167,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind168,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind169,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind170,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind171,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind172,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind173,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind174,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind175,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind176,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind177,./WGS_linked/stage1/./WGS_linked_tmp_EM_linked_haploid_file1_ind178  # Derived. Comma separated list of the stage 1 output files names.

###################
stage:2
## 

s2outputrootvec:./WGS_linked/stage2/./WGS_linked_tmp_mainrun.linked_haploid_file1_ind1-178  # Derived. Comma separated list of the stage 2 output files names.

###################
stage:3
## 

fsmcmcoutputvec:  # Derived. Comma separated list of the stage 3 output files names.
old_fsmcmcoutputvec:  # Derived. Comma separated list of the stage 3 output files names, if we need to continue a too-short MCMC run.

###################
stage:4
## 

fstreeoutputvec:  # Derived. Comma separated list of the stage 4 output files names.

section:fsprojectend
stage:2  # Derived. Don't mess with this! The internal measure of which stage of processing we've reached. Change it via -reset or -duplicate.
