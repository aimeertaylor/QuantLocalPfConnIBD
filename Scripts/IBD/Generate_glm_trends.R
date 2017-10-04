#####################################################################################
# Script to generate logistic regression results from WGS, Barcode and Barcode data 
# with no intra-clinic clones 
# Takes ~ 4 hr to run to completion on the server
#####################################################################################
rm(list = ls())
nrep <- 1000
air <- FALSE
if(air){path <- '../../RData/'}else{path <- './'} # Change to '~/TM_border' if server
if(air){source('../../FunctionFiles/simtests.R')}else{source(sprintf('%ssimtests.R', path))}
if(air){source('../../FunctionFiles/glm_trends.R')}else{source(sprintf('%sglm_trends.R', path))}
if(air){load('../../RData/formulas.RData')}else{load(sprintf('%sformulas.RData', path))}
if(air){load('../../RData/WGS_threshold.RData')}else{load(sprintf('%sWGS_threshold.RData', path))}
if(air){load('../../RData/Barcode_threshold.RData')}else{load(sprintf('%sBarcode_threshold.RData', path))}
if(air){load('../../RData/Barcode_nowithinclinicclones_threshold.RData')}else{load(sprintf('%sBarcode_nowithinclinicclones_threshold.RData', path))}


#====================================================================================
# WGS 
# 468 sec, nrep = 1000 and speed = TRUE on server
#====================================================================================
system.time( 
  for(glmformula in names(formulas)){
    glm_WGS_results <- glm_trends(X = WGS, 
                                  distances = "ProbIBD_tail",
                                  glmformula = formulas[[glmformula]],
                                  nrep = nrep, 
                                  years = 14, 
                                  response_specified = FALSE, 
                                  speed = TRUE)
    save(glm_WGS_results, file = sprintf('%s/GLM_WGS_%s.RData', path, glmformula))
  }
)


#============================================================================
# Barcode
# 7533.051 sec for nrep = 1000 with speed = TRUE on server
#============================================================================
system.time (
  for(glmformula in names(formulas)[1:2]){
    glm_barcode_results <- glm_trends(X = Barcode,
                                      distances = 'ProbIBD_93_tail',
                                      glmformula = formulas[[glmformula]],
                                      nrep = nrep,
                                      years = 8:10,
                                      response_specified = FALSE,
                                      speed = TRUE)
    save(glm_barcode_results, file = sprintf('%sRData/GLM_barcode_%s.RData', path, glmformula))
  }
)


#============================================================================
# Barcode no within-clinic clones
# 5172.417 sec with nrep = 1000 and speed = TRUE on server
#============================================================================
system.time(
  for(glmformula in names(formulas)[1:2]){
    glm_barcode_results <- glm_trends(X = Barcode_nowithinclinicclones,
                                      distance = 'ProbIBD_93_tail',
                                      formulas[[glmformula]],
                                      nrep = nrep,
                                      years = 8:10,
                                      response_specified = FALSE,
                                      speed = TRUE)
    save(glm_barcode_results, file = sprintf('%sRData/GLM_barcode_nowithinclinicclones_%s.RData', path, glmformula))
  }
)
