# -----------------------------------------------------------------------------#
# Preliminaries ################################################################
#------------------------------------------------------------------------------#

# Please manually set the processed data folder to be the working directory.

# Start clean:
rm(list=ls())

# Decide whether all the output messages should be shown:
showOutput4All <- FALSE

# Load the package for the beta regressions:
library("glmmTMB")
library("performance")

# Load the data:
path2BlinkPausesFile   <- './07_listener_blink_pause_proportions/listener_blink_pause_proportions.csv'
path2SynchedBlinksFile <- './09_synched_blink_proportions/synched_blink_proportions.csv'
blinkPausesData        <- read.csv(path2BlinkPausesFile, header = TRUE)
synchedBlinksData      <- read.csv(path2SynchedBlinksFile, header = TRUE)

# Determine where the results should be saved:
pathoutBetaResults     <- './beta_regression_results'

if (!dir.exists(pathoutBetaResults)){
  dir.create(pathoutBetaResults)
}

# Convert the conditions to factors so they can be used for regression coding:
blinkPausesData$Condition    <- as.factor(blinkPausesData$Condition)
synchedBlinksData$Condition  <- as.factor(synchedBlinksData$Condition)



# Create a loop to calculate the blink-pause statistics ------------------------

referenceConditions <- c("V","A","AV-nolips")
effectConditions    <- c("AV","AV","AV")
blinkPauseResults   <- data.frame()

for (i in 1:length(referenceConditions)) {
  
  # Subset the relevant conditions:
  conditionSubset     <- blinkPausesData$Condition == referenceConditions[i] | 
                         blinkPausesData$Condition == effectConditions[i]
  datSubset           <- blinkPausesData[conditionSubset, ]
  
  # Set the reference:
  datSubset$Condition <- relevel(datSubset$Condition, ref = referenceConditions[i])
  
  # Fit the model:
  betaRegModel        <- glmmTMB(Proportion ~ Condition  + (1 | Subject), 
                                 data = datSubset, 
                                 family = beta_family(link = "logit") )
  regModelSummary     <- summary(betaRegModel)
  
  # Adjust the p-value to be one sided
  slope               <- regModelSummary[["coefficients"]][["cond"]][2,1]
  if(slope > 0) {
    effectPVal        <- regModelSummary[["coefficients"]][["cond"]][2,4] / 2
  } else {
    effectPVal        <- 1 - regModelSummary[["coefficients"]][["cond"]][2,4] / 2
  }
  
  # Estimate the effect size based on Lorah (2018):
  fixedAndRandVars    <- r2(betaRegModel) # Uses Nakagawa's method
  fixedExplainedVar   <- as.numeric(fixedAndRandVars$R2_marginal)
  effectSizeEst       <- fixedExplainedVar / (1-fixedExplainedVar)
  
  # Extract other useful information:
  intercept           <- regModelSummary[["coefficients"]][["cond"]][1,1]
  totalExplainedVar   <- as.numeric(fixedAndRandVars$R2_conditional)
  intraClassCorrs     <- icc(betaRegModel)
  adjustedIcc         <- intraClassCorrs$ICC_adjusted
  
  # Display the results:
  if (showOutput4All){
    print(paste("Comparison between ", referenceConditions[i], " and ", effectConditions[i]))
    print(regModelSummary)
    print(paste("One sided p-value of the slope: ", effectPVal))
    print(fixedAndRandVars)
    print(intraClassCorrs)
  }
  
  
  # Store the results:
  currResults <- data.frame(Reference     = referenceConditions[i], 
                            Effect        = effectConditions[i],
                            Intercept     = intercept,
                            B1            = slope,
                            pVal_OneSided = effectPVal,
                            f2            = effectSizeEst,
                            ICC           = adjustedIcc,
                            R2Fixed       = fixedExplainedVar,
                            R2Total       = totalExplainedVar)
  blinkPauseResults <- rbind(blinkPauseResults,currResults)
  
} # End of the for loop

blinkPauseResults[,-(1:2)] <- signif(blinkPauseResults[,-(1:2)], digits = 3)
print("Results for the blink-pauses:")
print(blinkPauseResults)


# Create a loop to calculate the synched-blink statistics: ---------------------

referenceConditions  <- c("V","A")
effectConditions     <- c("AV","AV")
synchedBlinksResults <- data.frame()

for (i in 1:length(referenceConditions)) {
  # Subset the relevant conditions:
  conditionSubset     <- synchedBlinksData$Condition == referenceConditions[i] | 
                         synchedBlinksData$Condition == effectConditions[i]
  datSubset           <- synchedBlinksData[conditionSubset, ]
  
  # Set the reference:
  datSubset$Condition <- relevel(datSubset$Condition, ref = referenceConditions[i])
  
  # Fit the model:
  betaRegModel        <- glmmTMB(Proportion ~ Condition  + (1 | Subject), 
                                 data = datSubset, 
                                 family = beta_family(link = "logit") )
  regModelSummary     <- summary(betaRegModel)
  
  # Adjust the p-value to be one sided
  slope               <- regModelSummary[["coefficients"]][["cond"]][2,1]
  if(slope > 0) {
    effectPVal        <- regModelSummary[["coefficients"]][["cond"]][2,4] / 2
  } else {
    effectPVal        <- 1 - regModelSummary[["coefficients"]][["cond"]][2,4] / 2
  }
  
  # Estimate the effect size based on Lorah (2018):
  fixedAndRandVars    <- r2(betaRegModel)
  fixedExplainedVar   <- as.numeric(fixedAndRandVars$R2_marginal)
  effectSizeEst       <- fixedExplainedVar / (1-fixedExplainedVar)
  
  # Extract other useful information:
  intercept           <- regModelSummary[["coefficients"]][["cond"]][1,1]
  totalExplainedVar   <- as.numeric(fixedAndRandVars$R2_conditional)
  intraClassCorrs     <- icc(betaRegModel)
  adjustedIcc         <- intraClassCorrs$ICC_adjusted
  
  # Display the results:
  if (showOutput4All){
    print(paste("Comparison between ", referenceConditions[i], " and ", effectConditions[i]))
    print(regModelSummary)
    print(paste("One sided p-value of the slope: ", effectPVal))
    print(fixedAndRandVars)
    print(intraClassCorrs)
  }
  
  
  # Store the results:
  currResults <- data.frame(Reference     = referenceConditions[i], 
                            Effect        = effectConditions[i],
                            Intercept     = intercept,
                            B1            = slope,
                            pVal_OneSided = effectPVal,
                            f2            = effectSizeEst,
                            ICC           = adjustedIcc,
                            R2Fixed       = fixedExplainedVar,
                            R2Total       = totalExplainedVar)
  synchedBlinksResults <- rbind(synchedBlinksResults,currResults)
  
} # End of the for loop

synchedBlinksResults[,-(1:2)] <- signif(synchedBlinksResults[,-(1:2)], digits = 3)
print("Results for the synched-blinks")
print(synchedBlinksResults, digits = 3)



# Save the results -------------------------------------------------------------

# Save the blink-pause results:
write.csv(blinkPauseResults, 
          paste(pathoutBetaResults,'/blink_pause_results.csv', sep = ""), 
          col.names = TRUE)


# Save the synched-blink results:
write.csv(synchedBlinksResults, 
          paste(pathoutBetaResults, '/synched_blink_results.csv', sep = ""),
          col.names = TRUE)


