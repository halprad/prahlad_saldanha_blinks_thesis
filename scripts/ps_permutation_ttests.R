# -----------------------------------------------------------------------------#
# Preliminaries ################################################################
#------------------------------------------------------------------------------#

# Please manually set the processed data folder to be the working directory.

# Start clean:
rm(list=ls())

# Load the package for the permutation ttest:
library("MKinfer")

# Load the data:
path2BlinkRateFile <- './14_listener_blink_rates/listener_blink_rates.csv'
dat <- read.csv(path2BlinkRateFile)
dat$Condition <- as.factor(dat$Condition)

# Determine where the results should be saved:
pathoutPermTtestResults     <- './permutation_ttest_results'

if (!dir.exists(pathoutPermTtestResults)){
  dir.create(pathoutPermTtestResults)
}

# Do the ttests: ---------------------------------------------------------------


avRows  <- dat$Condition == "AV"
vRows   <- dat$Condition == "V"
aRows   <- dat$Condition == "A"
nlRows  <- dat$Condition == "AV-nolips"

vResults  <- perm.t.test(dat$Proportion[avRows],dat$Proportion[vRows],  alternative = "greater", paired = TRUE, nperm = 10000)
aResults  <- perm.t.test(dat$Proportion[avRows],dat$Proportion[aRows],  alternative = "greater", paired = TRUE, nperm = 10000)
nlResults <- perm.t.test(dat$Proportion[avRows],dat$Proportion[nlRows], alternative = "greater", paired = TRUE, nperm = 10000)


# Organize the results: --------------------------------------------------------

allResults <-data.frame( 
               Comparison      = rbind("AV > V", "AV > A", "AV > AV-nolips"),
               pValue_OneSided = rbind(vResults$perm.p.value, aResults$perm.p.value, nlResults$perm.p.value),
               meanDifference  = rbind(as.numeric(vResults$perm.estimate), 
                                       as.numeric(aResults$perm.estimate),
                                       as.numeric(nlResults$perm.estimate)),
               seDifference    = rbind(vResults$perm.stderr, aResults$perm.stderr, nlResults$perm.stderr)
               )

allResults[ ,2:4] <- signif(allResults[ ,2:4], digits = 3)


# Save: ------------------------------------------------------------------------

write.csv(allResults, 
          paste(pathoutPermTtestResults,'/permutation_ttest_results.csv', sep = ""), 
          col.names = TRUE, row.names = FALSE)
            
