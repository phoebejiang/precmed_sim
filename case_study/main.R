# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: main.R
#
# Purpose: Main function to run the stratified CV of PM method implementation 
#  to a subset of MarketScan data in batches of 5 (PM methods and CV batch can 
#  be run separately)
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(haven))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(pscl))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(mpath))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(fastDummies))
suppressPackageStartupMessages(library(listdtr))
suppressPackageStartupMessages(library(DTRreg))

# homepath <- "/home/pjiang/pmms/" # CHANGE THIS #
# setwd(homepath) 
source("./utility.R")
source("./01-preprocessing.R")
source("./02-propensityscore.R")
source("./03-dWOLS.R")
source("./04-regression-based.R")
source("./05-listdtr.R")
source("./06-LuScore.R")
source("./eachCV.R")

args <- commandArgs(trailingOnly = TRUE)

# Constants
n.fold <- 10     # number of folds in each CV iteration
n.cv <- 25       # total number of CV iterations 
base.seed <- 999 # randomization seed

# User-specified constants
method <- args[1]               # PM method, an argument from command line
batch <- as.numeric(args[2])    # the batch index, could be from 1 to n.cv/batch_size, an argument from command line
batch_size <- 5

## X variables to be included in each model
categoricalvars <- NULL
continuousvars <- c("female", "white", "prmsgr", "age", "weightbl", "diagyrs", "rlps1yr", "trelmos", 
                    "edssbl", "tm25zbl", "nhptzbl", "pasatzbl", "chrt2_5bl", "sf36pcsbl", "sf36mcsbl")


## Outcome depends on the method and user
candidates <- list(allDMF = "postrelapse_num",
                   allGA = "postrelapse_num",
                   linear = "logarr0001", # default = "logarr0001", lower is better
                   weightedLinear = "logarr0001", # default = "logarr0001", lower is better
                   weightedPoisson = "postrelapse_num", 
                   negBin = "postrelapse_num",
                   weightedNegBin = "postrelapse_num",
                   dWOLS = "mlogarr0001", # default = "mlogarr0001", add minus because higher is better
                   listDTR2 = "mlogarr0001", # default = "mlogarr0001", add minus because higher is better
                   listDTR3 = "mlogarr0001", # default = "mlogarr0001", add minus because higher is better
                   poisson = "postrelapse_num",
                   boosting = "postrelapse_num",
                   twoReg = "postrelapse_num",
                   contrastReg = "postrelapse_num") 
outcome <- candidates[[method]]
cat("\nOutcome is:", outcome, "\n")


#######################################################################
############################# Read in data ############################
#######################################################################

# Read in preprocessed CONFIRM case study data
ds <- readRDS(homepath, "intermediate_results/case_study/preprocess.RDS") 

ds %<>% 
  mutate(trt = ifelse(trt01p == "BG00012 240 mg BID", 1, 0)) %>% 
  dplyr::select(-trt01p) %>% 
  mutate(logarr0001 = log((inecp / (time_inecp / 365.25)) + 0.001), # log ARR = log (number of INEC relapses / years of follow up + 0.001)
         mlogarr0001 = -logarr0001, # -log ARR = -log (number of INEC relapses / years of follow up + 0.001)
         postrelapse_num = inecp,
         FUweight = time_inecp/sum(time_inecp),
         offset = log(time_inecp / 365.25)) 

# Add PS and IPTW
ds <- IPTWfun(data = ds, 
              PSmodel = trt ~ age + female + weightbl + white + diagyrs + prmsgr + rlps1yr + trelmos + edssbl + tm25zbl + nhptzbl + pasatzbl + chrt2_5bl + sf36pcsbl + sf36mcsbl)

# Format data: y = postrelapse_num, other outcomes supplied in xcontinuousvars
temp <- format.countdata(data = ds, yvar = outcome, timevar = "time_inecp", trtvar = "trt", 
                         # xcontinuousvars = c(continuousvars, "offset", "iptw", "ps", "FUweight", "postrelapse_num"),
                         xcontinuousvars = c(continuousvars), # no need to include offset, iptw, ps, FUweight, postrelapse_num (will be included in the function)
                         xcategoricalvars = categoricalvars, imputation.method = NULL, scale = FALSE)
input <- data.frame(y = temp$y, trt = temp$trt, time = log(temp$time), temp$x) 
rm(temp)
rm(ds)

# The Batch CV loop
start <- Sys.time()
batchcv <- vector("list", batch_size)
names(batchcv) <- paste0("batch.ind", 1:batch_size)

for (cv.i in 1:batch_size){
  
  seed = base.seed*10000 + cv.i*100 + batch
  cat("\n\nCV iteration =", cv.i, "out of", batch_size, "with seed", seed)
  
  # Create CV folds
  set.seed(seed)
  folds <- createFolds(input$trt, k = n.fold, list = TRUE) # Stratified CV
  
  # Run each CV with the given method
  batchcv[[paste0("batch.ind", cv.i)]] <- eachCV(data = input, 
                                                 method = method, 
                                                 outcome = "y", # candidates -> outcome + plus formatting = outcome is always y 
                                                 folds = folds, 
                                                 n.fold = n.fold, 
                                                 categoricalvars = categoricalvars, 
                                                 continuousvars = continuousvars, 
                                                 seed = seed,
                                                 RCT = T)
}

end <- Sys.time()
cat("\nTime elapsed: ", end - start, "s\n\n")


output_dir <- paste0("./intermediate_results/case_study/case_study_stratified", n.fold, "foldCV")
if (!dir.exists(output_dir)) dir.create(output_dir)
save(batchcv, paste0("./intermediate_results/case_study/case_study_stratified", n.fold, "foldCV/case_study_stratified", n.fold, "foldCV_", method, "_", outcome,  "_batch", batch, ".RData"))



