# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: simmain.R
# 
# Purpose: Main function to run the stratified CV of PM method implementation 
#  to simulated data in batches (PM methods and CV batch be run separately)
#  across different sample sizes and magnitudes of HTE (fixed 20x5 proportion)
# 
#  This was run on high-performance cluster via simmain.sh
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(haven))
suppressPackageStartupMessages(library(Hmisc))
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
source("./02-propensityscore.R")
source("./03-dWOLS.R")
source("./04-regression-based.R")
source("./05-listdtr.R")
source("./06-LuScore.R")
source("./eachCV.R")

args <- commandArgs(trailingOnly = TRUE)
print(args)

# Constants
n.fold <- 10     # number of folds in each CV iteration
n.cv <- 25       # total number of CV iterations 
base.seed <- 999 # randomization seed
RCT <- T         # randomized trial, if TRUE
big.n <- 1000000  # sample size of the large independent test set to get true value, e.g. 1million

# User-specified constants
method <- args[1]               # PM method, an argument from command line
yvar <- args[2]                 # y-variable specified by user, could be logarr0001 = log(arr+0.001) or logarr1 = log(arr+1) or log((NRL+0.01)/offsetbyyear)
batch <- as.numeric(args[3])    # the batch index, could be from 1 to n.cv/batch_size, an argument from command line
n <- as.numeric(args[4])        # sample size of the randomly generated simulated data
beta <- eval(parse(text = args[5]))  # level of heterogeneity 
cat("\nLevel of heterogeneity beta =", beta, "\n")
percentiles <- eval(parse(text = args[6])) # percentiles of subgroups
cat("\nSubgroup proportions =", percentiles, "\n")

batch_size <- ifelse(n %in% c(100, 250, 500, 1000), 5, 1)  # number of CV iterations in each batch, could be from 1 to n.cv

## X variables to be included in each model
categoricalvars <- c("female", "prevDMTefficacy")
continuousvars <- c("ageatindex_centered", "prerelapse_num", "premedicalcost")

## Outcome depends on the method and user
candidates <- list(allA1 = "postrelapse_num",
                   allA0 = "postrelapse_num",
                   linear = yvar, # default = "logarr0001", lower is better
                   negBin = "postrelapse_num",
                   dWOLS = paste("m", yvar, sep = ""), # default = "mlogarr0001", add minus because higher is better
                   listDTR2 = paste("m", yvar, sep = ""), # default = "mlogarr0001", add minus because higher is better
                   listDTR3 = paste("m", yvar, sep = ""), # default = "mlogarr0001", add minus because higher is better
                   poisson = "postrelapse_num",
                   boosting = "postrelapse_num",
                   twoReg = "postrelapse_num",
                   contrastReg = "postrelapse_num") 
outcome <- candidates[[method]]
cat("\nOutcome is:", outcome, "\n")

# Simulate a large test set 
sim.big <- simdata(n = big.n, RCT = RCT, beta = beta, seed = base.seed, percentiles = percentiles)
cat("\nA random sample is simulated with seed", base.seed, "with dimension: ", dim(sim.big$data), "as the one independent large test data for calculation of ture value function.\n")

# One time run to get V(d)
trueV <- getTrueOptimalValue(n = big.n, beta = beta, beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), RCT = RCT,
                             percentiles = percentiles, seed = 0)
cat("\n###################################
    \n##The empirical V(d) is", trueV, "##
    \n###################################\n")

#######################################################################
############################# Read in data ############################
#######################################################################

# X categorical variables (formatted) to be included in each model
formatted_categoricalvars <- c("female", "prevDMTefficacy_Medium.and.high.efficacy", "prevDMTefficacy_None")

# The Batch CV loop
start <- Sys.time()
batchcv <- vector("list", batch_size)
names(batchcv) <- paste0("batch.ind", 1:batch_size)
  
for (cv.i in 1:batch_size){
  
  cat("\n\nCV iteration =", cv.i, "out of", batch_size)
  seed = base.seed + cv.i + batch*10
  set.seed(seed)
  
  # Simulate a random sample 
  sim <- simdata(n = n, RCT = RCT, beta = beta, seed = seed, percentiles = percentiles)$data 
  
  # Format data
  temp <- format.countdata(data = sim, yvar = outcome, timevar = "finalpostdayscount", trtvar = "trt", 
                           xcontinuousvars = c(continuousvars, "postrelapse_num", "offset", "FUweight"), 
                           xcategoricalvars = categoricalvars, imputation.method = NULL)
  input <- data.frame(y = temp$y, trt = factor(temp$trt), time = log(temp$time), temp$x)
  rm(sim)
  cat("\nA random sample is simulated with seed", seed, "with dimension: ", dim(input), "for the current CV iteration.\n")
  
  # Create CV folds
  folds <- createFolds(input$trt, k = n.fold, list = TRUE) # Stratified CV
  
  batchcv[[paste0("batch.ind", cv.i)]] <- eachCV(data = input, 
                                                 method = method, 
                                                 outcome = "y", # candidates -> outcome + plus formatting = outcome is always y 
                                                 folds = folds, 
                                                 n.fold = n.fold, 
                                                 categoricalvars = formatted_categoricalvars, 
                                                 continuousvars = continuousvars,
                                                 RCT = RCT,
                                                 seed = seed,
                                                 sim.big = sim.big)
}

end <- Sys.time()
cat("\nTime elapsed: ", end - start, "s\n\n")


output_dir <- paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified", n.fold, "foldCV")
if (!dir.exists(output_dir)) dir.create(output_dir)
save(batchcv, paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simulations_n", n, "_stratified", n.fold, "foldCV_", method, "_", outcome,  "_batch", batch, "_beta", paste0(round(beta, 2), collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".RData"))



