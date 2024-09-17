# ------------------------------------------------------------------
# Product: A1 vs. A0
# Protocol: MarketScan simulated data
# Project: Precision Medicine MS
# 
# Program name: simmain.R
# Developer/Programmer: pj
# Date: 08MAR2020
# 
# Purpose: Main function to run the stratified CV of PM method implementation 
#  to simulated data in batches (PM methods and CV batch be run separately)
#    
# 
# Platform: Windows
# R Version: 4.0.1
# 
#   Modifications:
# 
# Date			  By			  Description
# --------		--------	-----------------------------
# 08MAR2021   pj        Start the main structure (based on main.R in MarketScan data)
# 11MAR2021   gs        Replaced estimated ps with RCT option
#                       No nuisance variables for now, trt as a factor to match
#                       Comment methods we won't be using for now
# 11MAR2021   pj        Calculate PS/IPTW after splitting into training/test set
# 26APR2021   pj        Add the V(d.hat) results to output; add V(d) calculation prior to the for loop
# 16MAY2021   pj        Round beta to 2 decimals
# 25OCT2021   pj        Explore proportion 
# 01AUG2023   pj        Adapt script to new cluster paths and change to 1:1 Rx allocation
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(haven))
# suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(pscl))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(mpath))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(fastDummies))
suppressPackageStartupMessages(library(listdtr))
suppressPackageStartupMessages(library(DTRreg))



# Define function
simmain <- function(main_path, method, n, magnitude, distribution, batch_ind, n.cv = 25, batch_size = 1){
  
  #' Main function to run the PM simulation for a given scenario specified by \code{method}, 
  #' \code{n}, \code{magnitude}, \code{distribution} via a certain number of CV iterations
  #' 
  #' @param method PM method to use; choose from {'allA1', 'allA0', 'linear', 'poisson','dWOLS', 'listDTR2', 'boosting', 'twoReg', 'contrastReg'}
  #' @param n number of total sample size; choose from {500, 1000, 2500, 5000, 10000}
  #' @param magnitude magnitude of the HTE; choose from {'no', 'low', 'medium', 'high'}
  #' @param distribution distribution of the HTE patient groups; 
  #'    choose from {'symm 20x5', 'symm 10-15-50-15-10', 'symm 10-30-20-30-10', 'asymm 55-30-15', 'asymm 55-15-15-15'}
  #' @param batch_ind the batch index for the current run, could be from 1 to \code{batch_size}
  #' @param n.cv number of total CV iterations
  #' @param batch_size number of CV iterations in each batch, could be from 1 to n.cv
  
  setwd(main_path)
  source("utility.R")
  source("eachCV.R")
  source("01-setup.R")
  source("02-propensityscore.R")
  source("03-dWOLS.R")
  source("04-regression-based.R")
  source("05-listdtr.R")
  source("06-LuScore.R")
  
  # Constants
  n.fold <- 10     # number of folds in each CV iteration
  base.seed <- 999 # randomization seed
  RCT <- T         # randomized trial, if TRUE
  big.n <- 1000000 # sample size of the large independent test set to get true value, e.g. 1million
  yvar = "logarr0001" # transformed y variable for PM method that cannot accept count outcomes
  stopifnot(batch_size <= n.cv)
  cat("### BATCH No.", batch_ind, "ouf of", n.cv / batch_size, "###")
  
  # Parameters
  params <- convertParameters(magnitude, distribution, verbose = T)
  beta <- params$beta
  percentiles <- params$percentiles
    
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
  cat("\n* Outcome =", outcome)
  
  # Simulate a large test set 
  sim.big <- simdata(n = big.n, RCT = RCT, beta = beta, seed = base.seed, percentiles = percentiles)
  cat("\nA random sample is simulated with seed", base.seed, "with dimension:", dim(sim.big$data), "as the one independent large test data for calculation of ture value function.")
  
  # One time run to get V(d)
  trueV <- getTrueOptimalValue(n = big.n, beta = beta, beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), RCT = RCT,
                               percentiles = percentiles, seed = 0)
  cat("\nThe empirical V(d) is", trueV, "\n")
  
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
    
    cat("\n## CV iteration =", cv.i, "out of", batch_size, "##")
    seed = base.seed + cv.i + batch_ind*n.fold
    set.seed(seed)
    
    # Simulate a random sample 
    sim <- simdata(n = n, RCT = RCT, beta = beta, seed = seed, percentiles = percentiles)$data 
    
    # Format data
    temp <- format.countdata(data = sim, yvar = outcome, timevar = "finalpostdayscount", trtvar = "trt", 
                             xcontinuousvars = c(continuousvars, "FUweight"), 
                             xcategoricalvars = categoricalvars, 
                             RCT = RCT, imputation.method = NULL)
    input <- data.frame(y = temp$y, trt = factor(temp$trt), time = log(temp$time), temp$x)
    rm(sim)
    cat("\nA random sample is simulated with seed", seed, "with dimension: ", dim(input), "for the current CV iteration.")
    
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
  
  output_dir <- paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV")
  if (!dir.exists(output_dir)) dir.create(output_dir)
  saveRDS(batchcv, paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simulations_n", n, "_stratified", n.fold, "foldCV_", method, "_", outcome,  "_batch", batch_ind, "_beta", paste0(round(beta, 2), collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".RData"))
}

# Example 
if (F){
  # User-specified constants
  args <- commandArgs(trailingOnly = TRUE)
  method <- args[1]               # PM method, an argument from command line
  yvar <- args[2]                 # y-variable specified by user, could be logarr0001 = log(arr+0.001) or logarr1 = log(arr+1) or log((NRL+0.01)/offsetbyyear)
  batch_ind <- as.numeric(args[3])    # the batch index, could be from 1 to n.cv/batch_size, an argument from command line
  n <- as.numeric(args[4])        # sample size of the randomly generated simulated data
  beta <- eval(parse(text = args[5]))  # level of heterogeneity 
  percentiles <- eval(parse(text = args[6])) # percentiles of subgroups
  
  simmain(method, yvar, batch_ind, n, beta, percentiles)
}

# Debugging
if (F){
  method = 'listDTR2'
  yvar = "logarr0001"
  method = "negBin"
  yvar = "postrelapse_num"
  batch_ind = 1
  n = 500
  beta = c(-0.2, -0.2, -0.2, -0.2, -0.2)
  percentiles = seq(0, 1, by = 0.2)
  cv.i = 1
  fold.i = 1
}
