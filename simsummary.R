# ------------------------------------------------------------------
# Product: A1 vs. A0
# Protocol: MarketScan simulated data
# Project: Precision Medicine MS
# 
# Program name: simsummary.R
# Developer/Programmer: pj
# Date: 01DEC2020
# 
# Purpose: Collect CV iteration batch results of the simulated data and 
#          summarize the mean and sd of value across PM methods.
#          Run this after simmain.R finishes for all methods.
# 
# Platform: Windows
# R Version: 4.0.0
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   15MAR2021 pj      Start the script (borrowed from summary.R in MarketScan)
#   03MAY2021 pj      Revise to incorporate arguments beta and percentiles
#   07MAY2021 pj      Adapt to new output structure (list)
#   10MAY2021 pj      Add betas and percentiles in all output data and filenames and summarize vs.dhat
#   16MAY2021 pj      dhats.big takes too much memory and processing power so not outputting
#   31MAY2021 gs      Calculate difference vhat.dhat - vs.dhat for each iteration/fold
#   14JUN2021 pj      Add absolute difference to V.hat for MAE calculation later
#   25JUN2021 gs      Add relative absolute difference
#   14JUL2021 pj      Minor change to sdVold (remove /n.fold)
#   24NOV2021 pj      Explore proportions
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))

setwd("C:/Users/xjiang1/OneDrive - Biogen/Documents/Innovation/PMMS/precmed_sim/") # CHANGE THIS #
source("utility.R")

# Define the summary function
simsummarize <- function(main_path, n, magnitude, distribution, n.cv = 25, batch_size = 1){
  
  setwd(main_path)
  
  # Parameters
  params <- convertParameters(magnitude, distribution, verbose = T)
  beta <- params$beta
  percentiles <- params$percentiles
  beta.text <- params$beta.text
  percentiles.text <- params$percentiles.text
  
  # Find a list of the CV batch results
  wd <- paste0('./outputs/simulations_n', n,'_stratified10foldCV/')
  
  results <- list.files(wd, pattern = paste0("^simulations.*_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".RData")) # select only files starting with "simulations"
  n <- str_extract(results, "(?<=simulations_n)[0-9]+(?=_)") %>% as.numeric() %>% unique()
  
  n.fold <- str_extract(results, "(?<=stratified)[0-9]+(?=foldCV)") %>% as.numeric() %>% unique()
  n.batch <- str_extract(results, "(?<=batch)[0-9]+(?=_beta)") %>% as.numeric() %>% max()
  n.cv <- n.batch * batch_size
  stopifnot(n.fold == 10)
  stopifnot(n.cv == n.cv)
  cat("\nn.cv and n.batch is:", n.cv, n.batch, "\n")
  
  cat("\nNumber of outputs found: ", length(results), "\n")
  cat("\n############################################################\n")
  
  # Read in each result in a loop
  vhats.dhat <- vs.dhat <- dhats <- NULL
  # dhats.big <- NULL
  for (result in results){
    
    method_outcome = str_extract(result, "(?<=foldCV\\_).*(?=\\_batch)")
    batch_index2 = as.numeric(str_extract(result, "(?<=batch)[0-9]+(?=\\_beta)"))
    
    if (!(method_outcome %in% c("listDTR3_mlogarr0001"))){
      # Read output file
      batchcv <- readRDS(paste0(wd, result))
      vhat.dhat <- v.dhat <- dhat <- dhat.big <- data.frame()
      
      for (name in names(batchcv)){
        # Get estimated values, vhat.dhat
        vhat.dhat  <- rbind(vhat.dhat, 
                            batchcv[[name]] %>% 
                              map_df(~bind_rows(names(.x) %>% str_detect("vhat") %>% keep(.x, .)), .id = "fold") %>% 
                              mutate(batch = name))
        # Get true values, v.dhat
        v.dhat  <- rbind(v.dhat, 
                         batchcv[[name]] %>% 
                           map_df(~bind_rows(names(.x) %>% str_detect("v.dhat") %>% keep(.x, .)), .id = "fold") %>% 
                           mutate(batch = name))
        # Get estimated rule from CV test fold, dhat
        dhat  <- rbind(dhat, 
                       batchcv[[name]] %>% 
                         map_df(~bind_rows(names(.x) %>% str_detect("^dhat$") %>% keep(.x, .)), .id = "fold") %>% 
                         mutate(batch = name))
        # # Get estimated rule from large test set, dhat.big
        # dhat.big  <- rbind(dhat.big, 
        #                    batchcv[[name]] %>% 
        #                      map_df(~bind_rows(names(.x) %>% str_detect("^dhat.big$") %>% keep(.x, .)), .id = "fold") %>% 
        #                      mutate(batch = name))
      }
      
      # Add other info to vhat.dhat
      vhat.dhat %<>% 
        mutate(method_outcome = method_outcome,
               batch_index2 = batch_index2,
               batch_size = batch_size,
               batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
               iteration = batch_index + (batch_index2 - 1) * batch_size) %>% 
        dplyr::select(-batch_index2, -batch_size, -batch)
      
      # Add other info to v.dhat
      v.dhat %<>% 
        mutate(method_outcome = method_outcome,
               batch_index2 = batch_index2,
               batch_size = batch_size,
               batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
               iteration = batch_index + (batch_index2 - 1) * batch_size) %>% 
        dplyr::select(-batch_index2, -batch_size, -batch)
      
      # Add other info to dhat
      dhat %<>% 
        mutate(method_outcome = method_outcome,
               batch_index2 = batch_index2,
               batch_size = batch_size,
               batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
               iteration = batch_index + (batch_index2 - 1) * batch_size) %>% 
        dplyr::select(-batch_index2, -batch_size, -batch)
      
      # # Add other info to dhat.big
      # dhat.big %<>% 
      #   mutate(method_outcome = method_outcome,
      #          batch_index2 = batch_index2,
      #          batch_size = batch_size,
      #          batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
      #          iteration = batch_index + (batch_index2 - 1) * batch_size) %>% 
      #   dplyr::select(-batch_index2, -batch_size, -batch)
      
      # Combine results over all methods
      vhats.dhat <- rbind(vhats.dhat, vhat.dhat)
      vs.dhat <- rbind(vs.dhat, v.dhat)
      dhats <- rbind(dhats, dhat)
      # dhats.big <- rbind(dhats.big, dhat.big)
    }  
  }
  
  dhats %<>% mutate(beta = beta.text, percentiles = percentiles.text)
  # dhats.big %<>% mutate(beta = beta.text, percentiles = percentiles.text)
  
  cat("\n")
  print(dim(vhats.dhat))
  print(head(vhats.dhat %>% arrange(method_outcome, iteration, fold)))
  table(vhats.dhat$iteration)
  table(vhats.dhat$batch_index)
  table(vhats.dhat$fold)
  table(vhats.dhat$method_outcome)
  cat("\n############################################################\n")
  cat("\n")
  print(head(vs.dhat))
  cat("\n############################################################\n")
  cat("\n")
  print(head(dhats))
  cat("\n############################################################\n")
  # cat("\n")
  # print(head(dhats.big))
  # cat("\n############################################################\n")
  
  #######################################################################
  ############################# Calculate Summary ############################
  #######################################################################
  
  vhats.dhat %<>%
    group_by(method_outcome) %>%
    dplyr::summarize(n.batches = n(),
              n.nonnaU = sum(!is.na(U)),
              n.nonnaW = sum(!is.na(W)),
              meanVold = mean(U/W, na.rm = T),
              meanV = sum(U, na.rm = T)/sum(W, na.rm = T),
              sdVold = sd(U/W, na.rm = T),
              sdV = sqrt(sum(sumRj2, na.rm = T) / (n.fold * (n.fold * n.cv  - 1))),
              .groups = "keep") %>%
    ungroup %>%
    arrange(desc(meanV)) %>%
    mutate(n = n,
           beta = beta.text,
           percentiles = percentiles.text)
  
  vs.dhat %<>% 
    group_by(method_outcome) %>% 
    dplyr::summarize(n.batches = n(),
              n.nonTrueV = sum(!is.na(v.dhat)),
              meanTrueV = mean(v.dhat, na.rm = T),
              sdTrueV = sd(v.dhat, na.rm = T),
              .groups = "keep") %>% 
    ungroup %>%
    arrange(desc(meanTrueV)) %>%
    mutate(n = n,
           beta = beta.text,
           percentiles = percentiles.text)
  
  print(vhats.dhat %>% dplyr::select(-contains("n.nonna")))
  print(vs.dhat)
  write_csv(vhats.dhat, paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_vhats.dhat_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
  write_csv(vs.dhat, paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_vs.dhat_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
  write_csv(dhats, paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_dhats_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
  # write_csv(dhats.big, paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_dhats.big_n", n, "_stratified", n.fold, "folds_beta", paste0(beta.text, collapse = "-"), "_perc", paste0(percentiles.text, collapse = "-"), ".csv"))
}


# Example
if (F){
  # Input arguments and constants
  args <- commandArgs(trailingOnly = TRUE) 
  n <- as.numeric(args[1])                                       # sample size of the simulated data (retrieved from user input) (based on what simmain.R has)
  batch_size <- ifelse(n %in% c(500, 1000), 5, 1)                # from simmain.R
  cat("\nSize of simulated data: ", n, "\n")
  beta <- eval(parse(text = args[2]))                            # level of heterogeneity 
  cat("\nLevel of heterogeneity beta =", beta, "\n")
  percentiles <- eval(parse(text = args[3]))                     # percentiles of subgroups
  cat("\nSubgroup proportions =", percentiles, "\n")
  
  # Run the summarize function with each parameter combination
  simsummarize(n, beta, percentiles)
}
