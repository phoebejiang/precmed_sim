# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: simsummary.R
# 
# Purpose: Collect CV iteration batch results of the simulated data and summarize the mean and sd of value.
#          Run this after simmain.R finishes. 
#
#          This was run on high-performance cluster via simsummary.sh
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))

# homepath <- "/home/pjiang/pmms/" # CHANGE THIS #
# setwd(homepath) 

#######################################################################
############################# Read in data ############################
#######################################################################

# Constants and setups
args <- commandArgs(trailingOnly = TRUE) 
n <- as.numeric(args[1])                                       # sample size of the simulated data (retrieved from user input) (based on what simmain.R has)
batch_size <- ifelse(n %in% c(500, 1000), 5, 1)                # from simmain.R
cat("\nSize of simulated data: ", n, "\n")
beta <- eval(parse(text = args[2]))                            # level of heterogeneity 
cat("\nLevel of heterogeneity beta =", beta, "\n")
percentiles <- eval(parse(text = args[3]))                     # percentiles of subgroups
cat("\nSubgroup proportions =", percentiles, "\n")

# Find a list of the CV batch results
wd <- paste0('./simulations/samplesize/outputs/simulations_n', n,'_stratified10foldCV/')

results <- list.files(wd, pattern = paste0("^simulations.*_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".RData")) # select only files starting with "simulations"
n <- str_extract(results, "(?<=simulations_n)[0-9]+(?=_)") %>% as.numeric() %>% unique()

n.fold <- str_extract(results, "(?<=stratified)[0-9]+(?=foldCV)") %>% as.numeric() %>% unique()
n.batch <- str_extract(results, "(?<=batch)[0-9]+(?=\\.RData)") %>% as.numeric() %>% max()
n.cv <- n.batch * batch_size


cat("\nNumber of outputs found: ", length(results), "\n")
cat("\n############################################################\n")

# Read in each result in a loop
vhats.dhat <- vs.dhat <- dhats <- dhats.big <- NULL
for (result in results){

  method_outcome = str_extract(result, "(?<=foldCV\\_).*(?=\\_batch)")
  batch_index2 = as.numeric(str_extract(result, "(?<=batch)[0-9]+(?=\\_beta)"))
  
  if (!(method_outcome %in% c("listDTR3_mlogarr0001"))){

    # Read output file
    batchcv <- load(paste0(wd, result))
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
    }
    
    # Retrieve slurm outputs
    slurmoutput <- readLines(paste0('./simulations/samplesize/code/slurmoutputs/slurm-simmain_n', 
                                    n,'_', n.fold, 'folds/slurm-main_n', n, '_', n.fold, 'folds_', 
                                    str_replace(method_outcome, "mlogarr", "logarr"), '_', 
                                    str_extract(batch_index2, "[0-9]+"), "_", args[2], "_", args[3], '.out'))
    
    realtime <- slurmoutput[str_detect(slurmoutput, "^real")]
    cat("\n", result, method_outcome, batch_index2, realtime)
    minute <- str_extract(realtime, "(?<=real).*(?=m)") %>% as.numeric()
    if (is.na(minute)){
      second <- str_extract(realtime, "(?<=real).*(?=s)") %>% as.numeric()
    } else {
      second <- str_extract(realtime, "(?<=m).*(?=s)") %>% as.numeric()
    }
    cat('\n', minute, second, sum(minute * 60, second, na.rm = T))
    elapsedtime <- slurmoutput[str_detect(slurmoutput, "Time elapsed")]
    elapsedTimeInSeconds <- str_extract(elapsedtime, "(?<=Time elapsed:).*(?=s)") %>% as.numeric()
    
    # Add other info to vhat.dhat
    vhat.dhat %<>% 
      mutate(method_outcome = method_outcome,
             batch_index2 = batch_index2,
             batch_size = batch_size,
             realTimeInSeconds = sum(minute * 60, second, na.rm = T),
             elapsedTimeInSeconds = elapsedTimeInSeconds, 
             batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
             iteration = batch_index + (batch_index2 - 1) * batch_size) %>% 
      select(-batch_index2, -batch_size, -batch)
    
    # Add other info to v.dhat
    v.dhat %<>% 
      mutate(method_outcome = method_outcome,
             batch_index2 = batch_index2,
             batch_size = batch_size,
             batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
             iteration = batch_index + (batch_index2 - 1) * batch_size) %>% 
      select(-batch_index2, -batch_size, -batch)
      
    # Add other info to dhat
    dhat %<>% 
      mutate(method_outcome = method_outcome,
             batch_index2 = batch_index2,
             batch_size = batch_size,
             batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
             iteration = batch_index + (batch_index2 - 1) * batch_size) %>% 
      select(-batch_index2, -batch_size, -batch)
        
    # Combine results over all methods
    vhats.dhat <- rbind(vhats.dhat, vhat.dhat)
    vs.dhat <- rbind(vs.dhat, v.dhat)
    dhats <- rbind(dhats, dhat)
    dhats.big <- rbind(dhats.big, dhat.big)
  }  
}

dhats %<>% 
  mutate(beta = args[2],
         percentiles = args[3])
dhats.big %<>% 
  mutate(beta = args[2],
         percentiles = args[3])


#######################################################################
############################# Calculate Summary ############################
#######################################################################
# Add difference Vhat.dhat - v.dhat
vhats.dhat$valueDiff <- vhats.dhat$U/vhats.dhat$W - vs.dhat$v.dhat
# Add absolute difference |Vhat.dhat - v.dhat|
vhats.dhat$valueAbsDiff <- abs(vhats.dhat$valueDiff)
# Add relative absolute difference |Vhat.dhat - v.dhat/vs.dhat$v.dhat|
vhats.dhat$valueAbsDiffRel <- abs(vhats.dhat$valueDiff/vs.dhat$v.dhat)

# Summarize by method
vhats.dhat %<>% 
  group_by(method_outcome) %>% 
  summarize(n.batches = n(),
            n.nonnaU = sum(!is.na(U)),
            n.nonnaW = sum(!is.na(W)),
            meanVold = mean(U/W, na.rm = T),
            meanV = sum(U, na.rm = T)/sum(W, na.rm = T),
            sdVold = sd(U/W, na.rm = T), 
            meanU = mean(U, na.rm = T),
            meanW = mean(W, na.rm = T),
            sdV = sum((U / meanW - meanU * W / ((meanW)^2))^2, na.rm = T) / (n.fold * (n.fold * n.cv  - 1)),
            meanVdiff = mean(valueDiff, na.rm = T),
            sdVdiff = sd(valueDiff, na.rm = T),
            meanVabsDiff = mean(valueAbsDiff, na.rm = T),
            sdVabsDiff = sd(valueAbsDiff, na.rm = T),
            meanVabsDiffRel = mean(valueAbsDiffRel, na.rm = T),
            sdVabsDiffRel = sd(valueAbsDiffRel, na.rm = T), 
            medianRealTimeInSeconds = median(realTimeInSeconds, na.rm = T),
            medianElapsedTimeInSeconds = median(elapsedTimeInSeconds, na.rm = T),
            minRealTimeInSeconds = min(realTimeInSeconds, na.rm = T),
            maxRealTimeInSeconds = max(realTimeInSeconds, na.rm = T),
            minElapsedTimeInSeconds = min(elapsedTimeInSeconds, na.rm = T),
            maxElapsedTimeInSeconds = max(elapsedTimeInSeconds, na.rm = T),
            .groups = "keep") %>% 
  select(-meanU, -meanW) %>%
  ungroup %>%
  arrange(desc(meanV)) %>%
  mutate(n = n,
         beta = args[2],
         percentiles = args[3])

vs.dhat %<>% 
  group_by(method_outcome) %>% 
  summarize(n.batches = n(),
            n.nonTrueV = sum(!is.na(v.dhat)),
            meanTrueV = mean(v.dhat, na.rm = T),
            sdTrueV = sd(v.dhat, na.rm = T),
            .groups = "keep") %>% 
  ungroup %>%
  arrange(desc(meanTrueV)) %>%
  mutate(n = n,
         beta = args[2],
         percentiles = args[3])

print(vhats.dhat %>% select(-contains("min"), -contains("max"), -meanV, -sdV))
print(vs.dhat)
write_csv(vhats.dhat, paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_vhats.dhat_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
write_csv(vs.dhat, paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_vs.dhat_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
write_csv(dhats, paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_dhats_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
write_csv(dhats.big, paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_dhats.big_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))

