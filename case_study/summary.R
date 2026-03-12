# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: summary.R
#
# Purpose: Collect CV iteration batch results and summarize the mean and sd of value.
#          Run this after main.R finishes. 
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(corrplot))

# homepath <- "/home/pjiang/pmms/" # CHANGE THIS #
# setwd(homepath) 

#######################################################################
############################# Read in data ############################
#######################################################################

# Define constants
n.fold <- 10
n.cv <- 25
batch_size <- 5

# Load CV outputs from main.R
wd <- paste0("./intermediate_results/case_study/case_study_stratified", n.fold, "foldCV/")
results <- list.files(wd, pattern = "^case_study.*") # select only files starting with "case_study"

cat("\nNumber of outputs found: ", length(results), "\n")
cat("\n############################################################\n")

# Read in each result in a loop
vhats.dhat <- NULL
for (result in results){

  method_outcome = str_extract(result, "(?<=foldCV\\_).*(?=\\_)")
  method = str_split(method_outcome, "_")[[1]][1]
  batch_index2 = as.numeric(str_extract(result, "(?<=batch)[0-9]+(?=.RData)"))
    
  slurmoutput <- readLines(paste0('./case_study/slurmoutputs/slurm-main_', n.fold, 'folds_', method, '_', batch_index2, '.out'))

  realtime <- slurmoutput[str_detect(slurmoutput, "^real")]
  cat("\n", result, method_outcome, batch_index2, realtime)
  minute <- str_extract(realtime, "(?<=real).*(?=m)") %>% as.numeric()
  if (is.na(minute)){
    second <- str_extract(realtime, "(?<=real).*(?=s)") %>% as.numeric()
  }
  else {
    second <- str_extract(realtime, "(?<=m).*(?=s)") %>% as.numeric()
  }
  cat('\n', minute, second, sum(minute * 60, second, na.rm = T))
  elapsedtime <- slurmoutput[str_detect(slurmoutput, "Time elapsed")]
  elapsedTimeInSeconds <- str_extract(elapsedtime, "(?<=Time elapsed:).*(?=s)") %>% as.numeric()
  
  batchcv <- load(paste0(wd, result))
  vhat.dhat <- data.frame()

  for (name in names(batchcv)){
    # Get estimated values, vhat.dhat
    vhat.dhat <- rbind(vhat.dhat, batchcv[[name]] %>% map_df(~bind_rows(names(.x) %>% str_detect("vhat") %>% keep(.x, .)), .id = "fold") %>% mutate(batch = name))
  }
  
  # Add other info to vhat.dhat
  vhat.dhat %<>% 
    mutate(method_outcome = method_outcome,
           batch_index2 = batch_index2,
           batch_index = as.numeric(str_extract(batch, "[0-9]+")), # Convert counter from batch specific iteration to total iteration
           iteration = batch_index + (batch_index2 - 1) * batch_size, 
           realTimeInSeconds = sum(minute * 60, second, na.rm = T),
           elapsedTimeInSeconds = elapsedTimeInSeconds)
  
  # Combine results over all methods
  vhats.dhat <- rbind(vhats.dhat, vhat.dhat)
}


write_csv(vhats.dhat, paste0("./intermediate_results/case_study/case_study_stratified", n.fold, "foldCV/main_CV_results_rawvhats.dhat_stratified", n.fold, "folds.csv"))

