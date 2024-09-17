# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: simsummary_alln.R
# Developer/Programmer: pj
# Date: 04FEB2021
# 
# Purpose: Collect CV iteration batch results and summarize the mean and sd of value, agreement and
#          accuracy across different sample sizes, levels of heterogeneity, and percentiles.
#          Run this after main.sh and summary.R finishes for all methods and n's.
# 
# Platform: Windows
# R Version: 4.0.3
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   04FEB2021 pj      Start the main structure
#   24MAR2021 pj      Allow different subset sizes in MarketScan and simulation
#   09APR2021 pj      Minor changes
#   27APR2021 pj      Minor changes
#   10MAY2021 pj      Adapt to other betas or percentiles and plot trueV
#   17MAY2021 pj      Revise line plots, add agree table, minor edits
#   31MAY2021 gs      Add accuracy by sample size, minor change to agreement plot
#                     Move method_vec after read data
#                     Add plot for difference Vhat.dhat - V.dhat 
#                     Agreement based on folds with results
#   01JUN2021 pj      Minor edits   
#   07JUN2021 gs      Reorganize sections to have data formatting, then plots
#                     Save important data sets in RData space with tags for HTE and percentile
#   14JUN2021 pj      Move to cluster (based on summary_sample_size.R but with HPC directories and args)
#   25JUN2021 gs      Add relative absolute difference
#   02DEC2021 pj      Start the main structure
#   13DEC2021 pj      Add equal symmetry 20x5 proportion
#   07JAN2022 pj      Update trueA in dhat.cvs2 to accommodate for asymmetric proportions
#   17JAN2022 gs      Add accuracy calculations removing neutral group
#   24JAN2022 pj      Adjust code to n=5000
#   09MAY2022 pj      Add n=500
#   14JUN2022 gs      Add n tag in datasets exported in RData at the end of summarize()
#   12SEP2023 pj        Convert to HPC paths
#   02APR2024 pj      Keep Vold (simple avg/sd) and V (adj for CV correlation)
# ------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggnewscale))

setwd("C:/Users/xjiang1/OneDrive - Biogen/Documents/Innovation/PMMS/precmed_sim/") # CHANGE THIS #
source("01-setup.R")
source("utility.R")


# Define function
simsummarize_alln <- function(main_path, ns, magnitude, distribution){
  
  setwd(main_path)
  
  # Parameters
  params <- convertParameters(magnitude, distribution, verbose = T)
  beta <- params$beta
  percentiles <- params$percentiles
  
  # Find neutral group
  if (0 %in% beta) {
    neutral <- which(beta == 0) # look for the neutral group which has beta = 0 (log(1) = 0)
  } else {
    neutral <- 3 # if beta = -0.2 x 5, neutral group is set to be the middle group
  }
  
  # Constants (should be the same as simmain.R and main.R)
  n.fold <- 10     # number of folds in each CV iteration
  n.cv <- 25       # total number of CV iterations 
  base.seed <- 999 # randomization seed
  RCT <- T         # randomized trial, if TRUE
  big.n <- 1000000  # sample size of the large independent test set to get true value, e.g. 1million
  
  # Formatting
  myblue <- rgb(37, 15, 186, maxColorValue = 255)
  mygrey <- rgb(124, 135, 142, maxColorValue = 255)
  
  arg.batch <- 25
  # if (ns %in% c(500, 1000)) arg.batch = 5
  # arg.batch <- c(5, 5, 25, 25, 25) # this is the batch argument from command line
  
  #######################################################################
  ############################# Read in data ############################
  #######################################################################
  
  cvs <- true.cvs <- dhat.cvs <- NULL
  
  for (n in ns){
    vhats.dhat <- read_csv(paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_vhats.dhat_n", n, "_stratified", n.fold, "folds_beta", 
                                  paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"), col_types = cols())
    cvs <- rbind(cvs, vhats.dhat %>% 
                   mutate(n.subset = n) %>% 
                   dplyr::select(n.subset, n.nonnaU, n.nonnaW, method_outcome, n.batches, meanVold, meanV, sdVold, sdV, beta, percentiles))
    
    vs.dhat <- read_csv(paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_vs.dhat_n", n, "_stratified", n.fold, "folds_beta", 
                               paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"), col_types = cols())
    true.cvs <- rbind(true.cvs, vs.dhat %>% 
                        mutate(n.subset = n) %>% 
                        dplyr::select(n.subset, method_outcome, n.batches, n.nonTrueV, meanTrueV, sdTrueV, beta, percentiles))
    
    dhats <- read_csv(paste0("./outputs/simulations_n", n, "_stratified", n.fold, "foldCV/simmain_CV_results_dhats_n", n, "_stratified", n.fold, "folds_beta", 
                             paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"), col_types = cols())
    dhat.cvs <- rbind(dhat.cvs, dhats %>% mutate(n.subset = n) %>% dplyr::select(-batch_index))
    
  }
  
  alldmf <- "All A1"
  allteri <- "All A0"
  
  
  method.vec <- c(allteri, alldmf,  
                  "Poisson", "Weighted\n Poisson", "Negative\n Binomial", "Weighted\n NegBin", 
                  "Linear", "Weighted\n Linear", "dWOLS", 
                  "Boosting", "Two\n Regressions", "Contrast\n Regression", "List DTR\n (2 nodes)", "List DTR\n (3 nodes)")
  
  # One time run to get V(d)
  trueV <- getTrueOptimalValue(n = big.n, beta = beta, beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), RCT = RCT,
                               percentiles = percentiles, seed = 0)
  cat("\n## The optimal empirical V(d) is", trueV, "##")
  
  
  trueWorstV <- getTrueWorstValue(n = big.n, beta = beta, beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), RCT = RCT,
                                  percentiles = percentiles, seed = 0)
  cat("\n## The worst empirical V(d) is", trueWorstV, "##")
  
  
  
  #######################################################################
  ############################# Format data ############################
  #######################################################################
  
  #### ---- Cross-validated V.hat and sd plot by sample size for each PM method ---- ####
  
  # Format outputs for plotting
  cvs2 <- cvs %>% 
    mutate(method = stringr::str_split(method_outcome, "_") %>% map_chr(., 1),
           outcome = stringr::str_split(method_outcome, "_") %>% map_chr(., 2)) %>% 
    filter(outcome %in% c("postrelapse", "mlogarr0001", "logarr0001")) %>% 
    mutate(method = case_when(
      method == "contrastReg" ~ "Contrast\n Regression",
      method == "allA0" ~ allteri,
      method == "twoReg" ~ "Two\n Regressions",
      method == "boosting" ~ "Boosting",
      method == "allA1" ~ alldmf,
      method == "listDTR3" ~ "List DTR\n (3 nodes)",
      method == "listDTR2" ~ "List DTR\n (2 nodes)",
      method == "dWOLS" ~ "dWOLS",
      method == "weightedPoisson" ~ "Weighted\n Poisson",
      method == "poisson" ~ "Poisson",
      method == "weightedNegBin" ~ "Weighted\n NegBin",
      method == "weightedLinear" ~ "Weighted\n Linear",
      method == "linear" ~ "Linear",
      method == "negBin" ~ "Negative\n Binomial"
    ),
    method = factor(method, 
                    levels = method.vec, 
                    labels = method.vec)
    ) %>% 
    filter( !(method %in% c("Weighted\n Linear", "Weighted\n Poisson", "Weighted\n NegBin")))
  
  true.cvs2 <- true.cvs %>% 
    mutate(method = stringr::str_split(method_outcome, "_") %>% map_chr(., 1),
           outcome = stringr::str_split(method_outcome, "_") %>% map_chr(., 2)) %>% 
    filter(outcome %in% c("postrelapse", "mlogarr0001", "logarr0001")) %>% 
    mutate(method = case_when(
      method == "contrastReg" ~ "Contrast\n Regression",
      method == "allA0" ~ allteri,
      method == "twoReg" ~ "Two\n Regressions",
      method == "boosting" ~ "Boosting",
      method == "allA1" ~ alldmf,
      method == "listDTR3" ~ "List DTR\n (3 nodes)",
      method == "listDTR2" ~ "List DTR\n (2 nodes)",
      method == "dWOLS" ~ "dWOLS",
      method == "weightedPoisson" ~ "Weighted\n Poisson",
      method == "poisson" ~ "Poisson",
      method == "weightedNegBin" ~ "Weighted\n NegBin",
      method == "weightedLinear" ~ "Weighted\n Linear",
      method == "linear" ~ "Linear",
      method == "negBin" ~ "Negative\n Binomial"
    ),
    method = factor(method, 
                    levels = method.vec, 
                    labels = method.vec)
    ) %>% 
    filter( !(method %in% c("Weighted\n Linear", "Weighted\n Poisson", "Weighted\n NegBin")))
  
  # Add trueV and trueWorstV to cvs2
  cvs2 <- cvs2 %>% mutate(trueV = trueV, trueWorstV = trueWorstV)
  
  ## Concatenate V(d.hat) and V.hat(d.hat)
  cvs3 <- cvs2 %>% 
    dplyr::select(n.subset, n.nonnaU, method_outcome, n.batches, meanVold, sdVold, beta, percentiles, method, outcome) %>% 
    rename(n.nonna = n.nonnaU,
           meanV = meanVold, # TODO: change to meanV
           sdV = sdVold) %>% # TODO: change to sdV
    mutate(whichV = "V.hat")
  true.cvs3 <- true.cvs2 %>% 
    dplyr::select(n.subset, n.nonTrueV, method_outcome, n.batches, meanTrueV, sdTrueV, beta, percentiles, method, outcome) %>% 
    rename(n.nonna = n.nonTrueV,
           meanV = meanTrueV,
           sdV = sdTrueV) %>% 
    mutate(whichV = "V")
  all.cvs3 <- rbind(cvs3, true.cvs3) %>%
    mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                           labels = c("No HTE", "Low HTE", "Medium HTE", "High HTE", "Medium HTE", "Medium HTE")),
           symm = factor(percentiles, levels = c("c(0,0.2,0.4,0.6,0.8,1)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                         labels = c("Symmetric", "Asymmetric", "Asymmetric", "Symmetric", "Symmetric")),
           proportion = factor(percentiles, levels = c("c(0,0.2,0.4,0.6,0.8,1)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                               labels = c("20%-20%-20%-20%-20%", "55%-15%-15%-15%", "55%-30%-15%", "10%-15%-50%-15%-10%", "10%-30%-20%-30%-10%"))
    )
  
  all.cvs3.wide <- all.cvs3 %>% 
    pivot_wider(names_from = whichV, values_from = meanV:sdV) %>% 
    mutate(trueWorstV = trueWorstV,
           trueV = trueV,
           VR_V.hat = (meanV_V.hat - trueWorstV) / (trueV - trueWorstV),
           VR_V = (meanV_V - trueWorstV) / (trueV - trueWorstV))
  
  
  #### ---- Cross-validated d.hat by sample size for each PM method ---- ####
  dhat.cvs2 <- dhat.cvs %>% 
    mutate(method = stringr::str_split(method_outcome, "_") %>% map_chr(., 1),
           outcome = stringr::str_split(method_outcome, "_") %>% map_chr(., 2)) %>% 
    filter(outcome %in% c("postrelapse", "mlogarr0001", "logarr0001")) %>% 
    mutate(method = case_when(
      method == "contrastReg" ~ "Contrast\n Regression",
      method == "allA0" ~ allteri,
      method == "twoReg" ~ "Two\n Regressions",
      method == "boosting" ~ "Boosting",
      method == "allA1" ~ alldmf,
      method == "listDTR3" ~ "List DTR\n (3 nodes)",
      method == "listDTR2" ~ "List DTR\n (2 nodes)",
      method == "dWOLS" ~ "dWOLS",
      method == "weightedPoisson" ~ "Weighted\n Poisson",
      method == "poisson" ~ "Poisson",
      method == "weightedNegBin" ~ "Weighted\n NegBin",
      method == "weightedLinear" ~ "Weighted\n Linear",
      method == "linear" ~ "Linear",
      method == "negBin" ~ "Negative\n Binomial"
    ),
    method = factor(method, 
                    levels = method.vec, 
                    labels = method.vec)
    ) %>% 
    filter( !(method %in% c("Weighted\n Linear", "Weighted\n Poisson", "Weighted\n NegBin")))
  
  #### ---- Agreement by sample size ---- ####
  dhat.concat <- dhat.cvs2 %>% 
    mutate(iteration.fold = (iteration - 1) * 10 + as.numeric(str_extract(fold, "[0-9]+"))) %>% # there should be 250 iteration.folds (= 25 iterations * 10 folds)
    dplyr::select(beta, percentiles, n.subset, method, iteration.fold, dhat) %>% 
    group_by(beta, percentiles, n.subset, method, iteration.fold) %>% 
    mutate(i = 1:n()) %>%  # give each test observation of an iteration.fold a unique index 
    ungroup %>% 
    group_by(beta, percentiles, n.subset, method) %>% 
    mutate(iteration.fold.i = 1:n()) %>% # unique index for each test observation across iteration.folds
    ungroup 
  
  dhat.agreement <- vector(mode = "list", length = length(ns))
  names(dhat.agreement) <- paste0("n", ns)
  methods <- levels(dhat.concat$method)
  methods <- methods[!str_detect(methods, "Weighted|3")]
  m <- length(methods)
  for(this.n in ns) {
    C <- matrix(nrow = m, ncol = m)
    colnames(C) <- methods
    rownames(C) <- methods
    for(k in seq_len(m)){
      for(j in seq(k, m)){
        data.k <- dhat.concat %>% filter(n.subset == this.n, method == methods[k])
        data.j <- dhat.concat %>% filter(n.subset == this.n, method == methods[j])
        data.jk <- data.k %>% full_join(data.j, by = c("beta", "percentiles", "n.subset", "iteration.fold", "i"))
        # C[k, j] <- C[j, k] <- sum(data.jk$dhat.x == data.jk$dhat.y, na.rm = T) / (250*(this.n/10))
        C[k, j] <- C[j, k] <- sum(data.jk$dhat.x == data.jk$dhat.y, na.rm = T) / sum(is.na(data.jk$dhat.x) == FALSE & is.na(data.jk$dhat.y) == FALSE)
        # cat("\n", k, j)
      }
    }
    dhat.agreement[[paste0("n", this.n)]] <- C
  }
  
  #### ---- Accuracy by sample size and method ---- ####
  
  ## X variables to be included in each model
  categoricalvars <- c("female", "prevDMTefficacy")
  continuousvars <- c("ageatindex_centered", "prerelapse_num", "premedicalcost")
  
  ## Retrieve true d 
  if (length(unique(beta)) == 1){ 
    # If no heterogeneity
    dhat.cvs2 <- dhat.cvs2 %>% mutate(d = 1) # for beta = -0.2, true optimal treatment is A1
  } else {
    # If heterogeneity: generate random sample with correct seed and retrieve true d
    base.seed <- 999 
    n.fold <- 10
    n.cv <- 25 
    
    dhat.cvs2$d <- dhat.cvs2$neutral <- rep(NA, nrow(dhat.cvs2))
    
    for(i in 1:length(ns)) {
      n <- ns[i]
      batch <- arg.batch[i]
      batch_size <- ifelse(n %in% c(500, 1000), 5, 1)
      it <- 0 # Keep track of iteration in dhat.cvs2
      
      for(b in 1:batch) { 
        # batch = 5, batch_size = 5, n.fold = 10, total 250
        # batch = 25, batch_size = 1, n.fold = 10, total 250
        for(cv.i in 1:batch_size) {
          it <- it + 1
          seed <- base.seed + cv.i + b*10
          set.seed(seed)
          # print(seed)
          
          # Simulate a random sample 
          sim <- simdata(n = n, RCT = RCT, beta = beta, seed = seed, percentiles = percentiles)$data 
          
          # Calculate the true optimal treatment 
          sim <- sim %>% 
            mutate(trueA = ifelse(as.numeric(Iscore) < neutral, 1, 0), # TODO: this definition of trueA would change if the definition score groups changes
                   trueA = ifelse(as.numeric(Iscore) == neutral, trt, trueA)) # true optimal A
          
          # Format data
          temp <- format.countdata(data = sim, yvar = "postrelapse_num", timevar = "finalpostdayscount", trtvar = "trt",
                                   xcontinuousvars = c(continuousvars, "postrelapse_num", "offset", "FUweight"),
                                   xcategoricalvars = categoricalvars, imputation.method = NULL)
          input <- data.frame(y = temp$y, trt = factor(temp$trt), time = log(temp$time), temp$x)
          cat("\nA random sample is simulated with seed", seed, "with dimension: ", dim(input), "for the current CV iteration.\n")
          
          # Create CV folds
          folds <- createFolds(input$trt, k = n.fold, list = TRUE) # Stratified CV, follow the same as the simmain.R where folds were created on input$trt instead of sim$trt
          # print(folds[[1]])
          # print(head(input))
          # print(lapply(folds, length))
          # stop()

          for (fold.i in 1:n.fold){
            testdata <- sim[folds[[fold.i]],]
            # number of methods which succeeded for the given fold/batch. The "is.na(dhat) == FALSE" is to remove methods that didn't produce results for that fold/batch
            # TODO: check that this works
            nr <- nrow(dhat.cvs2 %>% filter(n.subset == n & fold == paste0("fold", fold.i) & iteration == it & is.na(dhat) == FALSE))
            dhat.cvs2$d[which(dhat.cvs2$n.subset == n & dhat.cvs2$fold == paste0("fold", fold.i) & dhat.cvs2$iteration == it & is.na(dhat.cvs2$dhat) == FALSE)] <- rep(testdata$trueA, nr/nrow(testdata))
            # Add column to identify who belongs to neutral group
            dhat.cvs2$neutral[which(dhat.cvs2$n.subset == n & dhat.cvs2$fold == paste0("fold", fold.i) & dhat.cvs2$iteration == it & is.na(dhat.cvs2$dhat) == FALSE)] <- ifelse(testdata$Iscore == neutral, 1, 0)
          }
        } # end of all cv iterations for sample size n
      } # end of all sample sizes
    }
  }

  
  ## Calculate % accuracy I(dhat == d)/n for each iteration & summary statistics
  dhat.cvs2 %<>%
    mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                           labels = c("No HTE", "Low HTE", "Medium HTE", "High HTE", "Medium HTE", "Medium HTE")),
           symm = factor(percentiles, levels = c("c(0,0.2,0.4,0.6,0.8,1)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                         labels = c("Symmetric", "Asymmetric", "Asymmetric", "Symmetric", "Symmetric")),
           proportion = factor(percentiles, levels = c("c(0,0.2,0.4,0.6,0.8,1)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                               labels = c("20%-20%-20%-20%-20%", "55%-15%-15%-15%", "55%-30%-15%", "10%-15%-50%-15%-10%", "10%-30%-20%-30%-10%"))
    )
  dat.accuracy <- dhat.cvs2 %>% 
    group_by(beta, percentiles, n.subset, method, fold, iteration, levhte) %>% 
    dplyr::summarise(accuracy = sum(dhat == d)/n(), 
                     accuracy.wo.neutral = sum(dhat[neutral == 0] == d[neutral == 0])/sum(neutral == 0))
  dat.accuracy.summary <- dat.accuracy %>% 
    group_by(beta, percentiles, n.subset, method, levhte) %>% 
    dplyr::summarise(mean.acc = mean(accuracy, na.rm = TRUE), 
                     sd.acc = sd(accuracy, na.rm = TRUE),
                     q1.acc = quantile(accuracy, prob = 0.25, na.rm = TRUE),
                     q3.acc = quantile(accuracy, prob = 0.75, na.rm = TRUE),
                     mean.acc.wo.neutral = mean(accuracy.wo.neutral, na.rm = TRUE), 
                     sd.acc.wo.neutral = sd(accuracy.wo.neutral, na.rm = TRUE),
                     q1.acc.wo.neutral = quantile(accuracy.wo.neutral, prob = 0.25, na.rm = TRUE),
                     q3.acc.wo.neutral = quantile(accuracy.wo.neutral, prob = 0.75, na.rm = TRUE))
  
  
  #### ---- Save all formatted data sets in RData  ---- ####
  # Rename data sets to export with tags
  assign(paste0("cvs2_", magnitude, "_", distribution), cvs2)
  assign(paste0("true.cvs2_", magnitude, "_", distribution), true.cvs2)
  assign(paste0("all.cvs3_", magnitude, "_", distribution), all.cvs3)
  assign(paste0("all.cvs3.wide_", magnitude, "_", distribution), all.cvs3.wide)
  assign(paste0("dhat.agreement_", magnitude, "_", distribution), dhat.agreement)
  assign(paste0("dat.accuracy_", magnitude, "_", distribution), dat.accuracy)
  assign(paste0("dat.accuracy.summary_", magnitude, "_", distribution), dat.accuracy.summary)
  assign(paste0("dhat.cvs2_", magnitude, "_", distribution), dhat.cvs2)
  save(file = paste0("./outputs/simsummary_sample_size_dataplot_", magnitude, "_", distribution, ".RData"), list = ls()[grepl("_", ls())])
}


if (F){
  # Input argument = n 
  args <- commandArgs(trailingOnly = TRUE)
  
  # Sample size
  ns <- as.numeric(args[1]) # {500, 1000, 2500, 5000, 100000} # Sample size for the simulated data #
  
  percentile_group <- list(c(0,0.55,0.65,0.75,0.85,1), 
                           c(0, 0.55, 0.65, 0.7, 0.85, 1), 
                           c(0, 0.1, 0.4, 0.6, 0.9, 1), 
                           c(0, 0.1, 0.25, 0.75, 0.9, 1), 
                           c(0,0.2,0.4,0.6,0.8,1))
  beta_group <- list(c(-0.92,-0.69,-0.69,-0.69,0), 
                     c(-0.92,-0.69,-0.69,0,0.1), 
                     c(-0.92,-0.69,0,0.1,0.18), 
                     c(-0.92,-0.69,0,0.1,0.18), 
                     c(-0.92,-0.69,0,0.1,0.18))
  
  # Run the summary function with all parameters
  for (i in 1:5){
    print(beta_group[[i]])
    print(percentile_group[[i]])
    simsummarize_all(ns, beta_group[[i]], percentile_group[[i]])
  }
}


