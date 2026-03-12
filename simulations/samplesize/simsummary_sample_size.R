# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: simsummary_sample_size.R
#
# Purpose: Collect CV iteration batch results and summarize the mean and sd of value, agreement and
#          accuracy across different sample sizes. Generate plots by sample size.
#          Run this after simmain.sh and simsummary.R finishes for all sample sizes 
# 
#          This was run on high-performance cluster via simsummary_sample_size.sh
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(caret))

# homepath <- "/home/pjiang/pmms/" # CHANGE THIS #
# setwd(homepath) 

source("./utility.R")

# User specified
args <- commandArgs(trailingOnly = TRUE)
input <- "simulation" # {"MarketScan", "simulation"}
beta <- eval(parse(text = args[1]))  # level of heterogeneity  # {rep(-0.2, 5), c(-1.2,-0.69,0,0.1,0.41), c(-0.36,-0.29,0,0.05,0.1), c(-0.92, -0.69, 0, 0.1, 0.18)}
cat("\nLevel of heterogeneity beta =", beta, "\n")
percentiles <- eval(parse(text = args[2])) # percentiles of subgroups # {seq(0, 1, by = 0.2)}
cat("\nSubgroup proportions =", percentiles, "\n")

# Constants (should be the same as simmain.R and main.R)
n.fold <- 10     # number of folds in each CV iteration
n.cv <- 25       # total number of CV iterations 
base.seed <- 999 # randomization seed
RCT <- T         # randomized trial, if TRUE
big.n <- 1000000  # sample size of the large independent test set to get true value, e.g. 1million

if (F){
  trueOptimalValuePHB(n = 100000, beta = c(-2, -1, 0, 1, 2)) # high level of heterogeneity
  trueOptimalValuePHB(n = 100000, beta = c(-1, -0.5, 0, 0.5, 1)) # low level of heterogeneity
  trueOptimalValuePHB(n = 100000, beta = c(0, 0, 0, 0, 0)) # no heterogeneity
  
  trueOptimalValueGAB(n = 10000000, beta = c(-2, -1, 0, 1, 2))
  trueOptimalValueGAB(n = 1000000, beta = c(-1, -0.5, 0, 0.5, 1))
  trueOptimalValueGAB(n = 1000000, beta = c(0, 0, 0, 0, 0))
}

# Formatting
myblue <- rgb(37, 15, 186, maxColorValue = 255)
mygrey <- rgb(124, 135, 142, maxColorValue = 255)

#######################################################################
############################# Read in data ############################
#######################################################################

cvs <- true.cvs <- dhat.cvs <- NULL
if (input == "MarketScan"){
  
  ns <- c(800, 2000, 4000, 8599)
  for (n in ns){
    if (n != 8599){
      cv <- read_csv(paste0("./simulations/samplesize/outputs/main_CV_results_n", n, "_stratified10folds_resampled.csv")) 
    } else{
      cv <- read_csv(paste0("./simulations/samplesize/outputs/main_CV_results_n", n, "_stratified10folds.csv"))  
    }
    
    cvs <- rbind(cvs, cv %>% 
                   mutate(n.subset = n) %>% 
                   dplyr::select(n.subset, method_outcome, n.batches, meanVold, meanV, sdVold, sdV, medianRealTimeInSeconds))
  }
  alldmf <- "All DMF"
  allteri <- "All TERI"
  
} else if (input == "simulation"){
 
  ns <- c(500, 1000, 2500, 5000, 10000)
  for (n in ns){
    vhats.dhat <- read_csv(paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified10foldCV/simmain_CV_results_vhats.dhat_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
    cvs <- rbind(cvs, vhats.dhat %>% 
                   mutate(n.subset = n) %>% 
                   dplyr::select(n.subset, n.nonnaU, n.nonnaW, method_outcome, n.batches, meanVold, meanV, sdVold, sdV, 
                   meanVdiff, sdVdiff, meanVabsDiff, sdVabsDiff, meanVabsDiffRel, sdVabsDiffRel, medianRealTimeInSeconds, beta, percentiles))
    
    vs.dhat <- read_csv(paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified10foldCV/simmain_CV_results_vs.dhat_n", n, "_stratified", n.fold, "folds_beta", paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
    true.cvs <- rbind(true.cvs, vs.dhat %>% 
                        mutate(n.subset = n) %>% 
                        dplyr::select(n.subset, method_outcome, n.batches, n.nonTrueV, meanTrueV, sdTrueV, beta, percentiles))
    
    dhats <- read_csv(paste0("./simulations/samplesize/outputs/simulations_n", n, "_stratified10foldCV/simmain_CV_results_dhats_n", n, "_stratified", n.fold, "folds_beta", 
                             paste0(beta, collapse = "-"), "_perc", paste0(percentiles, collapse = "-"), ".csv"))
    
    dhat.cvs <- rbind(dhat.cvs, dhats %>% mutate(n.subset = n) %>% dplyr::select(-batch_index))
  }
  
  alldmf <- "All A1"
  allteri <- "All A0"
} 

method.vec <- c(allteri, alldmf,  
                "Poisson", "Weighted\n Poisson", "Negative\n Binomial", "Weighted\n NegBin", 
                "Linear", "Weighted\n Linear", "dWOLS", 
                "Boosting", "Two\n Regressions", "Contrast\n Regression", "List DTR\n (2 nodes)", "List DTR\n (3 nodes)")

# One time run to get V(d)
trueV <- getTrueOptimalValue(n = big.n, beta = beta, beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), RCT = RCT,
                             percentiles = percentiles, seed = 0)
cat("\n###################################
    \n##The optimal empirical V(d) is", trueV, "##
    \n###################################\n")


trueWorstV <- getTrueWorstValue(n = big.n, beta = beta, beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003), RCT = RCT,
                                percentiles = percentiles, seed = 0)
cat("\n###################################
    \n##The worst empirical V(d) is", trueWorstV, "##
    \n###################################\n")




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
  filter( !(method %in% c("Weighted\n Linear", "Weighted\n Poisson", "Weighted\n NegBin"))) %>%
  mutate(trueWorstV = trueWorstV, trueV = trueV)

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


## Concatenate V(d.hat) and V.hat(d.hat)
cvs3 <- cvs2 %>% 
  dplyr::select(n.subset, n.nonnaU, method_outcome, n.batches, meanVold, sdVold, beta, percentiles, method, outcome) %>% 
  rename(n.nonna = n.nonnaU,
         meanV = meanVold,
         sdV = sdVold) %>% 
  mutate(whichV = "V.hat")
true.cvs3 <- true.cvs2 %>% 
  dplyr::select(n.subset, n.nonTrueV, method_outcome, n.batches, meanTrueV, sdTrueV, beta, percentiles, method, outcome) %>% 
  rename(n.nonna = n.nonTrueV,
         meanV = meanTrueV,
         sdV = sdTrueV) %>% 
  mutate(whichV = "V")
all.cvs3 <- rbind(cvs3, true.cvs3)

all.cvs3.wide <- all.cvs3 %>% 
  pivot_wider(names_from = whichV, values_from = meanV:sdV) %>%
  mutate(trueWorstV = trueWorstV, trueV = trueV) %>% 
  mutate(VR_V.hat = (meanV_V.hat - trueWorstV) / (trueV - trueWorstV),
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
      C[k, j] <- C[j, k] <- sum(data.jk$dhat.x == data.jk$dhat.y, na.rm = T) / sum(is.na(data.jk$dhat.x) == FALSE & is.na(data.jk$dhat.y) == FALSE)
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
  arg.batch <- c(5, 5, 25, 25, 25) # this is the batch argument from command line
  
  dhat.cvs2$d <- rep(NA, nrow(dhat.cvs2))
  
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
        
        # Simulate a random sample 
        sim <- simdata(n = n, RCT = RCT, beta = beta, seed = seed, percentiles = percentiles)$data 
        
        # Calculate the true optimal treatment 
        sim <- sim %>% 
          mutate(trueA = ifelse(as.numeric(Iscore) <= 2, 1, 0), # TODO: this definition of trueA would change if the definition score groups changes
                 trueA = ifelse(Iscore == 3, trt, trueA)) 
                 
         # Format data
        temp <- format.countdata(data = sim, yvar = "postrelapse_num", timevar = "finalpostdayscount", trtvar = "trt",
                                 xcontinuousvars = c(continuousvars, "postrelapse_num", "offset", "FUweight"),
                                 xcategoricalvars = categoricalvars, imputation.method = NULL)
        input <- data.frame(y = temp$y, trt = factor(temp$trt), time = log(temp$time), temp$x)
        cat("\nA random sample is simulated with seed", seed, "with dimension: ", dim(input), "for the current CV iteration.\n")

        # Create CV folds
        folds <- createFolds(input$trt, k = n.fold, list = TRUE) # Stratified CV, follow the same as the simmain.R where folds were created on input$trt instead of sim$trt
        
        for (fold.i in 1:n.fold){
          testdata <- sim[folds[[fold.i]],]
          # number of methods which succeeded for the given fold/batch. The "is.na(dhat) == FALSE" is to remove methods that didn't produce results for that fold/batch
         
          nr <- nrow(dhat.cvs2 %>% filter(n.subset == n & fold == paste0("fold", fold.i) & iteration == it & is.na(dhat) == FALSE))
          dhat.cvs2$d[which(dhat.cvs2$n.subset == n & dhat.cvs2$fold == paste0("fold", fold.i) & dhat.cvs2$iteration == it & is.na(dhat.cvs2$dhat) == FALSE)] <- rep(testdata$trueA, nr/nrow(testdata))
          stopifnot(nr %% nrow(testdata) == 0)
        }
      } # end of all cv iterations for sample size n
    } # end of all sample sizes
  }
}

## Calculate % accuracy I(dhat == d)/n for each iteration & summary statistics
dat.accuracy <- dhat.cvs2 %>% group_by(beta, percentiles, n.subset, method, fold, iteration) %>% dplyr::summarise(accuracy = sum(dhat == d)/n())
dat.accuracy.summary <- dat.accuracy %>% group_by(beta, percentiles, n.subset, method) %>% dplyr::summarise(mean.acc = mean(accuracy, na.rm = TRUE), 
                                                                                                            sd.acc = sd(accuracy, na.rm = TRUE),
                                                                                                            q1.acc = quantile(accuracy, prob = 0.25, na.rm = TRUE),
                                                                                                            q3.acc = quantile(accuracy, prob = 0.75, na.rm = TRUE))

 
#### ---- Save all formatted data sets in RData  ---- ####
# Tag data sets with level of heterogeneity and percentile
if(all(beta == rep(-0.2, 5))) levhte <- "no"
if(all(beta == c(-1.2, -0.69, 0, 0.1, 0.41))) levhte <- "high"
if(all(beta == c(-0.92, -0.69, 0, 0.1, 0.18))) levhte <- "medium"
if(all(beta == c(-0.36, -0.29, 0, 0.05, 0.1))) levhte <- "low"

if(all(percentiles == seq(0, 1, by = 0.2))) pergr <- "20x5"

cat("\n levhte:", levhte, ", pergr:", pergr, "\n")

# Rename data sets to export with tags
assign(paste0("cvs2_", levhte, "_", pergr), cvs2)
assign(paste0("true.cvs2_", levhte, "_", pergr), true.cvs2)
assign(paste0("all.cvs3_", levhte, "_", pergr), all.cvs3)
assign(paste0("all.cvs3.wide_", levhte, "_", pergr), all.cvs3.wide)
assign(paste0("dhat.agreement_", levhte, "_", pergr), dhat.agreement)
assign(paste0("dat.accuracy_", levhte, "_", pergr), dat.accuracy)
assign(paste0("dat.accuracy.summary_", levhte, "_", pergr), dat.accuracy.summary)

save(file = paste0("./simulations/samplesize/outputs/simsummary_sample_size_dataplot_", levhte, "_", pergr, ".RData"), list = ls()[grepl("_", ls())])
# This output file need to be moved or copied to "./intermediate_results/simulations_samplesize/" to be run by "simulation_analysis.R"

