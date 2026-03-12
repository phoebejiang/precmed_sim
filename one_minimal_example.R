# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: one_minimal_example.R
# 
# Purpose: One minimal, non-Slurm representative example to excute simmain.R.
#          Specific configurations: method = Poisson, 
#                                   sample size = 500, 
#                                   batch index = 1, 
#                                   magnitude of heterogeneity = no
#                                   distribution of heterogeneity = equal symmetric
#
#  This script can be run locally without an HPC environment.
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
suppressPackageStartupMessages(library(mpath))
suppressPackageStartupMessages(library(fastDummies))

source("./utility.R")
source("./02-propensityscore.R")
source("./06-LuScore.R")
source("./eachCV.R")
source("./manuscript_func_final.R")

# Constants
n.fold <- 10     # number of folds in each CV iteration
n.cv <- 25       # total number of CV iterations 
base.seed <- 999 # randomization seed
RCT <- T         # randomized trial, if TRUE
big.n <- 1000000  # sample size of the large independent test set to get true value, e.g. 1million
colB <- brewer.pal(n = 9, name = "PuBu") # Blue for estimated value
colR <- brewer.pal(n = 9, name = "OrRd") # Red for value
colG <- brewer.pal(n = 9, name = "BuGn") # Green for level of heterogeneity
display_methods <- c("Poisson", "Boosting", 
                     "Contrast\n Regression", "List DTR\n (2 nodes)")

# User-specified constants
method <- "poisson"               # PM method, an argument from command line
yvar <- "logarr0001"              # y-variable specified by user, could be logarr0001 = log(arr+0.001) or logarr1 = log(arr+1) or log((NRL+0.01)/offsetbyyear)
# batch <- 1                        # the batch index, could be from 1 to n.cv/batch_size, an argument from command line
n <- 500                          # sample size of the randomly generated simulated data
beta <- c(-0.2,-0.2,-0.2,-0.2,-0.2)  # magnitude of heterogeneity is no
cat("\nLevel of heterogeneity beta =", beta, "\n")
percentiles <- seq(0,1,by=0.2)    # distribution of heterogeneity is equal symmetric, i.e., 20% for all subgroups
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

vhats.dhat <- vs.dhat <- dhats <- NULL
for (batch in 1:5){
  
  cat("\n##############################")
  cat("\nExecuting batch No.", batch)
  cat("\n##############################\n")
  
  # Run the CV (based on simmain.R)
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
  
  # Summary outputs (based on simsummary.R)
  vhat.dhat <- v.dhat <- dhat <- data.frame()
  
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
  
  # Combine results over all methods
  vhats.dhat <- rbind(vhats.dhat, vhat.dhat)
  vs.dhat <- rbind(vs.dhat, v.dhat)
  dhats <- rbind(dhats, dhat)

}

vhats.dhat %<>% 
  mutate(method_outcome = paste0(method, "_", outcome))

vs.dhat %<>% 
  mutate(method_outcome = paste0(method, "_", outcome))

dhats %<>% 
  mutate(method_outcome = paste0(method, "_", outcome))

end <- Sys.time()
cat("\nTime elapsed: ", end - start, "s\n\n")

#######################################################################
############################# Calculate Summary ############################
#######################################################################
# Summarize by method
vhats.dhat %<>% 
  group_by(method_outcome) %>% 
  dplyr::summarize(n.batches = n(),
            n.nonnaU = sum(!is.na(U)),
            n.nonnaW = sum(!is.na(W)),
            meanVold = mean(U/W, na.rm = T),
            meanV = sum(U, na.rm = T)/sum(W, na.rm = T),
            sdVold = sd(U/W, na.rm = T), 
            meanU = mean(U, na.rm = T),
            meanW = mean(W, na.rm = T),
            sdV = sum((U / meanW - meanU * W / ((meanW)^2))^2, na.rm = T) / (n.fold * (n.fold * n.cv  - 1)),
            .groups = "keep") %>% 
  dplyr::select(-meanU, -meanW) %>%
  ungroup %>%
  arrange(desc(meanV)) %>%
  mutate(n = n)

vs.dhat %<>% 
  group_by(method_outcome) %>% 
  dplyr::summarize(n.batches = n(),
            n.nonTrueV = sum(!is.na(v.dhat)),
            meanTrueV = mean(v.dhat, na.rm = T),
            sdTrueV = sd(v.dhat, na.rm = T),
            .groups = "keep") %>% 
  ungroup %>%
  arrange(desc(meanTrueV)) %>%
  mutate(n = n)

vhats.dhat %<>% 
  mutate(beta = paste0("c(", paste(beta, collapse = ", "), ")"),
         percentiles = paste0("c(", paste(percentiles, collapse = ", "), ")"))

vs.dhat %<>% 
  mutate(beta = paste0("c(", paste(beta, collapse = ", "), ")"),
         percentiles = paste0("c(", paste(percentiles, collapse = ", "), ")"))

dhats %<>% 
  mutate(beta = paste0("c(", paste(beta, collapse = ", "), ")"),
         percentiles = paste0("c(", paste(percentiles, collapse = ", "), ")"))

# Collect CV iteration batch results and summarize (based on simsummary_sample_size.R)
# Since we only have one sample size, this code chunk only needs to be run onces
cvs <- true.cvs <- dhat.cvs <- NULL
ns <- c(500)
for (n in ns){
  cvs <- rbind(cvs, vhats.dhat %>% 
               mutate(n.subset = n) %>% 
               dplyr::select(n.subset, n.nonnaU, n.nonnaW, method_outcome, 
                             n.batches, meanVold, meanV, sdVold, sdV, 
                             beta, percentiles))

  true.cvs <- rbind(true.cvs, vs.dhat %>% 
                      mutate(n.subset = n) %>% 
                      dplyr::select(n.subset, method_outcome, n.batches, 
                                    n.nonTrueV, meanTrueV, sdTrueV, beta, percentiles))
  
  dhat.cvs <- rbind(dhat.cvs, dhats %>% mutate(n.subset = n))
}

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


#### ---- Cross-validated V.hat and sd plot by sample size for each PM method ---- ####
alldmf <- "All DMF"
allteri <- "All TERI"

method.vec <- c(allteri, alldmf,  
                "Poisson", "Weighted\n Poisson", "Negative\n Binomial", "Weighted\n NegBin", 
                "Linear", "Weighted\n Linear", "dWOLS", 
                "Boosting", "Two\n Regressions", "Contrast\n Regression", "List DTR\n (2 nodes)", "List DTR\n (3 nodes)")


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
  mutate(trueWorstV = trueWorstV, 
         trueV = trueV,
         VR_V.hat = (meanV_V.hat - trueWorstV) / (trueV - trueWorstV),
         VR_V = (meanV_V - trueWorstV) / (trueV - trueWorstV),
         levhte = "No")

# Visualize the example output
# which will show as one dot in Figure 3
plotV(data = all.cvs3.wide, 
       thisV = "meanV_V", 
       display_methods = display_methods, 
       colorgroup = "levhte",
       output = F)
