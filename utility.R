# ------------------------------------------------------------------
# Product: DMF (Tecfidera) vs. Teriflunomid
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: F_utility.R
# Developer/Programmer: pj
# Date: 12NOV2020
# 
# Purpose: Utility functions for the MarketScan PM implementation
#    
# 
# Platform: Windows
# R Version: 4.0.1
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   12NOV2020 pj      Start the script
#   16NOV2020 gs      Added dWOLS
#                     Change default categorical and continuous vars in dummyXvars and eachCV
#                     Added 4 methods by Lu
#   16NOV2020 pj      Minor edits (e.g. data.matrix, split variables into conti and categ)
#   30NOV2020 pj      Keep postrelapse_num in data regardless of y (format.countdata) and add all_of to avoid warnings
#   23APR2021 pj      Revise getValue() 
#   26APR2021 pj      Add getTrueValue(), getTrueOptimalValue()
#   22NOV2021 pj      Move getvalue functions from setup.R
#   03DEC2021 pj      Update true optimal and true worst value function to accommodate for < 5 groups\
#   10APR2024 pj      Revise getValue() to add sumRi^2
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

#library(tidyverse)
#library(magrittr)
#library(haven)
#library(Hmisc)
#library(ggalt)
#library(caret)
#library(fastDummies)


readSAS <- function(pathname){
  
  #' Read SAS file with given path name and lower column names to lower cases
  #' @param pathname string, path name to the sas7bdat file
  #' @return data.frame 
  
  return(read_sas(pathname) %>% rename_all(tolower))
}

colLabels <- function(df){
  setNames(stack(lapply(df, label))[2:1], c("Varcode", "Variables"))
}

checkComplete <- function(data){
  
  #' Check completeness for all columns in  data.frame
  data %>% summarise_all(funs(100 * mean(!is.na(.)))) %>% t() 
}

propTable <- function(x){
  
  #' Proportion frequency table
  table(x) %>% prop.table()
}


getValue <- function(y, a, d.hat, p.hat, fu){
  
  #' Estimated value function of an estimated ITR, V.hat(d.hat), the weighted approach
  #' 
  #' @param y - Outcome
  #' @param a - Treatment received 0/1
  #' @param d.hat - Estimated optimal treatment 0/1
  #' @param p.hat - Estimated propensity score, P(A=1|X)
  #' @param fu - Follow-up time
  
  pa <- a*p.hat + (1-a)*(1-p.hat) # P(A=a|X)
  u <- y*(a == d.hat)/pa
  w <- (a == d.hat)*fu/pa
  sumU <- sum(y*(a == d.hat)/pa)
  sumW <- sum((a == d.hat)*fu/pa)
  meanU <- mean(y*(a == d.hat)/pa)
  meanW <- mean((a == d.hat)*fu/pa)
  sumRj2 <- sum((u / sumU - w * sumU / sumW / sumW)^2) # to calculate var(vhat)
  sumRj2.mean <- sum((u / meanU - w * meanU / meanW / meanW)^2) # to calculate var(vhat)
  return(list(U = sumU, W = sumW, sumRj2 = sumRj2, sumRj2.mean = sumRj2.mean))
}

randomSample <- function(data, n.subset = 500, seed){
  
  #' Generate a random subset of data
  set.seed(seed)
  subset <- data %>% slice_sample(n = n.subset)
}

dummyXvars <- function(data, 
                       categoricalvars = c("female", "cvd", "diabetes", "prevDMTefficacy", "insuranceplan", "premedicationcosts_cat"),
                       continuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "numSymptoms", "severityScore")){
  
  #' Convert the provided categorical X variables into dummy variables (with the most frequent category as reference group)
  #' and concatenate with continuous X variables and save as output
  #' 
  #' @param data Input data; data.frame
  #' @param categoricalvars Categorical X variables; vector of strings
  #' @param continuousvars Continuous X variables; vector of strings
  #' 
  #' @return X - a matrix of continuous X variables and dummy categorical X variables
  
  dum <- dummy_cols(data[, categoricalvars], remove_most_frequent_dummy = T) %>% dplyr::select(-one_of(categoricalvars))
  X <- data.matrix(cbind(data[, continuousvars], dum))
  
  return(X)
  
}

format.countdata <- function(data, yvar, timevar, trtvar, xcontinuousvars, xcategoricalvars, RCT = F, imputation.method = NULL, scale = FALSE){
  
  # TODO: Update this function with new addition in PMMS/R/utility.R
  
  #' Transform the input count data to the format that the package functions accept
  #'
  #' 1. Rename outcome, time, treatment, covariates to y, time, trt, x
  #' 2. Convert covariates to a matrix with categorical variables turned into dummy variables
  #' 3. Handles missing data in X with specified imputation method; cannot impute outcome or treatment data
  #'
  #' @param data The original input data; data.frame
  #' @param yvar The column name of the count outcome; string
  #' @param timevar The column name of the person-year follow-up time variable; string
  #' @param trtvar The column name of the treatment variable; string
  #' @param xcontinuousvars The column names of the continuous covariates (no formatting done); string or a vector of strings
  #' @param xcategoricalvars The column names of the categorical covariates (to be turned into dummies with the most frequent category as the reference group); string or a vector of strings
  #' @param imputation.method The imputation method could be:
  #'                NULL, if no imputation specified then complete case analysis
  #'                "missforest", imputation on the missing data with random forest
  #' @param scale An indicator of standardizing the continuous variables (not including ps or iptw) to mean 0 and SD 1; boolean
  #'
  #' @return a list of four objects with no missing values:
  #'           y - the formatted and/or imputed count response variable; vector
  #'           time - the formatted and/or imputed  person-year follow-up time variable; vector
  #'           trt - the formatted and/or imputed treatment variable; vector
  #'           x - the formatted and/or imputed covariate variables; matrix

  suppressPackageStartupMessages(library(fastDummies))
  suppressPackageStartupMessages(library(missForest))
  
  data2 <- data %>% rename(y = all_of(yvar), time = all_of(timevar), trt = all_of(trtvar))
  if (yvar == "postrelapse_num") data2$postrelapse_num = data2$y # Keep postrelapse_num in the data always because we need it to calculate getValue
  if (yvar == "inecp") data2$inecp = data2$y
 
  if (is.null(imputation.method)){ 
    if (is.null(xcategoricalvars)) {
      data2 = data2 %>% dplyr::select(y, time, trt, FUweight, offset, postrelapse_num, all_of(xcontinuousvars))  
      } else {
        data2 = data2[, unique(c("y", "time", "trt", xcontinuousvars, xcategoricalvars, "offset", "postrelapse_num"))]
      }
    
  } else if (imputation.method == "missforest"){
    x.imp <- missForest(data2[, c(xcontinuousvars, xcategoricalvars)], verbose = FALSE)$ximp
    data2 = cbind(data2[, c("y", "time", "trt")], x.imp)
  }
  data3 <- data2 %>% drop_na()

  if (!is.null(xcategoricalvars)){
    xcat <- dummy_cols(data3[, xcategoricalvars, drop = FALSE], remove_first_dummy = T, remove_selected_columns = T)
    x <- as.matrix(cbind(data3[, xcontinuousvars, drop = FALSE], xcat))
  } else {
    x <- as.matrix(data3[, xcontinuousvars, drop = FALSE])
  }
  if (scale){
    x <- scale(x %>% dplyr::select(-Fuweight, -offset, -postrelapse_num))
    # apply(x, 2, mean) # check okay, should be all 0
    # apply(x, 2, sd) # check okay, should be all 1 
  }

  x <- cbind(x, data3[, c("offset", "postrelapse_num")])
  
  return(list(y = data3$y, time = data3$time, trt = data3$trt, x = x))
}

convertParameters <- function(magnitude, distribution, verbose){
  
  # Current setting only allows distribution to vary for medium level of magnitude
  if (distribution != "symm 20x5") {
    stopifnot(magnitude == "medium")
  } 
  
  # Parameter Dictionary (mapping between label and value)
  betas <- c("no" = "c(-0.2,-0.2,-0.2,-0.2,-0.2)", 
             "low" = "c(-0.36,-0.29,0,0.05,0.1)", 
             "medium" = "c(-0.92,-0.69,0,0.1,0.18)", 
             "high" = "c(-1.2,-0.69,0,0.1,0.41)")
  percentiless <- c("symm 20x5" = "c(0,0.2,0.4,0.6,0.8,1)",
                    "symm 10-15-50-15-10" = "c(0,0.1,0.25,0.75,0.9,1)", 
                    "symm 10-30-20-30-10" = "c(0,0.1,0.4,0.6,0.9,1)", 
                    "asymm 55-30-15" = "c(0,0.55,0.65,0.75,0.85,1)", 
                    "asymm 55-15-15-15" = "c(0,0.55,0.65,0.7,0.85,1)") 
  
  # Assign parameter values according to their label
  if (distribution != "symm 20x5"){
    beta_group <- c("symm 20x5" = "c(-0.92,-0.69,0,0.1,0.18)", 
                    "symm 10-15-50-15-10" = "c(-0.92,-0.69,0,0.1,0.18)", 
                    "symm 10-30-20-30-10" = "c(-0.92,-0.69,0,0.1,0.18)",
                    "asymm 55-30-15" = "c(-0.92,-0.69,-0.69,-0.69,0)",
                    "asymm 55-15-15-15" = "c(-0.92,-0.69,-0.69,0,0.1)")
    beta <- beta_group[distribution]
    percentiles <- percentiless[distribution]
  } else{
    beta <- betas[magnitude]
    percentiles <- percentiless[distribution]
  }
  
  if (verbose){
    cat("\n* Magnitude of HTE =", beta)
    cat("\n* Distribution of HTE =", percentiles) 
  }
  
  return(list(beta = eval(parse(text = unname(beta))), 
              percentiles = eval(parse(text = unname(percentiles))),
              beta.text = unname(beta),
              percentiles.text = unname(percentiles)))
}
