# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: 05_listdtr.R
# Date: 06NOV2020
# 
# Purpose: Estimate treatment rule with list-based DTR
# 
# Platform: Windows
# R Version: 4.0.3
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   06NOV2020 pj      Start the script 
#   11NOV2020 pj      Investigate the results from the HPC (slow on all data)
#   16NOV2020 pj      Wrap method in a function
#   07DEC2020 gs      Update functions when test data != NULL
#   26APR2021 pj      Change outputs to single Vhat(dhat) to four items: d.hat test fold, Vhat(dhat), d.hat large test, V(dhat)
# ------------------------------------------------------------------

#remove(list = ls())
#library(tidyverse)
#library(magrittr)
#library(listdtr)
#library(fastDummies)

#setwd("/home/pjiang/PMMS/MarketScan/code/")
#setwd("C:/Users/xjiang1/OneDrive - Biogen/Documents/Innovation/Code/PMMS/TruvenMarketScan/")
#source("./01-preprocessing.R")
#source("./02-propensityscore.R")


itrLIST <- function(traindata, testdata, outcome, treatment, categoricalvars, continuousvars,
        maxlen = 2L, seed = seed + fold.i, sim.big = NULL){
  #' ITR with list-based DTR 
  #' Listdtr assumes higher outcome is better but our outcome is assumed to be better if lower
  #' 
  #' @param traindata - training data with both arms; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - outcome variable, assumed to be favorable if lower; string
  #' @param treatment - treatment variable (assume binary); string
  #' @param categoricalvars - Categorical X variables; vector of strings
  #' @param continuousvars - Continuous X variables; vector of strings
  #' @param maxlen - maximum number of nodes (i.e. if-else statements); integer (ends with "L")
  #' @param seed - randomization seed; integer
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
  #' 
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'      OR if testdata == NULL: a list with:
  #'              - score: dataframe with ID and optimal ITR (itr_listDTR)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`        
  
  output <- list()
  
  # Prepare X variables into a matrix (with dummy categorical variables)
  Xtrain <- as.matrix(traindata[, c(categoricalvars, continuousvars)])
  if(is.null(testdata)){
    Xtest <- Xtrain
  } else {
    Xtest <- as.matrix(testdata[, c(categoricalvars, continuousvars)])
  }
    
  # Fit the list-DTR model
  listmod <- listdtr(y = traindata[[outcome]], 
                     a = traindata[[treatment]],
                     x = Xtrain,
                     stage.x = rep(1, ncol(Xtrain)),
                     seed = seed, maxlen = maxlen)
  
  # Optimal ITR in test data
  itr <- predict(listmod, xnew = Xtest, stage = 1)
  
  if(is.null(testdata)){
    output <- list(coef = listmod, score = data.frame(ID = traindata$ID, itr_listDTR = itr))
    colnames(output$score)[2] <- paste0("itr_listDTR", maxlen, "_", outcome)
  } else {
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testdata[['ps']], fu = exp(testdata[['time']]) / 365.25)
  }
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, 
                                 yvar = "postrelapse_num", 
                                 timevar = "finalpostdayscount", 
                                 trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    Xtest.big <- as.matrix(testdata.big[, c(categoricalvars, continuousvars)])
    itr.big <- predict(listmod, xnew = Xtest.big, stage = 1)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }
  
  return(output)
}

