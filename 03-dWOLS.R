# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: 03_dWOLS.R
# Date: 05OCT2020
# 
# Purpose: Estimate treatment rule with dWOLS
#    
# 
# Platform: Windows
# R Version: 4.0.1
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   05OCT2020 gs      Start the script
#   06OCT2020 pj      Read in preprocessed MarketScan data from a function and clean up
#   28OCT2020 gs      Improve dWOLS application
#   09NOV2020 gs      dWOLS with IPTW
#   12NOV2020 gs      Wrap dWOLS into a function for CV
#   03DEC2020 gs      Update function when test data == NULL
#   11MAR2021 gs      Minor change to accommodate factor trt
#   15MAR2021 pj      Move factor->numeric trt change to eachCV.R
#   26APR2021 pj      Change outputs to single Vhat(dhat) to four items: d.hat test fold, Vhat(dhat), d.hat large test, V(dhat)
#   14MAY2021 pj      Set itr.big as numeric instead of factor bc model.matrix can't handle factors with < 2 levels
# ------------------------------------------------------------------

library(DTRreg)

#source("./utility.R") 

# ------------------------------------------------------------------ #
####                    dWOLS function for CV                     ####
# ------------------------------------------------------------------ #

itrDWOLS <- function(traindata, testdata, outcome, Xoutcome, Xinteraction, dWOLSweight = "IPTW", sim.big = NULL){
  #' ITR based on doubly robust dWOLS
  #' 
  #' @param traindata - training data with both arms; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - assumed to be favorable if lower, continuous for dWOLS; string
  #' @param Xoutcome - covariates to be included as main effects in the outcome model; a vector of strings
  #' @param Xinteraction - covariates with an interaction term with treatment, must be a subset of Xoutcome; a vector of strings
  #' @param weight - weights used in dWOLS; "IPTW" or "overlap"
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
  #' 
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with:
  #'              - fit: dataframe with coefficient of ITR with dWOLS (coef) and SE based on 200 bootstrap
  #'              - score: dataframe with ID, score (score_dwols), optimal ITR (itr_dwols)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`
  
  output <- list()
  
  # Save outcome input
  text.outcome <- outcome
  
  # Specify subset variables
  X.beta <- list(as.formula(paste0("~ ", paste0(Xoutcome, collapse = " + "))))
  X.psi <- list(as.formula(paste0("~ ", paste0(Xinteraction, collapse = " + "))))
  outcome <- as.numeric(traindata[,which(colnames(traindata) == outcome)])
  
  # Calculate SE if testdata == NULL
  var = "none"
  if(is.null(testdata)) var = "bootstrap"
  
  # Fit dWOLS
  if(dWOLSweight == "IPTW"){
    weightiptw <- function(w) 1/w
    mod <- DTRreg(outcome = outcome, blip.mod = X.psi, tf.mod = X.beta, treat.mod = list(trt ~ 1), method = "dwols", data = traindata, treat.mod.man = list(traindata$ps), weight = weightiptw, var.estim = var)
  } else {
    mod <- DTRreg(outcome = outcome, blip.mod = X.psi, tf.mod = X.beta, treat.mod = list(trt ~ 1), method = "dwols", data = traindata, treat.mod.man = list(traindata$ps), var.estim = var)
  }
  
  # Derive optimal treatment in test
  if(is.null(testdata)){
    scoredwols <- model.matrix(X.psi[[1]], data = traindata) %*% mod$psi[[1]]
    itr <- factor(ifelse(scoredwols > 0, 1, 0))
    
    # output score = -scoredwols to match other methods i.e. score > 0 means treatment 1 is better, < 0 means treatment 0 is better
    output <- list(fit = data.frame(coef = mod$psi[[1]], SE = sqrt(diag(mod$covmat[[1]]))), score = data.frame(ID = traindata$ID, score_dwols = -scoredwols, itr_dwols = itr))
    colnames(output$score)[2:3] <- c(paste0("score_dwols_", text.outcome), paste0("itr_dwols_", text.outcome))
  } else {
    itr <- factor(ifelse(model.matrix(X.psi[[1]], data = testdata) %*% mod$psi[[1]] > 0, 1, 0))
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testdata[['ps']], fu = exp(testdata[['time']]) / 365.25)
  }
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    itr.big <- ifelse(model.matrix(X.psi[[1]], data = testdata.big) %*% mod$psi[[1]] > 0, 1, 0)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }
  
  return(output)
}


