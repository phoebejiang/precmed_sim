# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: 04_regression-based.R
# Date: 06NOV2020
# 
# Purpose: Estimate treatment rule with regression-based models 
#   (Linear, Poisson, Negative Binomial, Zero-Inflated Negative Binomial)
#    
# 
# Platform: Windows
# R Version: 4.0.3
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   06NOV2020 pj      Start the script based on 04-diagnostics&regression.R: write Poisson into a function and apply CV
#   11NOV2020 pj      Write the Cv part into a main.R function; only keeping the functions
#   16NOV2020 pj      Add lasso penalization to linear and poisson with glmnet
#   24NOV2020 pj      Update data with formatted inputs
#   07DEC2020 gs      Update functions when test data == NULL 
#   29JAN2021 gs      Add LASSO option to itrLinear and itrPoisson
#   26APR2021 pj      Change outputs to single Vhat(dhat) to four items: d.hat test fold, Vhat(dhat), d.hat large test, V(dhat)
# ------------------------------------------------------------------

#remove(list = ls())
#library(tidyverse)
#library(magrittr)
#library(listdtr)
#library(fastDummies)
#library(MASS)
#library(pscl)
#library(caret)
#library(glmnet)
#library(mpath)

#source("./01-preprocessing.R")
#source("./02-propensityscore.R")
#source("./utility.R")

# TODO: add all variables with lasso penalization 

itrLinear <- function(traindata1, traindata0, testdata, outcome, categoricalvars, continuousvars, testps, LASSO = F, trainweight1 = NULL, trainweight0 = NULL, sim.big = NULL){
  #' Regression-based two-sample linear regression ITR
  #' add LASSO penalization if number of variables exceeds 20
  #' 
  #' @param traindata1 - training data but with only arm = 1; data.frame
  #' @param traindata0 - training data but with only arm = 0; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - outcome variable, assumed to be favorable if lower; string
  #' @param categoricalvars - categorical X variables to be included in the linear regression; vector of strings
  #' @param continuousvars - continuous X variables to be included in the linear regression; vector of strings
  #' @param testps - estimated propensity score for the test data; vector of size nrow(testdata)
  #' @param LASSO - whether LASSO penalization should be used; boolean
  #' @param trainweight1 - prior weight to be used in the fitting process of traindata1; vector of size nrow(traindata1), or NULL if no weighting
  #' @param trainweight0 - prior weight to be used in the fitting process of traindata0; vector of size nrow(traindata0), or NULL if no weighting
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
  #' 
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with:
  #'              - fit: dataframe with coefficients in trt=1 (coef1) and trt=0 (coef0)
  #'              - score: dataframe with ID, score (score_linear), optimal ITR (itr_linear)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`              
  
  output <- list()
    
  # Save traindata ID 
  ID <- c(traindata0$ID, traindata1$ID)
  
  variables = c(categoricalvars, continuousvars)
  if (LASSO == F){
    # Linear model with pre-specified variables
    mod1 <- lm(as.formula(paste0(outcome, " ~ ", paste0(variables, collapse = " + "))), 
                weight = trainweight1, data = traindata1)
    mod0 <- lm(as.formula(paste0(outcome, " ~ ", paste0(variables, collapse = " + "))), 
                weight = trainweight0, data = traindata0)
    
    # Predictions
    if(is.null(testdata)){
      traindata <- rbind(traindata0, traindata1)
      pred1 <- predict(mod1, traindata)
      pred0 <- predict(mod0, traindata)
      coef1 <- mod1$coefficients
      coef0 <- mod0$coefficients
    } else {
      pred1 <- predict(mod1, testdata)
      pred0 <- predict(mod0, testdata)
    }
  } else {
    # Prepare X variables into a matrix (with dummy categorical variables)
    X1 <- as.matrix(traindata1 %>% dplyr::select(all_of(variables)))
    X0 <- as.matrix(traindata0 %>% dplyr::select(all_of(variables)))

    # Find the best lambda with 10-fold CV
    cvfit1 <- cv.glmnet(x = X1, y = traindata1[[outcome]], family = "gaussian", weights = trainweight1, alpha = 1, data = traindata1)
    cvfit0 <- cv.glmnet(x = X0, y = traindata0[[outcome]], family = "gaussian", weights = trainweight0, alpha = 1, data = traindata0)
    
    # Train the model with tuned lambda
    mod1 <- glmnet(x = X1, y = traindata1[[outcome]], family = "gaussian", weights = trainweight1, alpha = 1, lambda = cvfit1$lambda.min, data = traindata1)
    mod0 <- glmnet(x = X0, y = traindata0[[outcome]], family = "gaussian", weights = trainweight0, alpha = 1, lambda = cvfit0$lambda.min, data = traindata0)
    
    # Predictions
    if(is.null(testdata)){
      X <- rbind(X0, X1)
      pred1 <- predict(mod1, X)
      pred0 <- predict(mod0, X)
      coef1 <- mod1$beta
      coef0 <- mod0$beta
      
    } else {
      Xtest <- as.matrix(testdata %>% dplyr::select(all_of(variables)))
      pred1 <- predict(mod1, Xtest)
      pred0 <- predict(mod0, Xtest)
    }
  }
  
  itr <- as.numeric(pred1 <= pred0)
  
  if (is.null(testdata)){
    # Output score such that score > 0 means treatment 0 is preferred, < 0 means treatment 1 is preferred
    score <- as.numeric(pred1 - pred0)
    output <- list(fit = data.frame(coef1, coef0), score = data.frame(ID = ID, score_linear = score, itr_linear = itr))
    colnames(output$score)[2:3] <- c(paste0("score_linear_", outcome), paste0("itr_linear_", outcome))
  } else{
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testps, fu = exp(testdata[['time']]) / 365.25)
  } # lower value is better

  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    if (LASSO == F) {
      pred1.big <- predict(mod1, testdata.big)
      pred0.big <- predict(mod0, testdata.big)
    } else {
      Xtest.big <- as.matrix(testdata.big %>% dplyr::select(all_of(variables)))
      pred1.big <- predict(mod1, Xtest.big)
      pred0.big <- predict(mod0, Xtest.big)
    }
    itr.big <- as.numeric(pred1.big <= pred0.big)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }

  return(output) 

}

itrPoisson <- function(traindata1, traindata0, testdata, outcome, offset, categoricalvars, continuousvars, testps, LASSO = F, trainweight1 = NULL, trainweight0 = NULL, sim.big = NULL){
  #' Regression-based two-sample Poisson ITR
  #' add LASSO penalization if number of variables exceeds 20
  #' 
  #' @param traindata1 - training data but with only arm = 1; data.frame
  #' @param traindata0 - training data but with only arm = 0; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - outcome variable, assumed to be favorable if lower; string
  #' @param offset - offset variable in Poisson regression; string
  #' @param categoricalvars - Categorical X variables to be included in the Poisson regression; vector of strings
  #' @param continuousvars - Continuous X variables to be included in the Poisson regression; vector of strings
  #' @param testps - estimated propensity score for the test data; vector of size nrow(testdata)
  #' @param LASSO - whether LASSO penalization should be used; boolean
  #' @param trainweight1 - prior weight to be used in the fitting process of traindata1; vector of size nrow(traindata1), or NULL if no weighting
  #' @param trainweight0 - prior weight to be used in the fitting process of traindata0; vector of size nrow(traindata0), or NULL if no weighting
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
  #' 
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with:
  #'              - fit: dataframe with coefficients in trt=1 (coef1) and trt=0 (coef0)
  #'              - score: dataframe with ID, score (score_poisson), optimal ITR (itr_poisson)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`              
  
  output <- list()
  
  # Save traindata ID 
  ID <- c(traindata0$ID, traindata1$ID)
  
  variables = c(categoricalvars, continuousvars)
  if (LASSO == F){
    # Poisson regresion with pre-specified variables
    mod1 <- glm(as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + "))), 
                 family = "poisson", weight = trainweight1, data = traindata1)
    mod0 <- glm(as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + "))), 
                 family = "poisson", weight = trainweight0, data = traindata0)
    
    # Predictions
    if(is.null(testdata)){
      traindata <- rbind(traindata0, traindata1)
      pred1 <- predict(mod1, traindata)
      pred0 <- predict(mod0, traindata)
      coef1 <- mod1$coefficients
      coef0 <- mod0$coefficients
    } else {
      pred1 <- predict(mod1, testdata)
      pred0 <- predict(mod0, testdata)
    }
  } else {
    # Prepare X variables into a matrix (with dummy categorical variables)
    X1 <- as.matrix(traindata1 %>% dplyr::select(all_of(variables)))
    X0 <- as.matrix(traindata0 %>% dplyr::select(all_of(variables)))
    
    # Find the best lambda with 10-fold CV
    cvfit1 <- cv.glmnet(x = X1, y = traindata1[[outcome]], family = "poisson", offset = traindata1[["offset"]], weights = trainweight1, alpha = 1, data = traindata1)
    cvfit0 <- cv.glmnet(x = X0, y = traindata0[[outcome]], family = "poisson", offset = traindata0[["offset"]], weights = trainweight0, alpha = 1, data = traindata0)
    
    # Train the model with tuned lambda
    mod1 <- glmnet(x = X1, y = traindata1[[outcome]], family = "poisson", offset = traindata1[["offset"]], weights = trainweight1, alpha = 1, lambda = cvfit1$lambda.min, data = traindata1)
    mod0 <- glmnet(x = X0, y = traindata0[[outcome]], family = "poisson", offset = traindata0[["offset"]], weights = trainweight0, alpha = 1, lambda = cvfit0$lambda.min, data = traindata0)
    
    # Predictions 
    if(is.null(testdata)){
      X <- rbind(X0, X1)
      pred1 <- predict(mod1, X, newoffset = traindata[['offset']])
      pred0 <- predict(mod0, X, newoffset = traindata[['offset']])
      coef1 <- mod1$beta
      coef0 <- mod0$beta
    } else {
      Xtest <- as.matrix(testdata %>% dplyr::select(all_of(variables)))
      pred1 <- predict(mod1, Xtest, newoffset = testdata[['offset']])
      pred0 <- predict(mod0, Xtest, newoffset = testdata[['offset']])
    }
  }

  itr <- as.numeric(pred1 <= pred0)
  
  if(is.null(testdata)){
    # Output score such that score > 0 means treatment 0 is preferred, < 0 means treatment 1 is preferred
    score <- as.numeric(pred1 - pred0)
    output <- list(fit = data.frame(coef1, coef0), score = data.frame(ID = ID, score_poisson = score, itr_poisson = itr))
  } else{
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testps, fu = exp(testdata[['time']]) / 365.25)
  } # lower value is better
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    if (LASSO == F) {
      pred1 <- predict(mod1, testdata.big)
      pred0 <- predict(mod0, testdata.big)
    } else {
      Xtest.big <- as.matrix(testdata.big %>% dplyr::select(all_of(variables)))
      pred1.big <- predict(mod1, Xtest.big, newoffset = testdata.big[['offset']])
      pred0.big <- predict(mod0, Xtest.big, newoffset = testdata.big[['offset']])
    }
    itr.big <- as.numeric(pred1.big <= pred0.big)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }
  
  return(output)
}

itrNegBin <- function(traindata1, traindata0, testdata, outcome, offset, categoricalvars, continuousvars, testps, LASSO = F, trainweight1 = NULL, trainweight0 = NULL, sim.big = NULL){
  #' Regression-based two-sample Negative Binomial ITR
  #' add LASSO penalization if number of variables exceeds 20
  #' 
  #' @param traindata1 - training data but with only arm = 1; data.frame
  #' @param traindata0 - training data but with only arm = 0; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - outcome variable, assumed to be favorable if lower; string
  #' @param offset - offset variable in Negative Binomial regression; string
  #' @param categoricalvars - Categorical X variables to be included in the NegBin regression; vector of strings
  #' @param continuousvars - Continuous X variables to be included in the NegBin regression; vector of strings
  #' @param testps - estimated propensity score for the test data; vector of size nrow(testdata)
  #' @param LASSO - whether LASSO penalization should be used; boolean
  #' @param trainweight1 - prior weight to be used in the fitting process of traindata1; vector of size nrow(traindata1), or NULL if no weighting
  #' @param trainweight0 - prior weight to be used in the fitting process of traindata0; vector of size nrow(traindata0), or NULL if no weighting
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
  #'
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with:
  #'              - fit: dataframe with coefficients in trt=1 (coef1) and trt=0 (coef0). If LASSO is used, coefficients for lambda=100 are returned (default)
  #'              - score: dataframe with ID, score (score_NB), optimal ITR (itr_NB)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`              
  
  output <- list()
  
  # Save traindata ID 
  ID <- c(traindata0$ID, traindata1$ID)
  
  variables = c(categoricalvars, continuousvars)
  if (LASSO == F){
    # Negative binomial regression with pre-specified variables
    mod1 <- glm.nb(as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + "))), 
                 weights = trainweight1, data = traindata1,
                 maxit = 500)
    mod0 <- glm.nb(as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + "))), 
                 weights = trainweight0, data = traindata0,
                 maxit = 500)
    
    # Predictions
    if(is.null(testdata)){
      traindata <- rbind(traindata0, traindata1)
      pred1 <- predict(mod1, traindata)
      pred0 <- predict(mod0, traindata)
      coef1 <- mod1$coefficients
      coef0 <- mod0$coefficients
    } else {
      pred1 <- predict(mod1, testdata)
      pred0 <- predict(mod0, testdata)
    }
  } else {                                                                                                                                                                 
    ## Default is LASSO and we try 100 lambda values (100 is default)
    formula <- as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + ")))
    mod1 <- glmregNB(formula = formula, data = traindata1, weights = trainweight1, nlambda = 100, alpha = 1)
    mod0 <- glmregNB(formula = formula, data = traindata0, weights = trainweight0, nlambda = 100, alpha = 1)
    
    # ## Find the best lambda with 10-fold CV
    cvfit1 <- cv.glmregNB(formula = formula, data = traindata1, weights = trainweight1, alpha = 1)
    cvfit0 <- cv.glmregNB(formula = formula, data = traindata0, weights = trainweight0, alpha = 1)
    # 
    # TODO: 
    # Error in { : 
    # task 70 failed - "NA/NaN/Inf in foreign function call (arg 1)"
    # In addition: Warning message:
    #   In model.matrix.default(Terms, mf, contrasts) :
    #   non-list contrasts argument ignored
    #
    # GS: got error message only for cvfit1, worked with warnings when formula specified as a function of continuousvars only (maybe problem is with low-prevalent binary variables?)
    
    # ## Make predictions on the test data
    # pred1 <- predict(object = mod1, newx = testdata, which = cvfit1$lambda.which) # use the lambda of the best model
    # pred0 <- predict(object = mod0, newx = testdata, which = cvfit0$lambda.which)
    
    if(is.null(testdata)){
      traindata <- rbind(traindata0, traindata1)
      pred1 <- predict(object = mod1, newx = traindata, newoffset = traindata[['offset']], which = which.min(BIC(mod1))) # use the lambda of the best model
      pred0 <- predict(object = mod0, newx = traindata, newoffset = traindata[['offset']], which = which.min(BIC(mod0)))
      # TODO: return coefficients from best lambda
      coef1 <- mod1$coefficients # save coefficients with lambda=100
      coef0 <- mod0$coefficients
    } else {
      pred1 <- predict(object = mod1, newx = testdata, newoffset = testdata[['offset']], which = which.min(BIC(mod1))) # use the lambda of the best model
      pred0 <- predict(object = mod0, newx = testdata, newoffset = testdata[['offset']], which = which.min(BIC(mod0)))
    }
  }
  
  itr <- as.numeric(pred1 <= pred0)
  
  if(is.null(testdata)){
    # Output score such that score > 0 means treatment 0 is preferred, < 0 means treatment 1 is preferred
    score <- as.numeric(pred1 - pred0)
    output <- list(fit = data.frame(coef1, coef0), score = data.frame(ID = ID, score_NB = score, itr_NB = itr))
  } else{
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testps, fu = exp(testdata[['time']]) / 365.25)
  } # lower value is better
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    if (LASSO == FALSE) {
      pred1.big <- predict(mod1, testdata.big)
      pred0.big <- predict(mod0, testdata.big)
    } else {
      pred1.big <- predict(mod1, newx = testdata.big, newoffset = testdata.big[['offset']], which = which.min(BIC(mod1)))
      pred0.big <- predict(mod0, newx = testdata.big, newoffset = testdata.big[['offset']], which = which.min(BIC(mod0)))
    }
    itr.big <- as.numeric(pred1.big <= pred0.big)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }
  
  return(output)
}

itrZeroInfl <- function(traindata1, traindata0, testdata, outcome, offset, categoricalvars, continuousvars, testps, trainweight1 = NULL, trainweight0 = NULL, sim.big = NULL){
  #' Regression-based two-sample Zero-Inflated Negative Binomial ITR
  #' add LASSO penalization if number of variables exceeds 20
  #' 
  #' @param traindata1 - training data but with only arm = 1; data.frame
  #' @param traindata0 - training data but with only arm = 0; data.frame
  #' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
  #' @param outcome - outcome variable, assumed to be favorable if lower; string
  #' @param offset - offset variable in Negative Binomial regression; string
  #' @param categoricalvars - Categorical X variables to be included in the zero-inflated regression; vector of strings
  #' @param continuousvars - Continuous X variables to be included in the zero-inflated regression; vector of strings
  #' @param testps - estimated propensity score for the test data; vector of size nrow(testdata)
  #' @param trainweight1 - prior weight to be used in the fitting process of traindata1; vector of size nrow(traindata1), or NULL if no weighting
  #' @param trainweight0 - prior weight to be used in the fitting process of traindata0; vector of size nrow(traindata0), or NULL if no weighting
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
  #'
  #' @return if testdata != NULL: a list with:
  #'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
  #'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
  #'     OR if testdata == NULL: a list with:
  #'              - score: dataframe with ID, score (score_zeroinf), optimal ITR (itr_zeroinf)
  #' @return In addition, if sim.big != NULL: the list also contains:
  #'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
  #'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`              
  
  output <- list()
  
  # Save traindata ID 
  ID <- c(traindata0$ID, traindata1$ID)
  
  variables = c(categoricalvars, continuousvars)
  if (length(variables) <= 20){
    # No penalization if fewer than 20 variables
    mod1 <- zeroinfl(as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + "))), 
                   dist = "negbin", weights = trainweight1, data = traindata1)
    mod0 <- zeroinfl(as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + "))), 
                   dist = "negbin", weights = trainweight0, data = traindata0)
      # The formula can be used to specify both components of the model: 
      # If a formula of type y ~ x1 + x2 is supplied, then the same regressors are employed in both components. 
      # This is equivalent to y ~ x1 + x2 | x1 + x2. 
      # A different set of regressors could be specified for the count and zero-inflation component, 
      # e.g., y ~ x1 + x2 | z1 + z2 + z3 giving the count data model y ~ x1 + x2 conditional on (|) the zero-inflation model y ~ z1 + z2 + z3.
    
    if(is.null(testdata)){
      traindata <- rbind(traindata0, traindata1)
      pred1 <- predict(mod1, traindata)
      pred0 <- predict(mod0, traindata)
    } else {
      pred1 <- predict(mod1, testdata)
      pred0 <- predict(mod0, testdata)
    }
  } else {
    # LASSO if more than 20 variables
    # Ref: https://cran.r-project.org/web/packages/mpath/vignettes/static_german.pdf
    ## Default is LASSO and we try 10 lambda values
    formula <- as.formula(paste0(outcome, " ~ offset(", offset, ") + ", paste0(variables, collapse = " + ")))
    mod1 <- zipath(formula = formula, data = traindata1, family = "negbin",  weights = trainweight1, nlambda = 100)
    mod0 <- zipath(formula = formula, data = traindata0, family = "negbin",  weights = trainweight0, nlambda = 100)
  
    # ## Find the best lambda with 10-fold CV
    # cvfit1 <- cv.zipath(formula = formula, data = traindata1, family = "negbin", weights = trainweight1, nlambda = 10)
    # cvfit0 <- cv.zipath(formula = formula, data = traindata0, family = "negbin", weights = trainweight0, nlambda = 10)
    # 
    # TODO: Error in cv.zipath(formula = formula, data = traindata0, family = "negbin",  : 
    #  argument "X" is missing, with no default
    # ## Make predictions on the test data
    # pred1 <- predict(object = mod1, newdata = testdata, which = cvfit1$lambda.which) # use the lambda of the best model
    # pred0 <- predict(object = mod0, newdata = testdata, which = cvfit0$lambda.which)
    
    if(is.null(testdata)){
      traindata <- rbind(traindata0, traindata1)
      pred1 <- predict(object = mod1, newdata = traindata, which = which.min(BIC(mod1))) # use the lambda of the best model
      pred0 <- predict(object = mod0, newdata = traindata, which = which.min(BIC(mod0)))
    } else {
      pred1 <- predict(object = mod1, newdata = testdata, which = which.min(BIC(mod1))) # use the lambda of the best model
      pred0 <- predict(object = mod0, newdata = testdata, which = which.min(BIC(mod0)))
    }
  }
  
  itr <- as.numeric(pred1 <= pred0)
  
  if(is.null(testdata)){
    # Output score such that score > 0 means treatment 0 is preferred, < 0 means treatment 1 is preferred
    score <- as.numeric(pred1 - pred0)
    # TODO (maybe?): return coefficients needed to estimate ITR
    output <- list(score = data.frame(ID = ID, score_zeroinf = score, itr_zeroinf = itr))
  } else{
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testps, fu = testdata[['finalpostdayscount']] / 365.25)
  } # lower value is better
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    pred1.big <- predict(mod1, newdata = testdata.big, which = which.min(BIC(mod1)))
    pred0.big <- predict(mod0, newdata = testdata.big, which = which.min(BIC(mod0)))
    itr.big <- as.numeric(pred1.big <= pred0.big)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }
  
  return(output) 
  
}

