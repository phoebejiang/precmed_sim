# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: 06_LuScore.R
# Date: 16NOV2020
# 
# Purpose: Estimate treatment rule with 4 methods in Yadlowsky (2020)
#    
# 
# Platform: Windows
# R Version: 4.0.1
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   16NOV2020 gs      Start the script
#   03DEC2020 gs      Update functions when test data==NULL
#   07DEC2020 gs      Add more outputs with testdata==NULL
#   11MAR2021 gs      Update DR function to accommodate RCT data and trt factor or numeric
#                     Add -time in boosting
#   15MAR2021 pj      Move factor->numeric trt change to eachCV.R
#   26APR2021 pj      Change outputs to single Vhat(dhat) to four items: d.hat test fold, Vhat(dhat), d.hat large test, V(dhat)
#   14MAY2021 pj      Set itr.big as numeric instead of factor bc model.matrix can't handle factors with < 2 levels
# ------------------------------------------------------------------

#library(gbm)

#### Lu's basic functions to be used in the ITR functions below
####  see file R package "PRECMED", CATE.R script for details
twoarmglmcount.dr <- function(y, x, time, trt, ps, f1.predictor, f0.predictor, error.control = 1e-3, max.iter = 150, tune = c(0.5, 2.0)){
  y <- y*exp(-time)
  x.aug <- cbind(1, x)
  p.aug <- length(x.aug[1, ])
  
  # Estimate coefficients with NR algorithm
  beta <- rep(0, p.aug)
  mae <- Inf
  iter <- 0
  epsilon.min <- Inf
  while(mae > error.control && iter <= max.iter && epsilon.min > 0){
    eta <- as.numeric(exp(x.aug %*% beta))
    error <- (trt*(y - eta*f0.predictor/2 - f1.predictor/2)*(1 - ps) - (1 - trt)*(eta*y - eta*f0.predictor/2 - f1.predictor/2)*ps)/(eta*ps + (1 - ps))
    score <- colSums(x.aug*error)
    slopewt <- (y + f0.predictor*(trt/ps - 1)/2 + f1.predictor*((1 - trt)/(1 - ps) - 1)/2)*eta*ps*(1 - ps)/(eta*ps + (1 - ps))^2
    slope <- t(x.aug*slopewt) %*% x.aug
    epsilon.min <- eigen(slope)$value[p.aug]
    if(iter == 0) epsilon0 <- epsilon.min + epsilon.min*(epsilon.min < 0) # fixed to all iterations
    beta <- beta + solve(slope + diag(tune[2]*abs(epsilon0), p.aug, p.aug)) %*% score*tune[1] # adding the diagonal matrix to slove potential singualrity issues
    mae <- sum(abs(score))
    iter <- iter + 1
  }
  
  converge1 <- 1*(mae <= error.control)
  converge2 <- 0
  
  # If NP did not converge, solve for minimizing the L2-norm (sum of squares) of the score function to avoid inversing the slope inverse (generally slower than NP)
  if(converge1 == 0){
    lossf <- function(beta){
      eta <- as.numeric(exp(x.aug %*% beta))
      eta.max <- max(eta)
      error <- (trt*(y - eta*f0.predictor/2 - f1.predictor/2)*(1 - ps) - (1 - trt)*(eta*y - eta*f0.predictor/2 - f1.predictor/2)*ps)/(eta*ps + (1 - ps))
      score <- colSums(x.aug*error)
      return(sum(score^2))
    }
    
    initial.value <- lossf(rep(0, p.aug))
    fit <- optim(rep(0, p.aug), fn = lossf) #, control=list(trace=T))
    beta <- fit$par
    converge2 <- 1*(abs(fit$value) < initial.value/100)
  }
  
  beta <- as.vector(beta)
  eta <- as.numeric(exp(x.aug %*% beta))
  error <- (trt*(y - eta*f0.predictor/2 - f1.predictor/2)*(1 - ps) - (1 - trt)*(eta*y - eta*f0.predictor/2 - f1.predictor/2)*ps)/(eta*ps + (1 - ps))
  score <- colSums(x.aug*error)
  slopewt <- (y + f0.predictor*(trt/ps - 1)/2 + f1.predictor*((1 - trt)/(1 - ps) - 1)/2)*eta*ps*(1 - ps)/(eta*ps + (1 - ps))^2
  slope <- t(x.aug*slopewt) %*% x.aug
  sigma <- solve(slope) %*% (t(x.aug*error^2) %*% x.aug) %*% solve(slope)
  
  return(list(coef = beta, vcov = sigma, converge = (converge1 + converge2 > 0)))
}

onearmglmcount.dr <- function(y, x, time, trt, ps, f.predictor){
  f.predictor <- as.vector(f.predictor)
  y <- y*exp(-time)
  
  withCallingHandlers({
    fit <- glm(y ~ log(f.predictor) + x, family = "poisson", weight = trt/ps)
    yhat <- exp(cbind(1, log(f.predictor), x) %*% fit$coef)
    fit2 <- glm(yhat ~ x, family = "poisson")
  },
  warning = function(w) { # don't change the = to <- in withCallingHandlers
    if (grepl("non-integer", conditionMessage(w)))
      invokeRestart("muffleWarning") # suppress warnings in glm(): "In dpois(y, mu, log = TRUE) : non-integer x = 0.557886."
  })
  
  return(fit2$coef)
}

##### ITR function with 4 methods Lu

#' ITR based on Poisson regression
#' 
#' @param traindata1 - training data but with only arm = 1; data.frame
#' @param traindata0 - training data but with only arm = 0; data.frame
#' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
#' @param categoricalvars - Categorical X variables to be included in the zero-inflated regression; vector of strings
#' @param continuousvars - Continuous X variables to be included in the zero-inflated regression; vector of string
#' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
#' 
#' @return if testdata != NULL: a list with:
#'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
#'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
#'     OR if testdata == NULL: a list with:
#'              - fit: dataframe with coefficient of Poisson fit trt=1 (coef1) and trt=0 (coef0), and corresponding SE (SE1 and SE0)
#'              - score: dataframe with ID, score (score_pois), optimal ITR (itr_pois)
#' @return In addition, if sim.big != NULL: the list also contains:
#'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
#'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`  

itrLuPoisson <- function(traindata1, traindata0, testdata, categoricalvars, continuousvars, sim.big = NULL){
  
  output <- list()
  
  # Save traindata ID 
  ID <- c(traindata0$ID, traindata1$ID)
  
  # Separate into training trt = 1, trt = 0, remove extra variables
  traindata1 <- traindata1 %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars), y, time)
  traindata0 <- traindata0 %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars), y, time)
  
  # Fit Poisson by treatment arm
  fit1 <- glm(y ~. - time + offset(time), family = "poisson", data = traindata1)
  fit0 <- glm(y ~. - time + offset(time), family = "poisson", data = traindata0)
  beta1.ini <- fit1$coef
  beta0.ini <- fit0$coef
  delta2 <- beta1.ini - beta0.ini
  
  # Predict score
  if(is.null(testdata)){
    traindata <- rbind(traindata0, traindata1)
    xtot <- traindata %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
  } else{
    xtot <- testdata %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
  }
  xtot <- as.matrix(cbind(1, xtot))
  scorepois <- as.numeric(as.matrix(xtot) %*% delta2)
  
  # Optimal treatment: scorepois > 0 means treatment 0 is preferred, pois < 0 means treatment 1 is prefered
  itr <- factor(ifelse(scorepois >= 0, 0, 1))
 
  if(is.null(testdata)){
    output <- list(fit = data.frame(coef1 = beta1.ini, SE1 = summary(fit1)$coef[,2], coef0 = beta0.ini, SE0 = summary(fit0)$coef[,2]), score = data.frame(ID = ID, score_lupois = scorepois, itr_lupois = itr))
  } else {
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testdata[['ps']], fu = exp(testdata[['time']]) / 365.25)
  }
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    xtot.big <- testdata.big %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
    xtot.big <- as.matrix(cbind(1, xtot.big))
    scorepois.big <- as.numeric(as.matrix(xtot.big) %*% delta2)
    itr.big <- ifelse(scorepois.big >= 0, 0, 1)
    output$dhat.big <- itr.big
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }
  
  return(output)
}

  
#' ITR based on boosting
#' 
#' @param traindata1 - training data but with only arm = 1; data.frame
#' @param traindata0 - training data but with only arm = 0; data.frame
#' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
#' @param categoricalvars - Categorical X variables to be included in the zero-inflated regression; vector of strings
#' @param continuousvars - Continuous X variables to be included in the zero-inflated regression; vector of strings
#' @param tree.depth Depth of individual trees in boosting (usually 2-3); integer
#' @param n.trees Maximum number of trees in boosting (usually 100-1000); integer
#' @param plot.gbmperf Plot the performance measures in the GBM method; boolean
#' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
#'
#' @return if testdata != NULL: a list with:
#'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
#'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
#'     OR if testdata == NULL: a list with:
#'              score: dataframe with ID, score (score_boost), optimal ITR (itr_boost)
#' @return In addition, if sim.big != NULL: the list also contains:
#'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
#'              - v.dhat: true value of the estimated ITR of every observation in `sim.big` 

itrLuBoosting <- function(traindata0, traindata1, testdata, categoricalvars, continuousvars, tree.depth = 2, n.trees = 200, plot.gbmperf = F, sim.big = NULL){
  
  output <- list()
  
  # Save traindata ID 
  ID <- c(traindata0$ID, traindata1$ID)
  
  # Separate into training trt = 1, trt = 0, remove extra variables
  traindata1 <- traindata1 %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars), y, time)
  traindata0 <- traindata0 %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars), y, time)
  
  # Apply boosting to training data
  set.seed(100)
  fit1.gbm <- gbm(y ~. - time + offset(time), data = traindata1, distribution = "poisson", interaction.depth = tree.depth, n.trees = n.trees, cv.folds = 5)
  best1.iter <- max(10, gbm.perf(fit1.gbm, method = "cv", plot.it = plot.gbmperf))

  set.seed(100)
  fit0.gbm <- gbm(y ~. - time + offset(time), data = traindata0, distribution = "poisson", interaction.depth = tree.depth, n.trees = n.trees, cv.folds = 5)
  best0.iter <- max(10, gbm.perf(fit0.gbm, method = "cv", plot.it = plot.gbmperf))
  
  # Predict score in test data
  if(is.null(testdata)){
    traindata <- rbind(traindata0, traindata1)
    withCallingHandlers({
      predict0 <- predict(object = fit0.gbm, newdata = traindata, n.trees = best0.iter)
      predict1 <- predict(object = fit1.gbm, newdata = traindata, n.trees = best1.iter)
    },
    warning=function(w) {
      if (grepl("does not add the offset", conditionMessage(w)))
        invokeRestart("muffleWarning")  
    })
  } else {
    withCallingHandlers({
      predict0 <- predict(object = fit0.gbm, newdata = testdata, n.trees = best0.iter)
      predict1 <- predict(object = fit1.gbm, newdata = testdata, n.trees = best1.iter)
    },
    warning=function(w) {
      if (grepl("does not add the offset", conditionMessage(w)))
        invokeRestart("muffleWarning")  
    })
  }
 
  scoreboost <- predict1 - predict0 # log(CATE1) - log(CATE0)
  
  # Optimal treatment: scoreboost > 0 means treatment 0 is preferred, scoreboost < 0 means treatment 1 is prefered
  itr <- factor(ifelse(scoreboost >= 0, 0, 1))
  
  if(is.null(testdata)){
    output <- list(score = data.frame(ID = ID, score_boost = scoreboost, itr_boost = itr))
  } else {
    output$dhat <- itr
    output$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr, p.hat = testdata[['ps']], fu = exp(testdata[['time']]) / 365.25)
  }
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    withCallingHandlers({
      predict0.big <- predict(object = fit0.gbm, newdata = testdata.big, n.trees = best0.iter)
      predict1.big <- predict(object = fit1.gbm, newdata = testdata.big, n.trees = best1.iter)
    },
    warning=function(w) {
      if (grepl("does not add the offset", conditionMessage(w)))
        invokeRestart("muffleWarning")  
    })
    scoreboost.big <- predict1.big - predict0.big
    itr.big <- ifelse(scoreboost.big >= 0, 0, 1)
    output$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr.big, betas = sim.big$betas)
  }
  
  return(output)
}


#' ITR based on two regressions and contrast regression
#' 
#' @param traindata - training data with both arms; data.frame
#' @param testdata - test data with both arms (assume two arms) and same variables as traindata; data.frame
#' @param categoricalvars - Categorical X variables to be included in the zero-inflated regression; vector of strings
#' @param continuousvars - Continuous X variables to be included in the zero-inflated regression; vector of strings
#' @param RCT Whether treatment is randomized. If RCT=T, the PS is the proportion of patients treated with DMF
#' @param tree.depth Depth of individual trees in boosting (usually 2-3); integer
#' @param n.trees Maximum number of trees in boosting (usually 100-1000); integer
#' @param Kfold Number of folds (parts) used in cross-fitting to partition the data; integer
#' @param B Number of time cross-fitting is repeated to reduce Monte Carlo variability; integer
#' @param seed.cf Randomization seed for cross-fitting partitions; integer
#' @param plot.gbmperf Plot the performance measures in the GBM method; boolean
#' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
#' 
#' @return if testdata != NULL: a list with two elements, tworeg and contrast reg, each with:
#'              - dhat: the estimated ITR of every observation in `testdata` (1 means arm 1 and 0 means arm 0)
#'              - vhat.dhat: estimated value of the estimated ITR of every observation in `testdata`
#'     OR if testdata == NULL: a list with two elements, tworeg and contrast reg, each with:
#'               - fit: coefficient of log(CATE) (coef) and SE (for contrast reg only)
#'               - score: dataframe with ID, score (score_tworeg or score_contrastreg), optimal ITR (itr_tworeg or itr_contrastreg)
#' @return In addition, if sim.big != NULL: a list with two elements, tworeg and contrast reg, each with:
#'              - dhat.big: the estimated ITR of every observation in `sim.big` (1 means arm 1 and 0 means arm 0)
#'              - v.dhat: true value of the estimated ITR of every observation in `sim.big`  
                       
itrLuDR <- function(traindata, testdata, categoricalvars, continuousvars, RCT = F, tree.depth = 2, n.trees = 200, Kfold = 6, B = 3, seed.cf = 3, plot.gbmperf = F, sim.big = NULL){
  
  output <- list()
  
  # Save traindata ID 
  ID <- traindata$ID
  
  # Prepare data
  if(RCT == F){
    traindata_ps <- traindata %>% dplyr::select(trt, ageatindex_centered, female_1, regioncensus_NORTHEAST, regioncensus_SOUTH, regioncensus_WEST, cci_1, cci_2, cci_..3, prevDMTefficacy_Medium.efficacy, prevDMTefficacy_High.efficacy, prevDMTefficacy_None, premedicalcost, prerelapse_num, severityScore, insuranceplan_High.Deductible, insuranceplan_Capitated, insuranceplan_other.unknown, hospitalization_1, premedicationcosts_cat_.7.829.30.274., premedicationcosts_cat_.30.274.68.053., premedicationcosts_cat_.68.053)
    # Formula for PS model in cross fitting
    fps <- as.formula("trt ~ ageatindex_centered + female_1 + regioncensus_NORTHEAST + regioncensus_SOUTH + regioncensus_WEST + cci_1 + cci_2 + cci_..3 + prevDMTefficacy_Medium.efficacy + prevDMTefficacy_High.efficacy + prevDMTefficacy_None + premedicalcost + prerelapse_num + severityScore + insuranceplan_High.Deductible + insuranceplan_Capitated + insuranceplan_other.unknown + hospitalization_1 + premedicationcosts_cat_.7.829.30.274. + premedicationcosts_cat_.30.274.68.053. + premedicationcosts_cat_.68.053")
  } 
  traindata <- traindata %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars), y, trt, time)
  
  # Prepare for cross-fitting
  N1 <- ifelse(is.numeric(traindata$trt), sum(traindata$trt), sum(as.numeric(traindata$trt) - 1))
  N0 <- nrow(traindata) - N1
  N <- N1 + N0
  p.aug <- ncol(traindata %>% dplyr::select(-y, -trt, -time)) + 1
  
  index1 <- rep(1:Kfold, floor(N1/Kfold))
  if(N1 > Kfold*floor(N1/Kfold)) index1 <- c(index1, 1:(N1 - Kfold*floor(N1/Kfold)))
  
  index0 <- rep(1:Kfold, floor(N0/Kfold))
  if(N0 > Kfold*floor(N0/Kfold)) index0 <- c(index0, Kfold + 1 - 1:(N0 - Kfold*floor(N0/Kfold)))
  
  delta3.mat <- delta4.mat <- matrix(NA, B, p.aug)
  sigma4.mat <- matrix(0, p.aug, p.aug)
  
  # Cross fitting
  converge <- rep(NA, B)
  for(bb in 1:B){
    cat("\n   Bootstrap:", bb, "out of", B, "\n")
    set.seed(bb + seed.cf)
    index1cv <- sample(index1, N1, F)
    index0cv <- sample(index0, N0, F)
    index <- rep(NA, N)
    index[traindata$trt == 1] <- index1cv
    index[traindata$trt == 0] <- index0cv
    
    f1.predictcv <- f0.predictcv <- pscv <- rep(NA, N)
    for(k in 1:Kfold){
      datatot_train <- traindata[index != k, ]
      if(RCT == F) datatot_train_ps <- traindata_ps[index != k, ]
      trt_train <- datatot_train$trt[index != k]
      
      datatot_test <- traindata[index == k, ]
      if(RCT == F) datatot_test_ps <- traindata_ps[index == k, ]
      x_test <- datatot_test %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
      
      data1 <- datatot_train %>% filter(trt == 1) %>% dplyr::select(-trt)
      set.seed(100)
      fit1.gbm <- gbm(y ~. - time + offset(time), data = data1, distribution = "poisson", interaction.depth = tree.depth, n.trees = n.trees, cv.folds = 5)
      best1.iter <- max(10, gbm.perf(fit1.gbm, method = "cv", plot.it = T))
      withCallingHandlers({
        f1.predictcv[index == k] <- predict(object = fit1.gbm, newdata = datatot_test, n.trees = best1.iter, type = "response")
      },
      warning=function(w) {
        if (grepl("does not add the offset", conditionMessage(w)))
          invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
      })
      
      data0 <- datatot_train %>% filter(trt == 0) %>% dplyr::select(-trt)
      set.seed(100)
      fit0.gbm <- gbm(y ~. - time + offset(time), data = data0, distribution = "poisson", interaction.depth = tree.depth, n.trees = n.trees, cv.folds = 5)
      best0.iter <- max(10, gbm.perf(fit0.gbm, method = "cv", plot.it = T))
      withCallingHandlers({
        f0.predictcv[index == k] <- predict(object = fit0.gbm, newdata = datatot_test, n.trees = best0.iter, type = "response")
      },
      warning=function(w) {
        if (grepl("does not add the offset", conditionMessage(w)))
          invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
      })
      
      if(RCT == F){
        pscv[index == k] <- IPTWfun(PSmodel = fps, data = datatot_train_ps, newdata = datatot_test_ps)$ps
      } else {
        pscv[index == k] <- sum(datatot_train$trt == 1)/nrow(datatot_train)
      }
    }
    
    ## bb-th cross fitting two regression estimator
    xb <- as.matrix(traindata %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars)))
    # trtb <- as.numeric(traindata$trt == 1) # DEPRECATED. input trt comes as numeric
    beta1.final <- onearmglmcount.dr(y = traindata$y, x = xb, time = traindata$time, trt = traindata$trt, ps = pscv, f.predictor = f1.predictcv)
    beta0.final <- onearmglmcount.dr(y = traindata$y, x = xb, time = traindata$time, trt = 1 - traindata$trt, ps = 1 - pscv, f.predictor = f0.predictcv)
    delta3.mat[bb, ] <- as.vector(beta1.final - beta0.final)
    
    ## bb-th cross fitting contrast regression estimator
    fit_two <- twoarmglmcount.dr(y = traindata$y, x = xb, time = traindata$time, trt = traindata$trt, ps = pscv, f1.predictor = f1.predictcv, f0.predictor = f0.predictcv)
    delta4.mat[bb, ] <- fit_two$coef
    converge[bb] <- fit_two$converge
    if(converge[bb] == T) sigma4.mat <- sigma4.mat + fit_two$vcov
    
    # Print intermediate output
#    cat(paste0("\n\n", bb, " out of ", B, " cross-fitting iterations\n"))
#    cat(paste0("Two regressions estimator (iteration ", bb, "):\n"))
#    print(delta3.mat[bb,])
#    cat(paste0("Contrast regression estimator (iteration ", bb, "):\n"))
#    print(delta4.mat[bb,])
  }
  
  # Final two regression estimator
  delta3 <- colMeans(delta3.mat)
  
  # Final contrast regression estimator
  converge4 <- (sum(converge) > 0)
  if(converge4 == T){
    delta4 <- colMeans(delta4.mat[converge == T, , drop = F])
    sigma4 <- sigma4.mat/sum(converge)
  } else {
    delta4 <- colMeans(delta4.mat)
    sigma4 <- sigma4.mat
  }
  
  # Predict score in test data
  if(is.null(testdata)){
    xtot <- traindata %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
  } else{
    xtot <- testdata %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
  }
  xtot <- as.matrix(cbind(1, xtot))
  scoretworeg <- as.numeric(as.matrix(xtot) %*% delta3)
  scorecontrastreg <- as.numeric(as.matrix(xtot) %*% delta4)
  
  # Optimal treatment: scorepois > 0 means treatment 0 is preferred, pois < 0 means treatment 1 is prefered
  itr_tworeg <- factor(ifelse(scoretworeg >= 0, 0, 1))
  itr_contrastreg <- factor(ifelse(scorecontrastreg >= 0, 0, 1))
  
  if(is.null(testdata)){
    output <- list(tworeg = list(fit = data.frame(coef = delta3), score = data.frame(ID = ID, score_tworeg = scoretworeg, itr_tworeg = itr_tworeg)),
                   contrastreg = list(fit = data.frame(coef = delta4, SE = sqrt(diag(sigma4.mat))), score = data.frame(ID = ID, score_contrastreg = scorecontrastreg, itr_contrastreg = itr_contrastreg)))
  } else {
    output$valueTwoReg$dhat <- itr_tworeg
    output$valueContrastReg$dhat <- itr_contrastreg
    output$valueTwoReg$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], a = testdata[["trt"]], d.hat = itr_tworeg, p.hat = testdata[['ps']], fu = exp(testdata[['time']]) / 365.25)
    output$valueContrastReg$vhat.dhat <- getValue(y = testdata[['postrelapse_num']], a = testdata[['trt']], d.hat = itr_contrastreg, p.hat = testdata[['ps']], fu = exp(testdata[['time']]) / 365.25)
  }
  
  # If large independent data are provided
  if (!is.null(sim.big)){
    temp.big <- format.countdata(data = sim.big$data, yvar = "mlogarr0001", timevar = "finalpostdayscount", trtvar = "trt", 
                                 xcontinuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "postrelapse_num", "offset", "FUweight"), 
                                 xcategoricalvars = c("female", "prevDMTefficacy"), imputation.method = NULL)
    testdata.big <- data.frame(y = temp.big$y, trt = temp.big$trt, time = log(temp.big$time), temp.big$x)
    xtot.big <- testdata.big %>% dplyr::select(all_of(categoricalvars), all_of(continuousvars))
    xtot.big <- as.matrix(cbind(1, xtot.big))
    scoretworeg.big <- as.numeric(as.matrix(xtot.big) %*% delta3)
    scorecontrastreg.big <- as.numeric(as.matrix(xtot.big) %*% delta4)
    itr_tworeg.big <- ifelse(scoretworeg.big >= 0, 0, 1)
    itr_contrastreg.big <- ifelse(scorecontrastreg.big >= 0, 0, 1)
    output$valueTwoReg$dhat.big <- itr_tworeg.big
    output$valueContrastReg$dhat.big <- itr_contrastreg.big
    output$valueTwoReg$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr_tworeg.big, betas = sim.big$betas)
    output$valueContrastReg$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = itr_contrastreg.big, betas = sim.big$betas)
  }
  
  return(output)
 }




