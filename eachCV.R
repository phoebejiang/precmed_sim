# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: eachCV.R
# Developer/Programmer: pj
# Date: 12NOV2020
# 
# Purpose: Performs one CV iteration for different methods. Called in main.R.
#    
# 
# Platform: Windows
# R Version: 4.0.1
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   12NOV2020 pj      Start the main structure
#   16NOV2020 gs      Add dWOLS, Lu's 4 approaches
#   23NOV2020 pj      Adapt code to HPC as a function
#   24NOV2020 gs      Format input data in the beginning instead of in each CV
#   01DEC2020 pj      Added three additional methods: allA0, allA1, listDTR3
#                     Change from outputting value to outputting the denom and num separately
#   02DEC2020 pj      Add tryCatch to twoReg 
#   29JAN2021 gs      Add linear and Poisson LASSO
#   11MAR2021 pj      Calculate PS/IPTW after splitting into training/test set
#   15MAR2021 pj      Convert factor trt to numeric for all methods
#                     Correct PS for test data (PS trained from training not testing)
#   20MAR2021 pj      Add more sample sizes in the beginning to see the shape
#   26APR2021 pj      Add code to get V(d.hat)
#   10MAY2021 pj      Wrap listDTR2 with tryCatch and modify output if error
#   22JUN2021 pj      Allow both "allA1" and "allDMF", "allA0" and "allTERI"
#   14MAR2022 pj      Add "allGA" to methods for CONFIRM data
# ------------------------------------------------------------------


#######################################################################
############################# Function for One CV Iteration ############################
#######################################################################

eachCV <- function(data, method, outcome, folds, n.fold = 10, 
                   categoricalvars = c("female_1", "cvd_1", "diabetes_1", "prevDMTefficacy_Medium.efficacy", "prevDMTefficacy_High.efficacy", "prevDMTefficacy_None", "insuranceplan_High.Deductible", "insuranceplan_Capitated", "insuranceplan_other.unknown", "premedicationcosts_cat_.7.829.30.274.", "premedicationcosts_cat_.30.274.68.053.", "premedicationcosts_cat_.68.053"), 
                   continuousvars = c("ageatindex_centered", "prerelapse_num", "premedicalcost", "numSymptoms", "severityScore"),
                   RCT = F,
                   sim.big = NULL,
                   seed = 999){
  
  #' One iteration of cross-validation
  #' Estimate PS in each training and testing dataset after each splitting
  #' 
  #' X variables can be specified by user but certain outcome (postNRL) and offset are fixed by the method
  #' 
  #' @param data Input data; data.frame
  #' @param method PM method to be implemented; string
  #'               Possible methods: 'allA1' or 'allDMF','allA0' or 'allTERI','linear','weightedLinear','linearLASSO','weightedPoisson','poissonLASSO','negBin','weightedNegBin','dWOLS','listDTR2','listDTR3','poisson','boosting','twoReg','contrastReg'
  #' @param outcome Response variable (whether higher outcome is better depends on method) specified by user; string
  #' @param folds Row indices of the input data in the same CV fold are saved in each element of list; list of length n.fold
  #' @param n.fold Number of CV folds; integer
  #' @param categoricalvars Categorical X variables; vector of strings
  #' @param continuousvars Continuous X variables; vector of strings
  #' @param RCT Whether treatment is randomized (T) or not. If RCT=T, the PS is taken as the proportion of patients treated with A1
  #' @param sim.big a large test set object (independent from data) to get true value function; list object with elements `data` and `betas` 
  #' @param seed randomization seed
  #' 
  #' @return eachcv.results - Estimated value for each CV fold; vector of size n.fold
  
  stopifnot(length(folds) == n.fold)
  
  set.seed(seed)
  eachcv <- vector("list", n.fold)
  names(eachcv) <- paste0("fold", 1:n.fold)
  
  # Cross-validation starts  
  for (fold.i in 1:n.fold){
    
    # cat("\n CV fold =", fold.i, "out of", n.fold)
    traindata <- data[-folds[[fold.i]],]
    testdata <- data[folds[[fold.i]],]

    # Transform trt to numeric
    traindata$trt <- as.numeric(traindata$trt == 1)
    if(is.null(testdata) == F) testdata$trt <- as.numeric(testdata$trt == 1)
    
    # Calculate PS/IPTW
    if (RCT == F){
      # Retrieve categorical variables (dummies) in each train/set and create formula for PS model
      # (in case some train/test set does not have all categories)
      vartr <- colnames(traindata)[str_detect(colnames(traindata), "region|cci|prevDMTefficacy|insuranceplan|premedicationcosts")]
      # varte <- colnames(testdata)[str_detect(colnames(testdata), "region|cci|prevDMTefficacy|insuranceplan|premedicationcosts")]
      fpstr <- as.formula(paste("trt ~ ageatindex_centered + female_1 + premedicalcost + prerelapse_num + severityScore + hospitalization_1 +", 
                                paste0(vartr, collapse = "+")))
      # fpste <- as.formula(paste("trt ~ ageatindex_centered + female_1 + premedicalcost + prerelapse_num + severityScore + hospitalization_1 +", 
      #                           paste0(varte, collapse = "+")))           
      traindata <- IPTWfun(data = traindata, PSmodel = fpstr)
      testdata <- IPTWfun(data = traindata, PSmodel = fpstr, newdata = testdata)
    } else {
      trainps <- mean(traindata$trt)
      traindata <- traindata %>% mutate(ps = trainps, iptw = ifelse(trt == 1, 1/ps, 1/(1 - ps)))
      testdata <- testdata %>% mutate(ps = trainps, iptw = ifelse(trt == 1, 1/ps, 1/(1 - ps)))
    }
 
    # TODO: think of a better way to call different methods than long if-elseif statements
    if (method %in% c("allA1")){
      
      eachcv[[paste0("fold", fold.i)]]$dhat <- rep(1, nrow(testdata))
      eachcv[[paste0("fold", fold.i)]]$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], 
                                                             a = testdata[["trt"]], 
                                                             d.hat = rep(1, nrow(testdata)), 
                                                             p.hat = testdata$ps, 
                                                             fu = exp(testdata[['time']]) / 365.25)
      if (!is.null(sim.big)){
            eachcv[[paste0("fold", fold.i)]]$dhat.big <- rep(1, nrow(sim.big$data))
            eachcv[[paste0("fold", fold.i)]]$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = rep(1, nrow(sim.big$data)), betas = sim.big$betas)
      } 
      
    } else if (method %in% c("allA0")){
      
      eachcv[[paste0("fold", fold.i)]]$dhat <- rep(0, nrow(testdata))
      eachcv[[paste0("fold", fold.i)]]$vhat.dhat <- getValue(y = testdata[["postrelapse_num"]], 
                                                             a = testdata[["trt"]], 
                                                             d.hat = rep(0, nrow(testdata)), 
                                                             p.hat = testdata$ps, 
                                                             fu = exp(testdata[['time']]) / 365.25)
      if (!is.null(sim.big)){
          eachcv[[paste0("fold", fold.i)]]$dhat.big <- rep(0, nrow(sim.big$data))
          eachcv[[paste0("fold", fold.i)]]$v.dhat <- getTrueValue(ds = sim.big$data, d.hat = rep(0, nrow(sim.big$data)), betas = sim.big$betas)
      }
      
    } else if (method == "linear"){
      
      traindata1 <- traindata %>% filter(trt == 1) # A1
      traindata0 <- traindata %>% filter(trt == 0) # A0

      # Weighted linear regression
      eachcv[[paste0("fold", fold.i)]] <- itrLinear(traindata1 = traindata1, 
                                                    traindata0 = traindata0, 
                                                    testdata = testdata, 
                                                    outcome = outcome,
                                                    categoricalvars = categoricalvars,
                                                    continuousvars = continuousvars,
                                                    testps = testdata$ps,
                                                    trainweight1 = traindata1$FUweight,
                                                    trainweight0 = traindata0$FUweight,
                                                    sim.big = sim.big)
      
    } else if (method == "negBin"){
      
      traindata1 <- traindata %>% filter(trt == 1) # A1
      traindata0 <- traindata %>% filter(trt == 0) # A0
      
      # Weighted negative binomial regression
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({itrNegBin(traindata1 = traindata1, 
                                                              traindata0 = traindata0, 
                                                              testdata = testdata, 
                                                              outcome = outcome, 
                                                              offset = "offset", 
                                                              categoricalvars = categoricalvars,
                                                              continuousvars = continuousvars,
                                                              testps = testdata$ps,
                                                              trainweight1 = NULL,
                                                              trainweight0 = NULL, 
                                                              sim.big = sim.big)
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat = NA)
      }
     )
      
    } else if (method == "dWOLS"){
      
      # dWOLS
      eachcv[[paste0("fold", fold.i)]] <- itrDWOLS(traindata = traindata,
                                                   testdata = testdata,
                                                   outcome = outcome,
                                                   Xoutcome = c(categoricalvars, continuousvars, "offset"),
                                                   Xinteraction = c(categoricalvars, continuousvars),
                                                   dWOLSweight = "IPTW", 
                                                   sim.big = sim.big)
      
    } else if (method == "listDTR2"){
      
      # List-based DTR
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({ itrLIST(traindata = traindata, 
                                                             testdata = testdata, 
                                                             outcome = outcome, # listdtr assumes higher better
                                                             treatment = "trt", 
                                                             categoricalvars = categoricalvars,
                                                             continuousvars = continuousvars,
                                                             maxlen = 2L,
                                                             seed = seed + fold.i, 
                                                             sim.big = sim.big)
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat = NA)
      }
      )
      
    } else if (method == "listDTR3"){
      
      # List-based DTR
      eachcv[[paste0("fold", fold.i)]] <- itrLIST(traindata = traindata, 
                                                        testdata = testdata, 
                                                        outcome = outcome, # listdtr assumes higher better
                                                        treatment = "trt", 
                                                        categoricalvars = categoricalvars,
                                                        continuousvars = continuousvars,
                                                        maxlen = 3L,
                                                        seed = seed + fold.i, 
                                                        sim.big = sim.big)
      
    } else if (method == "poisson"){ 
      
      traindata1 <- traindata %>% filter(trt == 1) # A1
      traindata0 <- traindata %>% filter(trt == 0) # A0
      
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({itrLuPoisson(traindata1 = traindata1,
                                                                 traindata0 = traindata0,
                                                                 testdata = testdata,
                                                                 categoricalvars = categoricalvars,
                                                                 continuousvars = continuousvars, 
                                                                 sim.big = sim.big)
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat= NA)
      }
      )
    } else if (method == "boosting"){ 
      
      traindata1 <- traindata %>% filter(trt == 1) # A1
      traindata0 <- traindata %>% filter(trt == 0) # A0
      
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({itrLuBoosting(traindata1 = traindata1,
                                                                  traindata0 = traindata0,
                                                                  testdata = testdata,
                                                                  categoricalvars = categoricalvars,
                                                                  continuousvars = continuousvars,
                                                                  tree.depth = 2, 
                                                                  n.trees = 100, 
                                                                  plot.gbmperf = F, 
                                                                  sim.big = sim.big)
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat = NA)
      }
     )

    } else if (method == "twoReg"){ 
      
      # TODO: get two regressions and contrast regression values as one method (runs the same code twice otherwise)
      
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({itrLuDR(traindata = traindata,
                                                            testdata = testdata,
                                                            categoricalvars = categoricalvars,
                                                            continuousvars = continuousvars,
                                                            RCT = RCT,            
                                                            tree.depth = 2, 
                                                            n.trees = 100, 
                                                            Kfold = 5, 
                                                            B = 3, 
                                                            seed.cf = 3, 
                                                            plot.gbmperf = F, 
                                                            sim.big = sim.big)$valueTwoReg
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat = NA)
      }
     )

    } else if (method == "contrastReg"){ 
      
      eachcv[[paste0("fold", fold.i)]] <- tryCatch({itrLuDR(traindata = traindata,
                                                            testdata = testdata, 
                                                            categoricalvars = categoricalvars,
                                                            continuousvars = continuousvars,
                                                            RCT = RCT,
                                                            tree.depth = 2, 
                                                            n.trees = 100, 
                                                            Kfold = 5, 
                                                            B = 3, 
                                                            seed.cf = 3, 
                                                            plot.gbmperf = F, 
                                                            sim.big = sim.big)$valueContrastReg
      }, error = function(error_message) {
        cat("\n")
        message(error_message)
        cat("\nSolution: Set U and W to NA and move on.\n")
        list(dhat = NA, vhat.dhat = list(U = NA, W = NA), dhat.big = NA, v.dhat = NA)
      }
     )
      
    } else{
      stop("Method not recognized! Pick one from ['allA1', 'allA0', 'allDMF', 'allGA', 'allTERI', 'linear', 'weightedLinear', 'linearLASSO', 'weightedPoisson', 'poissonLASSO', 'negBin', 'weightedNegBin', 'dWOLS', 'listDTR2', 'listDTR3', 'poisson', 'boosting', 'twoReg', 'contrastReg'].") 
    }

  }#end-of-loop fold.i
    
  return(eachcv)
}
