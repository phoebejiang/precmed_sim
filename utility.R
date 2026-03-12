# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: utility.R
# 
# Purpose: Utility functions for the PM implementation
# ------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(haven)
library(Hmisc)
library(fastDummies)
library(truncnorm)
library(reshape2)


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

simdata <- function(n, RCT = T, originaldata = NULL, seed = 999,
                    beta = c(-0.5, -0.25, 0, 0.25, 0.5),
                    beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003),
                    percentiles = seq(0, 1, by = 0.2)){
  
  #' Generate simulated data with settings based on original data
  #' Assume randomized treatment and independent covariates
  #' 
  #'   
  #' Here are the five different options of proportions (cumulative)
  #' seq(0, 1, by = 0.2)
  #' c(0, 0.1, 0.25, 0.75, 0.9, 1.0)
  #' c(0, 0.1, 0.4, 0.6, 0.9, 1.0)
  #' c(0, 0.55, 0.85, 1.0, 1.0)
  #' c(0, 0.55, 0.70, 0.85, 1.0, 1.0)
  #' 
  #' Here are the different beta values to determine different levels of heterogeneity
  #' c(-0.2,-0.2,-0.2,-0.2,-0.2) => no
  #' c(log(0.7), log(0.75), log(1), log(1.05), log(1.1) c(-0.36,-0.29,0,0.05,0.1) => low
  #' c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)) c(-0.92,-0.69,0,0.1,0.18) => medium
  #' c(log(0.3), log(0.5), log(1), log(1.1), log(1.5)) c(-1.2,-0.69,0,0.1,0.41) => high
  #'
  #'
  #' @param n sample size; integer
  #' @param RCT Whether treatment is randomized (T) or not. If RCT=T, treatment variable is randomized with binomial distribution.
  #' @param originaldata the dataset to compare the simulated data with; data.frame
  #' @param seed randomization seed; integer
  #' @param beta coefficients characterizing treatment effect heterogeneity; vector of length 5
  #'             beta[1]*trt*I(high responder to DMF) + 
  #'             beta[2]*trt*I(moderate responder to DMF) + 
  #'             beta[3]*trt*I(neutral) + 
  #'             beta[4]*trt*I(moderate responder to TERI) + 
  #'             beta[5]*trt*I(high responder to TERI)
  #'             In the absence of treatment effect heterogeneity, set all beta[1:5] with the same values
  #' @param beta.x coefficients for main effects of other covariates in the rate; vector of length 7
  #'               beta.x[1] (intercept)
  #'               beta.x[2]*ageatindex_centered
  #'               beta.x[3]*female
  #'               beta.x[4]*prerelapse_num
  #'               beta.x[5]*prevDMTefficacy== mediumhigh (reference low efficacy)
  #'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
  #'               beta.x[7]*premedicalcost
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @return data - a data frame of \code{n} rows and 8 columns of simulated data; data.frame
  #'         p1 - histogram plot comparing continuous variables between original data and simulated data 
  #'         p2 - bar plot comparing categorical variables between original data and simulated data

  myblue <- rgb(37, 15, 186, maxColorValue = 255)
  mygreen <- rgb(109, 173, 70, maxColorValue = 255)
  mygrey <- rgb(124, 135, 142, maxColorValue = 255)
  
  set.seed(seed)
  
  if (percentiles[1] != 0 | percentiles[length(percentiles)] != 1 | length(percentiles) != 6){message("Wrong values of percentiles!")}
  
  # Create an empty shell
  ds <- data.frame(matrix(NA, nrow = n, ncol = 10))
  colnames(ds) <- c("treatment", "ageatindex_centered", "female", "prerelapse_num",
                    "prevDMTefficacy", "premedicalcost", "numSymptoms",
                    "postrelapse_num", "finalpostdayscount", "group")
  
  # Define X, A, and time
  ds %<>% 
    mutate(trt =                 rbinom(n = n, size = 1, prob = 0.5), 
           treatment =           ifelse(trt == 1, "DMF", "TERI"),
           female =              rbinom(n = n, size = 1, prob = 0.75),
           ageatindex_centered = round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48, # rounded to integers
           prerelapse_num =      rpois(n = n, lambda = 0.44),
           prevDMTefficacy =     sample(x = c("None", "Low efficacy", "Medium and high efficacy"), 
                                        size = n, replace = TRUE, prob = c(0.45, 0.44, 0.11)),
           premedicalcost =      pmin(round(exp(rnorm(n = n, mean = 8.9, sd = 1.14)), 2), 600000), # rounded to 2 decimal points
           numSymptoms =         sample(x = c("0", "1", ">=2"), 
                                        size = n, replace = TRUE, prob = c(0.67, 0.24, 0.09)), # nuisance variable; do not include in the score
           finalpostdayscount =  ceiling(rgamma(n = n, shape = 0.9, scale = 500)), # rounded up to integers
           finalpostdayscount =  ifelse(finalpostdayscount > 2096, 2096, finalpostdayscount), # truncate at the max follow up day, 2096
           finalpostdayscount =  ifelse((finalpostdayscount > 2090) & (runif(1, 0, 1) < .5), 29, finalpostdayscount), # mimic the 1 month peak;  move roughly half of the large values to 29
           # finalpostdayscount = ifelse((finalpostdayscount < 29) & (runif(1, 0, 1) < .94), 29, finalpostdayscount), # mimic the 1 month peak # TODO: tune this between 0.94 and 0.95
           group =              "Simulated")
  
  # Define Y 
  xmat.score <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost, ds) %>% as.matrix()
  gamma <- matrix(c(-0.33, # Intercept
                    -0.001, # Age
                    0.05, # female
                    -0.002, # prerelapse_num
                    0.33, 0.02, # Medium/high and none DMT efficacy
                    -0.0000005), nrow = 7) # premedicalcost
  ds <- ds %>% mutate(score = exp(xmat.score %*% gamma),
                      Iscore = cut(score, quantile(score, percentiles), include.lowest = TRUE, labels = seq(1, 5)))
  
  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + 
                              Iscore + trt + trt * Iscore, ds) %>% as.matrix()
  
  # design matrix under all trt=0 and all trt=1
  ds_temp <- ds %>% mutate(trt1 = 1, trt0 = 0)
  xmat.rate0 <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + 
                               Iscore + trt0 + trt0 * Iscore, ds_temp) %>% as.matrix()
  
  xmat.rate1 <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + 
                               Iscore + trt1 + trt1 * Iscore, ds_temp) %>% as.matrix()
  
  betas <- matrix(c(beta.x[1], # Intercept
                    beta.x[2], # Age
                    beta.x[3], # female
                    beta.x[4], # prerelapse_num
                    beta.x[5], beta.x[6], # Medium/high and none DMT efficacy
                    beta.x[7], # premedicalcost
                    # Previous beta.x: -1.54 Intercept, -0.01 Age, 0.06 female, 0.47 prerelapse_num, 0.68, 0.13 Medium/high and none DMT efficacy, 0.000003 premedicalcost
                    0, 0, 0, 0, # Iscore categories 2:5 (reference group is Iscore1, high responders to DMF)
                    beta[1], beta[2] - beta[1], beta[3] - beta[1], beta[4] - beta[1], beta[5] - beta[1])) 
  # Treatment effect from high responders to DMF to high responders to TERI is: beta[1:5]
  # Here, the last five values of `betas` are:
  #    beta[1]*trt + 
  #    (beta[2]-beta[1])*I(moderate responder to DMF) +
  #    (beta[3]-beta[1])*I(neutral) +
  #    (beta[4]-beta[1])*I(moderate responder to TERI) +
  #    (beta[5]-beta[1])*I(high reponder to TERI) 
  #    e.g. -0.5 = beta1, -0.25 = beta1 + beta2, 0 = beta1 + beta3, 0.25 = beta1 + beta4, 0.5 = beta1 + beta5
  rate <-  exp(xmat.rate %*% betas)
  rate0 <- exp(xmat.rate0 %*% betas)
  rate1 <- exp(xmat.rate1 %*% betas)
  ds <- ds %>% mutate(rate = rate, rate1 = rate1, rate0 = rate0, postrelapse_num = rpois(n = n, lambda = rate * finalpostdayscount / 365.25))
  
  if (!is.null(originaldata)){
    # Combine with original data to generate comparison plots
    combined <- originaldata %>% 
      mutate(prevDMTefficacy = ifelse(prevDMTefficacy %in% c("High efficacy", "Medium efficacy"), "Medium and high efficacy", 
                                      as.character(prevDMTefficacy)),
             numSymptoms = ifelse(numSymptoms >= 2, ">=2", as.character(numSymptoms))) %>% 
      dplyr::select(trt, ageatindex_centered, female, prerelapse_num, prevDMTefficacy, premedicalcost, 
                    numSymptoms, postrelapse_num, finalpostdayscount) %>% 
      mutate(group = "MarketScan") %>% 
      rbind(ds %>% dplyr::select(-treatment, -score, - Iscore)) %>% 
      mutate(logpremedicalcost = log(premedicalcost + 1),
             logfinalpostdayscount = log(finalpostdayscount + 0.1)) %>% 
      dplyr::select(-premedicalcost, -finalpostdayscount)
    
    # Continuous variables
    melted1 <- melt(id.vars = "group", data = combined %>% dplyr::select(-prevDMTefficacy, -trt, -female, -numSymptoms)) %>% 
      mutate(value = as.numeric(value))
    p1 <- melted1 %>% 
      ggplot(aes(x = value, fill = group)) + 
      facet_wrap( ~ variable, scales = "free", nrow = 2) + 
      geom_histogram(aes(y = ..density.., group = group), binwidth = 1, position = "dodge") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = c(myblue, mygreen)) + 
      theme_classic() + 
      theme(legend.position = "top", 
            legend.title = element_blank(), 
            text = element_text(size = 20)) +
      labs(x = "Variable", y = "Proportion") 
    
    # Categorical variables
    melted2 <- melt(id.vars = "group", data = combined %>% dplyr::select(female, trt, prevDMTefficacy, numSymptoms, group))
    p2 <- melted2 %>% 
      ggplot(aes(x = value, group = group)) +
      facet_wrap( ~ variable, scales = "free") +
      geom_bar(aes(y = ..prop.., fill = group), stat = "count", position = 'dodge') + 
      scale_fill_manual(values = c(myblue, mygreen)) +
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            text = element_text(size = 20)) +
      labs(x = "Variable", y = "Proportion")
    
    ds <- ds %>% dplyr::select(-group)
    
  } else {
    ds <- ds %>% dplyr::select(-group)
  }
  
  # Create various forms of outcomes
  ds <- ds %>% mutate(arr = postrelapse_num / (finalpostdayscount / 365.25),
                      logarr0001 = log(arr + 0.001),
                      mlogarr0001 = -logarr0001,
                      logarr1 = log(arr + 1),
                      mlogarr1 = -logarr1,
                      logarr01r = log( (postrelapse_num + 0.1) / (finalpostdayscount / 365.25)),
                      mlogarr01r = -logarr01r,
                      offset = log(finalpostdayscount / 365.25),
                      FUweight = finalpostdayscount/sum(finalpostdayscount),
                      ID = paste0("ID", 1:n)) %>% 
    dplyr::select(-arr)
  
  if (!is.null(originaldata))  return(list(data = ds, betas = betas, percentiles = percentiles, plots = list(p1 = p1, p2 = p2)))
  if (is.null(originaldata))   return(list(data = ds, betas = betas, percentiles = percentiles))
}


format.countdata <- function(data, yvar, timevar, trtvar, xcontinuousvars, xcategoricalvars, RCT = F, imputation.method = NULL, scale = FALSE){
  
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
  }

  x <- cbind(x, data3[, c("FUweight", "offset", "postrelapse_num")])
  
  return(list(y = data3$y, time = data3$time, trt = data3$trt, x = x))
}

get_gamma_beta <- function(originaldata, x.score, x.rate, offset, percentiles = seq(0, 1, by = 0.2), seed = 999){
  
  #' Inform values of gamma (score) and beta (rate) based on original data
  #' 
  #' @param originaldata the dataset to compare the simulated data with; data.frame
  #' @param x.score covariates in the score; vector
  #' @param x.rate covariates in the Poisson outcome model (rate); vector
  #' @param offset offset variable in both Poisson outcome models; vector of data
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @param seed randomization seed; integer
  #' @return coef.score - proposed coefficients for the score
  #'         coef.rate - proposed coefficients for the Poisson rate
  
  set.seed(seed)
  
  # Get coefficients for the score
  form.score1 <- as.formula(paste0("postrelapse_num ~ ", paste0(x.score, collapse = " + ")))
  
  mod1 <- glm(form.score1, family = "poisson", data = originaldata, offset = offset, subset = trt == 1)
  mod0 <- glm(form.score1, family = "poisson", data = originaldata, offset = offset, subset = trt == 0)
  
  coef.score <- coef(mod1) - coef(mod0)
  
  # Add estimated score to the data
  originaldata <- originaldata %>% mutate(score = exp(predict(mod1, originaldata) - predict(mod0, originaldata)),
                                          Iscore = cut(score, quantile(score, percentiles), include.lowest = T))
  
  # Get coefficients for the rate
  form.rate <- as.formula(paste0("postrelapse_num ~ trt + Iscore + trt*Iscore + ", paste0(x.rate, collapse = " + ")))
  mod1 <- glm(form.rate, family = "poisson", offset = offset, data = originaldata)
  coef.rate <- coef(mod1)
  
  return(list(coef.score = coef.score, coef.rate = coef.rate))
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
  denom <- sum((a == d.hat)*fu/pa)
  num <- sum(y*(a == d.hat)/pa)
  return(list(U = num, W = denom))
}

getTrueValue <- function(ds, d.hat, betas){
  
  #' True value function of an estimated ITR, V(d.hat)
  #' 
  #' @param ds - A large test set, e.g. simdata$data
  #' @param d.hat - Estimated optimal treatment 0/1
  #' @param betas - Coefficients to generate `ds`, outputted from simdata
  
  ds$d.hat <- as.numeric(as.character(d.hat))
  xmat.rate <- model.matrix(as.formula(paste0("~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + 
                                              Iscore + d.hat + d.hat*Iscore")), ds) %>% as.matrix()  
  rate <- exp(xmat.rate %*% betas) # E^d.hat[Y] for each individual
  
  return(mean(rate))
}

getTrueOptimalValue <- function(n, beta, beta.x, RCT = T,
                                percentiles = seq(0, 1, by = 0.2), seed = 0){
  
  #' True value function of the true optimal ITR, V(d)
  #' An empirical way to calculate the true optimal value given the setup in simdata() function
  #' 
  #' @param n sample size (choose something large to get close to the theoretical truth); scalar
  #' @beta beta coefficients characterizing treatment effect heterogeneity; vector of length 5
  #'             beta[1]*trt*I(high responder to DMF) + 
  #'             beta[2]*trt*I(moderate responder to DMF) + 
  #'             beta[3]*trt*I(neutral) + 
  #'             beta[4]*trt*I(moderate responder to TERI) + 
  #'             beta[5]*trt*I(high responder to TERI)
  #'             In the absence of treatment effect heterogeneity, set all beta[1:5] with the same values
  #' @param beta.x coefficients for main effects of other covariates in the rate; vector of length 7
  #'               beta.x[1] (intercept)
  #'               beta.x[2]*ageatindex_centered
  #'               beta.x[3]*female
  #'               beta.x[4]*prerelapse_num
  #'               beta.x[5]*prevDMTefficacy== mediumhigh (reference low efficacy)
  #'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
  #'               beta.x[7]*premedicalcost
  #' @param RCT Whether treatment is randomized (T) or not. If RCT=T, treatment variable is randomized with binomial distribution.
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @return true optimal value function, $V = E[Y^(optA)]$
  
  sim <- simdata(n = n, RCT = RCT, beta = beta, beta.x = beta.x, seed = seed, percentiles = percentiles)
  
  if (length(unique(beta)) == 1){ 
    # If no heterogeneity
    ds <- sim$data %>% 
      mutate(trueA = ifelse(beta[1] > 0, 0, 1), # if beta > 0, A0 is better than A1; if beta < 0, A1 is better than A0
             trueA = ifelse(beta[1] == 0, trt, trueA)) # if beta = 0, A1 and A0 are the same
  } else {
    # If heterogeneity
    neutral <- which(beta == 0) # look for the neutral group which has beta = 0 (log(1) = 0)
    ds <- sim$data %>% 
      mutate(trueA = ifelse(as.numeric(Iscore) < neutral, 1, 0), # TODO: this definition of trueA would change if the definition score groups changes
             trueA = ifelse(as.numeric(Iscore) == neutral, trt, trueA)) # true optimal A
    # If all 5 groups
    # optimal A is DMF if in the score group 1 and 2 (high and moderate responder to DMF); 
    # optimal A is either DMF or TERI if in score group 3 (neutral group); 
    # optimal A is TERI if in score group 4 & 5 (high and moderate responder to TERI)
  }
  
  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + 
                              Iscore + trueA + trueA*Iscore, ds) %>% as.matrix()
  rate <- exp(xmat.rate %*% sim$betas) # Not FU, this is E[Y|time=1]
  
  return(mean(rate)) # empirical estimate over a large number of samples
}

getTrueWorstValue <- function(n, beta, beta.x, RCT = T,
                              percentiles = seq(0, 1, by = 0.2), seed = 0){
  
  #' True value function of the true worse ITR
  #' An empirical way to calculate the true worse value given the setup in simdata() function
  #' The worse value is the opposite of the true optimal value
  #' 
  #' @param n sample size (choose something large to get close to the theoretical truth); scalar
  #' @beta beta coefficients characterizing treatment effect heterogeneity; vector of length 5
  #'             beta[1]*trt*I(high responder to DMF) + 
  #'             beta[2]*trt*I(moderate responder to DMF) + 
  #'             beta[3]*trt*I(neutral) + 
  #'             beta[4]*trt*I(moderate responder to TERI) + 
  #'             beta[5]*trt*I(high responder to TERI)
  #'             In the absence of treatment effect heterogeneity, set all beta[1:5] with the same values
  #' @param beta.x coefficients for main effects of other covariates in the rate; vector of length 7
  #'               beta.x[1] (intercept)
  #'               beta.x[2]*ageatindex_centered
  #'               beta.x[3]*female
  #'               beta.x[4]*prerelapse_num
  #'               beta.x[5]*prevDMTefficacy== mediumhigh (reference low efficacy)
  #'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
  #'               beta.x[7]*premedicalcost
  #' @param RCT Whether treatment is randomized (T) or not. If RCT=T, treatment variable is randomized with binomial distribution.
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @return true optimal value function, $V = E[Y^(optA)]$
  
  sim <- simdata(n = n, RCT = RCT, beta = beta, beta.x = beta.x, seed = seed, percentiles = percentiles)
  
  if (length(unique(beta)) == 1){ 
    # If no heterogeneity
    ds <- sim$data %>% 
      mutate(trueWorstA = ifelse(beta[1] > 0, 1, 0), # if beta > 0, A0 is better than A1, so worst ITR yields A1; if beta < 0, A1 is better than A0, so worst ITR yield A0
             trueWorstA = ifelse(beta[1] == 0, trt, trueWorstA)) # if beta = 0, A1 and A0 are the same
  } else {
    # If heterogeneity
    neutral <- which(beta == 0) # look for the neutral group which has beta = 0 (log(1) = 0)
    ds <- sim$data %>% 
      mutate(trueWorstA = ifelse(as.numeric(Iscore) < neutral, 0, 1), # TODO: this definition of trueA would change if the definition score groups changes
             trueWorstA = ifelse(as.numeric(Iscore) == neutral, trt, trueWorstA)) # true optimal A
    # If all 5 groups
    # worst A is A0 if in the score group 1 and 2 (high and moderate responder to A1); 
    # worst A is either A0 or A1 if in score group 3 (neutral group); 
    # worst A is A1 if in score group 4 & 5 (high and moderate responder to A0)
  }
  
  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + 
                              Iscore + trueWorstA + trueWorstA*Iscore, ds) %>% as.matrix()
  rate <- exp(xmat.rate %*% sim$betas) # Not FU, this is E[Y|time=1]
  
  return(mean(rate)) # empirical estimate over a large number of samples
}