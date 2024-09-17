# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: 02_propensityscore.R
# Developer/Programmer: pj
# Date: 28SEP2020
# 
# Purpose: Estimate PS model in analytic dataset and output as an additional column
#    
# 
# Platform: Windows
# R Version: 4.0.1
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#   28SEP2020 gs      Start the script
#   06OCT2020 pj      Read in preprocessed MarketScan data from a function
#   07OCT2020 gs      Update PS model based on literature, output new plot
#   19OCT2020 gs      Explore alternatives to include premedication costs in PS model
#   17NOV2020 gs      Update IPTWfun
# ------------------------------------------------------------------

#library(tidyverse)

#' Estimate PS and append PS and IPTW to dataset
#'
#' @param PSmodel A formula with treatment variable of the left-hand side and covariates separated by "+" on the right-hand side
#' @param data Dataframe from which to fetch all variables in PSmodel
#' @param newdata Data for which we want to predict the PS. If NULL, dataframe supplied in \code{data} is used
#' 
#' @return Same dataframe as supplied by \code{newdata}, with 2 additional columns ps and itpw
#'
#' @export

IPTWfun <- function(PSmodel = trt ~ ageatindex_centered + female + regioncensus + cci + prevDMTefficacy + 
                      premedicalcost + prerelapse_num + severityScore + insuranceplan + hospitalization + 
                      premedicationcosts_cat, 
                    data, 
                    newdata = NULL){
  if(class(PSmodel) != "formula") stop("PSmodel must be a formula.")
  ps <- glm(PSmodel, family = "binomial", data = data)
  if(is.null(newdata)) newdata <- data
  # TODO: if one category absent in datatot_train, PS coefficient will be NA, need to make sure to keep only the columns in testdata for which coefficient is not NA
  newdata$ps <- predict(ps, newdata, type = "response")
  newdata <- newdata %>% mutate(iptw = ifelse(trt == 1, 1/ps, 1/(1 - ps)))
  return(newdata)
}


