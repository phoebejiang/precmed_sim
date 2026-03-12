# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: 01-preprocessing.R
#
# Purpose: First look at the IHD outputted analysis data of Truven MarketScan (see SAP for details)
# 
# ------------------------------------------------------------------


library(tidyverse)
library(magrittr)
library(Hmisc)

source("./utility.R")

preprocessMarketScanData <- function(wd, allvars = FALSE, NTDreplicate = FALSE, PScasestudy = FALSE){
  
  #' Read in and preprocess the raw MarketScan data outputted from IHD
  #' @param wd string, the working directory path
  #' @param allvars Whether all intermediate variables should be kept (if need more granular covariates, e.g. to apply ML methods); logical
  #' @param NTDreplicate Whether to mimic the exclusion criteria of the PM NTD analysis and define specific covariates; logical
  #'                     Exclusion: pre-treatment with any DMT other than GA or IFN 
  #'                     Covariates: Age, # prior DMTs, Prior GA, Prior interferon, # relapsed in the previous 24 months (we only have 12 months)
  #' @param PScasestudy Whether to apply the inclusion/exclusion criteria of the PS case study; logical
  #'                     Exclusion: pre-treatment with any DMT other than GA or IFN, or with no DMT
  #' @return data.frame, preprocessed MarketScan data of 8599 rows and 59 columns (78 columns if allvars = T)
  
  
  ##################################################
  #### Read in raw analysis data ####
  ##################################################
  
  # Set working directory
  setwd(wd)
  
  # Read in the MarketScan raw analysis data
  
  newvars <- readSAS("./analysisfile_newvariables.sas7bdat") # 11 additional variables to be used to computer MS severity score
  
  ms0 <- readSAS("./analysisfile.sas7bdat") %>% 
    mutate(treatment = ifelse(treatment == 1, "DMF", "TERI")) %>% # 8887
    left_join(newvars, by = "patientid") # merge with new MS severity score items
  
  ##################################################
  #### Check basic properties of raw data ####
  ##################################################
  
  if (FALSE) {
    checkComplete(ms0) # check completeness
    describe(ms0) # check all summary statisitics
    colLabels(ms0) # check column labels (if any)
    stopifnot(all(length(unique(ms0$patientid)) == nrow(ms0))) # one row per unique patient
    
    # Study periods
    summary(ms0$indexdate) # drug initiation date
    summary(ms0$finalenddate) # final end date; the minimum of the next three dates
    summary(ms0$discontinuationdate) # date of discontinuation of index drug
    summary(ms0$switchdmtdate) # date of switching to other DMT drugs
    summary(ms0$studyenddate) # last date of study period or last date of enrollment
    
    # Check finalenddate is the earliest date among study end date, discontinuation date, and DMT switch date
    tmp <- ms0 %>% transmute(earliest = pmin(studyenddate, discontinuationdate, switchdmtdate, na.rm = TRUE))
    stopifnot(all(tmp$earliest == ms0$finalenddate)) 
    
    # Age
    stopifnot(ms0$ageatindex >= 18 & ms0$ageatindex <= 64)
    hist(ms0$ageatindex) # skewed to the left
    
    # Check postrelapse_num mathces with the row summation of nonNA rdate's
    stopifnot(all(rowSums(!is.na(ms0 %>% select(contains("rdate")))) == ms0$postrelapse_num))
  }
  
  ##################################################
  #### Apply additional exclusion criteria (clean up rows) ####
  ##################################################
  # Reordered the exclusion criteria to more easily compare the 60 vs. 90 days cohorts
  
  # Exclude missing discontinuation (due to missing days of supply value) (n = 15)
  ms1 <- ms0 %>% filter(!is.na(discontinuation))
  
  # Exclude negative costs (n = 5) 
  ms2 <- ms1 %>% filter(premedicalcost >= 0, premedicationcosts >= 0) 
  if (FALSE) {
    summary(ms0$premedicalcost)
    table(ms0$premedicalcost < 0, ms0$premedicationcosts < 0)
  }
  
  # Exclude missing regions (n = 255)
  ms3 <- ms2 %>% filter(regioncensus != "")
  
  # # Exclude patients with exposure to index drug for < 90 days (n = 2209) # PAUSE FOR RIGHT NOW
  # ms90 <- ms3 %>% filter(finalpostdayscount >= 90)
  # # Check exclusion with exposure to index drug for < 60 days (n = 1606) 
  # ms60 <- ms3 %>% filter(finalpostdayscount >= 60)
  
  # Exclusion criteria to mimic NTD cohort: treatment-naive or prior GA or prior IFN
  if (NTDreplicate == TRUE) {
    # ms3 <- ms3 %>% filter(prega == 1 | preifn == 1) # 3797 patients
    ms3 <- ms3 %>% filter(prefingolimod == 0 & predaclizumab == 0 & prealemtuzumab == 0 & prenatalizumab == 0 & preocrelizumab == 0 & prerituximab == 0 & premitoxantrone == 0)
  }
  
  # Exclusion criteria to apply PS case study inclusion/exclusion: prior GA or prior IFN
  if (PScasestudy == TRUE) {
    # ms3 <- ms3 %>% filter(prega == 1 | preifn == 1) # 3797 patients
    ms3 <- ms3 %>% filter(prefingolimod == 0 & predaclizumab == 0 & prealemtuzumab == 0 & prenatalizumab == 0 & preocrelizumab == 0 & prerituximab == 0 & premitoxantrone == 0 & (preifn == 1 | prega == 1 | prepegifn == 1))
  }
  
  ##################################################
  #### Preprocess covariates (clean up columns) ####
  ##################################################
  
  # Merge and recategorize unbalanced covariates and add severity score
  ms4 <- ms3 %>%
    mutate(ageatindex_centered = as.numeric(scale(ageatindex, scale = FALSE)),
           female = ifelse(sex == "FEMALE", 1, 0),
           insuranceplan = ifelse(plan %in% c("PPO/EPO", "POS"), "POS/PPO/EPO", plan),
           employmentstatus = ifelse(employmentstatus %in% c("Early Retiree", "Retiree (status unknown)", "COBRA Continuee"), "Retiree/Cobra",
                                     ifelse(employmentstatus %in% c("Long Term Disability", "Medicare Eligible Retiree"), "Disability",
                                            ifelse(employmentstatus %in% c("Other/Unknown", "Surviving Spouse/Depend."), "Other/Unknown/Surviving Spouse/Depend.", employmentstatus))),
           employmentstatus_2cat = ifelse(employmentstatus == "Active Full Time", "Active Full Time", "Other"),
           cci = ifelse(charlsonscore >= 3, ">=3", charlsonscore),
           eci = ifelse(elixhauser >= 5, ">=5", elixhauser),
           # pre-index medication costs by quartile; -0.1 to accommodate 0 costs
           # premedicationcosts_cat = cut(premedicationcosts, breaks = c(-0.1, quantile(premedicationcosts, probs = c(1/4, 1/2, 3/4, 1))))
           premedicationcosts_cat = factor(cut(premedicationcosts, breaks = c(-0.1, 7829, 30274, 68053, max(premedicationcosts))), labels = c("[0-7,829]", "(7,829-30,274]", "(30,274-68,053]", ">68,053")),
           timetorelapse = ifelse(is.na(difftime(as.Date(rdate1), as.Date(indexdate), units = "days")), finalpostdayscount, difftime(as.Date(rdate1), as.Date(indexdate), units = "days")),
           relapseindicator = factor(ifelse(is.na(difftime(as.Date(rdate1), as.Date(indexdate), units = "days")), 0, 1)),
           discontinuationindicator = factor(ifelse(is.na(finalenddate == discontinuationdate) | finalenddate != discontinuationdate, 0, 1)),
           timetodiscontinuation = ifelse(discontinuationindicator == 0, finalpostdayscount, difftime(as.Date(discontinuationdate), as.Date(indexdate), units = "days")),
           prerelapse_cat = case_when(prerelapse_num == 0 ~ "0", prerelapse_num == 1 ~ "1", prerelapse_num > 1 ~ "2+")
    ) %>% 
    rowwise() %>% 
    mutate(numPrevDMT = sum(prefingolimod, preifn, prepegifn, prega, predaclizumab, prealemtuzumab, prenatalizumab, 
                            preocrelizumab, prerituximab, premitoxantrone, na.rm = T), 
           # previous DMT is rare for each DMT (except ifn or ga); consider summation or a weighted summation of previous DMT use (weight = % time on drug)
           # --> weighting is a great stats idea but may be losing clinical soundnesss
           # similar problem for symptoms and severity
           # other meds are balanced so no processing done
           numSymptoms = sum(demcns, opticnervevisual, dizzinessgiddiness, fatiguemalaise, neurobladder, othermyelitis, neuralgia, na.rm = T),
           numSeverity = sum(walkingaids, wheelchair, hospitalbeds, transportation, physicalmed, na.rm = T),
           # --> common to consider no previous DMT vs. most recent=injectable vs. most recent=infusion vs. most recent=oral
           prevDMTdate = pmax(prefingolimoddate, preifndate, prepegifndate, pregadate, predaclizumabdate, prealemtuzumabdate, prenatalizumabdate, 
                              preocrelizumabdate, prerituximabdate, premitoxantronedate, na.rm = T),
           prevDMT = ifelse(is.na(prevDMTdate), "None", NA), 
           prevDMT = ifelse(is.na(prevDMTdate) == FALSE & prevDMTdate == prefingolimoddate, "Oral", prevDMT),
           prevDMT = ifelse(is.na(prevDMTdate) == FALSE & prevDMTdate %in% c(preifndate, prepegifndate, pregadate, predaclizumabdate), "Injectable", prevDMT),
           prevDMT = ifelse(is.na(prevDMTdate) == FALSE & prevDMTdate %in% c(prealemtuzumabdate, prenatalizumabdate, preocrelizumabdate, prerituximabdate, premitoxantronedate), "Infusion", prevDMT),
           # common to categorize by efficacy of previous DMT
           prevDMTefficacy = ifelse(is.na(prevDMTdate), "None", NA), 
           prevDMTefficacy = ifelse(is.na(prevDMTdate) == FALSE & prevDMTdate %in% c(preifndate, prepegifndate, pregadate), "Low efficacy", prevDMTefficacy),
           prevDMTefficacy = ifelse(is.na(prevDMTdate) == FALSE & prevDMTdate %in% c(prefingolimoddate , predaclizumabdate), "Medium efficacy", prevDMTefficacy),
           prevDMTefficacy = ifelse(is.na(prevDMTdate) == FALSE & prevDMTdate %in% c(prealemtuzumabdate, prenatalizumabdate, preocrelizumabdate, prerituximabdate, premitoxantronedate), "High efficacy", prevDMTefficacy),
           # combine walking aids, wheelchair, hospital bel
           mobility = ifelse(sum(walkingaids,wheelchair,hospitalbeds) >= 1, 1, 0),
           severityScore = scorevisual + scorebrainstem + ((scorewalking + scorecerebellar) > 0) * 2 + scorepyramidal * 2 + 
             scoresensory + scorespeech + scorebladder * 2 + scorecerebral + scoregeneral + 
             ((walkingaids + wheelchair + hospitalbeds) > 0) + # which is = mobility
             (prevDMT == "Injectable") + (prevDMT == "Oral") + (prevDMT == "Infusion") * 3 + 
             (prerelapse_num > 0) * 3 + 
             hospitalization * 3 # add the MS severity score (see Nickolas 2017's algorithm)
    ) 
  
  if (allvars == FALSE) {
    ms4 <- ms4 %>%  
      dplyr::select(-plan, -starts_with("pre"), prerelapse_num, prerelapse_cat, prevDMTdate, prevDMT, prevDMTefficacy, premedicalcost, premedicationcosts, premedicationcosts_cat, prega, preifn,
                    -walkingaids, -wheelchair, -hospitalbeds, -transportation, -physicalmed,
                    -demcns, -opticnervevisual, -dizzinessgiddiness, -fatiguemalaise, -neurobladder, -othermyelitis, -neuralgia) %>% 
      as.data.frame()
  } else {
    ms4 <- ms4 %>% dplyr::select(-plan, -(ends_with("date") & starts_with("pre")), -starts_with("rdate"), prevDMTdate) %>% 
      as.data.frame()
  }
  
  
  # Remove covariates that are too rare (n stays the same) 
  ms5 <- ms4 %>% dplyr::select(-snf, -chf, -mi, -moderatetosevereliverdisease) %>% # skilled nursing facility admission (SNF), only 0.6% true 
    dplyr::select(-discontinuationdate, -switchdmtdate, -lastprerelapse, -starts_with("rdate"), -prevDMTdate) %>% # Remove covariates that have missing variables + not relevant for analysis
    filter(finalpostdayscount != 0) # delete patients with 0 finalpostdayscount (n=13)
  
  
  
  ##################################################
  #### Export analytic dataset ####
  ##################################################
  
  # Create outcomes
  msexport <- ms5 %>% mutate(arr = postrelapse_num / (finalpostdayscount / 365.25),
                             relapse1year = ifelse(timetorelapse <= 365 & relapseindicator == 1, 1, 0),
                             relapse1year = factor(ifelse(relapse1year == 0 & timetorelapse <= 365, NA, relapse1year)),
                             logarr0001 = log(arr + 0.001),
                             mlogarr0001 = -logarr0001,
                             logarr1 = log(arr + 1),
                             mlogarr1 = -logarr1,
                             logarr01r = log( (postrelapse_num + 0.1) / (finalpostdayscount / 365.25)),
                             mlogarr01r = -logarr01r,
                             offset = log(finalpostdayscount / 365.25),
                             FUweight = finalpostdayscount/sum(finalpostdayscount))
  
  # Additional pre-processing to have meaningful levels/level order for factor variable
  msexport <- msexport %>% mutate(trt = factor(ifelse(treatment == "DMF", 1, 0)), 
                                  treatment = as.factor(treatment),       
                                  sex = as.factor(sex),
                                  female = as.factor(female),
                                  regioncensus = as.factor(regioncensus),
                                  regioncensusdivision = as.factor(regioncensusdivision),
                                  employmentstatus_2cat = as.factor(employmentstatus_2cat),
                                  cvd = as.factor(cvd),
                                  diabetes = as.factor(diabetes),
                                  hypertension = as.factor(hypertension),
                                  pvd = as.factor(pvd),
                                  othermedsantidepressants = as.factor(othermedsantidepressants),
                                  othermedsbenzodiazepines = as.factor(othermedsbenzodiazepines),
                                  othermedsgabaanalogs = as.factor(othermedsgabaanalogs),
                                  othermedsnsaidscox2 = as.factor(othermedsnsaidscox2),
                                  othermedsnarcotics = as.factor(othermedsnarcotics),
                                  othermedssmr = as.factor(othermedssmr),
                                  othermedssteroidsoralinjection = as.factor(othermedssteroidsoralinjection),
                                  urinarycatheter = as.factor(urinarycatheter),
                                  hospitalization = as.factor(hospitalization),
                                  scorebladder = as.factor(scorebladder),
                                  scorebrainstem = as.factor(scorebrainstem),
                                  scorecerebellar = as.factor(scorecerebellar),
                                  scorecerebral = as.factor(scorecerebral),
                                  scoregeneral = as.factor(scoregeneral),
                                  scorepyramidal = as.factor(scorepyramidal),
                                  scoresensory = as.factor(scoresensory),
                                  scorespeech = as.factor(scorespeech),
                                  scorewalking = as.factor(scorewalking),
                                  scorevisual = as.factor(scorevisual),
                                  mobility = as.factor(mobility),
                                  cci = factor(cci, levels = c("0", "1", "2", ">=3"), labels = c("0", "1", "2", ">=3")),
                                  eci = factor(eci, levels = c("0", "1", "2", "3", "4", ">=5"), labels = c("0", "1", "2", "3", "4", ">=5")),
                                  prevDMTefficacy = factor(prevDMTefficacy, levels = c("Low efficacy", "Medium efficacy", "High efficacy", "None"), labels = c("Low efficacy", "Medium efficacy", "High efficacy", "None")),
                                  prevDMT = factor(prevDMT, levels = c("Injectable", "Oral", "Infusion", "None"), labels = c("Injectable", "Oral", "Infusion", "None")),
                                  insuranceplan = factor(insuranceplan, levels = c("POS/PPO/EPO", "High Deductible", "Capitated", "other/unknown"), labels = c("POS/PPO/EPO", "High Deductible", "Capitated", "other/unknown")),
                                  employmentstatus = factor(employmentstatus, levels = c("Active Full Time", "Active Part Time or Seasonal", "Disability", "Other/Unknown/Surviving Spouse/Depend.", "Retiree/Cobra")),
                                  preifn = as.factor(preifn),
                                  prega = as.factor(prega),
                                  prerelapse_cat = factor(prerelapse_cat, levels = c("0", "1", "2+"), labels = c("0", "1", "2+")))
  
  if(allvars == TRUE){
    msexport <- msexport %>% mutate(prefingolimod = as.factor(prefingolimod),
                                    prepegifn = as.factor(prepegifn),
                                    predaclizumab = as.factor(predaclizumab),
                                    prealemtuzumab = as.factor(prealemtuzumab),
                                    prenatalizumab = as.factor(prenatalizumab),
                                    preocrelizumab = as.factor(preocrelizumab),
                                    prerituximab = as.factor(prerituximab),
                                    premitoxantrone = as.factor(premitoxantrone),
                                    demcns = as.factor(demcns),
                                    opticnervevisual = as.factor(opticnervevisual),
                                    dizzinessgiddiness = as.factor(dizzinessgiddiness),
                                    fatiguemalaise = as.factor(fatiguemalaise),
                                    neurobladder = as.factor(neurobladder),
                                    othermyelitis = as.factor(othermyelitis),
                                    neuralgia = as.factor(neuralgia),
                                    walkingaids = as.factor(walkingaids),
                                    wheelchair = as.factor(wheelchair),
                                    hospitalbeds = as.factor(hospitalbeds),
                                    transportation = as.factor(transportation),
                                    physicalmed = as.factor(physicalmed))
  }
  
  return(msexport)
  
}



compareX <- function(xcontinuous, xbinary, xcategorical, trtcbp){
  
  #######################################################
  # Comparison of baseline characteristics by treatment #
  #######################################################
  
  #' @param xcontinuous matrix, of continuous variables
  #' @param xbinary matrix, of binary variables
  #' @param xcategorical matrix, of categorical variables (> 2 categories)
  #' @param trtcbp vector, of 1 or 0
  
  mux1=apply(xcontinuous[trtcbp==1,],2,mean)
  sdx1=apply(xcontinuous[trtcbp==1,],2,sd)
  mux0=apply(xcontinuous[trtcbp==0,],2,mean)
  sdx0=apply(xcontinuous[trtcbp==0,],2,sd)
  diffx=(mux1-mux0)/sqrt((sdx1^2+sdx0^2)/2)
  p.ttest=rep(NA, dim(xcontinuous)[2])
  for(b in 1:dim(xcontinuous)[2])
    p.ttest[b]=round(t.test(xcontinuous[,b]~trtcbp)$p.value, 3)
  comp.xcontinuous=cbind(mux1, sdx1, mux0, sdx0, diffx, p.ttest)
  rownames(comp.xcontinuous)=colnames(xcontinuous)  
  colnames(comp.xcontinuous) <- c("mean1", "sd1", "mean0", "sd0", "SMD", "p.value")
  
  propx1=apply(xbinary[trtcbp==1,],2,mean)
  countpx1=apply(xbinary[trtcbp==1,],2,sum)
  propx0=apply(xbinary[trtcbp==0,],2,mean)
  countpx0=apply(xbinary[trtcbp==0,],2,sum)
  diffpx=(propx1-propx0)/sqrt((propx1*(1-propx1)+propx0*(1-propx0))/2)
  p.ftest=rep(NA, dim(xbinary)[2])
  for(b in 1:dim(xbinary)[2])
  {D=matrix(0, 2, 2)
  D[1,1]=sum(xbinary[,b]==1 & trtcbp==1) 
  D[1,2]=sum(xbinary[,b]==1 & trtcbp==1)
  D[2,1]=sum(xbinary[,b]==0 & trtcbp==0)
  D[2,2]=sum(xbinary[,b]==0 & trtcbp==0)
  if (any(D < 5)){
    p.ftest[b] = round(fisher.test(D)$p.value, 3)
  } else {
    p.ftest[b]=round(chisq.test(D)$p.value, 3)
  }
  }
  comp.xbinary=cbind(countpx1, propx1, countpx0, propx0, diffpx, p.ftest)
  rownames(comp.xbinary)=colnames(xbinary)
  colnames(comp.xbinary) <- c("count1", "prop1", "count0", "prop0", "SMD", "p.value")
  
  comp.xcategorical = NULL
  ks = NULL
  tabx1 = apply(xcategorical[trtcbp == 1, ], 2, propTable) 
  tabx0 = apply(xcategorical[trtcbp == 0, ], 2, propTable)
  for (var in colnames(xcategorical)){
    p1 = as.matrix(tabx1[[var]])
    p0 = as.matrix(tabx0[[var]])
    diff = (p1-p0)[-1]  # T - C, for k-1 categories (removing the 1st category in alphabetical order)
    S = ((p1[-1] %*% t(p1[-1])) + (p0[-1] %*% t(p0[-1]))) / 2 # when k != l, []
    diag(S) = (p1[-1] * (1-p1[-1]) + p0[-1] * (1-p0[-1])) /2
    d = sqrt(t(diff) %*% solve(S) %*% diff)  # d = sqrt((T-C)' S^{-1} (T-C))
    # multivariate Mahalanobis distance method (Dalton 2008, A new standardized difference metric for multinomial samples)
    tab = table(xcategorical[,var], trt)
    if (any(tab < 5)){
      p.chitest = round(fisher.test(tab)$p.value, 4) # p-values come from chi-squared test
    } else {
      p.chitest = round(chisq.test(tab)$p.value, 4)
    }
    k = length(unique(xcategorical[,var]))- 1 # number of categories - 1
    ks = c(ks, k)
    comp.xcategorical = rbind(comp.xcategorical, cbind(rownames(tab)[-1], p1[-1], p0[-1], diff, rep(d, k), rep(p.chitest, k)))
  }
  rownames(comp.xcategorical) <- rep(colnames(xcategorical), ks)
  colnames(comp.xcategorical) <- c("category", "prop1", "prop0", "diff", "SMD", "p.value")
  
  return(list(comp.xcontinuous, comp.xbinary, comp.xcategorical))
}

