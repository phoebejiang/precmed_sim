# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: simulation_analysis_functions.R
# 
# Purpose: Define functions needed in "simulation_analysis.R"
# ------------------------------------------------------------------

library(tidyverse)
library(magrittr)

merge_samplesize_outputs <- function(noarg = "noarg"){
  # All sample size outputs have 20% percentiles as proportion
 
  # Read in RData spaces for level of heterogeneity simulations
  for (levhte in c("no", "low", "medium", "high")) {
    for (pergr in c("20x5")) {
      file_path <- paste0("./intermediate_results/simulations_samplesize/simsummary_sample_size_dataplot_", 
                          levhte, "_", pergr, ".RData")
          # Originally outputted here: "./simulations/samplesize/outputs/" from "simsummary_sample_size.R"
      try({
        if (file.exists(file_path)) {
          load(file_path)
          message("Loaded: ", file_path) 
        } else {
          message("File not found: ", file_path)
        }
      }, silent = TRUE)
    }
  }
  
  # Merge all results 
  all.cvs3 <- do.call("rbind",  mget(ls()[grepl("all.cvs3_", ls())]))
  all.cvs3.wide <- do.call("rbind",  mget(ls()[grepl("all.cvs3.wide_", ls())]))
  dat.accuracy <- do.call("rbind",  mget(ls()[grepl("dat.accuracy_", ls())]))
  dat.accuracy.summary <- do.call("rbind",  mget(ls()[grepl("dat.accuracy.summary_", ls())]))
  
  # Label levhte for plots
  all.cvs3 %<>%
    dplyr::mutate(whichV = factor(whichV, levels = c("V", "V.hat"), labels = c("V(dhat)", "Vhat(dhat)"))) %>% 
    dplyr::mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                           labels = c("No", "Low", "Medium", "High", "Medium", "Medium")),
           symm = factor(percentiles, levels = c("seq(0,1,by=0.2)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                         labels = c("Symmetric", "Asymmetric", "Asymmetric", "Symmetric", "Symmetric")),
           proportion = factor(percentiles, levels = c("seq(0,1,by=0.2)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                               labels = c("20%-20%-20%-20%-20%", "55%-15%-15%-15%", "55%-30%-15%", "10%-15%-50%-15%-10%", "10%-30%-20%-30%-10%"))
    )
  
  all.cvs3.wide %<>%
    dplyr::mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                           labels = c("No", "Low", "Medium", "High", "Medium", "Medium")),
           symm = factor(percentiles, levels = c("seq(0,1,by=0.2)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                         labels = c("Symmetric", "Asymmetric", "Asymmetric", "Symmetric", "Symmetric")),
           proportion = factor(percentiles, levels = c("seq(0,1,by=0.2)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                               labels = c("20%-20%-20%-20%-20%", "55%-15%-15%-15%", "55%-30%-15%", "10%-15%-50%-15%-10%", "10%-30%-20%-30%-10%"))
    )
  
  dat.accuracy %<>% mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                                           labels = c("No", "Low", "Medium", "High", "Medium", "Medium")),
                           symm = factor(percentiles, levels = c("seq(0,1,by=0.2)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                                         labels = c("Symmetric", "Asymmetric", "Asymmetric", "Symmetric", "Symmetric")),
                           pergr = factor(percentiles, levels = c("c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)", "seq(0,1,by=0.2)"),
                                          labels = c("asymm 55-15-15-15", "asymm 55-30-15", "symm 10-15-50-15-10", "symm 10-30-20-30-10", "symm 20-20-20-20-20"))
  )
  
  dat.accuracy.summary %<>% mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                                                   labels = c("No", "Low", "Medium", "High", "Medium", "Medium")),
                                   symm = factor(percentiles, levels = c("seq(0,1,by=0.2)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                                                 labels = c("Symmetric", "Asymmetric", "Asymmetric", "Symmetric", "Symmetric")),
                                   pergr = factor(percentiles, levels = c("c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)", "seq(0,1,by=0.2)"),
                                                  labels = c("asymm 55-15-15-15", "asymm 55-30-15", "symm 10-15-50-15-10", "symm 10-30-20-30-10", "symm 20-20-20-20-20"))
  )
  return(list(all.cvs3 = all.cvs3, all.cvs3.wide = all.cvs3.wide, 
              dat.accuracy = dat.accuracy, dat.accuracy.summary = dat.accuracy.summary))
}

merge_proportion_outputs <- function(selected_methods){
  
  #######################################################################
  ############################# Merge data ############################
  #######################################################################
  allns <- c(500, 1000, 2500, 5000, 10000) # c(500, 1000, 2500, 5000, 10000)
  percentile_group <- list(c(0,0.55,0.65,0.75,0.85,1), c(0, 0.55, 0.65, 0.7, 0.85, 1), c(0, 0.1, 0.4, 0.6, 0.9, 1), c(0, 0.1, 0.25, 0.75, 0.9, 1), c(0,0.2,0.4,0.6,0.8,1))
  beta_group <- list(c(-0.92,-0.69,-0.69,-0.69,0), c(-0.92,-0.69,-0.69,0,0.1), c(-0.92,-0.69,0,0.1,0.18), c(-0.92,-0.69,0,0.1,0.18), c(-0.92,-0.69,0,0.1,0.18))
  
  # Run the following code one by one with i = 1 - 5
  #### ---- Load all formatted data sets in RData  ---- ####
  cols1 <- c("n.subset", "n.nonna", "n.batches", "beta", "percentiles", "pergr", "method", "outcome", "VR_V.hat", "VR_V", "meanV_V", "sdV_V")
  cols2 <- c("n.subset", "beta", "percentiles", "method", "pergr", "levhte", "symm", "mean.acc", "sd.acc", "q1.acc", "q3.acc", "mean.acc.wo.neutral", "sd.acc.wo.neutral", "q1.acc.wo.neutral", "q3.acc.wo.neutral")
  output1 <- matrix(NA, nrow = length(allns)*length(percentile_group)*10, ncol = length(cols1)) %>% as_tibble()
  output2 <- matrix(NA, nrow = length(allns)*length(percentile_group)*10, ncol = length(cols2)) %>% as_tibble()
  colnames(output1) <- cols1
  colnames(output2) <- cols2

  for (j in 1:length(allns)){
    for (i in 1:length(percentile_group)){ # 1, 2 asymmetric, 3, 4 symmetric, 5 all equal 20% 
      ns <- allns[j]
      beta <- beta_group[[i]] # {rep(-0.2, 5), c(-1.2,-0.69,0,0.1,0.41), c(-0.36,-0.29,0,0.05,0.1), c(-0.92, -0.69, 0, 0.1, 0.18)}
      percentiles <- percentile_group[[i]] # {seq(0, 1, by = 0.2), }
      cat(ns, "||", i, "||", beta, "||", percentiles, "\n")
      
      levhte <- "medium"
      if(all(percentiles == c(0, 0.55, 0.65, 0.75, 0.85, 1))) {
        pergr <- "asymm 55-30-15"
      } else if(all(percentiles == c(0, 0.55, 0.65, 0.7, 0.85, 1))) {
        pergr <- "asymm 55-15-15-15"
      } else if(all(percentiles == c(0, 0.1, 0.4, 0.6, 0.9, 1))) {
        pergr <- "symm 10-30-20-30-10"
      } else if(all(percentiles == c(0, 0.1, 0.25, 0.75, 0.9, 1))) {
        pergr <- "symm 10-15-50-15-10"
      } else if(all(percentiles == c(0, 0.2, 0.4, 0.6, 0.8, 1))) {
        pergr <- "symm 20-20-20-20-20"
      }
      
      file_path <- paste0("./intermediate_results/simulations_proportions/simsummary_n", ns, "_dataplot_", 
                          levhte, "_", pergr, ".RData")
          # Originally outputted here: "./simulations/proportion/outputs/" from "simsummary_proportion.R" 
      try({
        if (file.exists(file_path)) {
          load(file_path)
          message("Loaded: ", file_path)
      
          assign("all.cvs3.wide", get(paste0("all.cvs3.wide_n", ns, "_", levhte, "_", pergr)))
          assign("dat.accuracy.summary", get(paste0("dat.accuracy.summary_n", ns, "_", levhte, "_", pergr)))
          
          out <- all.cvs3.wide %>% mutate(pergr = pergr) %>% dplyr::select(all_of(cols1)) %>% filter(method != "List DTR\n (3 nodes)")
          while (nrow(out) < 10) out %<>% rbind(NA) # in case there are methods that did not run successfully (e.g., listDTR) but we still need the remaining results
          start <- (j-1)*length(percentile_group)*10 + (i-1)*10 + 1
          end <- start + 10 - 1
          output1[start:end, ] <- out
          rm(all.cvs3.wide)
          
          out <- dat.accuracy.summary %>% 
            mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                                   labels = c("No", "Low", "Medium", "High", "Medium", "Medium")),
                   symm = factor(percentiles, levels = c("seq(0,1,by=0.2)", "c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)"),
                                 labels = c("Symmetric", "Asymmetric", "Asymmetric", "Symmetric", "Symmetric")),
                   pergr = factor(percentiles, levels = c("c(0,0.55,0.65,0.7,0.85,1)", "c(0,0.55,0.65,0.75,0.85,1)", "c(0,0.1,0.25,0.75,0.9,1)", "c(0,0.1,0.4,0.6,0.9,1)", "seq(0,1,by=0.2)"),
                                       labels = c("asymm 55-15-15-15", "asymm 55-30-15", "symm 10-15-50-15-10", "symm 10-30-20-30-10", "symm 20-20-20-20-20"))
            ) %>% dplyr::select(all_of(cols2)) %>% filter(method != "List DTR\n (3 nodes)")
          while (nrow(out) < 10) out %<>% rbind(NA) # in case there are methods that did not run successfully (e.g., listDTR) but we still need the remaining results
          start <- (j-1)*length(percentile_group)*10 + (i-1)*10 + 1
          end <- start + 10 - 1
          output2[start:end, ] <- out
          rm(dat.accuracy.summary)
          
        } else {
          message("File not found: ", file_path)
        }
      }, silent = TRUE)
    }
  }
  
  # Format percentile variables
  output1 %<>% 
    mutate(symmetry = ifelse(str_detect(pergr, "asymm"), "Asymmetric", "Symmetric"),
           perc_group = case_when(
             pergr == "asymm 55-30-15" ~ "55%, 30%, 15%",
             pergr == "asymm 55-15-15-15" ~ "55%, 15%, 15%, 15%",
             pergr == "symm 10-30-20-30-10" ~ "10%, 30%, 20%, 30%, 10%",
             pergr == "symm 10-15-50-15-10" ~ "10%, 15%, 50%, 15%, 10%",
             pergr == "symm 20-20-20-20-20" ~ "20%, 20%, 20%, 20%, 20%"
           ),
           neutral = case_when(
             pergr == "asymm 55-30-15" ~ "15%",
             pergr == "asymm 55-15-15-15" ~ "15%",
             pergr == "symm 10-30-20-30-10" ~ "20%",
             pergr == "symm 10-15-50-15-10" ~ "50%",
             pergr == "symm 20-20-20-20-20" ~ "20%"
           ),
           a1responder = case_when(
             pergr == "asymm 55-30-15" ~ "85%",
             pergr == "asymm 55-15-15-15" ~ "70%",
             pergr == "symm 10-30-20-30-10" ~ "40%",
             pergr == "symm 10-15-50-15-10" ~ "25%",
             pergr == "symm 20-20-20-20-20" ~ "40%"
           )) %>% 
    mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                           labels = c("No", "Low", "Medium", "High", "Medium", "Medium"))
    )
  
  output2 %<>% 
    mutate(symmetry = ifelse(str_detect(pergr, "asymm"), "Asymmetric", "Symmetric"),
           perc_group = case_when(
             pergr == "asymm 55-30-15" ~ "55%, 30%, 15%",
             pergr == "asymm 55-15-15-15" ~ "55%, 15%, 15%, 15%",
             pergr == "symm 10-30-20-30-10" ~ "10%, 30%, 20%, 30%, 10%",
             pergr == "symm 10-15-50-15-10" ~ "10%, 15%, 50%, 15%, 10%",
             pergr == "symm 20-20-20-20-20" ~ "20%, 20%, 20%, 20%, 20%"
           ),
           neutral = case_when(
             pergr == "asymm 55-30-15" ~ "15%",
             pergr == "asymm 55-15-15-15" ~ "15%",
             pergr == "symm 10-30-20-30-10" ~ "20%",
             pergr == "symm 10-15-50-15-10" ~ "50%",
             pergr == "symm 20-20-20-20-20" ~ "20%"
           ),
           a1responder = case_when(
             pergr == "asymm 55-30-15" ~ "85%",
             pergr == "asymm 55-15-15-15" ~ "70%",
             pergr == "symm 10-30-20-30-10" ~ "40%",
             pergr == "symm 10-15-50-15-10" ~ "25%",
             pergr == "symm 20-20-20-20-20" ~ "40%"
           )) %>% 
    mutate(levhte = factor(beta, levels = c("c(-0.2,-0.2,-0.2,-0.2,-0.2)", "c(-0.36,-0.29,0,0.05,0.1)", "c(-0.92,-0.69,0,0.1,0.18)", "c(-1.2,-0.69,0,0.1,0.41)", "c(-0.92,-0.69,-0.69,0,0.1)", "c(-0.92,-0.69,-0.69,-0.69,0)"),
                           labels = c("No", "Low", "Medium", "High", "Medium", "Medium"))
    )
  
  return(list(all.cvs3.wide_prop = output1, dat.accuracy.summary_prop = output2))
}


plotV <- function(data, 
                  thisV, 
                  display_methods, 
                  colorgroup,
                  outname = NULL,
                  output = F, 
                  manuscript_figure_folder = "./results/figures/"){
  
  # Which value are we presenting?
  if (thisV %in% c("meanV_V", "meanV_V.hat")){
    thisy <- thisV; thisVlab <- "Value Function"
  }
  if (thisV %in% c("VR_V", "VR_V.hat")) {
    thisy = thisV; thisVlab = "Value Ratio"
  }
  
  # Conditional color scale, depending on what the color group is
  color_scale <- if (colorgroup == "levhte") {
    scale_color_manual(values = colG[c(3, 5, 7, 9)])
  } else if (colorgroup == "pergr") {
    scale_color_manual(values = colB[c(3, 4, 5, 7, 9)])
  } 
  
  # Make the plot
  p <- data %>%
    filter(method %in% c(display_methods, "All A1"),
           !is.na(get(thisy))) %>%
    ggplot(aes(x = n.subset, y = get(thisy), color = get(colorgroup), group = get(colorgroup))) +
    geom_point(size = 3, shape = 16) + 
    geom_line(linewidth = 1) + 
    facet_grid( ~ method) + 
    color_scale +
    scale_x_continuous(breaks = ns, trans = "log") + 
    theme_bw() + 
    labs(x = "Sample size (log-transformed)", y = thisVlab) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(angle = 315, size = 12),
          strip.text.x = element_text(size = 13),
          legend.text = element_text(size = 13),
          legend.position = "top",
          legend.title = element_blank()) 
  
  # Add horizontal reference lines 0 and 1 if needed
  if (thisV %in% c("VR_V", "VR_V.hat") & colorgroup == "levhte"){
    p <- p +
      geom_hline(yintercept = 1, color = "gray", linetype = 2, linewidth = 1) +
      geom_hline(yintercept = 0, color = "gray", linetype = 2, linewidth = 1)
  }
  
  # Output the figure
  if (output) png(paste0(manuscript_figure_folder, outname), height = 8, width = 19)
  print(p)
  if (output) ggsave(paste0(manuscript_figure_folder, outname), height = 8, width = 19)
  if (output) dev.off()
  return(p)
}


calculateATE <- function(percentiles, 
                         beta.matrix) {
  
  #' Function to calculate ATE by subgroup and overall for different percentiles
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @param beta.matrix a matrix of beta coefficient to supply to argument beta in simdata() 
  #'                    each row is level of heterogeneity: no, low, medium, high
  #'                    each column is one of the 5 stratums: high responder to A1, moderate responder to A1, neutral, moderate responder to A0, high responder to A0
  
  # Outcome model to estimate the ATE
  if (sum(percentiles == c()) == 0) {
    cate.mod <- as.formula(postrelapse_num ~ trt + ageatindex_centered + female + prerelapse_num + premedicalcost + offset(offset))
  } else {
    cate.mod <- as.formula(postrelapse_num ~ trt + ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + offset(offset))
  }
  
  # Generate data for large sample size (n=1,000,000)
  # The constant 4 represents four magnitudes of HTE (no, low, medium, high)
  ATE <- matrix(NA, nrow = 4)
  ATE.var <- matrix(NA, nrow = 4)
  for(i in 1:4){
    data <- simdata(n = 1000000, RCT = T, originaldata = NULL, seed = 999,
                    beta = beta.matrix[i, ],
                    # beta.x = c(-1.54, -0.01, 0.06, 0.47, 0.68, 0.13, 0.0000003), 
                    beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003),
                    percentiles = percentiles)$data
  
    fit <- glm(cate.mod, family = "poisson", data = data)
    summ <- summary(fit)$coefficients[2,]
    ATE[i] <- summ["Estimate"] %>% exp() %>% unname() %>% round(2)
    ATE.var[i] <- var(log(data$rate))
  }
  
  return(list(ATE = ATE, ATE.var = ATE.var))
}



