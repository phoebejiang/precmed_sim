# ------------------------------------------------------------------
# Product: Tecfidera [A1] vs. Teriflunomid [A0]
# Protocol: MarketScan
# Project: Precision Medicine MS
# 
# Program name: simplot.R
# Developer/Programmer: pj
# Date: 04
# 
# Purpose: Collect CV iteration batch results and summarize the mean and sd of value, agreement and
#          accuracy across different sample sizes, levels of heterogeneity, and percentiles.
#          Run this after main.sh and summary.R finishes for all methods and n's.
# 
# Platform: Windows
# R Version: 4.0.3
# 
#   Modifications:
# 
#   Date			By			Description
# --------		--------	-----------------------------
#  22APR2024 pj      Start the main structure (based on manuscript_final.R)
#  
# ------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(haven))
suppressPackageStartupMessages(library(cowplot))


source("C:/Users/xjiang1/OneDrive - Biogen/Documents/Innovation/PMMS/TruvenMarketScan/manuscript_func_final.R")
setwd("C:/Users/xjiang1/OneDrive - Biogen/Documents/Innovation/MarketScan/Outputs/") # CHANGE THIS #

colB <- brewer.pal(n = 9, name = "PuBu") # Blue for estimated value
colR <- brewer.pal(n = 9, name = "OrRd") # Red for value
colG <- brewer.pal(n = 9, name = "BuGn") # Green for level of heterogeneity


# Sample size and HTE with Value Function and Value Ratio as Y

plotV <- function(data, thisV, display_methods, outname,
                  output = F, proportion = F,
                  figure_folder = "C:/Users/xjiang1/OneDrive - Biogen/Documents/Innovation/MarketScan/Manuscripts/Figures/"){
  if (thisV %in% c("V", "V.hat")){
    data <- data %>% filter(whichV == thisV) 
    thisy <- "meanV"
    thisVlab <- "Value Function"
  } 
  if (thisV %in% c("VR_V", "VR_V.hat")) {
    thisy = thisV
    thisVlab = "Value Ratio"
  }
  
  if (proportion == F){
    p <- data %>%
      filter(method %in% c(display_methods, "All A1")) %>%
      # filter(method %in% display_methods) %>%
      ggplot(aes(x = n.subset, y = get(thisy), color = levhte, group = levhte)) +
      geom_point(size = 3, shape = 16) + 
      #aes(shape = method)) + 
      geom_line(linewidth = 1) + 
      scale_color_manual(values = rev(colG)) +
      facet_grid( ~ method) + 
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
  } else {
    p <- data %>%
      filter(method %in% c(display_methods, "All A1")) %>%
      # filter(method %in% display_methods) %>%
      ggplot(aes(x = n.subset, y = get(thisy), color = levhte, group = levhte)) +
      geom_point(size = 3, shape = 16) + 
      #aes(shape = method)) + 
      geom_line(linewidth = 1) + 
      scale_color_manual(values = c(colG[c(3, 5, 7, 9)])) +
      facet_grid(proportion ~ method) + 
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
  }
  # Add 0 and 1 ref lines if VR
  if (thisV %in% c("VR_V", "VR_V.hat")){
    p <- p + 
      geom_hline(yintercept = 1, color = "gray", linetype = 2, linewidth = 1) + 
      geom_hline(yintercept = 0, color = "gray", linetype = 2, linewidth = 1)
  } 
  # Add error bars if V
  if (thisV %in% c("V(dhat)", "Vhat(dhat)")){
    p = p + geom_errorbar(aes(ymin = meanV - sdV/sqrt(n.nonna), 
                              ymax = meanV + sdV/sqrt(n.nonna)), width = 0.1, size = 1)
  }
  
  if (output) png(paste0(figure_folder, outname), 
                  height = 5, width = 12, units = "in", res = 300)
  print(p)
  if (output) dev.off()
  # if (output) ggsave(paste0(figure_folder, outname), res = 300)
  return(p)
}

mergeplot <- function(data, display_methods, outname, output = F){
  figa <- plotV(data = data, thisV = "V(dhat)", 
                display_methods = display_methods,
                outname = outname, output = output)
  
  figb <- plotV(data = data, thisV = "VR_V", 
                display_methods = display_methods,
                outname = outname, output = output)
  
  fig <- plot_grid(fig2, fig2 + theme(legend.position = "none"), nrow = 2)
  
}

plotVvsVhat <- function(data, display_methods, figure_folder, outname, output){
  p <- data %>%
    filter(method %in% c(display_methods, "All A1")) %>%
    ggplot(aes(x = n.subset, y = meanV, color = whichV)) +
    geom_point(size = 3) + 
    geom_line(aes(group = whichV)) +
    geom_errorbar(aes(ymin = meanV - sdV, ymax = meanV + sdV), linewidth = 0.5, width = 0.2) +
    facet_grid(levhte ~ method, scales = "free") +
    scale_x_continuous(breaks = ns, trans = "log") +
    theme_bw() + 
    scale_color_manual(values = c(colR[8], colR[5])) +
    ylim(c(0, 0.6)) +
    labs(x = "Sample size (log-transformed)",
         y = "Value Function +- SE",
         title = "Estimated vs true cross-validated value function of ITR estimated from selected PM methods") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(angle = 315, size = 12),
          strip.text.x = element_text(size = 13),
          legend.text = element_text(size = 13),
          legend.position = "bottom",
          legend.title = element_blank()) 
  
    if (output) png(paste0(figure_folder, outname), 
                    height = 5, width = 12, units = "in", res = 300)
    print(p)
    if (output) dev.off()
    return(p)
}

# mean accuracy HTE, sample size on x-axis, level of heterogeneity panel, color by method

plotAcc <- function(data, display_methods, figure_folder, outname, output){
  mygrey <- rgb(124, 135, 142, maxColorValue = 255)
  p <- data %>% 
    filter(method %in% c(display_methods, "All A1")) %>%
    ggplot(aes(x = n.subset, y = mean.acc, group = method, color = method)) +
    geom_point(size = 2, shape = 16) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 1, color = mygrey) +
    geom_hline(yintercept = 0.5, linetype = 3, linewidth = 1, color = mygrey) +
    scale_color_brewer(palette = "Paired") +
    scale_x_continuous(breaks = ns, trans = "log") + 
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~ levhte, nrow = 1) + 
    theme_bw() + 
    labs(x = "Sample Size (log-transformed)", y = "Mean accuracy") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(angle = 315, size = 12),
          strip.text.x = element_text(size = 13), 
          legend.position = "bottom",
          legend.title = element_blank()) 
  
  if (output) png(paste0(figure_folder, outname), 
                  height = 4, width = 8, units = "in", res = 300)
  print(p)
  if (output) dev.off()
  return(p)
}

plotAgree <- function(data, display_methods, figure_folder, outname, output){
  # Remove rows of all NA values because the methods were not selected
  dat2 <- data[rowSums(!is.finite(data)) != nrow(data), rowSums(!is.finite(data)) != nrow(data)]
  
  # Plot correlations
  p <- corrplot(dat2,
                method = "color", type = "lower", 
                addCoef.col = "orange",
                tl.cex = 1, cl.cex = 1, tl.col = "black", tl.srt = 30) # saved 565 (height) X 550
  
  if (output) png(paste0(figure_folder, outname), 
                  height = 5, width = 12, units = "in", res = 300)
  print(p)
  if (output) dev.off()
  return(p)
}
