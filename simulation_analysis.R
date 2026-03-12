# ------------------------------------------------------------------
# Project: Precision Medicine MS
# 
# Program name: simulation_analysis.R
#
# Purpose: Final version of all figures after running all simulation scripts
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
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(cowplot))

source("./utility.R")
source("./simulation_analysis_functions.R")

# User specified
percentiles <- seq(0, 1, by = 0.2) 
ns <- c(500, 1000, 2500, 5000, 10000)
nrepeats <- 25 # number of CV repetitions
display_methods <- c("Poisson", "Boosting", 
                     "Contrast\n Regression", "List DTR\n (2 nodes)")
base.seed = 999
manuscript_figure_folder <- "./results/figures/"

# Formatting
myblue <- rgb(37, 15, 186, maxColorValue = 255)
mygrey <- rgb(124, 135, 142, maxColorValue = 255)
mygrey <- "lightgrey"

colB <- brewer.pal(n = 9, name = "PuBu") # Blue for estimated value
colR <- brewer.pal(n = 9, name = "OrRd") # Red for value
colG <- brewer.pal(n = 9, name = "BuGn") # Green for level of heterogeneity


# Sample size and HTE magnitude outputs
samplesize_outputs <- merge_samplesize_outputs()
all.cvs3_samplesize <- samplesize_outputs$all.cvs3
all.cvs3.wide_samplesize <- samplesize_outputs$all.cvs3.wide
dat.accuracy_samplesize <- samplesize_outputs$dat.accuracy
dat.accuracy.summary_samplesize <- samplesize_outputs$dat.accuracy.summary

# Sample size and HTE distribution outputs
prop_outputs <- merge_proportion_outputs(selected_methods = display_methods)
all.cvs3.wide_prop <- prop_outputs$all.cvs3.wide_prop
dat.accuracy.summary_prop <- prop_outputs$dat.accuracy.summary_prop


#######################################################################
##########      Figure 1: Baseline characteristics   ##################
#######################################################################

# Simulate a large test set 
sim.big <- simdata(n = 100000, RCT = T, beta = c(-0.92,-0.69,0,0.1,0.18), seed = base.seed, 
                   percentiles = c(0,0.1,0.25,0.75,0.9,1))$data %>%
  mutate(ageatindex = ageatindex_centered + 48,
         prevDMTefficacy = str_replace(prevDMTefficacy, " efficacy", ""))
cat("\nA random sample is simulated with seed", base.seed, "with dimension: ", dim(sim.big), "as the one independent large test data for calculation of ture value function.\n")

# Figure 1
fig1a <- sim.big %>%
  ggplot(aes(x = ageatindex)) +
  geom_density(col = "navy", linewidth = 1) + 
  geom_vline(aes(xintercept = mean(ageatindex)), color = "darkorange", linetype = "dashed", linewidth = 1.5) +
  theme_minimal() + 
  labs(x = "Age at Index (year)", y = "Density")

fig1b <- sim.big %>%
  ggplot(aes(x = premedicalcost)) + 
  geom_density(col = "navy", linewidth = 1) + 
  geom_vline(aes(xintercept = mean(premedicalcost)), color = "darkorange", linetype = "dashed", linewidth = 1.5) +
               # , label = after_stat(paste(mean, "±", sd)))) + 
  scale_x_continuous(trans = "log10") + 
  theme_minimal() + 
  labs(x = "Pre-index medical cost ($), log-transformed", y = "Density")

pct_format = scales::percent_format(accuracy = 1)
fig1c <- ggplot(sim.big %>% 
                  dplyr::select(female, prerelapse_num,
                         prevDMTefficacy, numSymptoms, ID) %>% melt(id = "ID"), 
       aes(x = value)) +  
  facet_grid(~ variable, scales = "free", 
             labeller = labeller(variable = c("female" = "Female",
                                              "prerelapse_num" = "Number of relapses in pre-index",
                                              "prevDMTefficacy" = "Efficacy of previous DMT",
                                              "numSymptoms" = "Number of symptoms in pre-index"))) + 
  geom_bar(aes(y = after_stat(count / ave(count, PANEL, FUN = sum))), 
           fill = "darkorange", color = "navy") + 
  geom_text(
    aes(
      y = after_stat(count / ave(count, PANEL, FUN = sum)),
      label = sprintf(
        '%s',
        pct_format(after_stat(count / ave(count, PANEL, FUN = sum)))
      )
    ),
    stat = 'count', nudge_y = 0.03, colour = 'navy', fontface = "bold"
  ) + 
  theme_minimal() + 
  labs(x = "", y = "Proportion") +
  theme(legend.position = "None")

fig1 <- cowplot::plot_grid(fig1c, cowplot::plot_grid(fig1a, fig1b, nrow = 1), nrow = 2)
fig1
ggsave(paste0(manuscript_figure_folder, "Figure1.tiff"), height = 8, width = 10, dpi = 400)

#######################################################################
##########      Figure 2: Describe data generating mechanism     ######
#######################################################################

# Generate very large data set to get true ATE and ATE by subgroup
beta <- matrix(NA, nrow = 4, ncol = 5)
beta[1, ] <- rep(-0.2, 5)
beta[2, ] <- c(-0.36, -0.29, 0, 0.05, 0.1)
beta[3, ] <- c(-0.92, -0.69, 0, 0.1, 0.18) 
beta[4, ] <- c(-1.2, -0.69, 0, 0.1, 0.41) 

cate.mod <- as.formula(postrelapse_num ~ trt + ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost + offset(offset))
RR <- matrix(NA, nrow = 4, ncol = 6)

for(i in 1:4){
  data <- simdata(n = 1000000, RCT = TRUE, originaldata = NULL, seed = 999,
                  beta = beta[i, ],
                  beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003),
                  percentiles = seq(0, 1, by = 0.2))$data
  
  RR[i, 6] <- exp(coef(glm(cate.mod, family = "poisson", data = data))[2]) # ATE
  RR[i, 1] <- exp(coef(glm(cate.mod, family = "poisson", data = data, subset = Iscore == 1))[2]) # ATE1
  RR[i, 2] <- exp(coef(glm(cate.mod, family = "poisson", data = data, subset = Iscore == 2))[2]) # ATE2
  RR[i, 3] <- exp(coef(glm(cate.mod, family = "poisson", data = data, subset = Iscore == 3))[2]) # ATE3
  RR[i, 4] <- exp(coef(glm(cate.mod, family = "poisson", data = data, subset = Iscore == 4))[2]) # ATE4
  RR[i, 5] <- exp(coef(glm(cate.mod, family = "poisson", data = data, subset = Iscore == 5))[2]) # ATE5
  
  print(i)
}

magnitudeHTE <- data.frame(RR = as.vector(t(RR[, 1:5])), 
                       scenario = factor(rep(seq(1, 4), each = 5), 
                                         labels = c("No", "Low", "Medium", "High")), 
                       # level = factor(rep(c("No HTE", "Low HTE", "Medium HTE", "High HTE"), each = 5), levels = c("No HTE", "Low HTE", "Medium HTE", "High HTE")),
                       group = factor(rep(c("High A1", "Moderate A1", "Neutral", "Moderate A0", "High A0"), 4), 
                                      levels = c("High A1", "Moderate A1", "Neutral", "Moderate A0", "High A0")),
                       distribution = "Proportion of responder strata: equal symm 20-20-20-20-20")

ATE <- data.frame(RR = as.vector(t(RR[, 6])), 
                  scenario = factor(rep(seq(1, 4)), labels = c("No", "Low", "Medium", "High")), 
                  group = "ATE") # overall RR

print(ATE)

fig2a <- magnitudeHTE %>%
  ggplot(aes(y = RR, x = group, group = scenario, color = scenario)) + 
  geom_point(size = 4, shape = 22) +
  geom_line(aes(y = RR, x = group), size = 1.1) +
  geom_hline(yintercept = 1, size = 1, color = "darkgrey") +
  facet_wrap(~ distribution, nrow = 1) +
  scale_color_manual(values = c(colG[c(5, 7, 8, 9)])) + 
  theme_bw() +
  xlab("Stratum") +
  ylab("Rate ratio A1/A0") +
  labs(color = "Magnitude of HTE") + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(angle = 315, hjust = 0.95, vjust = 0.2),
        legend.position = "top")

pergr <- c("symm 20-20-20-20-20", 
           "symm 10-30-20-30-10",
           "symm 10-15-50-15-10",
           "asymm 55-30-15",
           "asymm 55-15-15-15")

distributionHTE <- data.frame(scenario = "Magnitude of HTE: Medium", 
                              group = rep(unique(magnitudeHTE$group), 
                                          times = length(pergr)),
                              percentiles = c(20, 20, 20, 20, 20,
                                              10, 30, 20, 30, 10,
                                              10, 15, 50, 15, 10,
                                              55, 30, 15, 0, 0, 
                                              55, 15, 15, 15, 0
                                              ),
                              pergr = rep(pergr, each = 5)) %>%
  mutate(pergr = case_when(
    pergr == "symm 20-20-20-20-20" ~ "equal symm 20-20-20-20-20",
    pergr == "symm 10-30-20-30-10" ~ "unequal symm2 10-30-20-30-10",
    pergr == "symm 10-15-50-15-10" ~ "unequal symm1 10-15-50-15-10",
    pergr == "asymm 55-30-15" ~ "asymm1 55-30-15",
    pergr == "asymm 55-15-15-15" ~ "asymm2 55-15-15-15"
  )) %>% 
  mutate(pergr = factor(pergr, levels = c("asymm2 55-15-15-15",
                                          "asymm1 55-30-15",
                                          "equal symm 20-20-20-20-20",
                                          "unequal symm2 10-30-20-30-10",
                                          "unequal symm1 10-15-50-15-10")))


fig2b <- distributionHTE %>% 
  ggplot(aes(x = pergr, y = percentiles, fill = group)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~ scenario, nrow = 1) + 
  scale_fill_manual(values = c(colR[9], colR[5], "grey", colB[5], colB[9])) +
  scale_y_reverse() +
  coord_flip() + 
  theme_bw() +
  labs(fill = "Stratum", x = "", y = "Proportion of patients") + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        legend.position = "top")
  
fig2 <- cowplot::plot_grid(fig2a, fig2b, nrow = 2, rel_heights = c(2, 1))
fig2
ggsave(paste0(manuscript_figure_folder, "Figure2.tiff"), height = 10, width = 10, dpi = 400)


#######################################################################
##########   Section 3.4 Results in Text    ###########################
#######################################################################

# Section 3.4 Results
# Results about Figure 2, n = 500
## All A1, no HTE
all1.500.no <- all.cvs3.wide_samplesize %>% filter(method == "All A1", n.subset == 500, levhte == "No") %>% 
  select(-method_outcome, -n.batches, -beta, -percentiles)
all1.500.no$meanV_V # 0.24 true V of dhat, estimated RR

## Poisson, no HTE
pois.500.no <- all.cvs3.wide_samplesize %>% filter(method == "Poisson", n.subset == 500, levhte == "No") %>% 
  select(-method_outcome, -n.batches, -beta, -percentiles)
pois.500.no$meanV_V # 0.26 true V of dhat, estimated RR
pois.500.no$VR_V # 0.67 VR of true V

## Poisson, high HTE
pois.500.high <- all.cvs3.wide_samplesize %>% filter(method == "Poisson", n.subset == 500, levhte == "High") %>% 
  select(-method_outcome, -n.batches, -beta, -percentiles)
pois.500.high$VR_V # 0.82

## All A1, high HTE
alla1.500.high <- all.cvs3.wide_samplesize %>% filter(method == "All A1", n.subset == 500, levhte == "High") %>% 
  select(-method_outcome, -n.batches, -beta, -percentiles)
alla1.500.high$VR_V # 0.58


#######################################################################
##########      Figure 3: Sample Size and HTE   #######################
#######################################################################

# Figure 3a. Sample size and HTE with Value Function as Y
fig3a <- plotV(data = all.cvs3.wide_samplesize, 
               thisV = "meanV_V", 
               display_methods = display_methods, 
               colorgroup = "levhte",
               output = F)

# Figure 3b. Sample size and HTE with Value Ratio as Y
fig3b <- plotV(data = all.cvs3.wide_samplesize, 
               thisV = "VR_V", 
               display_methods = display_methods, 
               colorgroup = "levhte",
               output = F)

# Figure 3
fig3 <- cowplot::plot_grid(fig3a, fig3b + theme(legend.position = "none"), nrow = 2)
fig3
ggsave(paste0(manuscript_figure_folder, "Figure3.tiff"), height = 7, width = 10, dpi = 400)


#######################################################################
##########     Figure 4: Proportions   ################################
#######################################################################

# Figure 4a. Proportions with Value Function as Y
fig4a <- plotV(data = all.cvs3.wide_prop, 
               thisV = "meanV_V",
               display_methods = display_methods, 
               colorgroup = "pergr",
               output = F)

# Figure 4b. Proportions with Value Ratio as Y
fig4b <- plotV(data = all.cvs3.wide_prop, 
               thisV = "VR_V",
               display_methods = display_methods, 
               colorgroup = "pergr",
               output = F)

# Figure 4
fig4 <- cowplot::plot_grid(fig4a, fig4b + theme(legend.position = "none"), nrow = 2)
fig4
ggsave(paste0(manuscript_figure_folder, "Figure4.tiff"), height = 8, width = 12, dpi = 400)


#######################################################################
##########      Figure 5: V.hat vs V     ##############################
#######################################################################

# Figure 5
levhte <- "no"
pergr <- "20x5"
load(paste0("./intermediate_results/simulations_samplesize/simsummary_sample_size_dataplot_", levhte, "_", pergr, ".RData"))

fig5data <- `all.cvs3_no_20x5` %>% 
  mutate(whichV = factor(whichV, levels = c("V", "V.hat"), 
                         labels = c("V(dhat)", "Vhat(dhat)")),
         levhte = "No HTE") %>% 
  filter(percentiles == "seq(0,1,by=0.2)", 
         method %in% display_methods
         ) 

# Separately plotted and overlaid 
fig5data %>%
  ggplot(aes(x = n.subset, y = meanV, color = whichV)) +
  geom_point(size = 3) + 
  geom_line(aes(group = whichV)) +
  geom_errorbar(aes(ymin = meanV - sdV, ymax = meanV + sdV),
                linewidth = 0.5, width = 0.2) +
  facet_grid(levhte ~ method, scales = "free") +
  scale_x_continuous(breaks = ns, trans = "log") +
  theme_bw() + 
  scale_color_manual(values = c(colR[8], colR[5])) +
  ylim(c(0.05, 0.45)) +
  labs(x = "Sample size (log-transformed)",
       y = "Value Function +- SE") +
       # title = "Estimated vs true cross-validated value function of ITR estimated from selected PM methods\n(with naive SE error bars)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 315, size = 12),
        strip.text.x = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_blank())  

ggsave(paste0(manuscript_figure_folder, "Figure5.tiff"), height = 5, width = 10, dpi = 400)


#######################################################################
########## SuppTable 1: Configuration of the simulation scenarios #####
#######################################################################

# Define HTE vector repeated
HTE.unique <- c("No", "Low", "Medium", "High")

# Define ARR Ratio (this is pre-specified for the simulation)
ARR_Ratio <- c(
  round(exp(rep(-0.2, 5)), 2), # No
  0.7, 0.75, 1, 1.05, 1.1,     # Low
  0.4, 0.5, 1, 1.1, 1.2,       # Medium
  0.3, 0.5, 1, 1.1, 1.5        # High
)

# Beta coefficients (log of ARR ratio)
# Each row of beta.matrix is No, Low, Medium, High HTE
beta.matrix <- matrix(round(log(ARR_Ratio), 2), ncol = 4) %>% t()

# Define proportions of HTE
Proportions <- c(
  0.1, 0.15, 0.50, 0.15, 0.1,  # Unequal symmetric 1: 10-15-50-15-10
  0.1, 0.3, 0.2, 0.3, 0.1,     # Unequal symmetric 2: 10-30-20-30-10
  0.55, 0.3, 0.15, 0, 0,       # Asymmetric 1: 55-30-15
  0.55, 0.15, 0.15, 0.15, 0    # Asymmetric 2: 55-15-15
)

# Each row of Proportions correspond to each element of Proportion.unique
Proportion.unique <- c("Unequal symmetric 1", "Unequal symmetric 2", "Asymmetric 1", "Asymmetric 2")

# Combine into tibble
supptab1.col <- tibble(
  HTE = factor(rep(HTE.unique, each = 5), levels = c("No", "Low", "Medium", "High")),
  Population = paste0("Stratum ", rep(1:5, times = length(unique(HTE)))),
  `Prop Patients` = rep(0.2, length(HTE)),
  `ARR Ratio` = ARR_Ratio
) %>%
  arrange(HTE, Population)

# Overall summary rows
overall <- calculateATE(percentiles = seq(0, 1, by = 0.2), beta.matrix = beta.matrix)
overall$ATE.var[1] <- 0
summary_rows <- tibble(
    HTE = factor(HTE.unique, levels = c("No", "Low", "Medium", "High")),
    Population = "Overall",
    `Prop Patients` = 1,
    `ARR Ratio` = sprintf("%.2f (%.2f)", 
                          overall$ATE, 
                          overall$ATE.var)
  )


# Bind summary rows on top of each group
supptab1.equal.symmetric <- supptab1.col %>%
  mutate(`ARR Ratio` = sprintf("%.2f", `ARR Ratio`)) %>%
  bind_rows(summary_rows) %>%
  dplyr::mutate(Population = factor(Population, levels = c("Overall", sort(unique(Population[Population != "Overall"])))),
                `Prop Patients` = sprintf("%d%%", round(`Prop Patients` * 100))) %>%
  arrange(HTE, Population)
         
# Combine into tibble
supptab1.row <- tibble(
  HTE = rep("Medium", 4 * 5),
  Population = paste0("Stratum ", rep(1:5, times = 4)),
  Scenario = rep(Proportion.unique, each = 5),
  `Prop Patients` = Proportions,
  `ARR Ratio` = rep(c(0.4, 0.5, 1, 1.1, 1.2), times = 4)
) %>%
  dplyr::mutate(Scenario = factor(Scenario, levels = Proportion.unique),
                `ARR Ratio` = ifelse(`Prop Patients` == 0, NA, `ARR Ratio`)) %>%
  arrange(HTE, Scenario, Population)

# Adjust ARR_Ratio and beta.matrix to calculate overall rows for different proportion scenarios
# For Asymmetric 1:   
# subgroup 1-2: 55% high A1, subgroup 3-4: 30% moderate A1, subgroup 5: 15% neutral
ARR_Ratio2 <- c(
  round(exp(rep(-0.2, 5)), 2), # No
  0.7, 0.7, 0.75, 0.75, 1,     # Low
  0.4, 0.4, 0.5, 0.5, 1,     # Medium
  0.3, 0.3, 0.5, 0.5, 1        # High
)
beta.matrix2 <- matrix(round(log(ARR_Ratio2), 2), ncol = 4) %>% t()
# For Asymmetric 2:
# subgroup 1: 55% high A1, subgroup 2: 15% moderate A1, subgroup 3-4: 15% neutral, subgroup 5: 15% moderate A0
ARR_Ratio3 <- c(
  round(exp(rep(-0.2, 5)), 2), # No
  0.7, 0.75, 1, 1, 1.05,     # Low
  0.4, 0.5, 1, 1, 1.1,     # Medium
  0.3, 0.5, 1, 1, 1.1        # High
)
beta.matrix3 <- matrix(round(log(ARR_Ratio3), 2), ncol = 4) %>% t()

# Overall summary rows
# Percentiles are cumulative proportions
# We only need the 3rd element of each "overall" because the 3rd row of beta.matrix is Medium HTE
overall1 <- calculateATE(percentiles = c(0, 0.1, 0.25, 0.75, 0.9, 1), beta.matrix = beta.matrix)    # Unequal symmetric 1: 10-15-50-15-10
overall2 <- calculateATE(percentiles = c(0, 0.1, 0.4, 0.6, 0.9, 1),   beta.matrix = beta.matrix)    # Unequal symmetric 2: 10-30-20-30-10
overall3 <- calculateATE(percentiles = c(0, 0.25, 0.55, 0.70, 0.85, 1), beta.matrix = beta.matrix2) # Asymmetric 1: 55-30-15    
overall4 <- calculateATE(percentiles = c(0, 0.55, 0.70, 0.80, 0.85, 1), beta.matrix = beta.matrix3) # Asymmetric 2: 55-15-15-15 

summary_rows2 <- tibble(
  HTE = rep("Medium", 4),
  Population = "Overall",
  Scenario = Proportion.unique,
  `Prop Patients` = 1,
  `ARR Ratio` = sprintf("%.2f (%.2f)", 
                        c(overall1$ATE[3], overall2$ATE[3], overall3$ATE[3], overall4$ATE[3]), # first element of each list is ARR ratio estimate
                        c(overall1$ATE.var[3], overall2$ATE.var[3], overall3$ATE.var[3], overall4$ATE.var[3]) # second element of each list is ARR ratio variance
  )  
) 

# Bind summary rows on top of each group
supptab1.medium.hte <- supptab1.row %>%
  mutate(`ARR Ratio` = sprintf("%.2f", `ARR Ratio`)) %>%
  bind_rows(summary_rows2) %>%
  dplyr::mutate(Population = factor(Population, levels = c("Overall", sort(unique(Population[Population != "Overall"])))),
                `Prop Patients` = sprintf("%d%%", round(`Prop Patients` * 100)),
                Scenario = factor(Scenario, levels = Proportion.unique)) %>%
  arrange(HTE, Scenario, Population)

# Create Supplementary Table 1
supptable1 <- data.frame(
  Scenario    = c("Scenario", "", "HTE"),
  Population  = c("", "", "Population"),
  
  Equal_sym   = c(
    "Equal symmetric",
    "sym-20-20-20-20-20",
    "Prop. patients | ARR ratio (σ²)"
  ),
  
  Unequal_1   = c(
    "Unequal symmetric 1",
    "sym-10-15-50-15-10",
    "Prop. patients | ARR ratio (σ²)"
  ),
  
  Unequal_2   = c(
    "Unequal symmetric 2",
    "sym-10-30-20-30-10",
    "Prop. patients | ARR ratio (σ²)"
  ),
  
  Asym_1      = c(
    "Asymmetric 1",
    "asym-55-30-15",
    "Prop. patients | ARR ratio (σ²)"
  ),
  
  Asym_2      = c(
    "Asymmetric 2",
    "asym-55-15-15-15",
    "Prop. patients | ARR ratio (σ²)"
  ),
  
  stringsAsFactors = FALSE
)

# Add the column of Equal Symmetric in Supplementary Table 1
add1 <- supptab1.equal.symmetric %>% 
  mutate(Equal_sym = paste(`Prop Patients`, `ARR Ratio`, sep = " | ")) %>% 
  rename(Scenario = HTE) %>%
  select(Scenario, Population, Equal_sym)
missing_cols <- setdiff(colnames(supptable1), colnames(add1))
add1[, missing_cols] <- NA  # add the missing columns
add1 <- add1[, colnames(supptable1)]
supptable1 <- rbind(supptable1, add1)

# Add the row of Medium HTE in Supplementary Table 1
add2 <- supptab1.medium.hte %>%
  mutate(merged = paste(`Prop Patients`, `ARR Ratio`, sep = " | ")) %>% 
  select(HTE, Population, Scenario, merged) %>%
  pivot_wider(names_from = Scenario, values_from = merged) %>%
  rename(Scenario = HTE,
         Unequal_1 = `Unequal symmetric 1`,
         Unequal_2 = `Unequal symmetric 2`,
         Asym_1 = `Asymmetric 1`,
         Asym_2 = `Asymmetric 2`)
supptable1 %<>%
  left_join(add2, by = "Population", suffix = c("", ".add2")) %>%
  mutate(
    Unequal_1 = if_else(Scenario == "Medium", Unequal_1.add2, Unequal_1), 
    Unequal_2 = if_else(Scenario == "Medium", Unequal_2.add2, Unequal_2),
    Asym_1 = if_else(Scenario == "Medium", Asym_1.add2, Asym_1),
    Asym_2 = if_else(Scenario == "Medium", Asym_2.add2, Asym_2)
  ) %>%
  select(-ends_with(".add2")) %>%
  filter(Scenario != "No" | (Scenario == "No" & Population == "Overall"))
View(supptable1)
write.csv(supptable1, "./results/SuppTable1.csv", row.names = FALSE)

#######################################################################
########## SuppFigure 1: Sample Size and HTE - accuracy   #############
#######################################################################

# mean accuracy HTE, line plot, sample size on x-axis, level of heterogeneity color
dat.accuracy.summary_samplesize %>% 
  filter(method %in% c(display_methods, "All A1")) %>%
  ggplot(aes(x = n.subset, y = mean.acc, group = factor(levhte), color = factor(levhte))) +
  geom_point(size = 3, shape = 16) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2, size = 1, color = mygrey) +
  geom_hline(yintercept = 0.5, linetype = 3, size = 1, color = mygrey) +
  scale_color_manual(values = colG[c(3,5,7,9)]) +
  scale_y_continuous(limits = c(0.45, 1)) +
  scale_x_continuous(breaks = ns, trans = "log") +
  facet_wrap(~ method, nrow = 1) + 
  theme_bw() + 
  labs(x = "Sample Size (log-transformed)", y = "Mean accuracy") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 315, size = 12),
        strip.text.x = element_text(size = 13), 
        legend.position = "bottom",
        legend.title = element_blank()) 

ggsave(paste0(manuscript_figure_folder, "SuppFigure1.tiff"), height = 7, width = 12)

#######################################################################
########## SuppFigure 2: Proportions - accuracy   #####################
#######################################################################

dat.accuracy.summary_prop %>% 
  filter(method %in% c(display_methods, "All A1")) %>%
  ggplot(aes(x = n.subset, y = mean.acc, group = pergr, color = pergr)) +
  geom_point(size = 3, shape = 16) + 
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2, size = 1, color = mygrey) +
    facet_wrap(~ method, nrow = 1) + 
  scale_color_manual(values = colB[c(3,4,5,7,9)]) +
  scale_y_continuous(limits = c(0.3, 1)) +
  scale_x_continuous(breaks = ns, trans = "log") + 
  theme_bw() + 
  theme_classic() + 
  labs(x = "Sample size (log-transformed)", y = "Mean accuracy") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 315, size = 12),
        strip.text.x = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "top",
        legend.title = element_blank()) 

ggsave(paste0(manuscript_figure_folder, "SuppFigure2.tiff"), height = 7, width = 12)

