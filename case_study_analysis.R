# ------------------------------------------------------------------
# Project: Precision Medicine MS
#
# Program name: case_study_analysis.R
#
# Purpose: Code to generate table and figure for the CONFIRM case study
# ------------------------------------------------------------------

##################################################
#### Set up wd, libraries, and functions ####
##################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tableone))

manuscript_figure_folder <- "./results/figures/"
# homepath <- "/home/pjiang/pmms/" # CHANGE THIS #
# setwd(homepath) 

#######################################################################
############################# Read in data ############################
#######################################################################

# Preprocessed the case study: CONFIRM trial data
# NOTE: This is pseudo-data generated to reproduce Table 1 exactly
# Values do not correspond to real trial participants
# Non-negative constraints are enforced for time variables,
# which may lead to minor numerical differences from published summaries.
ds <- readRDS("./intermediate_results/case_study/preprocess.RDS") 

# Estimated value functions
# High-level intermediate results can be shared
rawvhats.dhat <- read_csv("./intermediate_results/case_study/case_study_stratified10foldCV/main_CV_results_rawvhats.dhat_stratified10folds.csv")

# Define constants
covars <- c("age", "female", "weightbl", "white", "diagyrs", "prmsgr", "rlps1yr", "trelmos",
            "edssbl", "tm25zbl", "nhptzbl", "pasatzbl", "chrt2_5bl", "sf36pcsbl", "sf36mcsbl")

myblue <- rgb(37, 15, 186, maxColorValue = 255)
mygrey <- rgb(124, 135, 142, maxColorValue = 255)

method_vec <- c("All GA", "All DMF",
                "Poisson", "Negative\n Binomial",
                "Linear", "dWOLS",
                "Boosting", "Two\n Regressions", "Contrast\n Regression",
                "List DTR\n (2 nodes)")


#######################################################################
############################# Figure 6: CV results ############################
#######################################################################

# Cross-validated value function boxplot
figdata <- rawvhats.dhat %>%
  mutate(method = stringr::str_split(method_outcome, "_") %>% map_chr(., 1),
         outcome = stringr::str_split(method_outcome, "_") %>% map_chr(., 2)) %>%
  filter(!(method == "listDTR3")) %>%
  mutate(method = case_when(
    method == "contrastReg" ~ "Contrast\n Regression",
    method == "allGA" ~ "All GA",
    method == "twoReg" ~ "Two\n Regressions",
    method == "boosting" ~ "Boosting",
    method == "allDMF" ~ "All DMF",
    method == "listDTR2" ~ "List DTR\n (2 nodes)",
    method == "dWOLS" ~ "dWOLS",
    method == "weightedPoisson" ~ "Weighted\n Poisson",
    method == "poisson" ~ "Poisson",
    method == "weightedNegBin" ~ "Weighted\n NegBin",
    method == "weightedLinear" ~ "Weighted\n Linear",
    method == "linear" ~ "Linear",
    method == "negBin" ~ "Negative\n Binomial"
  ),
  method = factor(method,
                  levels = method_vec,
                  labels = method_vec),
  vhat = U/W
  )

fixeddata <- figdata %>% 
  filter(method %in% c("All GA", "All DMF", "Poisson", "Boosting", "Contrast\n Regression", "List DTR\n (2 nodes)")) %>%
  group_by(method) %>%
  summarise(medianVold = median(vhat))

# Figure 6
# Using median value of all GA and all DMF as ref lines
figdata %>%
  filter(method %in% c("Poisson", "Boosting", "Contrast\n Regression", "List DTR\n (2 nodes)")) %>%
  ggplot(aes(x = method, y = vhat)) +
  geom_boxplot(color = myblue) +
  geom_hline(yintercept = fixeddata$medianVold[which(fixeddata$method == "All DMF")], linetype = 2, size = 1, color = mygrey) +
  geom_hline(yintercept = fixeddata$medianVold[which(fixeddata$method == "All GA")], linetype = 3, size = 1, color = mygrey) +
  theme_bw() +
  labs(x = "Method", y = "Cross-validated value function") +
  theme(axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 0, size = 15),
        strip.text.x = element_text(size = 12))
ggsave(paste0(manuscript_figure_folder, "Figure6.tiff"), height = 4, width = 8, dpi = 400)


######################################################################
############################# Table 1: Baseline statistics ############################
#######################################################################

table1 <- tableone::CreateTableOne(vars = covars, strata = "trt", data = ds) %>% 
  print(smd = FALSE, quote = FALSE, noSpaces = TRUE, test = FALSE) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Variable")

write.csv(table1, "./results/tables/Table1.csv", row.names = FALSE)
