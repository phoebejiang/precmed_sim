library(tableone)
library(tidyverse)

# Define constants
covars <- c("age", "female", "weightbl", "white", "diagyrs", "prmsgr", "rlps1yr", "trelmos",
            "edssbl", "tm25zbl", "nhptzbl", "pasatzbl", "chrt2_5bl", "sf36pcsbl", "sf36mcsbl")

# Define functions
make_continuous <- function(n, mean, sd) {
  x <- seq_len(n) - mean(seq_len(n))   # centered
  x <- x / sd(x)                       # unit SD
  mean + sd * x                        # scale + shift
}

make_negbin <- function(n, mean, sd) {
  size <- mean^2 / (sd^2 - mean)
  rnbinom(n, size = size, mu = mean)
}

## Sample sizes
n_dmf <- 359
n_ga  <- 350

## Create the data frame with treatment
ds <- data.frame(
  trt = factor(c(rep("DMF", n_dmf), rep("GA", n_ga)))
)

## Add demographics

ds$age <- c(
  make_continuous(n_dmf, 37.8, 9.4),
  make_continuous(n_ga,  36.7, 9.1)
)

ds$female <- c(
  rep(1, 245), rep(0, n_dmf - 245),
  rep(1, 247), rep(0, n_ga  - 247)
)

ds$weightbl <- c(
  make_continuous(n_dmf, 71.9, 17.9),
  make_continuous(n_ga,  71.4, 19.1)
)

ds$white <- c(
  rep(1, 304), rep(0, n_dmf - 304),
  rep(1, 290), rep(0, n_ga  - 290)
)

## Add disease severity -

set.seed(123)
ds$diagyrs <- c(
  make_negbin(n_dmf, mean = 4.9, sd = 5.1),
  make_negbin(n_ga,  mean = 4.4, sd = 4.7)
)

ds$prmsgr <- c(
  rep(1, 101), rep(0, n_dmf - 101),
  rep(1, 103), rep(0, n_ga  - 103)
)

set.seed(123)
rlps1yr_dmf <- sample(
  x = c(0, 1, 2),
  size = n_dmf,
  replace = TRUE,
  prob = c(0.05, 0.60, 0.35)
)

rlps1yr_ga <- sample(
  x = c(0, 1, 2),
  size = n_ga,
  replace = TRUE,
  prob = c(0.02, 0.56, 0.42)
)

ds$rlps1yr <- c(rlps1yr_dmf, rlps1yr_ga)

set.seed(123)
trelmos_dmf <- rnbinom(
  n_dmf,
  mu   = 6.1,
  size = 6.1^2 / (4.2^2 - 6.1)
)

trelmos_ga <- rnbinom(
  n_ga,
  mu   = 6.3,
  size = 6.3^2 / (6.3^2 - 6.3)
)

ds$trelmos <- c(trelmos_dmf, trelmos_ga)

ds$edssbl <- c(
  make_continuous(n_dmf, 2.6, 1.2),
  make_continuous(n_ga,  2.6, 1.2)
)

## Add z-scores

ds$tm25zbl <- c(
  make_continuous(n_dmf,  0.01, 0.94),
  make_continuous(n_ga,  -0.05, 1.20)
)

ds$nhptzbl <- c(
  make_continuous(n_dmf,  0.00, 1.03),
  make_continuous(n_ga,  -0.01, 0.99)
)

ds$pasatzbl <- c(
  make_continuous(n_dmf,  0.02, 0.99),
  make_continuous(n_ga,   0.00, 1.01)
)

ds$chrt2_5bl <- c(
  make_continuous(n_dmf, 31.89, 12.44),
  make_continuous(n_ga,  31.94, 12.61)
)

ds$sf36pcsbl <- c(
  make_continuous(n_dmf, 43.07,  9.91),
  make_continuous(n_ga,  43.19, 10.16)
)

ds$sf36mcsbl <- c(
  make_continuous(n_dmf, 45.42, 11.37),
  make_continuous(n_ga,  44.77, 10.65)
)

# Formatting
ds <- ds %>% 
  mutate(female = factor(female, levels = c(0, 1)),
         white = factor(white, levels = c(0, 1)),
         prmsgr = factor(prmsgr, levels = c(0, 1)))

# Create table 1
table1 <- tableone::CreateTableOne(vars = covars, strata = "trt", data = ds) %>% 
  print(smd = FALSE, quote = FALSE, noSpaces = TRUE, test = FALSE) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Variable")
table1

# NOTE: This is pseudo-data generated to reproduce Table 1 exactly.
# Values do not correspond to real trial participants.
# Non-negative constraints are enforced for time variables,
# which may lead to minor numerical differences from published summaries.

saveRDS(ds, "./intermediate_results/case_study/preprocess.RDS")
