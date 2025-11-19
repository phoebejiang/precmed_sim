# precmed_sim


This repository contains R scripts and resources for performing statistical analyses and simulations related to the manuscript ***Impact of treatment effect heterogeneity on the estimation of individualized treatment rules for count outcomes***. It explores the performance of individualized treatment rules (ITRs) via various simulations scenarios such as sample size, mangitude of hetereogeneity of treatment effect (HTE), and the proportion or responder distribution of the HTE, as well as a real-world case study in multiple sclerosis. Informed by clinical insights, these learnings can hep us determine if HTE estimation is feasible, and if so, identify the most effective ITR.

## Folder Structure

### Root Directory
- **01-setup.R**: Script for setting up the environment and loading required libraries.
- **02-propensityscore.R**: Implements propensity score analysis.
- **03-dWOLS.R**: Contains code for dynamic Weighted Ordinary Least Squares (dWOLS) analysis.
- **04-regression-based.R**: Performs regression-based precision medicine statistical methods.
- **05-listdtr.R**: Implements list-based dynamic treatment regimes (DTR).
- **06-LuScore.R**: Contains code for the two novel roubly robust methods: two regressions and contrast regression.
- **eachCV.R**: Script for cross-validation procedures.
- **simmain.R**: Main script for running simulations.
- **simplot.R**: Generates plots for simulation results.
- **simsummary.R**: Summarizes simulation results.
- **simsummary_alln.R**: Summarizes simulation results across all sample sizes.
- **utility.R**: Contains utility functions used across the project.

### Outputs
- **outputs/**: Directory for storing output files, including:
  - `.Rapp.history`: RStudio history file.
  - `simsummary_sample_size_dataplot_no_symm 20x5.RData`: Example simulation output file for the scenario with no symmetry with varying sample sizes and mangitude of HTE.
  - `accuracy_prop1.png`: Accuracy plot for the comparison of the the estimated ITR and true ITR.
  - `agreement_n2500_prop1.png`: Agreement plot among all pairs of the estimated ITRs for sample size 2500.
  - `value_prop1.png`: Value function plot as the evaluation of estimated ITR in terms of patient outcome.
  - `v_vs_vhat_prop1.png`: Comparison plot of true value function (v) vs. estimated value function (v_hat).

### Documentation
- **README.md**: This file, providing an overview of the project.
- **references.bib**: Bibliography file for citations.
- **runme.qmd**: Quarto markdown file for generating reports.
- **runme.html**: HTML report generated from `runme.qmd`.

## Getting Started

### Prerequisites
- R (version 4.2.2 or higher)
- Required R packages: tidyverse, magrittr, haven, Hmisc, MASS, pscl, caret, glmnet, mpath, gbm, fastDummies, listdtr, DTRreg, stringr, ggrepel, corrplot, RColorBrewer, ggnewscale

### Running the Project
1. Set up the environment by running `requirement.yml`.
2. Execute the analysis scripts in the desired order:
   - Run `simmain.R` to run simulations.
   - Run `simsummary.R` or `simsummary_alln.R` to summarize results.
   - Run `simplot.R` to visualize results.

### Outputs
Simulation results and plots will be saved in the `outputs/` directory.

## References
Citations and references are included in the `references.bib` file.

## License
[Insert license information here.]

## Acknowledgments
[Include acknowledgments or credits here.]
