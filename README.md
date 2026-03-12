---
editor_options: 
  markdown: 
    wrap: 72
---

# Impact of treatment effect heterogeneity on the estimation of individualized treatment rules for count outcomes

This README file contains the instructions required to reproduce the tables 
and figures presented in the manuscript, [Impact of Treatment Effect Heterogeneity on the Estimation of Individualized Treatment Rules for Count Outcomes](https://onlinelibrary.wiley.com/doi/10.1002/bimj.70119)

The R scripts implement the simulation studies, case study analyses, and 
supporting computations described in the paper, enabling full reproducibility 
of the reported results.

## Contact Me

If you have any questions, please feel free to reach out to the corresponding author Xiaotong Jiang
([xiaotong.phoebe.jiang\@gmail.com](mailto:xiaotong.phoebe.jiang@gmail.com){.email}).

## Version Information

This project was developed using the following software and packages:

-   Operating System: Linux (simulations) and Windows 11 (other)
-   R Version: 4.1.0
-   Loaded R Packages: Below is the output of `sessionInfo()` after
    loading all required packages:

``` r
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] devtools_2.4.6     usethis_3.2.1      DTRreg_2.3         gbm_2.2.2          mpath_0.4-2.26     pscl_1.5.9        
 [7] MASS_7.3-65        glmnet_4.1-8       Matrix_1.7-3       codetools_0.2-20   Hmisc_5.2-4        fastDummies_1.7.5 
[13] reshape2_1.4.4     truncnorm_1.0-9    cowplot_1.2.0      haven_2.5.4        ggpubr_0.6.2       ggnewscale_0.5.2  
[19] RColorBrewer_1.1-3 caret_7.0-1        lattice_0.22-6     corrplot_0.95      ggrepel_0.9.6      magrittr_2.0.3    
[25] lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.4        readr_2.1.5       
[31] tidyr_1.3.1        tibble_3.2.1       ggplot2_4.0.1      tidyverse_2.0.0    listdtr_1.1
```

## Reproducing Figures and Tables

To reproduce the figures and tables in the manuscript and its
supplements, follow the steps below:

### Step 1: Run simulation scripts

Folder `./simulations/samplesize`:

Run the scripts on a high-performance computing (HPC) cluster in the
following order:

1.  `simmain.sh`
2.  `simsummary.sh`
3.  `simsummary_sample_size.sh`

Folder `./simulations/proportion`:

Run the scripts on a HPC cluster in the following order:

1.  `simmain.sh`

2.  `simsummary.sh`

3.  `simsummary_proportion.sh`

At the end of this step, intermediate summary results will be generated
(likely on the HPC cluster). Copy or move these files to the
`./intermediate_results` folder.

Note: To save time, pre-computed intermediate results are already
provided in the `./intermediate_results/simulations_samplesize` folder
and the `./intermediate_results/simulations_proportions` folder. You can
use files in these two folders to skip the time-consuming simulations in
Step 1 and proceed directly to running `simulation_analysis.R`.

### Step 2: Generate simulation results

Run `simulation_analysis.R` locally. The script includes comments
indicating which code sections produce each figure and table related to
the simulation study (Figures 1-5, Supplementary Table 1, Supplementary
Figures 1 and 2).

### Step 3: Generate the CONFIRM case study results

Folder `./case_study`:

Run the scripts on a HPC in the following order:

1.  `main.sh`
2.  `summary.sl`

At the end of this step, intermediate summary results will be generated
(likely on the HPC cluster). Copy or move these files to the
`./intermediate_results/case_study` folder. Run `case_study_analysis.R`
locally to generate Figure 6 and Table 1 related to the case study.

Note: To save time, pre-computed intermediate results are already
provided in the `./intermediate_results/case_study` folder, which are
`preprocess.RDS` (preprocessed case study data) and
./`case_study_stratified10foldCV/main_CV_results_rawhats.dhat_stratified10foldCVs.csv`
(intermediate precision medicine results of the preprocessed data). You
can use files in this folder to skip the calculations and proceed
directly to running `case_study_analysis.R`.

## Manual Alterations

No manual alterations to the code are required. All scripts are
self-contained and can be executed as-is. If any manual edits are
necessary, they will be explicitly documented here with file names, line
numbers, and the exact content of the edits.

Note that the scripts in the `./simulations/` folder should be run on a
HPC cluster, as computations for large sample sizes and certain
time-intensive methods, such as listDTR, two regressions, and contrast
regression, which can be very time-consuming.

## Reproducible Research

The simulation scripts in `./simulations/samplesize` and
`./simulations/proportion` are designed to be run on a HPC cluster.
However, to improve accessibility and facilitate research
reproducibility, we provide an example of a minimal simulation
configuration that can be executed locally within a reasonable amount of
time.

### How to run `simmain.R` locally

The example below demonstrates how to run `simmain.R` from the command
line with a specific set of arguments:

``` bash
time Rscript simmain.R poisson logarr0001 1 500 "c(-0.2, -0.2, -0.2, -0.2, -0.2)" "seq(0, 1, by = 0.2)"
```

This configuration corresponds to one simulation scenario where the
model is Poisson, the batch index is 1, the sample size is 500, no HTE
heterogeneity (`beta = c(-0.2, -0.2, -0.2, -0.2, -0.2)`), and equal
symmetric HTE responder distribution (`perc = seq(0, 1, by = 0.2)`).

This configuration was selected because it is relatively simple and can
be run locally without requiring an HPC environment. Different values
may be supplied for each of the six arguments to explore alternative
simulation settings but we recommend that users choose smaller faster
models such as Poisson instead of contrast regression or listDTR.

Here are other values of the arguments that you can try:

-   PM method: `allA1`, `allA0`, `linear`, `poisson`, `dWOLS`,
    `boosting`, `twoReg`, `contrastReg`, `listDTR2` (note that the last
    three methods are much more slower)

-   Outcome variable: `logarr0001` (only one value)

-   Batch index:

    -   1 to 5 (the first of five batches; each batch contains five
        cross-validation iterations) if sample size is 500 or 1000

    -   1 to 25 if sample size is 2500, 5000, or 10000

-   Sample size: `500`, `1000`, `2500`, `5000`, `10000`

-   Overall HTE magnitude `beta`:

    -   `"c(-0.2,-0.2,-0.2,-0.2,-0.2)"` (which corresponds to No HTE)

    -   `"c(-0.36,-0.29,0,0.05,0.1)"` (which corresponds to Low HTE)

    -   `"c(-0.92,-0.69,0,0.1,0.18)"` (which corresponds to Medium HTE)

    -   `"c(-1.2,-0.69,0,0.1,0.41)"` (which corresponds to High HTE)

-   HTE distribution `perc`:

    -   `"seq(0,1,by=0.2)"` (which corresponds to equal symmetric
        heterogeneity with 20% of patients in each of the five responder
        groups)

    -   `"c(0,0.1,0.25,0.75,0.9,1)"` (which corresponds to the symmetric
        heterogeneity with 10%-15%-50%-15%-10% of patients in each of
        the five responder groups)

    -   `"c(0,0.1,0.4,0.6,0.9,1)"` (which corresponds to the symmetric
        heterogeneity with 10%-30%-20%-30%-10% of patients in each of
        the five responder groups)

    -   `"c(0,0.55,0.65,0.75,0.85,1)"` (which corresponds to the
        asymmetric heterogeneity with 55%-15%-15%-15%-0% of patients in
        each of the five responder groups)

    -   `"c(0,0.55,0.65,0.7,0.85,1)"` (which corresponds to the
        asymmetric heterogeneity with 55%-30%-15%-0%-0% of patients in
        each of the five responder groups)

More detailed information about the arguments for `simmain.R` can be
found in `simmain.sh`.

### How to run `simsummary.R` locally

After `simmain.R` is run 5 times for all 5 batches (e.g., 1 to 5 for n =
500) and all PM models (e.g., Poisson), the example below demonstrates
how to run `simsummary.R` from the command line with a specific set of
arguments:

``` bash
time Rscript simsummary.R 500 "c(-0.2, -0.2, -0.2, -0.2, -0.2)" "seq(0, 1, by = 0.2)"
```

This configuration corresponds to one simulation scenario where the
sample size is 500, no HTE magnitude
(`beta = c(-0.2, -0.2, -0.2, -0.2, -0.2)`), and equal symmetric HTE
responder distribution (`perc = seq(0, 1, by = 0.2)`) across all batches
of the PM models that you run.

If you run `simmain.R` with other values of arguments, then you will
need to keep the same arguments of sample size, `beta`, and `perc`.

More detailed information about the arguments for `simsummary.R` can be
found in `simsummary.sh`.

### How to run `simsummary_samplesize.R` and `simsummary_proportion` locally

The example below demonstrates how to run `simsummary_samplesize.R` or
`simsummary_proportion.R`, which summarizes `simsummary.R` results
across all HTE magnitudes and distributions, from the command line with
a specific set of arguments:

``` bash
time Rscript simsummary_samplesize.R 500 
```

or

``` bash
time Rscript simsummary_proportion.R 500 
```

This configuration corresponds to one simulation scenario where the
sample size is 500 across all batches, PM models, HTE magnitudes, and
HTE distributions that you run.

Here are other values of the arguments that you can try if you run
`simmain.R` and `simsummary.R` with other sample sizes:

-   Sample size: `500`, `1000`, `2500`, `5000`, `10000`

Note:

-   `simsummary_samplesize.R` can be run after `simmain.R` and
    `simsummary.R` in the folder `./simulations/samplesize/` are run for
    all PM methods, batches, and HTE magnitudes and distributions.

-   `simsummary_proportion.R` can be run after `simmain.R` and
    `simsummary.R` in the folder `./simulations/proportion/` are run for
    all PM methods, batches, and HTE magnitudes and distributions.

More detailed information about the arguments for
`simsummary_samplesize.R` can be found in `simsummary_samplesize.sh` and
`simsummary_proportion.sh` for `simsummary_proportion.R`.

## File and Folder Structure

Below is a listing of the files and folders in the project, with brief
explanations of their content:

-   `01-preprocessing.R`: Preprocesses the input data.

-   `02-propensityscore.R`: Calculates propensity scores.

-   `03-dWOLS.R`: Implements the dWOLS method.

-   `04-regression-based.R`: Implements the regression-based methods:
    linear, Poisson, negative binomial, and zero-inflated negative
    binomial models.

-   `05-listdtr.R`: Implements the listDTR method.

-   `06-LuScore.R`: Implements the four methods: Poisson, boosting, two
    regressions, and contrast regression.

-   Folder `case_study`: Contains the CONFIRM case study scripts.

    -   `main.R`: Implements the precision medicine methods to the case
        study data, CONFIRM, using cross-validation.

    -   `main.sh`: A shell script to run the `main.R` for different
        methods and different batches.

    -   `main_sl.temp`: A template slurm script that the shell script
        above calls.

    -   `summary.R`: Summarizes results from `main.R` to intermediate
        results.

    -   `summary.sl`: A slurm script that runs the `summary.R` on
        cluster.

-   `case_study_analysis.R`: Final script that generates figure and
    table that are related to the CONFIRM case study (Figure 6 and Table
    1 presented in the manuscript).

-   `eachCV.R`: Generates cross-validation results using the precision
    medicine methods in `03-dWOLS.R` to `06-LuScore.R`.
    
-   `generate_pseudo_case_study_data.R`: Generate a privacy-preserving pseudo 
     dataset that mirrors the original case study data in size, structure, 
     and key summary characteristics, and outputs the `processed.RDS` file 
     for downstream analyses.

-   Folder `intermediate_results`: Contains intermediate results.

    -   `case_study`: Intermediate results related to the CONFIRM case
        study.

        -   `case_study_stratified10foldCV`: Intermediate summarized
            results of the case study after 10-fold stratified CV.

        -   `preprocess.RDS`: Preprocessed case study data derived from
            the CONFIRM trial. To safeguard participant privacy, this
            file contains a synthetic, anonymized dataset designed to
            mirror the original trial data in size, structure, and key
            summary statistics.

    -   `simulations_proportions`: Intermediate simulation results
        related to the different distributions of HTE and sample sizes
        for a fixed magnitude of HTE.

    -   `simulations_samplesize`: Intermediate simulation results
        related tot he different magnitudes of HTE and sample sizes for
        a fixed distribution of HTE.

-   Folder `results`: Contains all results reported in the manuscript.

    -   `figures`: All figures in the manuscript, from `Figure1.tiff` to
        `SuppFigure2.tiff`.

    -   `tables`: `SuppTable1.csv` and `Table 1.csv`.

-   `simulation_analysis_functions.R`: Functions used in the
    `simulation_analysis.R`.

-   `simulation_analysis.R`: Final script that generates figures and
    tables of the simulation analysis presented in the manuscript.

-   Folder `simulations`: Folder containing simulation scripts.

    -   `samplesize`:

        -   Contains simulation code for 4 different magnitudes of HTE
            (No, Low, Medium, High) across 5 sample sizes (n = 500,
            1000, 2500, 5000, 10000) for a given responder proportion
            (symmetric equal 20%-20%-20%-20%-20%).

        -   See Supplementary Table 1 for more details.

    -   `proportion`:

        -   Contains simulation code for 5 different proportions of HTE
            (symmetric 20%-20%-20%-20%-20%, symmetric
            10%-15%-50%-15%-10%, symmetric 10%-30%-20%-30%-10%,
            asymmetric 55%-30%-15%, asymmetric 55%-155-15%-15%) across 5
            sample sizes (n = 500, 1000, 2500, 5000, 10000) for a given
            HTE magnitude (Medium).

        -   See Supplementary Table 1 for more details.

-   `utility.R`: Utility functions used across scripts.

## Data Documentation

### Input Data Sets

-   Input simulation data are generated from the `simdata()`
    function in `utility.R`.

    -   They are used in the two `simmain.R` scripts in
        `./simulations/samplesize` and `./simulations/proportion`.

-   Case study input data `preprocess.RDS` in `./intermediate_results/case_study`
    
    - To safeguard participant privacy, this file contains a synthetic, 
      anonymized dataset designed to  mirror the original preprocssed CONFIRM trial 
      data in size, structure, and key summary statistics.

### Intermediate Results

-   Folder: `intermediate_results`

    -   Contains simulation results for different magnitudes of HTE and
        proportions of HTE; both of which vary across different sample
        sizes. Use datasets in this folder and run
        `simulation_analysis.R` and `case_study_analysis.R` directly to
        reproduce simulation figures and tables in the manuscript.
        
        
    -   Contains case study analysis results. Use dataset in folder 
        `./case_study_stratified10foldCV` and run `case_study_analysis.R` 
        directly to reproduce case study figures and tables in the manuscripts.

