# mpox-behaviour-analysis
Code to analyze change in sexual partner numbers among gay, bisexual, and other men who have sex with men during the 2022-2023 mpox outbreak in Canada.

## code structure
Briefly, the first number in the file prefix means the following
- `1`: preliminary,
- `2`: analyses of change in sexual partner numbers, and
- `3`: parameterization of a dynamic model of mpox transmission.

The second number indicates the sequence in code execution.

### preliminary
- [`1.1_make_data_regression.R`](1.1_make_data_regression.R) data cleaning and wrangling to prepare for the regression analysis of change in sexual partner numbers.
- [`1.2_make_data_contact_rate.R`](1.2_make_data_contact_rate.R) implementing respondent-driven sampling (RDS) and inverse probability of censoring (IPC) weights to prepare for the paramterization of contact rate for the dynamic model.

### analyses of change in sexual partner numbers
- [`2.1_descr_table1produces.R`](2.1_descr_table1produces.R) produces [__Table S5__] (Unadjusted and RDS-II adjusted distribution of variables among Engage Cohort Study visits in 2022).
- [`2.2_behavourial_nb.R`](2.2_behavourial_nb.R) runs the random-intercept negative binomial regression model for main and sensitivity analyses with alternative dates to define the mpox high transmission period and sexual activity level groupings. Convergence check is in [`behavourial-figures`](behavourial-figures) and posterior coefficients are in [`behavourial-outputs`](behavourial-outputs).
- [`2.3_behavourial_sop.R`](2.2_behavourial_sop.R) runs the random-intercept logistic regression model for sensitivity analyses with visits to sex-on-premises venues as the outcome. Convergence check is in [`behavourial-figures`](behavourial-figures) and posterior coefficients are in [`behavourial-outputs`](behavourial-outputs).

### parameterization of a dynamic model of mpox transmission
- [`3.1_fitted_distri_contact_rate.R`](3.1_fitted_distri_contact_rate.R) leverages a previously described the negative binomial regression model (https://github.com/pop-health-mod/mpox-engage-sex-networks) to perform post-stratification with RDS-II and IPCW, and obtain fitted distribution of the number of sexual partners. Convergence check is in [`parametrization-figures`](parametrization-figures).
- [`3.2_sa_cats_partition_contact_rate.R`](3.2_sa_cats_partition_contact_rate.R) partitions the fitted distribution of the number of sexual partners into 15 sexual activity groups and formats the age and HIV-serostatus mixing matrices based on the data and Milwid et al (2022) (https://link.springer.com/article/10.1186/s12879-022-07207-7).
- [`3.3_ve_meta_analysis.Rmd`](3.3_ve_meta_analysis.Rmd) performs meta analysis on the effectiveness of a single dose of mpox vaccination. Results are in [`meta-analysis-figures`](meta-analysis-figures).
