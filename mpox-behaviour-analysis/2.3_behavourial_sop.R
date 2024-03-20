# author: Fanyu Xiu
# library and data --------
library(dplyr)
library(readr)
library(tibble)
library(rstanarm)
library(bayestestR)
source("./src/helper_dates.R")
source("./src/helper_effectestim.R")

select <- dplyr::select
    
analysis <- "main" # or "sen_date7" or "sen_date6"
sa_grouping <- 7 # or 3 or 5
# outcome_var <- "groupsex_m" 
outcome_var <- "bath_m"

## define paths & prefixes based on the analysis being done
fig_path <- case_when(outcome_var == "groupsex_m" ~ sprintf("./behavourial-figures/%s/results-checks-groupsex", analysis),
                      outcome_var == "bath_m" ~ sprintf("./behavourial-figures/%s/results-checks-bathhouse", analysis))

out_distr_path <- case_when(outcome_var == "groupsex_m" ~ sprintf("./behavourial-outputs/%s/fitted-distr-groupsex", analysis),
                            outcome_var == "bath_m" ~ sprintf("./behavourial-outputs/%s/fitted-distr-bathhouse", analysis))
out_distr_pref <- case_when(outcome_var == "groupsex_m" ~ sprintf("groupsex-%s", sa_grouping),
                            outcome_var == "bath_m" ~ sprintf("bathhouse-%s", sa_grouping))

# load data and create covariates -----
ifelse(analysis == "main",
       data_22 <- read.csv("./data-3cities-feb-2023/behavourial_regression_3cities.csv"),
       ifelse(analysis == "sen_date6",
              data_22 <- read.csv("./data-3cities-feb-2023/behavourial_regression_3cities_sen_date6.csv"),
              data_22 <- read.csv("./data-3cities-feb-2023/behavourial_regression_3cities_sen_date7.csv")))
data_22 <- data_22 %>% 
  mutate(month_intv_centre = scale(month(date_intv_visit), scale = F),
         period_mpox_centre = scale(period_mpox, scale = F), ### centre the continuous variable to avoid colinearity issue
         mpox_visit = ifelse(period_mpox == 0, 0, 1))

### number of visits during year 2022
nrow(data_22) # 1957

### number of participants with visits during year 2022
length(unique(data_22$part_id)) # 1445

### summary of mpox outbreak coverage
summary(data_22$period_mpox)

### number of visits during mpox outbreak
nrow((filter(data_22, period_mpox != 0))) 

### number of participants with visits during mpox outbreak
length(unique(filter(data_22, period_mpox != 0)$part_id)) 

# save time by aggregating individuals
# before computing Pr(Y = y) for 0 to 300
nrow(data_22) # nb of individuals, could be the same person at different timepoints

## create dummy variables ----

data_22 <- data_22 %>%
  mutate(nb_part_anal = nb_part_anal_visit,
         nb_part_ttl = nb_part_ttl_visit,
         groupsex_m = groupsex_m_visit,
         groupsex_d = groupsex_d_visit,
         bath_m = bath_m_visit,
         bath_d = bath_d_visit,
         age_cats = age_cats_visit,
         hiv_cats = hiv_cats_visit,
         rel_status = rel_status_last) 


# check SPVs and groupsex
table(data_22$bath_d, useNA = "ifany")
table(data_22$groupsex_d, useNA = "ifany")

data_22 <- data_22 %>% 
  filter(groupsex_d != 1 & bath_d != 1)

# age; reference is 16-29
data_22 <- make_ind_age(data_22)

# partnership status
table(data_22$reg_partn_last, data_22$rel_status, useNA = "ifany")
data_22 <- make_ind_rel(data_22)

# hiv status 
table(data_22$hiv_cats, useNA = "ifany")

# activity groups
higher_grp <- paste0(">", sa_grouping)
lower_grp <- paste0("0-", sa_grouping)
data_22 <- data_22 %>% 
  mutate(sa_cats_ttl = ifelse(nb_part_ttl_last > sa_grouping, higher_grp, lower_grp),
         sa_cats_anal = ifelse(nb_part_anal_last > sa_grouping, higher_grp, lower_grp)) %>% 
  mutate(sa_cats_ttl = factor(sa_cats_ttl, levels = c(lower_grp, higher_grp)),
         sa_cats_anal = factor(sa_cats_anal, levels = c(lower_grp, higher_grp)))
table(data_22$sa_cats_ttl)
data_22 <- data_22 %>%
  mutate(sa_cats_ttl_higher = ifelse(sa_cats_ttl == higher_grp, 1, 0),
         period_mpox_centre_sa_cats_ttl_higher = period_mpox_centre * sa_cats_ttl_higher)
table(data_22$`sa_cats_ttl_higher`)

# ID
data_22$id <- match(data_22$part_id, unique(data_22$part_id))

# Fit Bayesian models ----

# variables to use
vars_model <- c("period_mpox_centre",
                "month_intv_centre",
                "age_30_39", "age_40_49", "age_50_59", "age_60_",
                "hiv_cats",
                "rel_y_excl", "rel_y_open", "rel_y_uncl")

formula_ttl <- as.formula(paste(paste(outcome_var,"~ (1|id) +"), paste(vars_model, collapse = "+")))

num_cores <- parallel::detectCores()
t0 <- Sys.time()
set.seed(77)
fit_ttl <- stan_glmer(
  formula_ttl,
  family = "binomial",
  iter = 4000,
  chains = 2,
  data = data_22,
  cores = num_cores,
  prior_intercept = normal(0, 10),
  prior = normal(0, 10), # note second argument is the sd
  prior_aux = cauchy(0, 5)
)
t1 <- Sys.time()
t1 - t0 
effective_sample(fit_ttl) 
prior_summary(fit_ttl)
saveRDS(fit_ttl, file = sprintf("%s/fit_bayes_ls_%s.rds", 
                                     out_distr_path,
                                     out_distr_pref))

## Inspect model convergence diagnostic ----
trace <- plot(fit_ttl, "trace", pars = c("(Intercept)", vars_model))
p_trace_plot <- trace + ggplot2::scale_color_discrete()
ggsave(sprintf("%s/model-checks-p6m-%s.png", fig_path, out_distr_pref),
       p_trace_plot, device = "png",
       height = 14, width = 30, units = "cm", dpi = 320)

# get model
cur_model <- tbl_coeff_fn(fit_ttl, alpha = 0.05, stat = "bayesian")
  
vars_model_full <- c("Mpox outbreak coverage",
                     "Month the visit took place",
                     "Age 30-39", "Age 40-49", "Age 50-59", "Age ≥60", 
                     "HIV Status",
                     "Exclusive Relationship", "Open Relationship", "Unclear Relationship")

cur_model <- add_column(cur_model,
                        name_full = c("exp(intercept)", vars_model_full),
                        .after = "name")

cur_model <- cur_model %>% 
  mutate(coeff = factor(name_full, levels = unique(name_full))) 

coeff_tbl <- cur_model %>% 
  dplyr::select(coeff, mean, se, conf_lb, conf_ub, n_eff, Rhat) %>%
  mutate(Mean = paste0(round(mean, 2), " (", round(conf_lb, 2), ", ", round(conf_ub, 2), ")"),
         SE = signif(se, 2)) %>%
  dplyr::select(coeff, Mean, SE, n_eff, Rhat)

# save coefficients 
write.csv(coeff_tbl, sprintf("%s/table_coef_%s.csv", out_distr_path, out_distr_pref))
