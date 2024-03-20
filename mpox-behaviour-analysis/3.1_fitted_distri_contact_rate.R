# # authors: Jorge Luis Flores Anato, Fanyu Xiu
# library and data --------
library(dplyr)
library(readr)
library(rstan)
library(ggplot2)
source("./src/helper_dates.R")
data_full <- read_csv("./data-3cities-feb-2023/engage_visits_3cities.csv")
select <- dplyr::select

## main analyses
outcome_var <- "nb_part_ttl"

## define paths & prefixes based on the analysis being done
fig_path <- case_when(outcome_var == "nb_part_ttl" ~ "./parametrization-figures/results-checks-ttl",
                      outcome_var == "nb_part_anal" ~ "./parametrization-figures/results-checks-anal"
)

out_distr_path <- case_when(outcome_var == "nb_part_ttl" ~ "./parametrization-outputs/fitted-distr-ttl",
                            outcome_var == "nb_part_anal" ~ "./parametrization-outputs/fitted-distr-anal")
out_distr_pref <- case_when(outcome_var == "nb_part_ttl" ~ "all",
                            outcome_var == "nb_part_anal" ~ "anal")

## load data
data_3cities_pre_covid <- read_csv("./data-3cities-feb-2023/survey_baseline_3cities.csv")
data_3cities_pre_mpox <- read_csv("./data-3cities-feb-2023/data_contact_rate_3cities.csv")

# create single dataset with all time periods & cities
data_3cities <- bind_rows(
  data_3cities_pre_covid,
  data_3cities_pre_mpox) %>%
  mutate(time_pt = factor(time_pt, 
                          levels = c("Pre-COVID", "2022-Pre-Mpox"))) %>%
  mutate(city = recode_factor(city, mtl = "Montreal", trt = "Toronto", van = "Vancouver"))

# create city marker
CITIES <- c("Montreal", "Toronto", "Vancouver")
TIMEPTS <- c("Pre-COVID", "2022-Pre-Mpox")

AGES <- sort(unique(data_3cities$age_cats))
HIV <- unique(data_3cities$hiv_cats)
AGEHIV <- paste(rep(AGES, each = 2), HIV, sep = ".")

table(data_3cities$time_pt, 
        data_3cities$city, 
        data_3cities$age_cats,
        data_3cities$hiv_cats,
        useNA = "ifany")

data_3cities <- data_3cities %>% 
    mutate(
      data_pt = factor(paste(city, time_pt, sep = "-"),
                       levels = paste(rep(CITIES, each = length(TIMEPTS)), TIMEPTS, sep = "-")),
      age_hiv = factor(paste(age_cats, hiv_cats, sep = "."),
                       levels = AGEHIV)
    )
  
CITIES_DATAPTS <- paste(
  rep(CITIES, each = length(TIMEPTS)), 
  rep(TIMEPTS, times = length(CITIES)), 
  sep = "-"
)

# save time by aggregating individuals
# before computing Pr(Y = y) for 0 to 300
nrow(data_3cities) # nb of individuals, could be the same person at different timepoints
data_aggrt <- data_3cities %>% 
    group_by(age_hiv) %>% 
    count(data_pt, rel_status,
          bath_m, bath_d,
          groupsex_m, groupsex_d, 
          apps_partn_m, apps_partn_d, 
          sex_work_m, sex_work_d)
count(data_aggrt, "age_hiv") # nb of combinations of covariates for each age-hiv group (across all city-datapoint)

# Fit Bayesian models ----

## create dummy variables
# age; reference is 16-29
data_3cities <- make_ind_age(data_3cities)

# partnership status
table(data_3cities$data_pt, data_3cities$rel_status, useNA = "ifany")
table(data_3cities$reg_partn, data_3cities$rel_status, useNA = "ifany")
data_3cities <- make_ind_rel(data_3cities)

# check SPVs and groupsex
table(data_3cities$data_pt, data_3cities$bath_d, useNA = "ifany")
table(data_3cities$data_pt, data_3cities$groupsex_d, useNA = "ifany")

# check apps variable
# NOTE: apps variable not present in FU visit
table(data_3cities$data_pt, data_3cities$apps_partn_m, useNA = "ifany")

# check sex work variable
table(data_3cities$data_pt, data_3cities$sex_work_d, useNA = "ifany")

# hiv status 
table(data_3cities$data_pt, data_3cities$hiv_cats, useNA = "ifany")

### fit model
negbin_model <- stan_model(file = "./src/helper_negbin_aggregate_agehiv.stan",
                           model_name = "negbin_partn")

# variables to use (apps_partn_m and apps_partn_d are removed from FU since it was not asked)
vars_model_base <- c("age_30_39", "age_40_49", "age_50_59", "age_60_",
                     "rel_y_excl", "rel_y_open", "rel_y_uncl", 
                     "hiv_cats",
                     "bath_m", "bath_d",
                     "groupsex_m", "groupsex_d", 
                     "apps_partn_m", "apps_partn_d",
                     "sex_work_m", "sex_work_d")
vars_model_fu <- setdiff(vars_model_base, 
                         c("apps_partn_m", "apps_partn_d"))

fit_bayes_ls <- create_city_list(CITIES_DATAPTS)

# prepare data_frame for aggregate data points
df_x_aggrt_ah <- data_3cities %>% 
    split(.$age_hiv)
data_x_aggrt_ah <- list()
  
for(cur_city in CITIES_DATAPTS){
  data_x_aggrt_ah_cur_city <- list()
  for(ah in AGEHIV){
      cur_city_index <- df_x_aggrt_ah[[ah]]$data_pt == cur_city
      cur_city_data <- df_x_aggrt_ah[[ah]][cur_city_index, ]
      data_x_aggrt_ah_cur_city[[ah]] <- cur_city_data %>%
        group_by(across(all_of(vars_model_base))) %>% 
        summarize(nb = n(), 
                  ipw_rds = sum(ipw_rds), 
                  .groups = "drop") %>% 
        select(nb, ipw_rds, all_of(vars_model_base))
      }
    
    cur_city_max_combo <- max(unname(unlist(lapply(data_x_aggrt_ah_cur_city, nrow))))
    
    data_x_aggrt_ah[[cur_city]] <- array(data = 0,
                                         dim = c(length(AGEHIV), cur_city_max_combo, ncol(data_x_aggrt_ah_cur_city[[ah]])),
                                         dimnames = list(AGEHIV, 1:cur_city_max_combo, colnames(data_x_aggrt_ah_cur_city[[ah]])))
    
    for(ah in AGEHIV){
      
      n_row <- nrow(data_x_aggrt_ah_cur_city[[ah]])
      n_col <- ncol(data_x_aggrt_ah_cur_city[[ah]])
      for(r in 1:n_row){
        for(c in 1:n_col){
          data_x_aggrt_ah[[cur_city]][ah, r, c] <- unname(unlist(data_x_aggrt_ah_cur_city[[ah]][r, c]))
        }
      }
    }} 

num_cores <- parallel::detectCores()
t0 <- Sys.time()

# for Pre-Pandemic
for(cur_city in CITIES_DATAPTS){
  # tracker
  if( grepl("-Pre-Pandemic", cur_city) ){
    print(gsub("-Pre-Pandemic", "", cur_city))
  }
  
  # choose Pre-Pandemic or follow-up variables
  if(grepl("-Pre-Pandemic", cur_city) ){
    vars_model <- vars_model_base
  } else {
    vars_model <- vars_model_fu
  }
  set.seed(777)
  fit_bayes_ls[[cur_city]] <- sampling(
      negbin_model,
      data = list(y = filter(data_3cities, data_pt == cur_city)[[outcome_var]],
                  # data on which model is fit
                  x = data_3cities[data_3cities$data_pt == cur_city, vars_model],
                  N = sum(data_3cities$data_pt == cur_city),
                  # data to compute predictions
                  n_ah = length(AGEHIV),
                  x_aggr_ah = data_x_aggrt_ah[[cur_city]][, , vars_model],
                  N_aggr_ah = dim(data_x_aggrt_ah[[cur_city]])[2],
                  K = length(vars_model),
                  x_end = 300,
                  ipc_rds_w_ah = data_x_aggrt_ah[[cur_city]][, , "ipw_rds"]),
      cores = num_cores,
      chains = 2, iter = 4000
    )
    }


t1 <- Sys.time()
t1 - t0 

## Inspect model convergence diagnostic ----
# convergence of model chains (traceplots)
for(cur_city in CITIES_DATAPTS){
  cur_p_trace_plot <- traceplot(fit_bayes_ls[[cur_city]], pars = c("alpha", "beta", "phi"))
  ggsave(sprintf("%s/model-checks-p6m-all-%s-%s.png", fig_path, which(cur_city == CITIES_DATAPTS), cur_city),
         cur_p_trace_plot, device = "png",
         height = 14, width = 30, units = "cm", dpi = 320)
}
rm(cur_p_trace_plot)

# r hat and effective sample size
ess_ls <- create_city_list(CITIES_DATAPTS)
for(cur_city in CITIES_DATAPTS){
  print(cur_city)
  # get model
  cur_model <- summary(fit_bayes_ls[[cur_city]])$summary
  
  # show only intercept and regression coefficients, ignore y_hat and y_pred
  row_param_names <- rownames(cur_model)
  row_param_names <- grep("alpha|beta|phi|shape|zi", row_param_names, value = T)
  
  # output
  print(round(cur_model[row_param_names, ], 3))
  
  # save
  ess_ls[[cur_city]] <- cur_model[row_param_names, c("mean", "se_mean", "2.5%", "50%", "97.5%", "n_eff")]
}
rm(cur_model, row_param_names)

## effective sample size
# save output for all coefficients
ess_ls_tbl <- vector("list", length(ess_ls))
for(i in 1:length(ess_ls)){
  ess_ls_tbl[[i]] <- as_tibble(ess_ls[[i]], rownames = "coeff")
  ess_ls_tbl[[i]] <- mutate(ess_ls_tbl[[i]], city.time = names(ess_ls)[i], .before = 1)
}
ess_tbl <- bind_rows(ess_ls_tbl)
rm(ess_ls_tbl)

write_csv(ess_tbl, sprintf("%s/stan_model_fit_ess-%s.csv", fig_path, out_distr_pref))

# save summary by city & time period
summarize_ess(ess_tbl, beta_only = F)
summarize_ess(ess_tbl, beta_only = T)

# Probability mass function ----
## density and PMF computations already performed in Stan

## Collapse PMF into single dataset ----
pmf_iter <- create_city_list(CITIES_DATAPTS)
pmf_wt_by_city <- create_city_list(CITIES_DATAPTS)

for(cur_city in CITIES_DATAPTS){
    
    for(ah in AGEHIV){
      # extract PMF iterations from stan
      index_ah <- which(AGEHIV == ah)
      pmf_tmp_ah <- extract(fit_bayes_ls[[cur_city]], pars = "pmf")$pmf[, index_ah, ]
      pmf_iter[[cur_city]][[ah]] <- pmf_tmp_ah
      
      # get credible intervals
      cred_l_ah <- vector("double", ncol(pmf_tmp_ah))
      cred_u_ah <- vector("double", ncol(pmf_tmp_ah))
      
      for(i in 1:ncol(pmf_tmp_ah)){
        cred_l_ah[i] <- quantile(pmf_tmp_ah[ ,i], .025)
        cred_u_ah[i] <- quantile(pmf_tmp_ah[ ,i], .975)
      }
      
      # create tibble with each city and time period
      pmf_wt_by_city[[cur_city]][[ah]] <- tibble(data_pt = cur_city,
                                                 age_hiv = ah,
                                                 y_pred = 0:300,
                                                 mean = colSums(pmf_tmp_ah) / nrow(pmf_tmp_ah),
                                                 cr.i_low = cred_l_ah,
                                                 cr.i_upp = cred_u_ah) %>%
        pivot_longer(cols = "age_hiv") %>%
        separate(value, into = c('age_cats', 'hiv_cats'), sep = "\\.")
      
      rm(cred_l_ah, cred_u_ah)}
    
    pmf_wt_by_city[[cur_city]] <- bind_rows(pmf_wt_by_city[[cur_city]])}
  
rm(pmf_tmp_ah)

# collapse into
pmf_wt_by_city <- bind_rows(pmf_wt_by_city)
pmf_wt_by_city$data_pt <- factor(pmf_wt_by_city$data_pt, levels = CITIES_DATAPTS)

## Verify PMF posterior distributions ----
# verify results by looking at the mean number of partners
group_var <- c("age_cats", "hiv_cats", "data_pt")

data_mean_nb_partn <- pmf_wt_by_city %>%
  group_by(across(all_of(group_var))) %>%
  summarize(mean_wt = sum(y_pred * mean),
            cr.i_low = sum(y_pred * cr.i_low),
            cr.i_upp = sum(y_pred * cr.i_upp),
            .groups = "drop") %>% 
  mutate(type = "neg. bin.")

data_mean_nb_partn

# verify that weights sum up to 1
pmf_wt_by_city %>% 
  group_by(across(all_of(group_var))) %>% 
  summarize(dens_ttl = sum(mean), 
            .groups = "drop") %>% 
  View()

# save full fitted pmf
write.csv(pmf_wt_by_city,
            sprintf("%s/pmf_weighted_%s.csv", out_distr_path, out_distr_pref),
            row.names = F)
