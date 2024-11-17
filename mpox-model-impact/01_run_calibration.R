# authors:  Fanyu Xiu, Mathieu Maheu-Giroux, Jorge Luis Flores

library(MpoxModelPack)
library(tidyverse)

run_cpp <- T # if run model in cpp or not
cal_cities <- c("mtl", "trt", "van")
# cal_analysis <- c("main", "contact_15", "contact_10", "standardized_vaccine_date", "VE_lb", "VE_ub", "VE_1", "prioritize_vaccine")
cal_analysis <- c(VE_1","prioritize_vaccine")

for(analysis in cal_analysis){
  
  t0 <- Sys.time()
  
  print(analysis)
  
  VE_analysis <- ifelse(analysis == "VE_lb", 0.358, ifelse(analysis == "VE_ub", 0.86, ifelse(analysis == "VE_1", 0.7187, 0.5149)))
  standardized_vaccine_date_analysis <- ifelse(analysis == "standardized_vaccine_date", T, F)
  contact_prop_analysis <- ifelse(analysis == "contact_15", 0.15, ifelse(analysis == "contact_10", 0.10, 0.20))
  contact_dur_analysis <- 2 
  prioritize_vaccine_analysis <- (analysis == "prioritize_vaccine")

  ## set starting value as the median of the prior
  set.seed(7)
  theta_ls <- par_ls <- list(mtl = NULL,
                             trt = NULL,
                             van = NULL)
  theta0 <- vector()
  theta0[1:3] <- median(rnorm(10000, mean = qlogis(0.5), sd = 0.5))
  theta0[4] <- median(rnorm(10000, mean = qlogis(0.87), sd = 1))
  theta0[8:10] <- median(rnorm(10000, mean = qlogis(5 / 100), sd = 1))
  theta0[5] <- median(rnorm(10000, mean = qlogis(0.67 / 0.80), sd = 5))
  theta0[6] <- median(rnorm(10000, mean = qlogis((0.80 - 0.47) / (1 - 0.47)), sd = 1.5))
  theta0[7] <- median(rnorm(10000, mean = qlogis((5 - 3) / (15 - 3)), sd = 1))

  # trial run
  llk_all(theta = theta0,
    cal_type = 1,
    cal_cpp = run_cpp,
    VE_llk = VE_analysis,
    contact_prop_llk = contact_prop_analysis,
    contact_dur_llk = contact_dur_analysis,
    if_prioritized_vac_llk = prioritize_vaccine_analysis)
  
  # Nelder-Mead better to get good starting values
  opt_nelder_mead <- optim(theta0, 
                         llk_all, 
                         cal_type = 1,
                         cal_cpp = run_cpp,
                         VE_llk = VE_analysis,
                         contact_prop_llk = contact_prop_analysis,
                         contact_dur_llk = contact_dur_analysis,
                         if_prioritized_vac_llk = prioritize_vaccine_analysis,
                         method = "Nelder-Mead", 
                         control = list(fnscale = - 1, trace = 0, maxit = 250), 
                         hessian = FALSE)
  
  # BFGS better to get the mode if supplied with good starting values
  opt <- optim(opt_nelder_mead$par, 
             llk_all, 
             cal_type = 1, 
             cal_cpp = run_cpp,
             VE_llk = VE_analysis,
             contact_prop_llk = contact_prop_analysis,
             contact_dur_llk = contact_dur_analysis,
             if_prioritized_vac_llk = prioritize_vaccine_analysis,
             method = "BFGS", 
             control = list(fnscale = -1, trace = 4, REPORT = 1, maxit = 250), 
             hessian = TRUE)
  
  theta <- opt$par # return MAP (maximum a posteriori = mode)
  
  for(cty in cal_cities){
  print(cty)
  index_city <- ifelse(cty == "mtl", 1, ifelse(cty == "trt", 2, 3))
  imported_low <- ifelse(cty == "van", 1, ifelse(cty == "mtl", 2, 3))
  imported_upp <- ifelse(cty == "van", 6, ifelse(cty == "mtl", 8, 8)) 
  fit_par <- data.frame(names = c("imported cases (tau)",
                                "transmission parameter (beta)", 
                                "assortativity (omega)", 
                                "rate ratio in high risk group (RR_H)", 
                                "rate ratio in low risk group (RR_L)", 
                                "duration infectiousness (1/gamma1)"),
                      value = round(c(imported_low + plogis(theta[index_city]) * (imported_upp - imported_low), 
                                      plogis(theta[4]),
                                      plogis(theta[index_city + 7]) * 100, 
                                      (0.7 + plogis(theta[5]) * (1 - 0.7)) * (0.47 + plogis(theta[6]) * (1 - 0.47)), 
                                      0.47 + plogis(theta[6]) * (1 - 0.47), 
                                      3 + plogis(theta[7]) * (15 - 3)),
                                    2)); print(fit_par)
  
  ## store fitted results
  # optimized parameters
  theta_ls[[cty]] <- theta # in the logit scale
  par_ls[[cty]] <- fit_par # transform back to regular scale
  
  init.pop.fn(cty, 1)
  load.params.fn(VE_analysis,
                 contact_prop_analysis,
                 standardized_vaccine_date_analysis,
                 prioritize_vaccine_analysis)
  period_mod_cty <- ifelse(cty == "mtl", 150, ifelse(cty == "trt", 170, 160))
  map_fit <- fn_model(city = cty,
                           import_cases_city = imported_low + plogis(theta[index_city]) * (imported_upp - imported_low),
                           bbeta_city = plogis(theta[4]),
                           omega_city = plogis(theta[index_city + 7]) * 100,
                           period_city = period_mod_cty,
                           RR_H_city = (0.7 + plogis(theta[5]) * (1 - 0.7)) * (0.47 + plogis(theta[6]) * (1 - 0.47)),
                           RR_L_city = 0.47 + plogis(theta[6]) * (1 - 0.47),
                           gamma1_city = 1 / (3 + plogis(theta[7]) * (15 - 3)),
                           VACCINATING = 1,
                           TRACING = 1,
                           cpp = run_cpp)
  province = ifelse(cty == "van", "BC", ifelse(cty == "trt", "ON", "QC"))
  plot_max = ifelse(cty == "trt", 20, 10)
  cases <- map_fit$cases[map_fit$time > days_imported[[cty]]]
  time <- seq(0, (length(cases) - 1) * 0.25 , 0.25)
  plot(cases ~ time,
       main = cty,
       ylim = c(0, plot_max))
  points(case_data[case_data$prov == province, ]$city_cases)
  }
  
  t1 <- Sys.time()
  print(t1 - t0)

  set.seed(77)
  ci_ls <- simul_fun_all(cal_type = 1,
                            hessian = opt$hessian, 
                            thetas = theta, 
                            sim = 1000,
                            SIR = TRUE, 
                            nsir = 5000, 
                            with_replacement = TRUE,
                            parallel = TRUE,
                            cal_cpp = run_cpp,
                            VE_sim = VE_analysis,
                            contact_prop_sim = contact_prop_analysis,
                            contact_dur_sim = contact_dur_analysis,
                            standardized_vaccine_date_sim = standardized_vaccine_date_analysis,
                            if_prioritized_vac_sim = prioritize_vaccine_analysis)
  print(ci_ls$par_ci)
  print(ci_ls$AF_ci)
  
  # add date to fitted model
  for(cty in cal_cities){
    init.pop.fn(cty, 1)
    load.params.fn(VE_analysis,
                   contact_prop_analysis,
                   standardized_vaccine_date_analysis,
                   prioritize_vaccine_analysis)
    ci_ls$result[[cty]]$date <- ci_ls$result[[cty]]$time - days_imported[cty] + min(case_city[[cty]]$date)
  }
  
  data_results <- list()
  data_mod_fit <- list()
  data_behavourial_fit <- list()
  data_vaccine_fit <- list()
  data_tracing_fit <- list()
  data_bt_fit <- list()
  data_bv_fit <- list()
  data_vt_fit <- list()
  data_nothing_fit <- list()
  data_AF <- list()
  data_xval_prop_case_age <- data_xval_prop_vaccine_age <- data_xval_prop_case_hiv <- data_xval_cum_vaccine <- list()
  
  # compile results into dataframes
  for(cty in cal_cities){
  city_name <- case_when(cty == "mtl" ~ "Montréal",
                        cty == "trt" ~ "Toronto",
                        cty == "van" ~ "Vancouver")
  
  data_results[[cty]] <- data.frame(estimate = c(ci_ls$par_ci[[cty]]$med, ci_ls$AF_ci[[cty]]$med),
                                    name = c(par_ls[[cty]]$names, ci_ls$AF_ci[[cty]]$names),
                                    lci = c(ci_ls$par_ci[[cty]]$lci, ci_ls$AF_ci[[cty]]$lci),
                                    uci =  c(ci_ls$par_ci[[cty]]$uci, ci_ls$AF_ci[[cty]]$uci),
                                    city_name = city_name)
  data_mod_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$med,
                                    cases_lci = ci_ls$result[[cty]]$lci,
                                    cases_uci = ci_ls$result[[cty]]$uci,
                                    time = ci_ls$result[[cty]]$time,
                                    date = ci_ls$result[[cty]]$date, 
                                    city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_behavourial_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$behavourial_med,
                                            cases_lci = ci_ls$result[[cty]]$behavourial_lci,
                                            cases_uci = ci_ls$result[[cty]]$behavourial_uci,
                                            time = ci_ls$result[[cty]]$time,
                                            date = ci_ls$result[[cty]]$date, 
                                            city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_vaccine_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$vaccine_med,
                                        cases_lci = ci_ls$result[[cty]]$vaccine_lci,
                                        cases_uci = ci_ls$result[[cty]]$vaccine_uci,
                                        time = ci_ls$result[[cty]]$time,
                                        date = ci_ls$result[[cty]]$date, 
                                        city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_tracing_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$tracing_med,
                                        cases_lci = ci_ls$result[[cty]]$tracing_lci,
                                        cases_uci = ci_ls$result[[cty]]$tracing_uci,
                                        time = ci_ls$result[[cty]]$time,
                                        date = ci_ls$result[[cty]]$date, 
                                        city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_bt_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$bt_med,
                                        cases_lci = ci_ls$result[[cty]]$bt_lci,
                                        cases_uci = ci_ls$result[[cty]]$bt_uci,
                                        time = ci_ls$result[[cty]]$time,
                                        date = ci_ls$result[[cty]]$date, 
                                        city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_bv_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$bv_med,
                                        cases_lci = ci_ls$result[[cty]]$bv_lci,
                                        cases_uci = ci_ls$result[[cty]]$bv_uci,
                                        time = ci_ls$result[[cty]]$time,
                                        date = ci_ls$result[[cty]]$date, 
                                        city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_vt_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$vt_med,
                                        cases_lci = ci_ls$result[[cty]]$vt_lci,
                                        cases_uci = ci_ls$result[[cty]]$vt_uci,
                                        time = ci_ls$result[[cty]]$time,
                                        date = ci_ls$result[[cty]]$date, 
                                        city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_nothing_fit[[cty]] <- data.frame(cases = ci_ls$result[[cty]]$nothing_med,
                                        cases_lci = ci_ls$result[[cty]]$nothing_lci,
                                        cases_uci = ci_ls$result[[cty]]$nothing_uci,
                                        time = ci_ls$result[[cty]]$time,
                                        date = ci_ls$result[[cty]]$date, 
                                        city_name = city_name) %>% 
    mutate(fit_target = (time > days_imported[cty]))
  data_AF[[cty]] <- data.frame(names = ci_ls$AF_ci[[cty]]$names, 
                              estimate = ci_ls$AF_ci[[cty]]$med, 
                              lci = ci_ls$AF_ci[[cty]]$lci,
                              uci = ci_ls$AF_ci[[cty]]$uci,
                              city_name = city_name)
  data_xval_prop_case_age[[cty]] <- data.frame(grp = names_age_cats,
                                               estimate = ci_ls$xval_prop_case_age[[cty]]$med,
                                         lci = ci_ls$xval_prop_case_age[[cty]]$lci,
                                         uci = ci_ls$xval_prop_case_age[[cty]]$uci,
                                         city_name = city_name)
  data_xval_prop_vaccine_age[[cty]] <- data.frame(grp = names_age_cats,
                                                  estimate = ci_ls$xval_prop_vaccine_age[[cty]]$med,
                                                  lci = ci_ls$xval_prop_vaccine_age[[cty]]$lci,
                                                  uci = ci_ls$xval_prop_vaccine_age[[cty]]$uci,
                                                  city_name = city_name)
  data_xval_prop_case_hiv[[cty]] <- data.frame(grp = names_hiv_cats,
                                               estimate = ci_ls$xval_prop_case_hiv[[cty]]$med,
                                               lci = ci_ls$xval_prop_case_hiv[[cty]]$lci,
                                               uci = ci_ls$xval_prop_case_hiv[[cty]]$uci,
                                               city_name = city_name)
  
  data_xval_cum_vaccine[[cty]] <- data.frame(grp = "vaccine",
                                             estimate = ci_ls$xval_cum_vaccine[[cty]]$med,
                                                lci = ci_ls$xval_cum_vaccine[[cty]]$lci,
                                                uci = ci_ls$xval_cum_vaccine[[cty]]$uci,
                                                city_name = city_name)
  
  }
  # list to df
data_results <- do.call(rbind.data.frame, c(data_results, make.row.names = FALSE))
data_mod_fit <- do.call(rbind.data.frame, c(data_mod_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0) 
data_behavourial_fit <- do.call(rbind.data.frame, c(data_behavourial_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0)
data_vaccine_fit <- do.call(rbind.data.frame, c(data_vaccine_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0)
data_tracing_fit <- do.call(rbind.data.frame, c(data_tracing_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0)
data_bt_fit <- do.call(rbind.data.frame, c(data_bt_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0)
data_bv_fit <- do.call(rbind.data.frame, c(data_bv_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0)
data_vt_fit <- do.call(rbind.data.frame, c(data_vt_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0)
data_nothing_fit <- do.call(rbind.data.frame, c(data_nothing_fit, make.row.names = FALSE)) %>% 
  subset(fit_target & time %% 1 == 0)
data_AF <- do.call(rbind.data.frame, c(data_AF, make.row.names = FALSE)) 
data_xval_prop_case_age <- do.call(rbind.data.frame, c(data_xval_prop_case_age, make.row.names = FALSE)) 
data_xval_prop_vaccine_age <- do.call(rbind.data.frame, c(data_xval_prop_vaccine_age, make.row.names = FALSE)) 
data_xval_prop_case_hiv <- do.call(rbind.data.frame, c(data_xval_prop_case_hiv, make.row.names = FALSE)) 
data_xval_cum_vaccine <- do.call(rbind.data.frame, c(data_xval_cum_vaccine, make.row.names = FALSE)) 

saveRDS(data_results, sprintf("./out/%s/data_results.rds", analysis))
saveRDS(data_mod_fit, sprintf("./out/%s/data_mod_fit.rds", analysis))
saveRDS(data_behavourial_fit, sprintf("./out/%s/data_behavourial_fit.rds", analysis))
saveRDS(data_vaccine_fit, sprintf("./out/%s/data_vaccine_fit.rds", analysis))
saveRDS(data_tracing_fit, sprintf("./out/%s/data_tracing_fit.rds", analysis))
saveRDS(data_bt_fit, sprintf("./out/%s/data_bt_fit.rds", analysis))
saveRDS(data_bv_fit, sprintf("./out/%s/data_bv_fit.rds", analysis))
saveRDS(data_vt_fit, sprintf("./out/%s/data_vt_fit.rds", analysis))
saveRDS(data_nothing_fit, sprintf("./out/%s/data_nothing_fit.rds", analysis))
saveRDS(data_AF, sprintf("./out/%s/data_AF.rds", analysis))
saveRDS(data_xval_prop_case_age, sprintf("./out/%s/data_xval_prop_case_age.rds", analysis))
saveRDS(data_xval_prop_vaccine_age, sprintf("./out/%s/data_xval_prop_vaccine_age.rds", analysis))
saveRDS(data_xval_prop_case_hiv, sprintf("./out/%s/data_xval_prop_case_hiv.rds", analysis))
saveRDS(data_xval_cum_vaccine, sprintf("./out/%s/data_xval_cum_vaccine.rds", analysis))

# wgt_resample <- ci_ls$wgt_resample
# hist(wgt_resample, breaks = 100, xlim = c(0, 0.05))

t2 <- Sys.time()
print(analysis)
print(t2 - t1)
} # for analysis
