#' @title prior_dens_fix
#' @author: Mathieu Maheu-Giroux, Jorge Luis Flores Anato, Fanyu Xiu
#' @description sum of log prior density for a given set of parameter values 
#' @param theta value of a vector of parameter in logit scale 
#' @return a value equals to log prior density
#' @rdname prior_dens_fix
#' @importFrom foreach %dopar%
#' @export 
prior_dens_fix <- function(theta) {
  log_prior <- dnorm(theta[4], mean = qlogis(0.87), sd = 1, log = TRUE) + # beta (transmission parameter)
    dnorm(theta[5], mean = log(1.5), sd = 1, log = TRUE) + # omega (mixing)
    dnorm(theta[6], mean = qlogis(0.67 / 0.80), sd = 5, log = TRUE) + # RR_multiplier
    dnorm(theta[7], mean = qlogis((0.80 - 0.47) / (1 - 0.47)), sd = 1.5, log = TRUE) +
    dnorm(theta[8], mean = qlogis((5 - 3) / (15 - 3)), sd = 1, log = TRUE)  # duration of infectiousness (median ~ 5 days) (using transformed qlogis to bound within reasonable values)
    return(log_prior)
}

#' @title llk_all
#' @description function to output posterior density of a parameter set using negative binomial likelihood 
#' @param theta value of a vector of parameter in logit scale 
#' @param cal_type which type of sexual act to be selected (all-type:1 or anal:2), Default: 1
#' @param cal_cpp whether to run calibration with rcpp (1) or R (0)
#' @param VE_llk vaccine effectiveness to run the calibration (for sensitivity analysis)
#' @param contact_prop_llk proportion contact traced and isolated to run the calibration (for sensitivity analysis)
#' @param contact_dur_llk days it takes for contact tracing (for sensitivity analysis)
#' @return a value equals to the log posterior density
#' @rdname llk_all
#' @export 
llk_all <- function(theta, 
                cal_type = 1, 
                cal_cpp,
                VE_llk,
                contact_prop_llk,
                contact_dur_llk) {
  
  log_prior_fix <- prior_dens_fix(theta)
  
  post_llk_sum <- 0
  
  for(cal_city in c("mtl", "trt", "van")){
    
    index_city <- ifelse(cal_city == "mtl", 1, ifelse(cal_city == "trt", 2, 3))
    log_prior_vary <- dnorm(theta[index_city], mean = qlogis(0.5), sd = 0.5, log = TRUE)
    
    imported_low <- ifelse(cal_city == "van", 1, ifelse(cal_city == "mtl", 2, 3))
    imported_upp <- ifelse(cal_city == "van", 3, ifelse(cal_city == "mtl", 6, 6))
    
    # getting predictions
    # load population in the selected cities
    init.pop.fn(city = cal_city, 
                type = cal_type)
    
    # load model parameter of the selected cities
    load.params.fn(VE_llk, 
                   contact_prop_llk,
                   contact_dur_llk)
    
    # change parameters to be calibrated in the function input
    model_output <- fn_model(city = cal_city,
                             import_cases_city = imported_low + plogis(theta[index_city]) * (imported_upp - imported_low),
                             bbeta_city = plogis(theta[4]),
                             omega_city = exp(theta[5]),
                             RR_H_city = (0.7 + plogis(theta[6]) * (1 - 0.7)) * (0.47 + plogis(theta[7]) * (1 - 0.47)),
                             RR_L_city = 0.47 + plogis(theta[7]) * (1 - 0.47),
                             gamma1_city = 1 / (3 + plogis(theta[8]) * (15 - 3)),
                             period_city = 150,
                             VACCINATING = 1,
                             TRACING = 1,
                             cpp = cal_cpp)
    
    # only look at modeled cases after imported cases
    # i.e. only look at modeled cases since the first reported case
    mod_inc <- model_output$cases[model_output$time %in% case_city[[cal_city]]$time_intro]

    # log-likelihood
    # log_lik_city <- sum(dpois(x = case_city[[cal_city]]$incidence,
    #                      lambda = mod_inc, 
    #                      log = TRUE))
    log_lik_city <- sum(dnbinom(x = case_city[[cal_city]]$incidence,
                              mu = mod_inc, 
                              size = 0.1,
                              log = TRUE))
    # posterior in log scale
    post_llk_sum <- log_lik_city + log_prior_vary + post_llk_sum
  }
  
  post_llk_sum <- post_llk_sum + log_prior_fix
  
  return(post_llk_sum) ## return the log posterior
}

#' @title getci
#' @description Function for obtaining 95% credible interval of the estimates
#' @param df a dataframe object
#' @return a credible interval (vector) of the estimates
#' @rdname getci
#' @export 
getci <- function(df) {
  ci <- apply(X = df, 
              MARGIN = 2,
              FUN = quantile, 
              probs = c(0.025, 0.5, 0.975), 
              na.rm = TRUE)
  lower <- ci[1, ]
  med <- ci[2, ]
  upper <- ci[3, ]
  df1 <- data.frame(lower, med, upper)
  return(df1)
}

#' @title simul_fun_all
#' @description simulation to obtain credible interval of all parameter sets (CrI)
#' @param cal_type which type of sexual act to be selected (all-type:1 or anal:2), Default: 1
#' @param hessian hessian matrix returned by the optim function
#' @param thetas maximum a posterior (mode) of parameters to be calibrated
#' @param sim number of resamples, Default: 1000
#' @param SIR use sampling importance resampling algorithm, Default: TRUE
#' @param nsir number of draws from the importance function, Default: 5000
#' @param with_replacement resample with replacement, Default: TRUE
#' @param track_sims tracking simulations , Default: FALSE
#' @param parallel to fasten the loops, Default: TRUE
#' @param cal_cpp whether to run calibration with rcpp (1) or R (0)
#' @param VE_sim vaccine effectiveness to run the calibration (for sensitivity analysis)
#' @param contact_prop_sim proportion contact traced and isolated to run the calibration (for sensitivity analysis)
#' @param contact_dur_sim days it takes for contact tracing to run the calibration (for sensitivity analysis)
#' @param standardized_vaccine_date_sim whether or not standardizing vaccination start dates across the cities (for sensitivity analysis)
#' @return a list of result (model fit for each parameter sets), par_ci (posterior CrI for parameters), AF_ci (posterior CrI of averted fractions), samples (all resampled parameter sets)
#' @seealso 
#'  \code{\link[Matrix]{solve-methods}}, \code{\link[Matrix]{nearPD}}
#'  \code{\link[mvtnorm]{Mvnorm}}, \code{\link[mvtnorm]{Mvt}}
#'  \code{\link[parallel]{makeCluster}}, \code{\link[parallel]{detectCores}}, \code{\link[parallel]{clusterApply}}
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[Rcpp]{sourceCpp}}
#'  \code{\link[foreach]{foreach}}
#' @rdname simul_fun_all
#' @export 
#' @importFrom Matrix solve nearPD
#' @importFrom mvtnorm rmvnorm rmvt dmvt
#' @importFrom parallel makeCluster detectCores clusterCall clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach
simul_fun_all <- function(cal_type,
                      hessian, 
                      thetas, 
                      sim = 1000, 
                      SIR = TRUE, 
                      nsir = 5000, 
                      with_replacement = TRUE, 
                      track_sims = FALSE,
                      parallel = TRUE, 
                      cal_cpp,
                      VE_sim,
                      contact_prop_sim,
                      contact_dur_sim,
                      standardized_vaccine_date_sim
                      ) { 
  
  # From the hessian, simulate the model
  vcova <- Matrix::solve(-hessian)
  
  # Test if positive semi-definite
  eS <- eigen(vcova, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -1e-06 * abs(ev[1L]))) {
    vcova <- Matrix::nearPD(vcova, corr = FALSE)$mat
    vcova <- matrix(vcova@x, 
                    nrow = vcova@Dim[1], 
                    ncol = vcova@Dim[2])
  }
  
  # Laplace approximation only (assuming posterior follows a Gaussian)
  if (SIR == FALSE) {
    samp <- mvtnorm::rmvnorm(n = sim, thetas, vcova)
  }
  # SIR method, without constraint from the Gaussian assumption
  if (SIR == TRUE) {
    ## write the function for normalizing the sample
    ## to prevent overflow of log-exp transformation
    ## by removing the maximum x when calculating exp
    log_sum_exp <- function(x) {
      xmax <- which.max(x)
      log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax] }
    
    # Sample parameters from multivariate t-distribution (the importance function in this case)
    # to get thicker tails (by using df = 2)
    par_sir <- mvtnorm::rmvt(n = nsir, 
                             delta = thetas, 
                             sigma = as.matrix(vcova), 
                             df = 2)
    # add the mode to the resampled values (to be sure that it is included in the CI)
    par_sir <- rbind(par_sir, thetas)
    # calculate density of the sampled points from the importance function
    prp_dens <- mvtnorm::dmvt(par_sir, 
                              thetas, 
                              as.matrix(vcova), 
                              df = 2, 
                              log = TRUE)
    
    if (parallel == FALSE)  {
      # You might get some warnings, that's OK, the function will return NA
      # and we replace them with -Inf in the llk. We suppress them.
     # it's normal to get -inf when value far away
      loglikelihood_sir <- suppressWarnings(
        apply(
          ### for each parameter to be calibrated
          par_sir, MARGIN = 1, FUN = function(x) {
            ## calculate the log posterior of its mode
            llk_all(theta = thetas,
                    cal_type = cal_type,
                    cal_cpp = cal_cpp,
                    VE_sim,
                    contact_prop_sim,
                    contact_dur_sim)}
        )
      )
    } else {
      cl <- parallel::makeCluster(parallel::detectCores()) # use all cores
      doParallel::registerDoParallel(cl)
      parallel::clusterCall(cl, function() library(magrittr))
      parallel::clusterExport(cl, c("contact_rate_prop", "case_data", "A_matrix", "llk_all", "prior_dens_fix", "init.pop.fn", "load.params.fn", "fn_model"))
      loglikelihood_sir <- foreach::foreach(i = 1:nrow(par_sir), .combine = "c") %dopar% {
        par_sir_i <- par_sir[i, ]
        llk_all(theta = par_sir_i,
            cal_type = cal_type,
            cal_cpp = cal_cpp,
            VE_sim,
            contact_prop_sim,
            contact_dur_sim)}
      parallel::stopCluster(cl)
    }
    
    loglikelihood_sir[!is.finite(loglikelihood_sir)] <- -Inf
    dens_ratio <- loglikelihood_sir - prp_dens ### difference between the log posterior and proposal density
    wgt <- exp(dens_ratio - log_sum_exp(dens_ratio))
    
    ### resample proportionally to the weight
    resampleid <- sample(nrow(par_sir), # sample from all sets of the parameters drawn from the importance function
                         size = sim, # size of resampling
                         replace = with_replacement, 
                         prob = wgt)
    wgt_resample <- wgt[resampleid]
    samp <- par_sir[resampleid, , drop = FALSE]
    ### value with higher proportion may be sampled many times
    ### so only count unique sample values
    nunique <- length(unique(resampleid))
    max_wgt <- max(wgt, na.rm = TRUE)
    if (max_wgt > 0.1) { print(paste('Caution, maximum weight is ',
                                     round(max_wgt, 2),
                                     " consider increasing sim or sampling without replacement (latter option is not the best)")) }
    print(paste(nunique, 'unique parameter sets resampled'))
  }     
  
  posterior_par_ci <- list()
  posterior_AF_ci <- list()
  result <- list()
  xval_prop_case_age <- xval_prop_case_hiv <- xval_cum_vaccine <- xval_prop_vaccine_age <- list()
  
  for(cty in c("mtl", "trt", "van")){
    
    index_city <- ifelse(cty == "mtl", 1, ifelse(cty == "trt", 2, 3))
    
    imported_low <- ifelse(cty == "van", 1, ifelse(cty == "mtl", 2, 3))
    imported_upp <- ifelse(cty == "van", 3, ifelse(cty == "mtl", 6, 6))
    
    ### make an empty dataset for storing modeled cases 
    ### using each set of resampled parameters
    init.pop.fn(cty, 1)
    load.params.fn(VE_sim,
                   contact_prop_sim,
                   contact_dur_sim,
                   standardized_vaccine_date_sim) 
    
    # reruning the model covering longer period
    period_mod_cty <- ifelse(cty == "mtl", 150, ifelse(cty == "trt", 170, 160))
    end_mod_cty <- days_imported[cty] + period_mod_cty
    val_case <- matrix(data = NA, 
                       nrow = sim, 
                       ncol = (end_mod_cty - sstart) / ddt + 1)

    val_bv_case <- val_bt_case <- val_vt_case <- val_behavourial_case <- val_vaccine_case <- val_tracing_case <- val_nothing_case <- val_inc <- val_case
    
    val_cum_vaccine <- rep(NA, sim)
    val_prop_vaccine_age <- matrix(data = NA, 
                                 nrow = sim, 
                                 ncol = n_age_cats,
                                 dimnames = list(NULL, names_age_cats))
    val_prop_case_age <- val_prop_vaccine_age
    val_prop_case_hiv <- matrix(data = NA, 
                                nrow = sim, 
                                ncol = n_hiv_cats,
                                dimnames = list(NULL, names_hiv_cats))
    AF_behavourial <- rep(NA, sim)
    AF_bt <- AF_vt <- AF_bv <- AF_vaccine <- AF_tracing <- AF_all <- AF_behavourial
    
    for (s in 1:sim) {
      
      if(track_sims & s %% 100 == 0){
        print(sprintf("%s out of %s sims", s, sim))}
      
      tmp <- fn_model(city = cty,
                      import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                      bbeta_city = plogis(samp[s, 4]),
                      omega_city = exp(samp[s, 5]),
                      period_city = period_mod_cty,
                      RR_H_city = (0.7 + plogis(samp[s, 6]) * (1 - 0.7)) * (0.47 + plogis(samp[s, 7]) * (1 - 0.47)),
                      RR_L_city = 0.47 + plogis(samp[s, 7]) * (1 - 0.47),
                      gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                      VACCINATING = 1,
                      TRACING = 1,
                      cpp = cal_cpp)
      tmp_behavourial <- fn_model(city = cty,
                                  import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                                  bbeta_city = plogis(samp[s, 4]),
                                  omega_city = exp(samp[s, 5]),
                                  period_city = period_mod_cty,
                                  RR_H_city = (0.7 + plogis(samp[s, 6]) * (1 - 0.7)) * (0.47 + plogis(samp[s, 7]) * (1 - 0.47)),
                                  RR_L_city = 0.47 + plogis(samp[s, 7]) * (1 - 0.47),
                                  gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                                  VACCINATING = 0,
                                  TRACING = 0,
                                  cpp = cal_cpp)
      tmp_vaccine <- fn_model(city = cty,
                              import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                              bbeta_city = plogis(samp[s, 4]),
                              omega_city = exp(samp[s, 5]),
                              period_city = period_mod_cty,
                              RR_H_city = 1,  
                              RR_L_city = 1,
                              gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                              VACCINATING = 1,
                              TRACING = 0,
                              cpp = cal_cpp)
      tmp_tracing <- fn_model(city = cty,
                              import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                              bbeta_city = plogis(samp[s, 4]),
                              omega_city = exp(samp[s, 5]),
                              period_city = period_mod_cty,
                              RR_H_city = 1,
                              RR_L_city = 1,
                              gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                              VACCINATING = 0,
                              TRACING = 1,
                              cpp = cal_cpp)
      tmp_bv <- fn_model(city = cty,
                                  import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                                  bbeta_city = plogis(samp[s, 4]),
                                  omega_city = exp(samp[s, 5]),
                                  period_city = period_mod_cty,
                                  RR_H_city = (0.7 + plogis(samp[s, 6]) * (1 - 0.7)) * (0.47 + plogis(samp[s, 7]) * (1 - 0.47)),
                                  RR_L_city = 0.47 + plogis(samp[s, 7]) * (1 - 0.47),
                                  gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                                  VACCINATING = 1,
                                  TRACING = 0,
                                  cpp = cal_cpp)
      tmp_bt <- fn_model(city = cty,
                         import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                         bbeta_city = plogis(samp[s, 4]),
                         omega_city = exp(samp[s, 5]),
                         period_city = period_mod_cty,
                         RR_H_city = (0.7 + plogis(samp[s, 6]) * (1 - 0.7)) * (0.47 + plogis(samp[s, 7]) * (1 - 0.47)),
                         RR_L_city = 0.47 + plogis(samp[s, 7]) * (1 - 0.47),
                         gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                         VACCINATING = 0,
                         TRACING = 1,
                         cpp = cal_cpp)
      tmp_vt <- fn_model(city = cty,
                         import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                         bbeta_city = plogis(samp[s, 4]),
                         omega_city = exp(samp[s, 5]),
                         period_city = period_mod_cty,
                         RR_H_city = 1,
                         RR_L_city = 1,
                         gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                         VACCINATING = 1,
                         TRACING = 1,
                         cpp = cal_cpp)
      tmp_nothing <- fn_model(city = cty,
                              import_cases_city = imported_low + plogis(samp[s, index_city]) * (imported_upp - imported_low),
                              bbeta_city = plogis(samp[s, 4]),
                              omega_city = exp(samp[s, 5]),
                              period_city = period_mod_cty,
                              RR_H_city = 1,
                              RR_L_city = 1,
                              gamma1_city = 1 / (3 + plogis(samp[s, 8]) * (15 - 3)),
                              VACCINATING = 0,
                              TRACING = 0,
                              cpp = cal_cpp)
      
      val_case[s, ] <- tmp$cases

      val_cum_vaccine[s] <- sum(tmp$X["V", niter_city, , , ])
      
      for(a in names_age_cats){
        val_prop_vaccine_age[s, a] <- sum(tmp$X["V", niter_city, a, , ]) / sum(tmp$X["V", niter_city, , , ])
        val_prop_case_age[s, a] <- sum(tmp$X["CS", niter_city, a, , ]) / sum(tmp$X["CS", niter_city, , , ])
      }
      for(h in names_hiv_cats){
        val_prop_case_hiv[s, h] <- sum(tmp$X["CS", niter_city, , , h]) / sum(tmp$X["CS", niter_city, , , ])
      }
      
      val_behavourial_case[s, ] <- tmp_behavourial$cases
      val_vaccine_case[s, ] <- tmp_vaccine$cases
      val_tracing_case[s, ] <- tmp_tracing$cases
      val_bt_case[s, ] <- tmp_bt$cases
      val_bv_case[s, ] <- tmp_bv$cases
      val_vt_case[s, ] <- tmp_vt$cases
      val_nothing_case[s, ] <- tmp_nothing$cases
      
      val_inc[s, ] <- tmp$inc 
      AF_behavourial[s] <- sum(tmp_nothing$inc - tmp_behavourial$inc) / sum(tmp_nothing$inc)
      AF_vaccine[s] <- sum(tmp_nothing$inc - tmp_vaccine$inc) / sum(tmp_nothing$inc)
      AF_tracing[s] <- sum(tmp_nothing$inc - tmp_tracing$inc) / sum(tmp_nothing$inc)
      AF_bv[s] <- sum(tmp_nothing$inc - tmp_bv$inc) / sum(tmp_nothing$inc)
      AF_bt[s] <- sum(tmp_nothing$inc - tmp_bt$inc) / sum(tmp_nothing$inc)
      AF_vt[s] <- sum(tmp_nothing$inc - tmp_vt$inc) / sum(tmp_nothing$inc)
      AF_all[s] <- sum(tmp_nothing$inc - tmp$inc) / sum(tmp_nothing$inc)
      }

    res <- getci(val_case)
    
    res_prop_case_age <- getci(val_prop_case_age)
    res_prop_case_hiv <- getci(val_prop_case_hiv)
    res_cum_vaccine <- quantile(val_cum_vaccine, c(0.025, 0.5, 0.975))
    res_prop_vaccine_age <- getci(val_prop_vaccine_age)
    
    res_behavourial_case <- getci(val_behavourial_case)
    res_vaccine_case <- getci(val_vaccine_case)
    res_tracing_case <- getci(val_tracing_case)
    res_nothing_case <- getci(val_nothing_case)
    
    res_bv_case <- getci(val_bv_case)
    res_bt_case <- getci(val_bt_case)
    res_vt_case <- getci(val_vt_case)
    
    res_behavourial <- quantile(AF_behavourial, c(0.025, 0.5, 0.975))
    res_vaccine <- quantile(AF_vaccine, c(0.025, 0.5, 0.975))
    res_tracing <- quantile(AF_tracing, c(0.025, 0.5, 0.975))
    res_bv <- quantile(AF_bv, c(0.025, 0.5, 0.975))
    res_bt <- quantile(AF_bt, c(0.025, 0.5, 0.975))
    res_vt <- quantile(AF_vt, c(0.025, 0.5, 0.975))
    res_all <- quantile(AF_all, c(0.025, 0.5, 0.975))
    
    result[[cty]] <- data.frame(time = seq(sstart, end_mod_cty, ddt),
                         lci = res$lower,
                         med = res$med,
                         uci = res$upper,
                         behavourial_lci = res_behavourial_case$lower,
                         behavourial_med = res_behavourial_case$med,
                         behavourial_uci = res_behavourial_case$upper,
                         vaccine_lci = res_vaccine_case$lower,
                         vaccine_med = res_vaccine_case$med,
                         vaccine_uci = res_vaccine_case$upper,
                         tracing_lci = res_tracing_case$lower,
                         tracing_med = res_tracing_case$med,
                         tracing_uci = res_tracing_case$upper,
                         bt_lci = res_bt_case$lower,
                         bt_med = res_bt_case$med,
                         bt_uci = res_bt_case$upper,
                         bv_lci = res_bv_case$lower,
                         bv_med = res_bv_case$med,
                         bv_uci = res_bv_case$upper,
                         vt_lci = res_vt_case$lower,
                         vt_med = res_vt_case$med,
                         vt_uci = res_vt_case$upper,
                         nothing_lci = res_nothing_case$lower,
                         nothing_med = res_nothing_case$med,
                         nothing_uci = res_nothing_case$upper)
    
    par_ci <- getci(samp)
    
    posterior_par_ci[[cty]] <- data.frame(names = c("imported cases (tau)",
                                                    "transmission parameter (beta)",
                                                    "assortativity (omega)", 
                                                    "rate ratio in high risk group (RR_H)", 
                                                    "rate ratio in low risk group (RR_L)", 
                                                    "duration infectiousness (1/gamma1)"),
                              lci = c(imported_low + plogis(par_ci$lower[index_city]) * (imported_upp - imported_low),
                                      plogis(par_ci$lower[4]),
                                      exp(par_ci$lower[5]), 
                                      (0.7 + plogis(par_ci$lower[6]) * (1 - 0.7)) * (0.47 + plogis(par_ci$lower[7]) * (1 - 0.47)), 
                                      0.47 + plogis(par_ci$lower[7]) * (1 - 0.47), 
                                      3 + plogis(par_ci$lower[8]) * (15 - 3)),
                              med = c(imported_low + plogis(par_ci$med[index_city]) * (imported_upp - imported_low),
                                      plogis(par_ci$med[4]),
                                      exp(par_ci$med[5]), 
                                      (0.7 + plogis(par_ci$med[6]) * (1 - 0.7)) * (0.47 + plogis(par_ci$med[7]) * (1 - 0.47)), 
                                      0.47 + plogis(par_ci$med[7]) * (1 - 0.47), 
                                      3 + plogis(par_ci$med[8]) * (15 - 3)),
                              uci = c(imported_low + plogis(par_ci$upper[index_city]) * (imported_upp - imported_low),
                                      plogis(par_ci$upper[4]),
                                      exp(par_ci$upper[5]), 
                                      (0.7 + plogis(par_ci$upper[6]) * (1 - 0.7)) * (0.47 + plogis(par_ci$upper[7]) * (1 - 0.47)), 
                                      0.47 + plogis(par_ci$upper[7]) * (1 - 0.47), 
                                      3 + plogis(par_ci$upper[8]) * (15 - 3)))
    posterior_AF_ci[[cty]] <-  data.frame(names = c("behavourial change", 
                                             "first-dose vaccination", 
                                             "contact tracing",
                                             "behavourial change and contact tracing",
                                             "behavourial change and first-dose vaccination",
                                             "contact tracing and first-dose vaccination", 
                                             "all three combined"),
                                   lci = c(res_behavourial["2.5%"],
                                           res_vaccine["2.5%"],
                                           res_tracing["2.5%"],
                                           res_bt["2.5%"],
                                           res_bv["2.5%"],
                                           res_vt["2.5%"],
                                           res_all["2.5%"]),
                                   med = c(res_behavourial["50%"],
                                           res_vaccine["50%"],
                                           res_tracing["50%"],
                                           res_bt["50%"],
                                           res_bv["50%"],
                                           res_vt["50%"],
                                           res_all["50%"]),
                                   uci = c(res_behavourial["97.5%"],
                                           res_vaccine["97.5%"],
                                           res_tracing["97.5%"],
                                           res_bt["97.5%"],
                                           res_bv["97.5%"],
                                           res_vt["97.5%"],
                                           res_all["97.5%"]))
    
    xval_prop_case_age[[cty]] <- data.frame(lci = res_prop_case_age$lower,
                                            med = res_prop_case_age$med,
                                            uci = res_prop_case_age$upper)
      
    xval_prop_case_hiv[[cty]] <- data.frame(lci = res_prop_case_hiv$lower,
                                            med = res_prop_case_hiv$med,
                                            uci = res_prop_case_hiv$upper)
    
    xval_cum_vaccine[[cty]] <- data.frame(lci = res_cum_vaccine["2.5%"],
                                          med = res_cum_vaccine["50%"],
                                          uci = res_cum_vaccine["97.5%"])
    
    xval_prop_vaccine_age[[cty]]<- data.frame(lci = res_prop_vaccine_age$lower,
                                              med = res_prop_vaccine_age$med,
                                              uci = res_prop_vaccine_age$upper)
    
  } # for city
   
  return(list(wgt_resample = wgt_resample,
              result = result,
              par_ci = posterior_par_ci,
              AF_ci = posterior_AF_ci,
              xval_prop_case_age = xval_prop_case_age,
              xval_prop_vaccine_age = xval_prop_vaccine_age,
              xval_prop_case_hiv = xval_prop_case_hiv,
              xval_cum_vaccine = xval_cum_vaccine,
              samples = samp))
}