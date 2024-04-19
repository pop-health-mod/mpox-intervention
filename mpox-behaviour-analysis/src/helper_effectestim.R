# author: Fanyu Xiu
## produce OR/RR with CIs from models -----
library(tidyverse)

### transform model summary into data frame form
tbl_coeff_fn <- function(model, 
                         alpha = 0.05,
                         stat = "bayesian"){
  
  z = qnorm(1 - alpha/2)
  
  if(stat == "bayesian"){
    tbl_coeff <- summary(model, 
            probs = c(alpha/2, 1 - alpha/2)) %>% 
      as.data.frame()
    tbl_coeff <- tbl_coeff[!grepl("id|PPD|log", rownames(tbl_coeff)), ]
    tbl_coeff <- tbl_coeff %>% 
      rename(se = sd,
             mean_log = mean,
             conf_ub = paste0((1-alpha/2)*100, "%"),
             conf_lb = paste0(alpha/2*100, "%"))
    
      name = rownames(tbl_coeff)
      
      tbl_coeff <- tbl_coeff %>%
        mutate(name = name,
             mean = exp(mean_log),
             conf_ub = exp(conf_ub),
             conf_lb = exp(conf_lb),
             .before = 1)
  }else{
  tbl_coeff <- data.frame(summary(model)$coefficients) %>%
      dplyr::select(estimate = which(grepl("Estimate", colnames(.), ignore.case = T)),
                    se = which(grepl("Std", colnames(.), ignore.case = T)),
                    pvalue = which(grepl("p", colnames(.), ignore.case = T))) %>% 
      mutate(name = rownames(summary(model)$coefficients),
             mean = exp(estimate),
             conf_lb = exp(estimate - z * se),
             conf_ub = exp(estimate + z * se),
             .before = 1) %>%
      rename(mean_log = estimate)
    
  }
  rownames(tbl_coeff) <- NULL
  return(tbl_coeff)
}

### an exposure's OR/RR when there is an interaction term
#### output is a dataframe showing RR/OR for each level of the covariate relative to the exposure
interaction_fn <- function(tbl_coeff, # data frame output from the function above
                           matrix_cov, # covariance matrix from the function above
                           # a continuous exposure and a categorical covariate
                           exposure_name = "period_mpox_centre", 
                           covariate_name, 
                           covariate_level,
                           covariate_ref,
                           alpha = 0.05,
                           stat = "bayesian"){
  
  z = qnorm(1 - alpha/2)
  
  RR_covariate <- data.frame(name = covariate_level,
                               mean = rep(NA, length(covariate_level)),
                               conf_lb = rep(NA, length(covariate_level)),
                               conf_ub = rep(NA, length(covariate_level)),
                               mean_log = rep(NA, length(covariate_level)),
                               se = rep(NA, length(covariate_level)))
  
  for(covariate_level_name in covariate_level){
    if(covariate_level_name == covariate_ref){
      RR_covariate[RR_covariate$name == covariate_level_name,  2:6] <- tbl_coeff[tbl_coeff$name == exposure_name, c("mean", "conf_lb", "conf_ub", "mean_log", "se")]
      next()
    }
    
    effect_name = paste0(covariate_name, covariate_level_name)
    
    interaction_name = ifelse(stat == "bayesian",
                              paste0(exposure_name,"_", effect_name),
                              paste0(exposure_name, ":", effect_name))
    
    RR_covariate[RR_covariate$name == covariate_level_name, "mean"] <- exp(tbl_coeff[tbl_coeff$name == exposure_name, "mean_log"] + tbl_coeff[tbl_coeff$name == interaction_name, "mean_log"])
    SE <- sqrt(matrix_cov[effect_name, effect_name] + matrix_cov[exposure_name, exposure_name] +  matrix_cov[exposure_name, exposure_name])
    
    RR_covariate[RR_covariate$name == covariate_level_name, "conf_lb"] <- RR_covariate[RR_covariate$name == covariate_level_name, "mean"] * exp(- z * SE)
    RR_covariate[RR_covariate$name == covariate_level_name, "conf_ub"] <- RR_covariate[RR_covariate$name == covariate_level_name, "mean"] * exp(z * SE)
    RR_covariate[RR_covariate$name == covariate_level_name, "mean_log"] <- tbl_coeff[tbl_coeff$name == exposure_name, "mean_log"] + tbl_coeff[tbl_coeff$name == interaction_name, "mean_log"]
    RR_covariate[RR_covariate$name == covariate_level_name, "se"] <- SE
    }
  RR_covariate$name <- paste0(covariate_name, RR_covariate$name)
  return(RR_covariate)
}

### plot coefficient 
plot_coeff_fn <- function(tbl_coeff, outcome_name, effect_measure = "RR"){
effect_name = ifelse(effect_measure == "RR", "Rate", "Odds")

tbl_coeff$name = factor(tbl_coeff$name, levels = unique(tbl_coeff$name)) # ensure that order doesn't change in plot

ggplot(filter(tbl_coeff, name != "(Intercept)"), aes(x = name)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(aes(y = mean, 
                      ymin = conf_lb, 
                      ymax = conf_ub),
                  position = position_dodge(width = 0.4)) +
  coord_cartesian(ylim = c(0, 3)) +
  labs(x = "Covariate", 
       y = paste0(effect_name, 
                  " ratio for ", 
                  outcome_name, 
                  "\n during the mpox outbreak compared to the rest of 2022")) +
  scale_colour_viridis_d(option = "C", end = 0.8) + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1))
}


### plot coefficient comparing weighted vs unweighted
plot_compare_wt_fn <- function(tbl_coeff_unwt,
                          tbl_coeff_wt, 
                          outcome_name, 
                          effect_measure = "RR"){
  tbl_coeff <- rbind(mutate(tbl_coeff_unwt, weight = "Unweighted"),
                     mutate(tbl_coeff_wt, weight = "Weighted"))
  
  effect_name = ifelse(effect_measure == "RR", "Rate", "Odds")
  
  tbl_coeff$name = factor(tbl_coeff$name, levels = unique(tbl_coeff$name)) # ensure that order doesn't change in plot
  
  ggplot(filter(tbl_coeff, name != "(Intercept)"), aes(x = name, col = weight)) +
    geom_hline(yintercept = 1) +
    geom_pointrange(aes(y = mean, 
                        ymin = conf_lb, 
                        ymax = conf_ub,
                        col = weight),
                    position = position_dodge(width = 0.4)) +
    coord_cartesian(ylim = c(0, max(tbl_coeff$conf_ub))) +
    labs(x = "Covariate", 
         y = paste0(effect_name, 
                    " ratio for ", 
                    outcome_name, 
                    "\n during the mpox outbreak compared to the rest of 2022")) +
    scale_colour_viridis_d(option = "C", end = 0.8) + 
    theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1))
}
