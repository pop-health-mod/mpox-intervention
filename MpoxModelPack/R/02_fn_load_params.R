#' @title load.params.fn
#' @author: Fanyu Xiu, Carla Doyle
#' @description Loads the model's fixed parameter values
#' @param VE vaccine effectiveness from meta-analysis, Default: "0.5149"
#' @param contact_prop proportion of contact traced and isolated after first reported case, Default: 0.20
#' @param contact_dur days it takes for contact tracing, Default: 2
#' @param standardized_vaccine_date whether or not standardizing vaccination start dates across the cities (for sensitivity analysis)
#' @return parameters for the model will be added to the global environment
#' @details age and HIV mixing matrices from Milwid et al, 2022 (https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-022-07207-7#Sec18)
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  load.params.fn(); mix_odds_ah
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}
#' @rdname load.params.fn
#' @export 

load.params.fn <- function(VE = 0.5149,
                           contact_prop = 0.20,
                           contact_dur = 2,
                           standardized_vaccine_date = F,
                           prioritized_vaccine = F){
  pars <- list()
  
  # group names 
  pars$names_age_cats <- sort(as.character(unique(contact_rate_prop$age_cats)))
  pars$names_hiv_cats <- sort(as.character(unique(contact_rate_prop$hiv_cats)))
  pars$names_sa_cats <- sort(as.character(unique(contact_rate_prop$sa_cats)))
  
  # group number of levels 
  pars$n_age_cats <- length( pars$names_age_cats )
  pars$n_hiv_cats <- length( pars$names_hiv_cats )
  pars$n_sa_cats <- length( pars$names_sa_cats )
  
  # initial compartment prevalence 
  subgroup <- array(data = NA, 
                    dim = c(pars$n_age_cats, pars$n_sa_cats, pars$n_hiv_cats),
                    dimnames = list(pars$names_age_cats, pars$names_sa_cats, pars$names_hiv_cats))
  sublist <- list(S = subgroup,
                  V = subgroup,
                 E = subgroup,
                 I = subgroup,
                 J = subgroup,
                 R = subgroup,
                 CS = subgroup,
                 CR = subgroup)
  pars$names_comp <- names(sublist)
  pars$n_comp <- length(pars$names_comp)
  
  init_prev <- list(mtl = sublist, 
                   trt = sublist, 
                   van = sublist)
  
  CITIES <- c("mtl", "trt", "van")
  
  # report fraction of cases in each city
  pars$report_frac <- vector()
  pars$report_frac["mtl"] <- 0.82
  pars$report_frac["trt"] <- 0.86
  pars$report_frac["van"] <- 0.77
  
  pars$days_imported <- vector()
  pars$days_imported["mtl"] <- 21
  pars$days_imported["trt"] <- 21
  pars$days_imported["van"] <- 14
  
  # days between the first reported cases and mass PrEP vaccination
  pars$days_RR <- vector() 
  pars$days_RR["mtl"] <- as.numeric(as.Date("2022-08-14") - as.Date("2022-04-28")) # 06-14
  pars$days_RR["trt"] <- as.numeric(as.Date("2022-08-14") - as.Date("2022-05-13")) # 06-20 https://www.publichealthontario.ca/-/media/Documents/M/2022/mpx-immunization-post-ontario-cases.pdf?rev=d5914e5a55a249d59a27e31c09a49305&la=fr
  pars$days_RR["van"] <- as.numeric(as.Date("2022-08-14") - as.Date("2022-05-25")) # 07-11 http://www.bccdc.ca/Health-Info-Site/Documents/Monkeypox/Epidemiological_Summary/Mpox_Surveillance_20230109.pdf

  ## dt for the model
  pars$ddt <- 0.25
  pars$sstart <- 0
  
  # vaccination threshold in terms of partner numbers in each city
  pars$bar <- vector()
  ifelse(prioritized_vaccine,
         pars$bar <- c("mtl"=2, "trt"=3, "van"=3) / 180,
         pars$bar <- c("mtl"=0, "trt"=0, "van"=0) / 180)
  
  ### first dose vaccinations and cases data
  for(city in CITIES){

    if(is.null(pop_dat[[city]])){next} # if didn't select the city, then skip
    
    # case data (provincial)
    province <- ifelse(city == "mtl", "QC", ifelse(city == "trt", "ON", "BC"))
    pars$case_city[[city]] <- subset(case_data, prov == province)
    pars$case_city[[city]]$time_intro <- pars$case_city[[city]]$time_conti + pars$days_imported[city] # x days after the imported cases
    pars$case_city[[city]]$incidence <- pars$case_city[[city]]$city_cases # city
    
    ## format date 
    pars$psi_before_first_reported_case <- 0
    pars$case_city[[city]]$date <- as.Date(pars$case_city[[city]]$date)
    ifelse(standardized_vaccine_date,
           pars$psi_t[[city]] <- pars$case_city[[city]]$first_doses_standardized, 
           pars$psi_t[[city]] <- pars$case_city[[city]]$first_doses)
    
    # proportion of isolated on a specific day
    pars$upsilon[[city]] <- c(rep(0, pars$days_imported[city]),
                              rep(contact_prop * pars$report_frac[city] * pexp(q = 2 + contact_dur, rate = 1/5.1, lower.tail = F), 
                                  400))
    
    # add city MSM population (accounted for ceiling when dividing into strata)
    # and contact rate (used to estimate force of infection, lambda)
    pars$N[[city]] <- sum(pop_dat[[city]]$n_pop) 
    pars$c_ash[[city]] <- pop_dat[[city]]$c_ash
    pars$N_ash[[city]] <- pop_dat[[city]]$n_pop
    
    # calculate the initial prevalence
    not_S_comp_names <- pars$names_comp[pars$names_comp != "S"]
    for(comp in not_S_comp_names){
      init_prev[[city]][[comp]] <- array(data = 0, 
          dim = c(pars$n_age_cats, pars$n_sa_cats, pars$n_hiv_cats),
          dimnames = list(pars$names_age_cats, pars$names_sa_cats, pars$names_hiv_cats))}
    
    # S compartment initial population (before case importation)
    init_prev[[city]][["S"]] <- pop_dat[[city]]$n_pop # equals to the population at the beginning
    
    # age-hiv odds, unique to each city
    OR5 = matrix(0, pars$n_age_cats, pars$n_age_cats) # age
    OR5[upper.tri(OR5,TRUE)] = mix_odds[[city]][2:16]
    OR5 = OR5 + t(OR5) - diag(diag(OR5))
    OR = matrix(0, pars$n_age_cats * pars$n_hiv_cats, pars$n_age_cats * pars$n_hiv_cats) # age:hiv
    OR[1: 5, 1: 5] = OR5 + mix_odds[[city]][1]
    OR[6:10, 1: 5] = OR5
    OR[1: 5, 6:10] = OR5
    OR[6:10, 6:10] = OR5 + mix_odds[[city]][1]
    pars$mix_odds_ah[[city]] <- OR
    
    # check total number of people stay the same after partitioning into strata
    if(pars$N[[city]] != sum(pop_dat[[city]]$n_pop)){
      stop("Error: total number of people doesn't stay the same after partitioning into strata")
      }
        
    
    # proportion of vaccination by age across cities
    ifelse(city == "mtl",
           vartheta_age <- c(`16-29` = 0.215,
                             `30-39` = 0.265,
                             `40-49` = 0.2,
                             `50-59` = 0.16,
                             `60+` = 0.16),
           ifelse(city == "trt", 
                  vartheta_age <- c(`16-29` = 0.221,
                                    `30-39` = 0.334,
                                    `40-49` = 0.182,
                                    `50-59` = 0.1315,
                                    `60+` = 0.1315),
                  vartheta_age <- c(`16-29` = 0.218,
                                    `30-39` = 0.2995,
                                    `40-49` = 0.191,
                                    `50-59` = 0.14575,
                                    `60+` = 0.14575)
           ))
    
    ### if select prioritized vaccination, this variable will be assigned 0 if the contact rate of an ash group <= median in each city (bar[city])
    ### this is when assuming no change in RR
    ### when there is change in RR, proportion will be based on the variable vartheta_ash_RR_city as defined in the function file
    pars$vartheta_ash[[city]] <- array(rep(vartheta_age, pars$n_sa_cats * pars$n_age_cats * pars$n_hiv_cats), 
                                       dim = c(pars$n_age_cats, pars$n_sa_cats, pars$n_hiv_cats),
                                       dimnames = list(pars$names_age_cats, pars$names_sa_cats, pars$names_hiv_cats)
    )
    
   
      for(a in pars$names_age_cats){
        for(s in pars$names_sa_cats){
          for(h in pars$names_hiv_cats){
            if(pars$c_ash[[city]][a, s, h] <= pars$bar[city]){
              # only vaccinated to groups with partner numbers > bar
              pars$vartheta_ash[[city]][a, s, h] <- 0} 
            
          }}}
    
    
    
    
    
    
    
    } # for city
  
  # initial population size in each compartment 
  pars$init_prev <- init_prev
  
  # natural history parameters
  pars$iota <- 1 - VE # 1 - vaccine effectiveness 1 dose
  pars$alpha <- 1/5.1 # 1 / latent period
  pars$gamma2 <- 1/14 # rate of recovery among infected individuals who are quarantined/isolated
  
  # report delay
  pars$report_delay <- 1/2
  
  
  
  
  # assign all to the global environment
  invisible(lapply( 1:length( pars ), function( x )
    assign( as.character( names( pars )[ x ] ),
            pars[[ x ]], 
            envir = .GlobalEnv ) ))
}
