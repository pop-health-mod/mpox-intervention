  #' @title fn_model
  #' @author: Fanyu Xiu, Carla Doyle
  #' @description Loads the model function 
  #' @param city model one city at a time: "mtl" = Montreal, "trt" = Toronto, "van" = Vancouver
  #' @param import_cases_city number of imported cases
  #' @param bbeta_city transmission probability per effective contact
  #' @param omega_city mixing parameters
  #' @param RR_H_city rate ratio of change in sexual partner numbers among GBM with > 7 sexual partners before 2022
  #' @param RR_L_city rate ratio of change in sexual partner numbers among GBM with <= 7 sexual partners before 2022
  #' @param gamma1_city 1/duration of effective infectious period
  #' @param period_city days from the first reported cases to the last modeled case
  #' @param TRACING 1: contact tracing and isolating GBM; 0: no contact tracing and isolation
  #' @param VACCINATING 1: vaccinating GBM with 1-dose Imvamune vaccines; 0: no vaccination
  #' @param cpp 1: use Rcpp; 0: use R
  #' @return output_t a list of cases (daily reported cases across time), time (days since first reported case), X (array of size in each compartment across time), inc (daily incidence across time)
  #' @rdname fn_model
  #' @export 
  fn_model <- function(city,
                       import_cases_city,
                       bbeta_city, 
                       omega_city, 
                       RR_H_city, 
                       RR_L_city,
                       gamma1_city, 
                       period_city = 150,
                       TRACING = 1,
                       VACCINATING = 1,
                       cpp 
                       ){ 
    if(is.null(S_matrix_prop[[city]])){
      stop(paste(city, "data wasn't loaded"))
    }
      
    ## make city-specific parameters
    # int
    days_imported_city <- days_imported[city]
    days_RR_city <- days_RR[city]
    end_city <- days_imported_city + period_city
    niter_city <- (end_city - sstart) / ddt + 1;
    
    # double
    report_frac_city <- report_frac[[city]]
    N_city <- N[[city]]
    
    # NumericVector
    N_s_city <- N_s[[city]]
    c_s_city <- c_s[[city]]
    N_a_city <- N_a[[city]]
    N_h_city <- N_h[[city]]
    N_ash_city <- N_ash[[city]]
    c_ash_city <- c_ash[[city]]
    psi_t_city <- psi_t[[city]]
    g_city <- g[[city]]
    upsilon_city <- upsilon[[city]]
    vartheta_city <- vartheta[[city]]
    
    # NumericMatrix
    A_matrix_city <- A_matrix[[city]]
    S_matrix_prop_city <- S_matrix_prop[[city]]
    
    # array
    init_prev_city <- init_prev[[city]]
    
    ## import cases as per fraction of only certain age groups (based on case data)
    and <- which(names_age_cats %in% c("30-39"))
    ## in the 5% highest activity group
    ind <- unname(which(cumsum(N_s_city / N_city) >= 0.95))
    ### import equally to I and E compartments
    for(comp in c("I", "E")){
      for(a in names_age_cats[and]){
        for(s in names_sa_cats[ind]){
          for(h in names_hiv_cats){
            ### divide by 2 since allocating the cases to E and I compartments equally
            init_prev_city[[comp]][a, s, h] <- import_cases_city / 2 * N_ash_city[a, s, h] / sum(N_ash_city[and, ind, ])
          }}}}
    ### S compartment initial population after case impotation
    init_prev_city[["S"]] <- N_ash_city - init_prev_city[["I"]] - init_prev_city[["E"]]
    
 if(cpp){
   # if code in cpp, convert all arrays into vector/matrix
   c_ash_city <- as.data.frame.table(c_ash_city)
   colnames(c_ash_city) <- c("age", "sa", "hiv", "c_ash")
   N_ash_city <- as.data.frame.table(N_ash_city)
   colnames(N_ash_city) <- c("age", "sa", "hiv", "N_ash")
   c_ash_city <- merge(c_ash_city, 
                       N_ash_city)
   for(comp in names_comp){
     init_prev_city[[comp]] <- as.data.frame.table(init_prev_city[[comp]])
     colnames(init_prev_city[[comp]]) <- c("age", "sa", "hiv", comp)
     c_ash_city <- merge(c_ash_city, init_prev_city[[comp]])
     init_prev_city[[comp]] <- NULL
   }
   init_prev_city <- c_ash_city
   c_ash_city <- init_prev_city$c_ash
   N_ash_city <- init_prev_city$N_ash
   init_prev_city <- as.matrix(init_prev_city[, names_comp])
   # import to global environment
   global_pars <- list(days_imported_city = days_imported_city,
                       days_RR_city = days_RR_city,
                       end_city = end_city,
                       niter_city = niter_city,
                       report_frac_city = report_frac_city,
                       N_city =  N_city,
                       N_s_city = N_s_city,
                       c_s_city = c_s_city,
                       N_a_city =  N_a_city,
                       N_h_city = N_h_city,
                       N_ash_city =  N_ash_city,
                       c_ash_city = c_ash_city,
                       psi_t_city =  psi_t_city,
                       g_city =  g_city,
                       upsilon_city = upsilon_city,
                       vartheta_city = vartheta_city,
                       A_matrix_city = A_matrix_city,
                       S_matrix_prop_city = S_matrix_prop_city,
                       init_prev_city = init_prev_city
     )
   invisible(lapply( 1:length( global_pars ), function( x )
       assign( as.character( names( global_pars )[ x ] ),
               global_pars[[ x ]], 
               envir = .GlobalEnv ) ) )
     # run model with rcpp using set of variable parameters
   output_t <- fn_model_cpp(bbeta_city = bbeta_city, 
                            omega_city = omega_city,
                            RR_H_city = RR_H_city,
                            RR_L_city = RR_L_city,
                            gamma1_city = gamma1_city,
                            TRACING = TRACING,
                            VACCINATING = VACCINATING)
   
   output_t$X <- array( output_t$X, 
                        dim = c( n_comp,
                                 niter_city,
                                 n_age_cats,
                                 n_sa_cats, 
                                 n_hiv_cats),
                        dimnames = list( names_comp,
                                         paste0("iter", 1:niter_city),
                                         names_age_cats, 
                                         names_sa_cats,
                                         names_hiv_cats))
      return(output_t)
 }else{
      ## time
      time <- seq(sstart,
                  end_city,
                  by = ddt)
      
      ## initial value for time-varying parameters
      ### city's number of people belonging to each sexual activity group at time = 0 (iteration = 1, a vector)
      N_t_s <- list()
      N_t_s[[1]] <- N_s_city
      ### city's probability of having contact to a sexual activity group at time = 0 (iteration = 1, a vector)
      g_t_s <- list()
      g_t_s[[1]] <- g_city
      ### city's matrix S at time = 0 (iteration = 1)
      S_matrix_t <- list()
      S_matrix_t[[1]] <- (1 - omega_city) * S_matrix_prop_city + omega_city * S_matrix_assor 
      
      ## output matrix dimension names
      strata_name = list(names_comp,
                       paste0("iter", 1:niter_city), 
                       names_age_cats, 
                       names_sa_cats,
                       names_hiv_cats)
      
      output_t <- rep(0, niter_city)
      inc_t = output_t
      
      ## assign initial values for each stratified compartment in output matrix
      X <- array(data = NA, 
               dim = c(n_comp, niter_city, n_age_cats, n_sa_cats, n_hiv_cats),
               dimnames = strata_name)
      for(comp in names_comp){
        for(a in names_age_cats){
          for(s in names_sa_cats){
            for(h in names_hiv_cats){
            ## recall that each compartment's initial pop is a list of array, and has dim n_age_cats *   n_sa_cats * n_hiv_cats
            X[comp, 1, a, s, h] <- init_prev_city[[comp]][a, s, h]
          }}}}
    
    # iterate model over time sequence (using Euler algorithm)
    for (t in 2:niter_city){
      ## extract % vaccinated people per time step depending on the day
      ## start from 0, corresponding to time_conti in PHAC case data
      day_index <- floor(ddt * (t - 1)) 
      max_day_data <- length(psi_t_city)
      psi_t <- ifelse(VACCINATING, 
                      ifelse(day_index >= max_day_data,
                             psi_t_city[max_day_data],
                             psi_t_city[day_index + 1]),
                      0)

      # extract % isolated people per time step depending on the day
      upsilon_t <- ifelse(TRACING,
                          upsilon_city[day_index + 1],
                          0)
      # compute the time-varying sexual activity mixing matrix
      N_t_s[[t]] <- vector()
      g_t_s[[t]] <- vector()

      for(s in names_sa_cats){
        ### number of participants who are mixing in each sexual activity group at iteration t (N_t_s, a vector)
        ### minus people in J since they are isolated so unable to mix with others
        N_t_s[[t]][s] <- N_t_s[[1]][s] - sum(X["J", t - 1, , s, ])
      }
      
      ### compute g at time t
      for(s in names_sa_cats){
        g_t_s[[t]][s] <- c_s_city[s] * N_t_s[[t]][s] / sum(c_s_city * N_t_s[[t]])}
      ### use g at time t to construct the proportionate mixing matrix at time t
      S_matrix_prop_t <- matrix(rep(g_t_s[[t]], each = n_sa_cats),
                                  nrow = n_sa_cats,
                                  dimnames = list(names_sa_cats, names_sa_cats))
     
      if(any(rowSums(S_matrix_prop_t) < 0.99) | any(rowSums(S_matrix_prop_t) > 1.01)){
        print("proportionate mixing matrix row sum not equal to 1")
        return(list(city = city,
                    iteration = t,
                    g_t_s,
                    N_t_s,
                    S_matrix_prop_t = S_matrix_prop_t))
        break}
      
      S_matrix_t[[t]] <- (1 - omega_city) * S_matrix_prop_t + omega_city * S_matrix_assor
      
      ### calculate lambda_t_ash through looping        
      for(a in names_age_cats){
          for(s in names_sa_cats){
            for(h in names_hiv_cats){
              
              summation_t_ash <- 0
              
              for(ap in names_age_cats){ 
                for(sp in names_sa_cats){
                  for(hp in names_hiv_cats){
                    #### probability part
                    p_t_ash_apsphp <- as.numeric(A_matrix_city[a, ap]) * as.numeric(S_matrix_t[[t]][s, sp]) * as.numeric(H_matrix[h, hp])
                    #### prevalence part
                    prev_t_ash_apsphp <- X["I", t - 1, ap, sp, hp] / (N_ash_city[ap, sp, hp] - X["J", t - 1, ap, sp, hp])         
                    summation_t_ash <- summation_t_ash + p_t_ash_apsphp * prev_t_ash_apsphp
                  }}}
      # if(t < 5){print(paste(t, a, s, h, summation_t_ash))}
      RR_t <- ifelse(c_ash_city[a, s, h] > 7 / 180 & day_index >= days_imported[city] & (day_index <= days_imported[city] + days_RR[city]),
                     RR_H_city, 
                     ifelse(day_index >= days_imported[city] & (day_index <= days_imported[city] + days_RR[city]),
                            RR_L_city,
                            1))
      lambda_t_ash <- c_ash_city[a, s, h] * bbeta_city * summation_t_ash * RR_t
      # disease natural history compartments
      X["S", t, a, s, h] <- X["S", t - 1, a, s, h] - ddt * (lambda_t_ash + psi_t * vartheta_city[a] / sum(X["S", t - 1, a, , ])) * X["S", t - 1, a, s, h] 
      X["V", t, a, s, h] <- X["V", t - 1, a, s, h] + ddt * (psi_t * vartheta_city[a] * X["S", t - 1, a, s, h] / sum(X["S", t - 1, a, , ]) - iota * lambda_t_ash * X["V", t - 1, a, s, h])
      X["E", t, a, s, h] <- X["E", t - 1, a, s, h] + ddt * (lambda_t_ash * X["S", t - 1, a, s, h] + iota * lambda_t_ash * X["V", t - 1, a, s, h] - alpha * X["E", t - 1, a, s, h])
      X["I", t, a, s, h] <- X["I", t - 1, a, s, h] + ddt * ((1 - upsilon_t) * alpha * X["E", t - 1, a, s, h] - gamma1_city * X["I", t - 1, a, s, h])
      X["J", t, a, s, h] <- X["J", t - 1, a, s, h] + ddt * (upsilon_t * alpha * X["E", t - 1, a, s, h] - gamma2 * X["J", t - 1, a, s, h])
      X["R", t, a, s, h] <- X["R", t - 1, a, s, h] + ddt * (gamma1_city * X["I", t - 1, a, s, h] + gamma2 * X["J", t - 1, a, s, h])

      # cases that become infectious (symptomatic or not)
      X["CS", t, a, s, h] <-  X["CS", t - 1, a, s, h] + ddt * (alpha * X["E", t - 1, a, s, h] - report_delay * X["CS", t - 1, a, s, h]) 
      # cases that are reported
      X["CR", t, a, s, h] <-  X["CR", t - 1, a, s, h] + ddt * (report_frac_city * report_delay * X["CS", t - 1, a, s, h]) 
      
      ## incidence is calculated as what entered into the E compartment
      inc_t[t] <- inc_t[t] + lambda_t_ash * X["S", t, a, s, h] + iota * lambda_t_ash * X["V", t, a, s, h]
      
      # debug for negative compartment size
      if (any(X[, t, a, s, h] < 0)){ 
        print("negative compartment")
        return(list(city = city,
                    iteration = t, 
                    age = a,
                    sex_act_grp = s,
                    hiv = h,
                    N_1_ash = sum(X[, 1, a, s, h]), #### initial number of people in the stratum ash
                    N_t_ash = sum(X[, t, a, s, h]), #### number of people in the stratum ash at iteration t
                    X_t_ash = X[, 1:t, a, s, h]))
        break}
      
      # debug for close population
      if (abs(sum(X[names_comp[1:6], t, a, s, h]) - sum(X[names_comp[1:6], 1, a, s, h])) > 0.00001) { 
        print("leaky stratum")
        return(list(city = city,
                    iteration = t, 
                    age = a,
                    sex_act_grp = s,
                    hiv = h,
                    N_1_ash = sum(X[, 1, a, s, h]),
                    N_t_ash = sum(X[, t, a, s, h]),
                    X_t_ash = X[, 1:t, a, s, h]))
        break}
      
      } ### a
      } ### s
      } ### h

    output_t[t] <- sum(report_frac_city * report_delay * X["CS", t, , , ])

    } ### for iteration
    
    output_t <- list(cases = output_t,
                       time = time,
                       X = X,
                       inc = inc_t)
      
    return(output_t)
      
 } ### for else (modeling with R)
} ### for the model function