// author: Fanyu Xiu
#include <Rcpp.h>
#include <boost/multi_array.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List fn_model_cpp(double bbeta_city, 
                            double omega_city,
                            double RR_H_city,
                            double RR_L_city,
                            double gamma1_city,
                            bool TRACING = 1,
                            bool VACCINATING = 1) {
  
  // ---- Set-Up ----
  
  Environment env = Environment::global_env();

  // read parameter values from global environment
  int sstart = env["sstart"];
  int n_comp = env["n_comp"];
  int n_age_cats = env["n_age_cats"];
  int n_sa_cats = env["n_sa_cats"];
  int n_hiv_cats = env["n_hiv_cats"];
  int niter_city = env["niter_city"];
  int days_imported_city = env["days_imported_city"];
  int days_RR_city = env["days_RR_city"];
  int N_city = env["N_city"];
  double report_frac_city = env["report_frac_city"];
  double ddt = env["ddt"];
  double iota = env["iota"];
  double alpha = env["alpha"];
  double gamma2 = env["gamma2"];
  double report_delay = env["report_delay"];
  Rcpp::NumericVector N_s_city = env["N_s_city"];
  Rcpp::NumericVector c_s_city = env["c_s_city"];
  Rcpp::NumericVector N_a_city = env["N_a_city"];
  Rcpp::NumericVector N_h_city = env["N_h_city"];
  Rcpp::NumericVector c_ash_city = env["c_ash_city"];
  Rcpp::NumericVector N_ash_city = env["N_ash_city"];
  Rcpp::NumericVector psi_t_city = env["psi_t_city"];
  Rcpp::NumericVector g_city = env["g_city"];
  Rcpp::NumericVector upsilon_city = env["upsilon_city"];
  Rcpp::NumericVector vartheta_city = env["vartheta_city"];
  Rcpp::NumericMatrix A_matrix_city = env["A_matrix_city"];
  Rcpp::NumericMatrix H_matrix = env["H_matrix"];
  Rcpp::NumericMatrix S_matrix_prop_city = env["S_matrix_prop_city"];
  Rcpp::NumericMatrix S_matrix_assor = env["S_matrix_assor"];
  Rcpp::NumericMatrix init_prev_city = env["init_prev_city"];

  // Indexing for compartments
  int S = 0;
  int V = 1;
  int E = 2;
  int I = 3;
  int J = 4;
  int R = 5;
  int CS = 6;
  int CR = 7;
  
  // number of iterations
  NumericVector time(niter_city);
  time[0] = sstart;

  // initial value for time-varying parameters
  Rcpp::NumericMatrix N_t_s(niter_city, n_sa_cats);
  N_t_s(0, _) = N_s_city;
  Rcpp::NumericMatrix g_t_s(niter_city, n_sa_cats);
  g_t_s(0, _) = g_city;
  List S_matrix_t(niter_city);
  NumericMatrix S_matrix_at0(n_sa_cats, n_sa_cats);
  for(int s = 0; s != n_sa_cats; ++ s){
    for(int sp = 0; sp != n_sa_cats; ++ sp){
      S_matrix_at0(s, sp) = (1 - omega_city) * S_matrix_prop_city(s, sp) + omega_city * S_matrix_assor(s, sp); 
    }
  }
  S_matrix_t[0] = S_matrix_at0;
  
  // to store output data
  NumericVector cases_t(niter_city);
  cases_t[0] = 0;
  NumericVector inc_t(niter_city);
  inc_t[0] = 0;
  
  // assign initial values for each stratified compartment in output matrix
  typedef boost::multi_array<double, 5> array_type_5d;
  array_type_5d X(boost::extents[n_comp][niter_city][n_age_cats][n_sa_cats][n_hiv_cats]); 
  for(int c = 0; c != n_comp; ++ c){ 
    for(int a = 0; a != n_age_cats; ++ a){
      for(int s = 0; s != n_sa_cats; ++ s){
        for(int h = 0; h != n_hiv_cats; ++ h){
          NumericVector comp_init = init_prev_city(_, c); // extract the compartment from the matrix
          X[c][0][a][s][h] = comp_init[a * n_sa_cats * n_hiv_cats + s * n_hiv_cats + h];
        }}}}
  
  
  // iterate model over time sequence (using Euler algorithm)
  for( int t = 1; t != niter_city; ++t ){
    
    // extract % vaccinated people per time step depending on the day
    // start from 0, corresponding to time_conti in PHAC case data
    int day_index = std::floor(ddt * t); 
    int max_day_data = psi_t_city.length();
    double psi_t;
    if(!VACCINATING){
      psi_t = 0;
    }else if(day_index >= max_day_data){
      psi_t = psi_t_city[max_day_data - 1];
    }else{
      psi_t = psi_t_city[day_index];
    }

    // extract % isolated people per time step depending on the day
    double upsilon_t;
    if(!TRACING){
      upsilon_t = 0;
    }else{
      upsilon_t = upsilon_city[day_index];
    }

    // compute the time-varying sexual activity mixing matrix
    for(int s = 0; s != n_sa_cats; ++s) {
    // number of participants who are mixing in each sexual activity group at iteration t (N_t_s, a vector)
    // minus people in J since they are isolated so unable to mix with others
      double J_t_s = 0;
      for(int a = 0; a != n_age_cats; ++a){
          for(int h = 0; h != n_hiv_cats; ++h){
            J_t_s += X[J][t - 1][a][s][h];
      }}
      N_t_s(t, s) = N_t_s(0, s) - J_t_s;
    }
    
    // compute g at time t
    for(int s = 0; s != n_sa_cats; ++s) {
      g_t_s(t, s) = c_s_city[s] * N_t_s(t, s) / sum(c_s_city * N_t_s(t, _));}
    /// use g at time t to construct the proportionate mixing matrix at time t
    NumericVector S_v_prop_t = rep_each(g_t_s(t, _), n_sa_cats);
    S_v_prop_t.attr("dim") = Dimension(n_sa_cats, n_sa_cats);
    Rcpp::NumericMatrix S_matrix_prop_t = as<NumericMatrix>(S_v_prop_t);
    
    NumericMatrix S_matrix_att(n_sa_cats, n_sa_cats);
    for(int s = 0; s != n_sa_cats; ++ s){
      for(int sp = 0; sp != n_sa_cats; ++ sp){
        S_matrix_att(s, sp) = (1 - omega_city) * S_matrix_prop_t(s, sp) + omega_city * S_matrix_assor(s, sp); 
      }
    }
    S_matrix_t[t] = S_matrix_att;
    
    // calculate lambda_t_ash and sum(X["S", t - 1, a, , ]) at t through looping
    for(int a = 0; a != n_age_cats; ++a){
      double S_a_tm1 = 0;
      
      for(int sp = 0; sp != n_sa_cats; ++sp){
        for(int hp = 0; hp != n_hiv_cats; ++hp){
          S_a_tm1 += X[S][t - 1][a][sp][hp];}}
      
      for(int s = 0; s != n_sa_cats; ++s){
        for(int h = 0; h != n_hiv_cats; ++h){
          
          double summation_t_ash = 0; 
          
          for(int ap = 0; ap != n_age_cats; ++ ap){
            for(int sp = 0; sp != n_sa_cats; ++ sp){
              for(int hp = 0; hp != n_hiv_cats; ++ hp){
                // probability part
                double p_t_ash_apsphp = A_matrix_city(a, ap) * S_matrix_att(s, sp) * H_matrix(h, hp);
                
                // prevalence part
                double prev_t_ash_apsphp = X[I][t - 1][ap][sp][hp] / (N_ash_city[ap * n_sa_cats * n_hiv_cats + sp * n_hiv_cats + hp] - X[J][t - 1][ap][sp][hp]);         
     
               summation_t_ash += p_t_ash_apsphp * prev_t_ash_apsphp;
               
              }}}

          // incorporate changing of contact rate during mpox (max 3 weeks)
          double RR_t;
          double c_ash = c_ash_city[a * n_sa_cats * n_hiv_cats + s * n_hiv_cats + h];
          if(c_ash > 7.0/180.0 && day_index >= days_imported_city && day_index <= days_imported_city + days_RR_city){ 
            RR_t = RR_H_city;
          }else if(day_index >= days_imported_city && day_index <= days_imported_city + days_RR_city){
            RR_t = RR_L_city;
          }else{
            RR_t = 1;
          }
          
          double lambda_t_ash = c_ash * bbeta_city * summation_t_ash * RR_t;
          // disease natural history compartments
          X[S][t][a][s][h] = X[S][t - 1][a][s][h] - ddt * (lambda_t_ash + psi_t * vartheta_city[a] / S_a_tm1) * X[S][t - 1][a][s][h];
          X[V][t][a][s][h] = X[V][t - 1][a][s][h] + ddt * (psi_t * vartheta_city[a] * X[S][t - 1][a][s][h] / S_a_tm1 - iota * lambda_t_ash * X[V][t - 1][a][s][h]);
          X[E][t][a][s][h] = X[E][t - 1][a][s][h] + ddt * (lambda_t_ash * X[S][t - 1][a][s][h] + iota * lambda_t_ash * X[V][t - 1][a][s][h] - alpha * X[E][t - 1][a][s][h]);
          X[I][t][a][s][h] = X[I][t - 1][a][s][h] + ddt * ((1 - upsilon_t) * alpha * X[E][t - 1][a][s][h] - gamma1_city * X[I][t - 1][a][s][h]);
          X[J][t][a][s][h] = X[J][t - 1][a][s][h] + ddt * (upsilon_t * alpha * X[E][t - 1][a][s][h] - gamma2 * X[J][t - 1][a][s][h]);
          X[R][t][a][s][h] = X[R][t - 1][a][s][h] + ddt * (gamma1_city * X[I][t - 1][a][s][h] + gamma2 * X[J][t - 1][a][s][h]);
          
          // cases that become infectious (symptomatic or not)
          X[CS][t][a][s][h] = X[CS][t - 1][a][s][h] + ddt * (alpha * X[E][t - 1][a][s][h] - report_delay * X[CS][t - 1][a][s][h]);
          // cases that are reported
          X[CR][t][a][s][h] = X[CR][t - 1][a][s][h] + ddt * (report_frac_city * report_delay * X[CS][t - 1][a][s][h]);
          cases_t[t] += report_frac_city * report_delay * X[CS][t][a][s][h]; // number of newly reported cases per day
          inc_t[t] += lambda_t_ash * X[S][t][a][s][h] + iota * lambda_t_ash * X[V][t][a][s][h]; // number of newly acquired incidence per day
          time[t] = time[t - 1] + ddt;
  }}}}
  
  // convert each array to a vector for return
  std::vector<double> X_vec;
  X_vec.reserve(n_comp * niter_city * n_age_cats * n_sa_cats * n_hiv_cats);
  for(int h = 0; h != n_hiv_cats; ++h){
    for(int s = 0; s != n_sa_cats; ++s){
      for(int a = 0; a != n_age_cats; ++a){
        for (int t = 0; t != niter_city; ++t) {
          for (int c = 0; c != n_comp; ++c) {
            X_vec.push_back(X[c][t][a][s][h]);
          }}}}}
  
  // store all vectors of results in a list
  List output_t = List::create(Named("cases") = cases_t,
                               _["time"] = time,
                               _["X"] = X_vec,
                               _["inc"] = inc_t);

  return(output_t);
  
}