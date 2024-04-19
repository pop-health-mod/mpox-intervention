// author: Fanyu Xiu
#include <Rcpp.h>
#include <boost/multi_array.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_mult_cpp(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2) {
  // function author: https://stackoverflow.com/a/19920823
  Rcpp::NumericVector multMatrix = m1 * m2;
  multMatrix.attr("dim") = Dimension(m1.nrow(), m1.ncol());
  Rcpp::NumericMatrix m = as<NumericMatrix>(multMatrix); 
  return m;
}

// [[Rcpp::export]]
NumericMatrix matrix_round_cpp(NumericMatrix mat) {
  // Get number of rows and columns in the matrix
  int nRows = mat.nrow();
  int nCols = mat.ncol();
  
  // Create a new matrix to store the rounded values
  NumericMatrix roundedMat(nRows, nCols);
  
  // Loop through each element of the matrix
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      // Round each element to two decimal places
      roundedMat(i, j) = round(mat(i, j) * 100.0) / 100.0;
    }
  }
  
  // Return the rounded matrix
  return roundedMat;
}

  
// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_exp_cpp(Rcpp::NumericMatrix m1) {
  Rcpp::NumericVector multMatrix = exp(m1);
  multMatrix.attr("dim") = Dimension(m1.nrow(), m1.ncol());
  Rcpp::NumericMatrix m = as<NumericMatrix>(multMatrix); 
  return m;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix outer_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  int n = x.size();
  int m = y.size();
  
  Rcpp::NumericMatrix result(n, m);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      result(i, j) = x[i] * y[j];
    }
  }
  
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix ipf_cpp(Rcpp::NumericMatrix M, Rcpp::NumericVector m1, Rcpp::NumericVector m2, double tol = 1e-12){
  // function author: Jesse Knight
  Rcpp::NumericVector r1, r2;
  int i, k;
  for (k = 0; k < 100; k++){
    r1 = m1 / rowSums(M);
    for (i = 0; i < M.nrow(); i++ ){ M(i,_) = M(i,_) * r1(i); }
    r2 = m2 / colSums(M);
    for (i = 0; i < M.ncol(); i++ ){ M(_,i) = M(_,i) * r2(i); }
    if (is_true(all( abs(r1-1) < tol & abs(r2-1) < tol))){ break; }
  }
  return M;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix apply_mix_odds_cpp(Rcpp::NumericMatrix M0,
                                   Rcpp::NumericMatrix OR,
                                   double tol = 1e-12){
  // function author: Jesse Knight, Fanyu Xiu
  M0 = M0 + tol; // fix NaN issues
  NumericMatrix M = ipf_cpp(matrix_mult_cpp(M0, matrix_exp_cpp(OR)), 
                            rowSums(M0), 
                            colSums(M0)); 
  return M;
  }

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
  double report_frac_city = env["report_frac_city"];
  double ddt = env["ddt"];
  double iota = env["iota"];
  double alpha = env["alpha"];
  double gamma2 = env["gamma2"];
  double report_delay = env["report_delay"];
  Rcpp::NumericVector RR_ash_city = env["RR_ash_city"];
  Rcpp::NumericVector c_ash_city = env["c_ash_city"];
  Rcpp::NumericVector N_ash_city = env["N_ash_city"];
  Rcpp::NumericVector psi_t_city = env["psi_t_city"];
  Rcpp::NumericVector upsilon_city = env["upsilon_city"];
  Rcpp::NumericVector vartheta_city = env["vartheta_city"];
  Rcpp::NumericMatrix mix_odds_ah_city = env["mix_odds_ah_city"];
  Rcpp::NumericMatrix mix_odds_s_city = env["mix_odds_s_city"];
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
  
  // mixing probability by {a,h,ap,hp} change with t
  typedef boost::multi_array<double, 4> array_type_4d;
  array_type_4d mix_ah4p(boost::extents[n_age_cats][n_hiv_cats][n_age_cats][n_hiv_cats]); 
  
  // mixing *contacts per-person* by {a,s,h,ap,sp,hp} change with t
  typedef boost::multi_array<double, 6> array_type_6d;
  array_type_6d mix_ash6c(boost::extents[n_age_cats][n_sa_cats][n_hiv_cats][n_age_cats][n_sa_cats][n_hiv_cats]); 
  
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
 
    // compute the time-varying mixing matrix
    
    /// number of participants who are mixing in each ash group at iteration t (n_ash_t, a 5*15*2 array)
    typedef boost::multi_array<double, 3> array_type_3d;
    array_type_3d n_ash_t(boost::extents[n_age_cats][n_sa_cats][n_hiv_cats]); 
      
    /// total contacts by {a,s,h}
    typedef boost::multi_array<double, 3> array_type_3d; 
    array_type_3d x_ash_t(boost::extents[n_age_cats][n_sa_cats][n_hiv_cats]); /// total contacts
    
    for(int a = 0; a != n_age_cats; ++ a){
      for(int s = 0; s != n_sa_cats; ++ s){
        for(int h = 0; h != n_hiv_cats; ++ h){
          /// number of participants who are mixing in each ash group at t
          /// minus people in J since they are isolated so unable to mix with others
          n_ash_t[a][s][h] = N_ash_city[a * n_sa_cats * n_hiv_cats + s * n_hiv_cats + h] - X[J][t - 1][a][s][h];
          
          /// contact rate at t for a,s,h
          double c_ash_city_t;
          c_ash_city_t = c_ash_city[a * n_sa_cats * n_hiv_cats + s * n_hiv_cats + h];
          if(day_index >= days_imported_city && day_index <= days_imported_city + days_RR_city){ 
              c_ash_city_t =  c_ash_city_t * RR_ash_city[a * n_sa_cats * n_hiv_cats + s * n_hiv_cats + h];
            }
          
          /// total contacts
            x_ash_t[a][s][h] = n_ash_t[a][s][h] * c_ash_city_t;
          }}}
      
    /// total contacts by {a,h} only
    NumericVector x_ah_t(n_age_cats * n_hiv_cats);
    for(int h = 0; h != n_hiv_cats; ++ h){
      for(int a = 0; a != n_age_cats; ++ a){
            x_ah_t[a + h * n_age_cats] = 0;
            for(int s = 0; s != n_sa_cats; ++ s){
              x_ah_t[a + h * n_age_cats] += x_ash_t[a][s][h];
            }}}
      
    // random mixing by {a,h}
    NumericMatrix mix_ah2_rand = outer_cpp(x_ah_t, x_ah_t) / sum(x_ah_t);
  
    // apply preferences by {a,h}
    NumericMatrix mix_ah2 = apply_mix_odds_cpp(mix_ah2_rand, mix_odds_ah_city);
    
    // make probability & reshape (10,10) -> (5,2,5,2)
    NumericVector row_sum_mix_ah2 = rowSums(mix_ah2);
    for(int a = 0; a != n_age_cats; ++ a){
      for(int h = 0; h != n_hiv_cats; ++ h){
        for(int ap = 0; ap != n_age_cats; ++ ap){
          for(int hp = 0; hp != n_hiv_cats; ++ hp){
            mix_ah4p[a][h][ap][hp] = mix_ah2(a + h * n_age_cats, ap + hp * n_age_cats) / row_sum_mix_ah2[a + h * n_age_cats];
          }}}}
  
      NumericVector x_s(n_sa_cats);
      NumericVector x_sp(n_sa_cats);
        
      for(int a = 0; a != n_age_cats; ++ a){
        for(int ap = 0; ap != n_age_cats; ++ ap){
          for(int h = 0; h != n_hiv_cats; ++ h){
            for(int hp = 0; hp != n_hiv_cats; ++ hp){
              
              for(int s = 0; s != n_sa_cats; ++ s){
                 // total contacts by {s} among {a,h}
                  x_s[s]  = x_ash_t[a][s][h] * mix_ah4p[a][h][ap][hp];
                 // total contacts by {sp} among {ap,hp}
                  x_sp[s] = x_ash_t[ap][s][hp] * mix_ah4p[ap][hp][a][h];
              }
              
              // random mixing by {s,sp}
              NumericMatrix mix_s2_rand = outer_cpp(x_s, x_sp) / (.5 * sum(x_s) + .5 * sum(x_sp));
              
              // contacts * mixing per-person
              NumericMatrix temp_mat = apply_mix_odds_cpp(mix_s2_rand, 
                                                          mix_odds_s_city);
              
              
              // if(t == 1 && a == 1 && h == 1 && ap == 1 && hp == 1 ){
              //   Rcout << "a = " << a << " ";
              //   Rcout << "h = " << h << " ";
              //   Rcout << "ap = " << ap << " ";
              //   Rcout << "hp = " << hp << "\n";
              //   // Rcout << mix_ah4p[a][h][ap][hp] << "\n";
              //   // Rcout << x_s << "\n";
              //   // Rcout << x_sp << "\n";
              //   // Rcout << n_ash_t[a][14][h]<< "\n";
              //   // Rcout << c_ash_city[a * n_sa_cats * n_hiv_cats + 14 * n_hiv_cats + h] << "\n";
              //   // Rcout << RR_ash_city[a * n_sa_cats * n_hiv_cats + 14 * n_hiv_cats + h] << "\n";
              //   // Rcout << x_ash_t[a][14][h]<< "\n";
              //   Rcout << temp_mat<< "\n";
              // }
              
              for(int s = 0; s != n_sa_cats; ++ s){
                for(int sp = 0; sp != n_sa_cats; ++ sp){
                  mix_ash6c[a][s][h][ap][sp][hp]  =  temp_mat(s, sp) / n_ash_t[a][s][h];
                }}
            
            
            
            }}}}
      
   
    for(int a = 0; a != n_age_cats; ++a){ 
      
      // sum(X["S", t - 1, a, , ]) at t through looping
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
               summation_t_ash += mix_ash6c[a][s][h][ap][sp][hp] * X[I][t - 1][ap][sp][hp] / n_ash_t[ap][sp][hp];
              }}}

          // calculate lambda_t_ash
          double lambda_t_ash = bbeta_city * summation_t_ash;
          
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