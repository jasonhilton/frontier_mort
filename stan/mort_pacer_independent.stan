data {
  int n_ages; // number of ages in data
  int N;
  int n_countries;
  int n_years;
  int n_forecast_years;

  int n_basis;
  int n_basis_re;
  int n_period_constraint;
  
  int<lower=0> age_index[N];
  int<lower=0> time_index[N];
  int<lower=0> country_index[N];


  int<lower=0> deaths[N];
  vector<lower=0>[N] expos;

  matrix[n_ages, n_basis] age_basis;
  matrix[n_ages, n_basis_re] age_basis_re;

  matrix[n_basis, n_basis] penalty_age;
  matrix[n_basis, n_basis] penalty_age_null;

  matrix[n_basis_re, n_basis_re] penalty_age_re;
  matrix[n_basis_re, n_basis_re] penalty_age_re_null;

  matrix[n_years + n_forecast_years, n_years + n_forecast_years - n_period_constraint] inv_constraint;
  matrix[n_years + n_forecast_years - n_period_constraint, n_years +n_forecast_years - n_period_constraint] L_k;

}

transformed data { 
  vector[N] log_expos;
  vector[n_basis] zeroes;
  vector[n_years] time_std;

  for (t in 1:n_years){
    time_std[t] = t;
  }
  time_std = (time_std - mean(time_std)) / sd(time_std);

  for (i in 1:n_basis){
    zeroes[i]=0;
  }
  log_expos = log(expos);
}

parameters { 
  
  matrix[n_basis_re, n_countries] beta_mu_countries;
  vector<lower=0>[n_countries] sigma_mu_countries;
  vector<lower=0>[n_countries] sigma_mu_null_countries;
  

  matrix[n_basis_re, n_countries] beta_imp_countries;
  vector<lower=0>[n_countries] sigma_improvement_countries;
  vector<lower=0>[n_countries] sigma_improvement_countries_null;
  
  matrix[n_countries, n_years + n_forecast_years - n_period_constraint] dkk_countries_z;
  vector<lower=0>[n_countries] sigma_kk;  

  real log_dispersion;
}

transformed parameters {
  
  real dispersion;
  dispersion = exp(log_dispersion);

}

model {

  // prior penalties
  matrix[n_basis_re, n_basis_re] full_penalty_mu_re;
  matrix[n_basis_re, n_basis_re] full_penalty_imp_re;
  


  // country specific
  matrix[n_ages, n_countries] mu_countries;
  matrix[n_ages, n_countries] improvement_countries;
  matrix[n_countries, n_years + n_forecast_years] kk_countries;
  matrix[n_years + n_forecast_years - n_period_constraint ,n_countries] dkk_countries;

  // final likelihood
  vector[N] log_rate;


  // priors --------------------------------------------------------------------
  // construct penalties

  // priors on country-specific parameters
  for (c in 1:n_countries){
    full_penalty_mu_re = ((1 / pow(sigma_mu_countries[c],2)) * penalty_age_re + 
                     (1 / pow(sigma_mu_null_countries[c], 2)) * penalty_age_re_null);

    full_penalty_imp_re = ((1 / pow(sigma_improvement_countries[c],2)) * penalty_age_re + 
                     (1 / pow(sigma_improvement_countries_null[c], 2)) * penalty_age_re_null);

    beta_mu_countries[1:n_basis_re, c] ~ multi_normal_prec(zeroes[1:n_basis_re], 
                                                           full_penalty_mu_re);  
    beta_imp_countries[1:n_basis_re, c] ~ multi_normal_prec(zeroes[1:n_basis_re], 
                                                           full_penalty_imp_re);
  }
  

  // priors on variance parameters

  sigma_mu_countries ~ normal(0, 5);
  sigma_mu_null_countries ~ normal(0, 5);
  sigma_improvement_countries ~ normal(0, 5);
  sigma_improvement_countries_null ~ normal(0, 5);

  sigma_kk ~ normal(0,3);

  // priors on period effects
  for (c in 1:n_countries){
    dkk_countries_z[c] ~ normal(0,1);
  }
  
  // data model ----------------------------------------------------------------
  

  // country-specific smooths
  mu_countries = age_basis_re * beta_mu_countries;
  improvement_countries = age_basis_re * beta_imp_countries;

  // Period effects
  for (i in 1:n_countries){
    dkk_countries[1:(n_years + n_forecast_years - n_period_constraint), i] = (L_k * 
      (dkk_countries_z[i] * 0.1 * sigma_kk[i])');
    kk_countries[i] = cumulative_sum(inv_constraint * 
      dkk_countries[1:(n_years + n_forecast_years - n_period_constraint), i])';
  }

  // combine for log rates
  for (i in 1:N){
    log_rate[i] = (mu_countries[age_index[i], country_index[i]]  +
                   improvement_countries[age_index[i], country_index[i]] *
                      time_std[time_index[i]] +
                   kk_countries[country_index[i], time_index[i]]); 
  }

  // Likelihood
  deaths ~ neg_binomial_2_log(log_expos + log_rate, dispersion);
}
