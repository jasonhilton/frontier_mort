data {
  int n_ages; // number of ages in data
  int N;
  int n_countries;
  int n_years;
  int n_forecast_years;

  int n_period_constraint;

  int n_basis;
  int n_basis_re;
  
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
  matrix[n_years + n_forecast_years - n_period_constraint, n_years + n_forecast_years - n_period_constraint] L_k;

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
  vector[n_basis] beta_mu;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_mu_null;


  vector[n_basis] beta_improvement;
  real<lower=0> sigma_improvement;
  real<lower=0> sigma_improvement_null;
  
  
  matrix<lower=0>[n_basis_re, n_countries] beta_mu_countries;
  real<lower=0> sigma_mu_countries;
  real<lower=0> sigma_mu_null_countries;
  

  matrix[n_basis_re, n_countries] beta_imp_countries;
  real<lower=0> sigma_improvement_countries;
  real<lower=0> sigma_improvement_countries_null;


  matrix[n_basis_re, n_countries] beta_quad_countries;
  real<lower=0> sigma_quad_countries;
  real<lower=0> sigma_quad_countries_null;
  
  
  vector<lower=0>[n_basis_re] sigma_mu_re;
  vector<lower=0>[n_basis_re] sigma_improvement_re;
  
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

  matrix[n_basis, n_basis] full_penalty_mu;
  matrix[n_basis, n_basis] full_penalty_improvement;
  matrix[n_basis, n_basis] full_penalty_quad;
  matrix[n_basis_re, n_basis_re] full_penalty_mu_re;
  matrix[n_basis_re, n_basis_re] full_penalty_imp_re;
  matrix[n_basis_re, n_basis_re] full_penalty_quad_re;
  
  // frontier 
  vector[n_ages] mu;
  vector[n_ages] improvement;

  // country specific
  matrix[n_ages, n_countries] mu_countries;
  matrix[n_ages, n_countries] improvement_countries;
  matrix[n_ages, n_countries] quad_countries;
  matrix[n_countries, n_years + n_forecast_years] kk_countries;
  matrix[n_years + n_forecast_years - n_period_constraint ,n_countries] dkk_countries;

  // final likelihood
  vector[N] log_rate;


  // priors --------------------------------------------------------------------
  // construct penalties

  full_penalty_mu = ((1 / pow(sigma_mu,2)) * penalty_age + 
                     (1 / pow(sigma_mu_null, 2)) * penalty_age_null);

  full_penalty_improvement = ((1 / pow(sigma_improvement,2)) * penalty_age + 
                     (1 / pow(sigma_improvement_null, 2)) * penalty_age_null);

 
  full_penalty_mu_re = ((1 / pow(sigma_mu_countries,2)) * penalty_age_re + 
                     (1 / pow(sigma_mu_null_countries, 2)) * penalty_age_re_null);

  full_penalty_imp_re = ((1 / pow(sigma_improvement_countries,2)) * penalty_age_re + 
                     (1 / pow(sigma_improvement_countries_null, 2)) * penalty_age_re_null);

  full_penalty_quad_re = ((1 / pow(sigma_quad_countries,2)) * penalty_age_re + 
                     (1 / pow(sigma_quad_countries_null, 2)) * penalty_age_re_null);

  // priors on frontier parameters
  beta_mu ~ multi_normal_prec(zeroes[1:n_basis], full_penalty_mu);

  beta_improvement ~ multi_normal_prec(zeroes[1:n_basis], full_penalty_improvement);

  // priors on country-specific parameters
  for (c in 1:n_countries){
    beta_mu_countries[1:n_basis_re, c] ~ multi_normal_prec(zeroes[1:n_basis_re], 
                                                           full_penalty_mu_re);  
    beta_imp_countries[1:n_basis_re, c] ~ multi_normal_prec(zeroes[1:n_basis_re], 
                                                           full_penalty_imp_re);
    beta_quad_countries[1:n_basis_re, c] ~ multi_normal_prec(zeroes[1:n_basis_re], 
                                                           full_penalty_quad_re);
  }
  
  // borrowing strength across countries
  for (i in 1:n_basis_re){
    beta_mu_countries[i] ~ normal(0, sigma_mu_re[i]);
    beta_imp_countries[i] ~ normal(0, sigma_improvement_re[i]);
  }

  // priors on variance parameters
  sigma_mu ~ normal(0, 5);
  sigma_mu_null ~ normal(0, 5);
  sigma_mu_countries ~ normal(0, 5);
  sigma_mu_null_countries ~ normal(0, 5);
  sigma_improvement ~ normal(0, 5);
  sigma_improvement_null ~ normal(0, 5);
  sigma_improvement_countries ~ normal(0, 5);
  sigma_quad_countries ~ normal(0, 5);
  sigma_quad_countries_null ~ normal(0, 5);
  sigma_improvement_countries_null ~ normal(0, 5);
  sigma_mu_re ~ exponential(0.2);
  sigma_improvement_re ~ exponential(0.2);
  sigma_kk ~ normal(0,3);

  // priors on period effects
  for (c in 1:n_countries){
    dkk_countries_z[c] ~ normal(0,1);
  }
  
  // data model ----------------------------------------------------------------
  
  // frontier smooths
  mu = age_basis * beta_mu;
  improvement = age_basis * beta_improvement * 0.1;

  // country-specific smooths
  mu_countries = age_basis_re * beta_mu_countries;
  improvement_countries = age_basis_re * beta_imp_countries;
  quad_countries = age_basis_re * beta_quad_countries;

  // Period effects
  for (i in 1:n_countries){
    dkk_countries[1:(n_years + n_forecast_years - n_period_constraint), i] = (L_k * 
      (dkk_countries_z[i] * 0.1 * sigma_kk[i])');
    kk_countries[i] = cumulative_sum(inv_constraint * 
      dkk_countries[1:(n_years + n_forecast_years - n_period_constraint), i])';
  }

  // combine for log rates
  for (i in 1:N){
    log_rate[i] = (mu[age_index[i]] + 
                   improvement[age_index[i]] * time_std[time_index[i]] +
                   mu_countries[age_index[i], country_index[i]] * 
                    exp(improvement_countries[age_index[i], country_index[i]] * 
                        time_std[time_index[i]] +
                        quad_countries[age_index[i], country_index[i]] * 
                        pow(time_std[time_index[i]], 2)) +
                   kk_countries[country_index[i], time_index[i]]
                   ); 
  }

  // Likelihood
  deaths ~ neg_binomial_2_log(log_expos + log_rate, dispersion);
}
