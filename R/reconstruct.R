

#' Get a dataframe of frontier mortality by year and age
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param quantilise Should the results be as quantiles or raw samples?
#' 
#' @return A dataframe of posterior estimates of frontier mortality by year and 
#' age, either raw samples or quantiles of these
get_frontier_mortality <- function(stan_data, mort_fit, quantilise=T){
  n_total_years <- stan_data$n_years + stan_data$n_forecast_years
  
  log_rate <- get_frontier_array(stan_data, mort_fit)
  
  log_rate_df <- get_log_rate_df_from_array(stan_data, log_rate, quantilise)
    
  return(log_rate_df)
}

# Get the Age x Year x sample array of frontier mortality
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' 
#' @return  A Age x Year x Sample array of posterior frontier mortality
get_frontier_array <- function(stan_data, mort_fit){
  n_total_years <- stan_data$n_years + stan_data$n_forecast_years
  
  mu_front <- as.matrix(mort_fit, "beta_mu")
  mu_front <- stan_data$age_basis %*% t(mu_front)
  
  imp_front <- as.matrix(mort_fit, "beta_improvement") * 0.1
  imp_front <- stan_data$age_basis %*% t(imp_front)
  
  n_total_years <- stan_data$n_years + stan_data$n_forecast_years
  time_ind <- 1:stan_data$n_years
  time_ind_fore <- 1:n_total_years
  time_std <- (time_ind_fore - mean(time_ind))/ sd(time_ind)
  
  imp_surf <- aperm(time_std %o% imp_front, c(2,1,3))
  
  mu_surf <- aperm(rep(1, n_total_years) %o%
                     mu_front, c(2,1,3))
  
  log_rate <- mu_surf + imp_surf
  return(log_rate)
}


#' Get the Age x Year x sample array of the deviations of a specific country
#' from the frontier (on the log scale)
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' 
#' @return Age x Year x Sample array of posterior smaples of country-specific
#' deviations
get_country_deviation_array <- function(stan_data, mort_fit, country){
  n_total_years <- stan_data$n_years + stan_data$n_forecast_years
  time_ind <- 1:stan_data$n_years
  time_ind_fore <- 1:n_total_years
  time_ind <- (time_ind_fore - mean(time_ind))/ sd(time_ind)
  
  mu_country <- get_country_specific_array(stan_data, mort_fit, country, 
                                           "beta_mu_countries")
  imp_country <- get_country_specific_array(stan_data, mort_fit, country,
                                            "beta_imp_countries")
  
  mu_surf <- aperm(rep(1, n_total_years) %o% mu_country, 
                   c(2,1,3))
  
  imp_surf <- aperm(time_ind %o% imp_country, 
                   c(2,1,3))
  
  if(stringi::stri_detect(stan_data$config$model, regex="quad")){
    quad_country <- get_country_specific_array(stan_data, mort_fit, country,
                                               "beta_quad_countries")
    quad_surf <- aperm((time_ind**2) %o% quad_country, 
                       c(2,1,3))
    imp_surf <- imp_surf + quad_surf
  }
  
  deviation <- mu_surf * exp(imp_surf)
  return(deviation)
}

# Get the array of a posterior estimates of a country specific parameter
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' @param par_name Name of the parameter to be extracted
#' 
#' @return Age x sample array of a posterior estimates of a 
#' country-specific smooth function
get_country_specific_array <- function(stan_data, mort_fit, country, 
                                       par_name){
  country_index <- stan_data$country_decode %>% 
    filter(Country==country) %>% select(country_f) %>%
    unlist()
  n_basis_re <- stan_data$n_basis_re
  country_subset <- (country_index - 1) * n_basis_re + 1:n_basis_re
  
  beta_mu_country <- as.matrix(mort_fit, par_name)[,country_subset]
  mu_country <- stan_data$age_basis_re %*% t(beta_mu_country)
  return(mu_country)
}

#' Get a dataframe containing the set of country specific deviations from the 
#' frontier at the intercept of the time index
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' 
#' @return A dataframe of posterior estimates of the smooth country-specific
#' deviation intercept
get_country_intercept <-function(stan_data, mort_fit, country){
  mu_country <- get_country_specific_array(stan_data, mort_fit, country,
                                           "beta_mu_countries")
  mu_country_df <- mu_country %>% apply(1, quantile, probs=1:99/100) %>%
    t() %>%
    as_tibble() %>%
    mutate(Age=1:100) %>%
    gather(Quantile, Deviation, -Age) %>%
    mutate(Quantile=as.numeric(gsub("%", "", Quantile))/100)
  return(mu_country_df)
}



#' Get the dataframe containing the linear part of the term describing
#' change in the country-specific deviations from the frontier, for a particular
#' country
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' 
#' @param country String with the HMD code of the country required
#' 
#' @return A dataframe of posterior estimates of the smooth country-specific
#' deviation improvements
get_country_improvement <-function(stan_data, mort_fit, country){
  imp_country <- get_country_specific_array(stan_data, mort_fit, 
                                            country, "beta_imp_countries")
  imp_country_df <- imp_country %>% apply(1, quantile, probs=1:99/100) %>%
    t() %>%
    as_tibble() %>%
    mutate(Age=1:100) %>%
    gather(Quantile, Improvement, -Age) %>%
    mutate(Quantile=as.numeric(gsub("%", "", Quantile))/100)
  return(imp_country_df)
}

#' Get the dataframe containing the quadratic part of the term describing
#' change in the country-specific deviations from the frontier, for a particular
#' country
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' 
#' @return A dataframe of posterior estimates of the smooth country-specific
#' deviation quadratic term
get_country_quadratic <-function(stan_data, mort_fit, country){
  imp_country <- get_country_specific_array(stan_data,mort_fit, 
                                            country, "beta_quad_countries")
  imp_country_df <- imp_country %>% apply(1, quantile, probs=1:99/100) %>%
    t() %>%
    as_tibble() %>%
    mutate(Age=1:100) %>%
    gather(Quantile, Quadratic, -Age) %>%
    mutate(Quantile=as.numeric(gsub("%", "", Quantile))/100)
  return(imp_country_df)
}


#' Get an Year X Age X Sample array with country-specific log-rates, net of 
#' period effect
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param frontier_array The array of posterior estimates of the frontier 
#' @param country String with the HMD code of the country required
#' 
#' @return The Year X Age X Sample array of posterior estimates of the  
#' country-specific log-rates, net of the period effect
get_country_central_array <- function(stan_data, mort_fit, frontier_array,
                                   country){

  deviation_array <- get_country_deviation_array(stan_data,mort_fit, country)
  country_rate_array <- frontier_array + deviation_array
  
  return(country_rate_array)
}

#' Get a Year X Age X Sample array with country-specific log-rates
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' @param period Does the mode include a period effect?
#' 
#' @return A Year X Age X Sample array with posterior samples of  
#' country-specific log-rates
get_country_log_rate_array  <- function(stan_data, mort_fit,
                                        country,
                                        period=T){
  
  model <- stan_data$config$model
  if(stringi::stri_detect(model, regex="independent")){
    log_rate <- get_log_rate_independent(stan_data,
                                         mort_fit, 
                                         country)
    return(log_rate)
  } 

  frontier_array <- get_frontier_array(stan_data, mort_fit)
  country_array <- get_country_central_array(stan_data, mort_fit, frontier_array,
                                          country)

  if (period==T){
    kk_country <- get_period_effect_matrix(stan_data, mort_fit, 
                                           country, diff=F)
    kk_country_array <- rep(1,stan_data$n_ages) %o% t(kk_country)
    log_rate <- country_array + kk_country_array
  }  else {
    log_rate <- country_array
  }

  return(log_rate)
  
}

#' Get data frame with posterior samples of log rates for a particular country
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' @param quantilise Should the results be as quantiles or raw samples?
#' @param period Does the model have a period effect
#' 
#' @return Dataframe countaining posterior estimates of country-specific 
#' log-rates
get_country_log_rate_df <- function(stan_data, mort_fit,
                                           country, quantilise=T,
                                           period=T){
  
  log_rate <- get_country_log_rate_array(stan_data, mort_fit,
                                         country,
                                         period=T)
  
  log_rate_df <- get_log_rate_df_from_array(stan_data, log_rate, quantilise)
  
  log_rate_df %<>% mutate(Country=country)
  return(log_rate_df)
}

#' Given an array of posterior log rates for a particular country, 
#' covert to a dataframe
#' 
#' @param stan_data The list of input data to the stan model
#' @param log_rate_array The Year x Age x Sample array of posterior samples of 
#' log rates
#' @param quantilise Should the results be as quantiles or raw samples?
#' 
#' @return A data frame containing either samples or quantiles of the log-rates 
#' in the input `log_rate_array` 
get_log_rate_df_from_array <- function(stan_data, log_rate_array, quantilise=T){
  if (quantilise==T){
    probs <- (1:99) / 100
    log_rate_array <- apply(log_rate_array, c(1,2), quantile, probs) %>% 
      aperm(c(2,3,1))
    col_name <- "Quantile"
    divisor <- 100
  } else {
    col_name <- "Sim"
    divisor <- 1
  }
  #n_total_years <- stan_data$n_years + stan_data$n_forecast_years
  n_total_years <- dim(log_rate_array)[2]
  start_year <- stan_data$config$start_year
  log_rate_df <- log_rate_array %>% matrix(ncol=dim(log_rate_array)[3]) %>%
    as_tibble() %>%
    mutate(Age=rep(1:stan_data$n_ages, n_total_years),
           Year=rep(start_year + 1:n_total_years - 1,
                    rep(stan_data$n_ages, n_total_years))) %>%
    gather(!!col_name, Log_rate, -Age, -Year) %>%
    mutate(!!col_name := as.numeric(gsub("V","", !!as.name((col_name)))) / 
             divisor)
  return(log_rate_df)
}


#' Constuct the period effect for a particular country.
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' @param diff Should the differenced period effects by calculated, or the 
#' natural ones?
#' 
#' @return The matrix of posterior estimates of country-specific period effects
get_period_effect_matrix <- function(stan_data, mort_fit, country,diff=T){
  n_total_years <- stan_data$n_years + stan_data$n_forecast_years
  
  
  n_period_constraint <- -diff(dim(stan_data$inv_constraint))
  
  dkk_z <- as.matrix(mort_fit, "dkk_countries_z")
  sigma_k <- as.matrix(mort_fit, "sigma_kk")
  country_index <- stan_data$country_decode %>% 
    filter(Country==country) %>% select(country_f) %>%
    unlist()
  n_countries <- dim(stan_data$country_decode)[1]
  period_index <-  seq(country_index, 
                       (n_total_years - n_period_constraint) * n_countries, 
                       n_countries)
  
  dkk_z <- dkk_z[,period_index]
  sigma_k <- sigma_k[,country_index]
  
  for(i in 1:dim(dkk_z)[1]){
    dkk_z[i,]  <- sigma_k[i] * dkk_z[i,] * 0.1
  }
  
  kk <- stan_data$L_k %*% t(dkk_z)
  
  kk <- stan_data$inv_constraint %*% kk
  
  if(diff==F){
    for(i in 1:dim(kk)[2]){
      kk[,i] <- cumsum(kk[,i])
    }
  }
  kk <- t(kk)
  
  return(kk)
}


#' Get the data frame with country-specific period effects
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' @param quantilise Should the results be as quantiles or raw samples?
#' 
#' @return A data frame containing posterior estimates of country specific 
#' period effects
get_period_effect <- function(stan_data, mort_fit, country, 
                              quantilise=T, diff=F){
  kk_country <- get_period_effect_matrix(stan_data, mort_fit, country,
                                         diff=diff)
  n_total_years <- dim(kk_country)[2]
  if (quantilise==T){
    probs <- (1:99) / 100
    kk_country <- apply(kk_country, 2, quantile, probs) %>% t()
    col_name <- "Quantile"
    divisor <- 100
  } else {
    col_name <- "Sim"
    divisor <- 1
  }
  
  if(diff==T){
    effect_name <- "Period_difference"
  } else {
    effect_name <- "Period_effect"
  }
  
  kk_country %<>% as_tibble() %>%
    mutate(Year=stan_data$config$start_year + 0:(n_total_years - 1)) %>%
    gather(!!col_name, !!effect_name, -Year) %>%
    mutate(!!col_name := as.numeric(gsub("V|%","", !!as.name((col_name)))) / 
             divisor)
  return(kk_country)
}


#' Get the data frame showing the posterior linear improvements in the mortality
#' frontier
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param quantilise Should the results be as quantiles or raw samples?
#' 
#' @return Data frame with posterior estimates of linear improvements in the 
#' frontier
#' 
get_frontier_improvement <- function(stan_data, mort_fit, quantilise=T){
  imp_front <- as.matrix(mort_fit, "beta_improvement") * 0.1
  imp_front <- stan_data$age_basis %*% t(imp_front)
  
  if (quantilise==T){
    probs <- (1:99) / 100
    imp_front <- apply(imp_front, 1, quantile, probs) %>% t()
    col_name <- "Quantile"
    divisor <- 100
  } else {
    col_name <- "Sim"
    divisor <- 1
  }
  
  imp_front %<>% as_tibble() %>% mutate(Age=1:100) %>%
    gather(!!col_name, Improvement, -Age) %>%
    mutate(!!col_name := as.numeric(gsub("V|%","", !!as.name((col_name)))) / 
             divisor)
  
  return(imp_front)
}

#' Get the data frame showing the posterior of the mortality
#' frontier at the intercept of the time index
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param quantilise Should the results be as quantiles or raw samples?
#' 
#' @return Data frame with posterior estimates of the frontier intercept
#' 
get_frontier_intercept <- function(stan_data, mort_fit, quantilise=T){
  imp_front <- as.matrix(mort_fit, "beta_mu")
  imp_front <- stan_data$age_basis %*% t(imp_front)
  
  if (quantilise==T){
    probs <- (1:99) / 100
    imp_front <- apply(imp_front, 1, quantile, probs) %>% t()
    col_name <- "Quantile"
    divisor <- 100
  } else {
    col_name <- "Sim"
    divisor <- 1
  }
  
  
  imp_front %<>% as_tibble() %>% mutate(Age=1:100) %>%
    gather(!!col_name, Intercept, -Age) %>%
    mutate(!!col_name := as.numeric(gsub("V|%","", !!as.name((col_name)))) / 
             divisor)
  
  return(imp_front)
}

#' Get the Year x Age x Sample array of log rates for the independence model
#' 
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param country String with the HMD code of the country required
#' @return Array of posterior estimates of log-rates for a particular country
#' 
#' 
get_log_rate_independent <- function(stan_data, mort_fit, country){
  n_total_years <- stan_data$n_years + stan_data$n_forecast_years
  time_ind <- 1:stan_data$n_years
  time_ind_fore <- 1:n_total_years
  time_ind <- (time_ind_fore - mean(time_ind))/ sd(time_ind)
  
  mu_country <- get_country_specific_array(stan_data, mort_fit, country, 
                                           "beta_mu_countries")
  imp_country <- get_country_specific_array(stan_data, mort_fit, country,
                                            "beta_imp_countries")
  
  mu_surf <- aperm(rep(1, n_total_years) %o% mu_country, 
                   c(2,1,3))
  
  imp_surf <- aperm(time_ind %o% imp_country, 
                    c(2,1,3))
  
  kk <- get_period_effect_matrix(stan_data, mort_fit, country,
                                 diff=F)
  kk_country_array <- rep(1,stan_data$n_ages) %o% t(kk)
  log_rate <- mu_surf + imp_surf + kk_country_array
  return(log_rate)
}


#' Get the data frame of log rates including negative binomial uncertainty
#' @param stan_data The list of input data to the stan model
#' @param mort_fit Object of class stan_fit with the fitted frontier model
#' @param emp_rate_df Empirical array for a specific country that includes
#' exposures
#' @param country String with the HMD code of the country required
#' @param quantilise Should the results be as quantiles or raw samples?
#' 
#' 
#' @return Data frame of posterior estimates of log rates for a specific country
#' , including negative binomial uncertainty
#' 
get_nb_log_rate_df <- function(stan_data, mort_fit, emp_rate_df, 
                               country, quantilise=T){
  log_rate <- get_country_log_rate_array(stan_data, mort_fit,
                                         country,period=T)
  dispersion <- as.matrix(mort_fit, "dispersion")
  start_year <- stan_data$config$start_year
  forecast_joy <- stan_data$config$joy
  emp_rate_df %<>% filter(Country==country)
  total_years <- stan_data$n_years + stan_data$config$n_forecast_years
  last_forecast_year <- start_year + total_years - 1
  last_data_year <- max(emp_rate_df$Year)
  end_year <- min(last_data_year, last_forecast_year)
  nb_rate <- get_nb_rate_array(log_rate,
                               emp_rate_df,
                               dispersion,
                               start_year,
                               end_year)
  
  log_rate_df<- get_log_rate_df_from_array(stan_data, log(nb_rate),
                                           quantilise=quantilise)
  return(log_rate_df)
}


#' Get the Year x Age x Sample array of log rates including negative 
#' binomial uncertainty for a specific country
#' @param log_rate_array Year x Age x Sample array of posterior log rates
#' @param emp_rate_df Empirical array for a specific country that includes
#' exposures
#' @param dispersion Matrix, with second dimension 1, with posterior samples
#' of the dispersion parameter
#' @param start_year Start of the period for which rates are to be calculated
#' @param end_year Start of the period for which rates are to be calculated
#' 
#' @return Year x Age x Sample array of posterior log rates including neg. bin.
#' uncertainty
#' 
get_nb_rate_array <- function(log_rate_array, emp_rate_df, dispersion,
                              start_year, end_year){
  exposure_mat <- emp_rate_df %>% select(Year, Age, Exposure) %>% 
    filter(Year >= start_year,Year <=end_year) %>%
    filter(Age>0, Age<=100) %>%
    spread(Year, Exposure) %>%
    select(-Age) %>%
    as.matrix()
  rate_array <- exp(log_rate_array)
  
  n_ages <- dim(log_rate_array)[1]
  n_years <- dim(exposure_mat)[2]
  n_iter <- dim(rate_array)[3]
  
  E_deaths <- vapply(1:n_iter, function(i){
     return(rate_array[,1:n_years,i] * exposure_mat)
  }, FUN.VALUE = exposure_mat)
  PP_rate <- vapply(1:n_iter, function(i){
    vapply(1:n_years, function(j){
      return(rnbinom(n_ages,
                     mu=E_deaths[,j,i],
                     size=dispersion[i])/ exposure_mat[,j])
    }, FUN.VALUE = exposure_mat[,1])
  }, FUN.VALUE = exposure_mat)
  return(PP_rate)
}