#' 
#' Find the Root Mean Squared Error and bias of a posterior rate estimate.
#' 
#' @param emp_rate_df Data frame of empirical mortality rates by year and age 
#' for a specific country
#' @param log_rate_df Data frame of posterior samples or quantiles of rates for 
#' the same country
#' 
#' Takes the median of the posterior samples of rates at each age and year
#' (or the median if that dataframe contains only the quantiles of samples), and
#' calculates the error with respect to the equivalent empirical quantity. 
#' The RMSE and bias are calculated from these errors
#' 
#' @return A one-row data frame with the average RMSE and bias over all years 
#' and ages
get_rmse <- function(log_rate_df, emp_rate_df){
  if("Quantile" %in% names(log_rate_df)){
    log_rate_df %<>% filter((0.5 - Quantile)<1e-4)
  } else {
    log_rate_df %<>% group_by(Year, Age) %>%
      summarise(Rate= mean(exp(Log_rate)))
  }
  log_rate_df %<>% ungroup() %>%
    left_join(emp_rate_df %>% rename(Emp_Rate=Rate) %>% select(-Log_rate))
  
  log_rate_df %<>% filter(Exposure > 0) %>%
    mutate(Error = Rate - Emp_Rate) %>%
    summarise(RMSE=sqrt(mean(Error**2,na.rm=T)),
              Bias=mean(Error,na.rm=T))
  return(log_rate_df)
}


#' Is an empirical quantile within a particular interval
#' @param emp_quantiles Vector of empirical quantiles
#' @param interval A specified scalar interval in [0,1]
#' 
#' @return A vector of logicals indicating which observations are inside the 
#' intervals
get_indicator <- function(emp_quantiles, interval){
  quantiles <- c(0.50- rev(interval/2), 0.50 + interval/2)
  Indicator <- (emp_quantiles > quantiles[1]) & (emp_quantiles < quantiles[2])
  return(Indicator)
}


#' Find the empirical coverage of a set of posterior forecasts
#' @param stan_data The data used to produce a set of posterior rates in a list
#' @param nb_rate_array A 3-d array with posterior rate samples including 
#' negative binomial uncertainty. The dimensions are Age X Year X Sample
#' @param country String with the HMD code of the country for which coverage is 
#' required
#' @param forecast_joy Integer indicating the forecast jump-off year (joy). This
#' is defined as the final year included in the fitting process.
#' @param emp_rate_df Dataframe of empirical rates including log rates for the 
#' country for which coverage is being estimated
#' @param intervals Vector of probability intervals for which coverage is to be 
#' calculated.
#' 
#' @return A data from containing the proportion of observations with the 
#' specified empirical intervals
get_coverage <- function(stan_data, nb_rate_array, country, 
                         forecast_joy, emp_rate_df, 
                         intervals=c(0.5,0.9)){
  nb_log_rate_df <- get_log_rate_df_from_array(stan_data,
                                               log(nb_rate_array), 
                                               quantilise=T)
  nb_log_rate_df %<>% left_join(emp_rate_df %>% filter(Country==country) %>% 
                                  rename(Emp_log_rate=Log_rate))
  
  coverage_df <- nb_log_rate_df %>% 
    filter(Year > forecast_joy) %>%
    mutate(Indicator=Log_rate>Emp_log_rate) %>% 
    filter(Indicator==T) %>% 
    group_by(Year,Age) %>%
    summarise(Quantile=min(Quantile))

  
  coverage <- map(intervals, function(int, cov_df){
    intname <-paste0("I_", int * 100)
    cov_df %<>% mutate(!!intname:=get_indicator(Quantile, int))
    covname <- paste0("C_", int * 100)
    cov_df %<>% ungroup() %>%
      summarise(!!covname := sum(!!rlang::sym(intname))/n())
  }, cov_df=coverage_df)
  
  return(do.call(cbind,coverage)) 
}
  
  
#' Get the data frame with assessment metrics for posterior forecasts for 
#' a given country
#' @param stan_data The data used to produce a set of posterior rates in a list
#' @param mort_fit stan_fit object containing the fitted frontier model
#' @param country String with the HMD code of the country for which assessment
#' is required
#' @param emp_rate_df Dataframe of empirical rates including log rates  and 
#' exposures for the country for which assessment is required
#' @param period Logical. Does the model in question include a period effect?
#' 
#' @return A dataframe containing the relevant forecast assessment summaries for
#' the specified country
get_assessment_df <- function(stan_data, mort_fit, country, emp_rate_df,
                              period=T){
  
  log_rate <- get_country_log_rate_array(stan_data, mort_fit, 
                                         country,period)
  
  dispersion <- as.matrix(mort_fit, "dispersion")
  start_year <- stan_data$config$start_year
  forecast_joy <- stan_data$config$joy
  end_year <- forecast_joy + stan_data$config$n_forecast_years
  nb_rate <- get_nb_rate_array(log_rate,
                               emp_rate_df = emp_rate_df %>% filter(Country==country),
                               dispersion,
                               start_year,
                               end_year)
  
  coverage <- get_coverage(stan_data, nb_rate,
                           country, forecast_joy,
                           emp_rate_df)
  ## RMSE
  log_rate_df<- get_log_rate_df_from_array(stan_data, log_rate,
                                           quantilise=F) 
  rmse <- get_rmse(log_rate_df, emp_rate_df %>% filter(Country==country))
  assess_df <- cbind(rmse, coverage) %>% mutate(Country=country)
  return(assess_df)
}
  
  
  
  
  
  

