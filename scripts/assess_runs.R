library(rstan)
library(magrittr)
library(ggplot2)
library(tidyr)
library(ggfan)
library(purrr)


library(dplyr)

source("R/data_functions.R")
source("R/assess.R")
source("R/reconstruct.R")

args <- commandArgs(trailingOnly = T)
run_stamp <- args[1]
results_base <- args[2]

if(is.na(run_stamp)){
  run_stamp <- tail(list.files("results"), 1)
}

if(is.na(results_base)){
  results_base <- "results"
}

result_path <- file.path("results", run_stamp)

rate_df <- readRDS("data/rate.rds")

res_file <- list.files(result_path, pattern="^mort_pacer")


result <- readRDS(file.path(result_path, res_file))

stan_data <- result$stan_data
countries <- stan_data$country_decode$Country
mort_fit <- result$mort_fit


#all countries 

assess_df <- map_df(countries, 
                    function(country) {
                      get_assessment_df(stan_data, mort_fit,
                                        country, rate_df)
                      }
                    )


saveRDS(assess_df, file.path(result_path, "assess.rds"))
