library(tibble)
library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
library(rstan)
library(yaml)
library(git2r)

args <- commandArgs(trailingOnly = T)
config <- args[1]

if (is.na(config)){
  config <- "base_local"
  message("This is a test configuration and is unlikely to produce enough samples
          for reliable posterior inference")
}

config <- read_yaml(file.path("config",
                              paste0(config, ".yaml")))


model <- config$model


source("R/data_functions.R")
rate_df <- readRDS("data/rate.rds")

if (length(config$countries)>0){
  if(config$countries=="restricted"){
    rate_df %<>% filter(Country %in% c("GBRTENW", "AUT", 
                                       "FRATNP", "JPN",
                                       "NLD"))
  }
}

rate_df <- rate_df %>%
  filter(Year>=config$start_year, 
         Year<=config$joy, Age <=100, Age>0)

stan_data <- make_stan_data(rate_df, config)

stan_data$config <- config
time_model <- stan_model(file.path("stan",
                                   paste0(model, ".stan")))

mort_fit <- sampling(time_model, iter=config$iter,
                     chains=config$chains, 
                     cores=config$cores, 
                     thin=config$thin,
                     data=stan_data)
# save the git commit for replicability.
repo <- ifelse(git2r::in_repository(),repository(),NA)
result <- list(stan_data=stan_data,
               mort_fit=mort_fit,
               git_commit= ifelse(git2r::in_repository(),
                                  commits(repo, n=1),
                                  NA))

date_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- file.path("results", date_now)
dir.create(out_path, recursive = T)
res_file_name <- paste0(model, 
                        "_",
                        config$joy,
                        ".rds")
saveRDS(result, file.path(out_path, 
                          res_file_name))

if (file.exists("backup_copy_file.sh")){
  # copy result_file
  command <- paste0("backup_copy_file.sh ", out_path, 
                    " ", res_file_name)
  system(command)
}

