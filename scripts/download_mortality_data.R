library(HMDHFDplus)
library(curl)
library(XLConnect)
library(tibble)
library(magrittr)
library(dplyr)
library(purrr)
library(here)

HMD_username <- Sys.getenv("HFD_user")
HMD_password <- Sys.getenv("HFD_pass")

## download HMD data --------------------------------------------------------------------
#  single year of age and single year exposure data

countries <- c("AUS", "AUT", "BEL", "DNK", "FIN", "FRATNP", "DEUTW",
               "ESP", "IRL", "JPN", "NLD", "NZL_NM","NOR", "PRT",
               "SWE", "CHE", "GBRTENW", "GBR_SCO", "USA")

get_data <- function(cntry,HMD_username,HMD_password){
  exp_hmd <- readHMDweb(CNTRY = cntry, item = "Exposures_1x1", fixup = TRUE,
                        username = HMD_username,
                        password = HMD_password)
  
  # UK death data single year of age and single year exposure data
  deaths_hmd  <- readHMDweb(CNTRY = cntry, item = "Deaths_1x1", fixup = TRUE,
                            username = HMD_username,
                            password = HMD_password)
  rate_hmd <- left_join(exp_hmd %>% select(-Male, -Total) %>% 
                          rename(Exposure=Female),
                        deaths_hmd %>% select(-Male, -Total) %>% 
                          rename(Deaths=Female)) %>% 
    mutate(Rate=Deaths / Exposure,
           Log_rate=log(Rate),
           Country=cntry)
  return(rate_hmd)
}

rate_df <- map_df(countries, get_data,HMD_username,HMD_password)

dir.create(file.path(here(), "data"), recursive=T)
saveRDS(rate_df, file.path(here(), "data","rate.rds"))

get_e0 <- function(cntry,HMD_username,HMD_password){
  readHMDweb(CNTRY = cntry, item = "E0per", fixup = TRUE,
                        username = HMD_username,
                        password = HMD_password) %>%
    mutate(Country=cntry)
}


e0_df <- map_df(countries, get_e0, HMD_username, HMD_password) %>% as_tibble()
saveRDS(e0_df , file.path(here(), "data", "e0.rds"))

