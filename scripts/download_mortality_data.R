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

#countries <- c("GBRTENW", "AUT", "NLD", "JPN", "FRATNP","ITA")

countries <- c("AUS", "AUT", "BEL", "DNK", "FIN", "FRATNP", "DEUTW",
               "ESP","ITA", "IRL", "JPN", "NLD", "NZL_NM","NOR", "PRT",
               "SWE", "CHE", "GBRTENW", "GBR_SCO", "USA")

get_data <- function(cntry){
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

rate_df <- map_df(countries, get_data)

dir.create(file.path(here(), "data"), recursive=T)
saveRDS(rate_df, file.path(here(), "data","rate.rds"))


