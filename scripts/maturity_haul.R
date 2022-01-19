# notes ----
#Add a column for maturity in PCR_2014_2017.csv using clutch codes/chela heights
  #Distribution based cutlines for males published in 2019 tech memo
#Add in haul data (bottom temp, depth, lat/long, day sampled)

# Author: Erin Fedewa
# last updated: 2022/1/18

# load ----
library(tidyverse)

#############################

#Append maturity 
pcr_master <- read.csv("./data/PCR_2014_2017.csv")
  head(pcr_master)
#Determine male maturity via distribution-based cutline method/clutch codes
pcr_master %>%
  rename(GIS_STATION=STATIONID) %>% 
  mutate(maturity = case_when((Sex == 2 & Clutch > 0) ~ 1,
                              (Sex == 2 & Clutch == 0) ~ 0,
                              (grepl("opilio", Species_Name) & Sex == 1 & (log(Chela) < -2.20640 + 1.13523 * log(Size)))| (grepl("opilio", Species_Name) & Sex == 1 & Size < 50) ~ 0,
                              (grepl("bairdi", Species_Name) & Sex == 1 & (log(Chela) < -2.67411 + 1.18884 * log(Size)))| (grepl("bairdi", Species_Name) & Sex == 1 & Size < 60) ~ 0, 
                              (grepl("opilio", Species_Name) & Sex == 1 & (log(Chela) > -2.20640 + 1.13523 * log(Size))) ~ 1,
                              (grepl("bairdi", Species_Name) & Sex == 1 & (log(Chela) > -2.67411 + 1.18884  * log(Size))) ~ 1)) -> pcr_mat

#############################

#Append haul data 
tanner_haul <- read.csv("./data/haul_tanner.csv")
  head(tanner_haul)
snow_haul <- read.csv("./data/haul_opilio.csv", skip=5)
  head(snow_haul)

tanner_haul %>%
  bind_rows(snow_haul) %>% 
  filter(AKFIN_SURVEY_YEAR %in% c(2014, 2015, 2016, 2017),
         HAUL_TYPE==3) %>%
  select(VESSEL, CRUISE, START_DATE, HAUL, MID_LATITUDE, MID_LONGITUDE,GIS_STATION,
         BOTTOM_DEPTH,GEAR_TEMPERATURE) %>%
  distinct() ->tanner_snow
  
glimpse(tanner_snow)
  
#Join haul and PCR datasets 
pcr_mat %>% 
  as_tibble() %>%
  left_join(tanner_snow, by = c("CRUISE", "VESSEL", "HAUL", "GIS_STATION")) %>%
  rename_with(tolower) -> mat_haul

#################################
  
#Write new master csv                               
 write_csv(mat_haul, file="./data/pcr_haul_master.csv")
  
  
  