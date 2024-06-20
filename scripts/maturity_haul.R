# notes ----
#GOAL: Create "pcr_haul_master.csv" dataset for modeling/analyses
#1) Add a column for maturity in PCR_2014_2017.csv using clutch codes/chela heights
     #Distribution based cutlines for males published in 2019 tech memo
#2) Add in haul data (bottom temp, depth, lat/long, day sampled)
#3) Calculate CPUE by species at each station
    #_cpue: total tanner/snow CPUE by station 
    #70under_cpue: <70mm CW tanner/snow CPUE by station - portion of pop w/ highest BCS prevalence 
    #imm_cpue: Using Pam's BCS protocol target sizes for 2015-2017 collections
        #NOTE- separate target sizes for Tanner E and W specified in protocol so Tanner E sizes used ONLY 

# Author: Erin Fedewa
# last updated: 2022/1/20

# load ----
library(tidyverse)

#############################
#Append maturity 
pcr_master <- read.csv("./data/PCR_2014_2017.csv")

#Determine male maturity via distribution-based cutline method/clutch codes
pcr_master %>%
  mutate(maturity = case_when((Sex == 2 & Clutch > 0) ~ 1,
                              (Sex == 2 & Clutch == 0) ~ 0,
                              (grepl("opilio", Species_Name) & Sex == 1 & (log(Chela) < -2.20640 + 1.13523 * log(Size)))| (grepl("opilio", Species_Name) & Sex == 1 & Size < 50) ~ 0,
                              (grepl("bairdi", Species_Name) & Sex == 1 & (log(Chela) < -2.67411 + 1.18884 * log(Size)))| (grepl("bairdi", Species_Name) & Sex == 1 & Size < 60) ~ 0, 
                              (grepl("opilio", Species_Name) & Sex == 1 & (log(Chela) >= -2.20640 + 1.13523 * log(Size))) ~ 1,
                              (grepl("bairdi", Species_Name) & Sex == 1 & (log(Chela) >= -2.67411 + 1.18884  * log(Size))) ~ 1)) %>%
  select(-X) -> pcr_mat

#############################
#Append EBS haul data 
tanner_haul <- read.csv("./data/haul_bairdi.csv")
snow_haul <- read.csv("./data/haul_opilio.csv")

#Combine Tanner and snow haul files 
tanner_haul %>%
  bind_rows(snow_haul) %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(YEAR %in% c(2014, 2015, 2016, 2017),
         HAUL_TYPE==3) %>%
  select(VESSEL, CRUISE, START_DATE, HAUL, MID_LATITUDE, MID_LONGITUDE,GIS_STATION,
         BOTTOM_DEPTH,GEAR_TEMPERATURE) %>%
  distinct() ->tanner_snow

#Join haul and PCR datasets 
pcr_mat %>% 
  filter(Index_Site %in% c(1:6)) %>% #remove NBS index site for join 
  as_tibble() %>%
  left_join(tanner_snow, by = c("CRUISE", "VESSEL", "HAUL", "STATIONID"="GIS_STATION")) %>%
  rename_with(tolower) -> mat_haul

#################################
#Add in CPUE data for each station 
  #NOTE: There are minor differences between the tanner & snow haul files (e.g. digits in area swept) that will result in duplicate entries 
  #after summarizing if trying to merge haul files and then summarizing together, hence species by species approach here 

#tanner crab
tanner_haul %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(YEAR %in% c(2014, 2015, 2016, 2017),
         HAUL_TYPE==3) %>%
  rename_with(tolower) %>%
  group_by(year, gis_station, area_swept) %>% 
  summarise(tanner_cpue = sum(sampling_factor[species_code == 68560], na.rm = T) / mean(area_swept),
            tanner70under_cpue = sum(sampling_factor[species_code == 68560 & width  <= 70], na.rm = T) / mean(area_swept),
            tannerimm_cpue = sum(sampling_factor[species_code == 68560 & width  <= 112 & sex==1 |
                                                   species_code == 68560 & width <=84 & sex==2], na.rm = T) / mean(area_swept)) %>%
  right_join(mat_haul, by = c("year"="year", "gis_station"="stationid")) -> tanner_cpue

#snow crab
snow_haul %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(YEAR %in% c(2014, 2015, 2016, 2017),
         HAUL_TYPE==3) %>%
  rename_with(tolower) %>%
  group_by(year, gis_station, area_swept) %>% 
  summarise(snow_cpue = sum(sampling_factor[species_code == 68580], na.rm = T) / mean(area_swept),
            snow70under_cpue = sum(sampling_factor[species_code == 68580 & width  <= 70], na.rm = T) / mean(area_swept), 
            snowimm_cpue = sum(sampling_factor[species_code == 68580 & width  <= 94 & sex==1 |
                                                 species_code == 68580 & width <=50 & sex==2], na.rm = T) / mean(area_swept))                                                                                 %>%
  right_join(tanner_cpue, by = c("year", "gis_station")) %>%
  select(-area_swept.x, -area_swept.y) -> cpue

##################################################

#Write new master csv                               
write_csv(cpue, file="./data/pcr_haul_master.csv")
  
  

