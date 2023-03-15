#Exploratory analysis: calculate the ratio of small (<60mm) to large tanner/snow 
  #crab in BCS study years to explore as a mechanism for annual increase in BCS
  #prevalence given host size susceptibility 

# Author: Erin Fedewa
# last updated: 2022/3/20

# load ----
library(tidyverse)

#Append EBS haul data 
tanner_catch <- read.csv("./data/haul_tanner.csv")
sc_catch <- read.csv("./data/haul_opilio.csv")

#Append EBS strata data 
strata_sc <- read.csv("./data/strata_opilio.csv")

#############################
#snow crab first
sc_catch %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}")),
         size_cat = ifelse(WIDTH_1MM <= 60, "small", "large")) %>%
  filter(HAUL_TYPE == 3,
         !is.na(size_cat),
         YEAR %in% c(2014, 2015, 2016,2017)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, size_cat) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue_cnt = ncrab / AREA_SWEPT) %>%
  # join to hauls that didn't catch crab 
  right_join(sc_catch %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3,
                      YEAR %in% c(2014, 2015, 2016,2017)) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(CPUE = 0)) %>%
  #join to stratum
  left_join(strata_sc %>%
              select(STATION_ID, SURVEY_YEAR, STRATUM, TOTAL_AREA) %>%
              filter(SURVEY_YEAR%in% c(2014, 2015, 2016,2017)) %>%
              rename_all(~c("GIS_STATION", "YEAR",
                            "STRATUM", "TOTAL_AREA"))) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA, size_cat) %>%
  summarise(MEAN_CPUE = mean(cpue_cnt , na.rm = T),
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA))) %>%
  group_by(YEAR, size_cat) %>%
  #Sum across strata
  summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6) -> abundance

#Plot
ggplot(abundance, aes(y=ABUNDANCE_MIL, x=YEAR)) +
  geom_point() +
  geom_line()