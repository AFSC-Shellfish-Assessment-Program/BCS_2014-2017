# notes ----
#Data exploration: BCS infection dynamics in C. bairdi 

# Author: Mike Litzow & Erin Fedewa
# last updated: 2022/11/7

# load ----
library(tidyverse)
library(lubridate)

#load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

###########################################################
# data wrangling
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(species_name == "Chionoecetes bairdi",
         index_site %in% c(1, 2, 3),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth, 
         gear_temperature, tanner70under_cpue, tannerimm_cpue) %>%
  rename(pcr = pcr_result,
                station = gis_station,
                longitude = mid_longitude,
                depth = bottom_depth,
                temperature = gear_temperature,
                index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station),
         depth = as.numeric(depth)) -> tanner.dat 

#examine <70mm Tanner cpue distribution to decide on transformation                                    
ggplot(tanner.dat, aes(tanner70under_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

  #fourth root
ggplot(tanner.dat, aes(tanner70under_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

  #log
ggplot(tanner.dat, aes(log(tanner70under_cpue))) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#examine protocol target size Tanner cpue distribution to decide on transf.
ggplot(tanner.dat, aes(tannerimm_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(tanner.dat, aes(tannerimm_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(tanner.dat, aes(log(tannerimm_cpue))) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#Transform both CPUE indices with 4th root
tanner.dat %>%
  mutate(fourth.root.cpue70 = tanner70under_cpue^0.25,
         fouth.root.cpueimm = tannerimm_cpue^0.25) -> tanner.dat 

###############################################
#Data exploration

nrow(tanner.dat) # 1285 samples!

#Sample sizes
tanner.dat %>%
  group_by(year, index, station) %>%
  count() %>%
  print(n=100) #Only one crab sampled at  some stations so attributing
#variance to station level will be difficult 

# BCS+/- occurrence by index/year
tanner.dat %>%
  select(year, index, pcr) %>%
  group_by(year, index) %>%
  summarise(PCR_0 = sum(pcr == 0),
            PCR_1 = sum(pcr == 1)) %>%
  pivot_longer(cols = c(-year, -index)) -> plot2

#% +/- stacked barplot
ggplot(plot2, aes(fill=name, y=value, x=year)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ index)

#Size composition sampled by index site/yr
tanner.dat %>%
  group_by(year, index) %>%
  ggplot() +
    geom_histogram(aes(x=size), position = "stack") +
  facet_grid(year~index)

#Size-frequency distribution of uninfected vrs infected
tanner.dat %>%
  ggplot(aes(size, fill=as.factor(pcr), color=as.factor(pcr))) +  
  geom_histogram(position="identity",alpha=0.5) +
  theme_bw()

#Size-frequency distribution of maturity status- only those w/ chela measurements! 
dat %>%
  filter(species_name == "Chionoecetes bairdi",
         index_site %in% c(1, 2, 3),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         maturity %in% c(0,1),
         pcr_result %in% c(1, 0)) %>%
  ggplot(aes(size, fill=as.factor(maturity), color=as.factor(maturity))) +  
  geom_histogram(position="identity",alpha=0.5) +
  theme_bw()

#Percent prevalence by size bin
tanner.dat %>% 
  mutate(size_bin = cut(size, breaks=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140))) %>%
  group_by(size_bin) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) %>%
  filter(size_bin != "NA") -> size
ggplot(size, aes(as.factor(size_bin), Prevalance)) +
  geom_col()

#This is something to keep in mind when interpreting changes in prev. across
#site and year- in 2016 and 2017 prevalence was very low, but is likely due to 
#most samples being taken from mature males (despite protocol specifying imm).
#Without pulling them, it's difficult to interpret changes in prevalence as 
#a true measure of the disease prev, vrs. just a factor of the subset of population 
#being sampled. 

###########################################
#Covariate data exploration

#Plot range of observed data by year/index site 
tanner.dat %>%
  group_by(year, index, station) %>%
  summarise(temperature = mean(temperature),
            depth = mean(depth),
            julian = mean(julian)) -> plot 

#Temperature
ggplot(plot, aes(temperature)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

#Depth
ggplot(plot, aes(depth)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

#Julian Day 
ggplot(plot, aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index) #differences among index areas

#Julian date vrs temperature 
ggplot(plot, aes(julian, temperature)) +
  geom_point() 

#Plot explanatory variables as predictors of proportion BCS+ by year/station 
tanner.dat %>%
  group_by(year, station) %>%
  summarise(size = mean(size), 
            julian = mean(julian),
            temperature = mean(temperature),
            CPUE70 = mean(fourth.root.cpue70),
            CPUEimm = mean(fouth.root.cpueimm),
            proportion.positive = sum(pcr == 1) / n()) -> plot3

# Julian day vrs  %positive plot
ggplot(plot3, aes(julian, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") #Large seasonal effect!

#Mean size-at-station vrs %positive plot 
ggplot(plot3, aes(size, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") #Large size effect!

#Temp-at-station vrs %positive plot 
ggplot(plot3, aes(temperature, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam")

# <70mm CPUE-at-station vrs %positive plot 
ggplot(plot3, aes(CPUE70, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam")

# Immature CPUE-at-station vrs %positive plot 
ggplot(plot3, aes(CPUEimm, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") 
#Very similar to above plot- CPUE metrics likely very similar 



















