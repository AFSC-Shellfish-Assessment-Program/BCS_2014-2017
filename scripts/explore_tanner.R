# notes ----
#Data exploration: BCS infection dynamics in C. bairdi 

# Author: Mike Litzow & Erin Fedewa
# last updated: 2022/9/7

# load ----
library(tidyverse)
library(lubridate)

#load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

# set plot theme
theme_set(theme_bw())

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

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
  print(n=100) #Only one crab sampled at  some stations

#Plot range of observed data by year/index site 
tanner.dat %>%
  group_by(year, index, station) %>%
  summarise(temperature = mean(temperature),
            depth = mean(depth),
            julian = mean(julian)) -> plot 

# BCS+/- occurrence by index/year
tanner.dat %>%
  select(year, index, pcr) %>%
  group_by(year, index) %>%
  summarise(PCR_0 = sum(pcr == 0),
            PCR_1 = sum(pcr == 1)) %>%
  pivot_longer(cols = c(-year, -index)) -> plot2

#% +/- stacked barplot
ggplot(plot2, aes(fill=name, y=value, x=year)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(~ index)
ggsave("./figs/tanner_index_year_incidence.png", width = 6, height = 6, units = "in")


###########################################
#Covariate data exploration

#Temperature
ggplot(plot, aes(temperature)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)
ggsave("./figs/tanner_bottom_temp_by_index.png", width = 5, height = 4, units = "in")

#Depth
ggplot(plot, aes(depth)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)
ggsave("./figs/tanner_depth_by_index.png", width = 5, height = 4, units = "in")

#Julian Day 
ggplot(plot, aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index) # big differences among index areas
ggsave("./figs/tanner_date_by_index.png", width = 5, height = 4, units = "in")

#Julian date vrs temperature 
ggplot(plot, aes(julian, temperature)) +
  geom_point() # julian day is no longer collinear with temperature once 2017 is included
ggsave("./figs/tanner_date_vs_bottom_temp.png", width = 6, height = 4, units = "in")


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
  geom_smooth(method = "gam")
ggsave("./figs/tanner_date_vs_percent_positive.png", width = 6, height = 4, units = "in")

#Mean size-at-station vrs %positive plot 
ggplot(plot3, aes(size, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam")
ggsave("./figs/tanner_size_vs_percent_positive.png", width = 6, height = 4, units = "in")

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
  geom_smooth(method = "gam") #Very similar to above plot 



















