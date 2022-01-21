# explore and analyze infection dynamics in C. bairdi

library(tidyverse)
library(plyr)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


# set plot theme
theme_set(theme_bw())

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

# add julian date
dat$julian = lubridate::yday(lubridate::parse_date_time(x = dat$start_date, orders="mdy", tz="US/Alaska"))                                

# examine Snow crab cpue distribution                                    
ggplot(dat, aes(snow70under_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(dat, aes(snow70under_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

# separate Tanner data
opilio.dat <- dat %>%
  dplyr::filter(species_name == "Chionoecetes opilio",
                index_site %in% c(4, 5, 6),
                year %in% c(2015:2017),
                sex %in% c(1, 2),
                pcr_result %in% c(1, 0)) %>%
  dplyr::select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth, gear_temperature, snow70under_cpue) %>%
  dplyr::rename(pcr = pcr_result,
                station = gis_station,
                longitude = mid_longitude,
                depth = bottom_depth,
                temperature = gear_temperature,
                index = index_site,
                fourth.root.cpue70 = snow70under_cpue) %>%
  dplyr::mutate(year = as.factor(year),
                sex = as.factor(sex),
                index = as.factor(index),
                station = as.factor(station),
                fourth.root.cpue70 = fourth.root.cpue70^0.25) # transforming cpue here

nrow(opilio.dat) # 1511 samples!

# examine sample size by index / year

# now the same for index/year
plot <- opilio.dat %>%
  dplyr::select(year, index, pcr) %>%
  dplyr::group_by(year, index) %>%
  dplyr::summarise(PCR_0 = sum(pcr == 0),
                   PCR_1 = sum(pcr == 1)) %>%
  tidyr::pivot_longer(cols = c(-year, -index))

ggplot(plot, aes(name, value)) +
  geom_bar(stat = "identity") +
  facet_grid(year ~ index)

ggsave("./figs/opilio_index_year_incidence.png", width = 6, height = 6, units = "in")


plot <- opilio.dat %>%
  dplyr::group_by(year, index, station) %>%
  dplyr::summarise(temperature = mean(temperature),
                   depth = mean(depth),
                   julian = mean(julian)) 

ggplot(plot, aes(temperature)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

ggsave("./figs/opilio_bottom_temp_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(depth)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

ggsave("./figs/opilio_depth_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index) # big differences among index areas

ggsave("./figs/opilio_date_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(julian, temperature)) +
  geom_point() # julian day is no longer collinear with temperature once 2017 is included

ggsave("./figs/opilio_date_vs_bottom_temp.png", width = 6, height = 4, units = "in")


# look at Julian day as a predictor of % positive
plot <- opilio.dat %>%
  dplyr::group_by(year, station) %>%
  dplyr::summarise(julian = mean(julian), 
                   n = dplyr::n(),
                   proportion.positive = sum(pcr == 1) / dplyr::n())

ggplot(plot, aes(julian, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam", se = F)

ggsave("./figs/opilio_date_vs_percent_positive.png", width = 6, height = 4, units = "in")


# and size as a predictor of % positive
plot <- opilio.dat %>%
  dplyr::group_by(year, station) %>%
  dplyr::summarise(size = mean(size), 
                   n = dplyr::n(),
                   proportion.positive = sum(pcr == 1) / dplyr::n())

ggplot(plot, aes(size, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam")

ggsave("./figs/opilio_size_vs_percent_positive.png", width = 6, height = 4, units = "in")
