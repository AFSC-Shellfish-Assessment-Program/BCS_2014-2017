## Objective 2: Annual Hematodinium prevalence estimated via a binomial model
##NOTE: This approach was used instead of hurdle models in "prevalance_hurdle.R"
  #for manuscript 

#Author: M. Malick

library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

# load PCR data
dat = read.csv("./data/pcr_haul_master.csv")

## Snow crab -----------------------------------------------
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date
  filter(species_name == "Chionoecetes opilio",
         index_site %in% c(4, 5, 6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         size > 0,
         pcr_result %in% c(1, 0)) %>%
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_latitude, bottom_depth,
         gear_temperature, snow70under_cpue, snowimm_cpue) %>%
  rename(pcr = pcr_result,
         station = gis_station,
         latitude = mid_latitude,
         depth = bottom_depth,
         temperature = gear_temperature,
         index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station),
         fourth.root.cpue70 = snow70under_cpue^0.25,
         fouth.root.cpueimm = snowimm_cpue^0.25) %>%
  group_by(year, station, fourth.root.cpue70, julian, depth, sex, index) %>%
  summarise(n_pos = sum(pcr),
            n_total = n(),
            Prevalance = (sum(pcr)/n())*100,
            avg_size = mean(size)) -> prev.snow

range(prev.snow$n_pos)
range(prev.snow$n_total)
hist(prev.snow$Prevalance)
hist(prev.snow$n_pos)

prev.snow %>%
    group_by(year) %>%
    summarise(mean_prev = mean(Prevalance),
              median_prev = median(Prevalance),
              sum_prev = sum(Prevalance))

g = ggplot(prev.snow) +
    geom_histogram(aes(x = Prevalance), bins = 10) +
    facet_wrap( ~ year) +
    theme_bw()
print(g)

#Snow Crab Binomial prevalence model using proportions: n positive/n total
bin_snow = brm(n_pos | trials(n_total) ~ year, #each row i.e. station is a single trial
               data = prev.snow,
               family = binomial(link = "logit"),
               cores = 4, chains = 4, 
               warmup = 1500, iter = 6000,
               save_pars = save_pars(all = TRUE))
saveRDS(bin_snow, file = "./output/bin_snow.rds")

bin_snow = readRDS("./output/bin_snow.rds")
summary(bin_snow)
pp_check(bin_snow, type = "dens_overlay", ndraws = 100)
check_hmc_diagnostics(bin_snow$fit)
neff_lowest(bin_snow$fit)
rhat_highest(bin_snow$fit)
plot(bin_snow, ask = FALSE)
bayes_R2(bin_snow)

## Look at conditional effects
ce_snow = brms::conditional_effects(bin_snow, effect = "year",)
plot(ce_snow)
ce_snow$year

## Generate conditional effects by hand for comparison
nd = data.frame(year = 2015:2017, n_total = c(1, 1, 1))
epred_snow = posterior_epred(bin_snow, newdata = nd)
apply(epred_snow, 2, quantile, prob = c(0.025, 0.5, 0.975))

## Compare to mgcv::gam -- nearly identical!
mat_snow = cbind(prev.snow$n_pos, prev.snow$n_total - prev.snow$n_pos)
mgcv_snow = mgcv::gam(mat_snow ~ year, family = binomial(link = "logit"), data = prev.snow)
summary(mgcv_snow)
predict(mgcv_snow, newdata = nd, type = "response")

## Fit a binomial prevalence model using 0/1 PCR status (all data)
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date
  filter(species_name == "Chionoecetes opilio",
         index_site %in% c(4, 5, 6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         size > 0,
         pcr_result %in% c(1, 0)) %>%
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_latitude, bottom_depth,
         gear_temperature, snow70under_cpue, snowimm_cpue) %>%
  rename(pcr = pcr_result,
         station = gis_station,
         latitude = mid_latitude,
         depth = bottom_depth,
         temperature = gear_temperature,
         index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station),
         fourth.root.cpue70 = snow70under_cpue^0.25,
         fouth.root.cpueimm = snowimm_cpue^0.25) -> opilio.dat

bin_snow_pcr = brm(pcr ~ year,
                   data = opilio.dat,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, 
                   warmup = 1500, iter = 6000)

saveRDS(bin_snow_pcr, file = "./output/bin_snow_pcr.rds")
bin_snow_pcr = readRDS("./output/bin_snow_pcr.rds")

pp_check(bin_snow_pcr, type = "dens_overlay", ndraws = 100)
ce_snow_pcr = brms::conditional_effects(bin_snow_pcr, effect = "year")
ce_snow_pcr$year  ## almost identical to previous model!!


## Compare estimated and raw prevalence from each dataset
library(data.table)
prev.snow.dt = as.data.table(prev.snow)
prev1 = prev.snow.dt[ , .(prev = sum(n_pos) / sum(n_total)), by = .(year)]
prev1_est = ce_snow$year

opilio.dat.dt = as.data.table(opilio.dat)
prev2 = opilio.dat.dt[ , .(prev = sum(pcr) / .N), by = .(year)]
prev2_est = ce_snow_pcr$year

data.frame(year = 2015:2017,
           prev_prop_raw = prev1$prev,
           prev_prop_est = prev1_est$estimate__,
           prev_cnt_raw = prev2$prev,
           prev_cnt_est = prev2_est$estimate__)


## Tanner crab ---------------------------------------------
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date
  filter(species_name == "Chionoecetes bairdi",
         index_site %in% c(1, 2, 3),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  dplyr::select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth,
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
         depth = as.numeric(depth)) %>%
  # group_by(year, station) %>%
  group_by(year, station, julian, depth, sex, index) %>%
  summarise(n_pos = sum(pcr),
            n_total = n(),
            Prevalance = (sum(pcr)/n())*100,
            avg_size = mean(size)) -> prev.tanner


range(prev.tanner$n_pos)
range(prev.tanner$n_total)
hist(prev.tanner$Prevalance)
hist(prev.tanner$n_pos)

prev.tanner %>%
    group_by(year) %>%
    summarise(mean_prev = mean(Prevalance),
              median_prev = median(Prevalance),
              sum_prev = sum(Prevalance))

g = ggplot(prev.tanner) +
    geom_histogram(aes(x = Prevalance), bins = 10) +
    facet_wrap( ~ year) +
    theme_bw()
print(g)

#Tanner Crab Binomial prevalence model using proportions: n positive/n total
bin_tanner = brm(n_pos | trials(n_total) ~ year,
                 data = prev.tanner,
                 family = binomial(link = "logit"),
                 cores = 4, chains = 4, 
                 warmup = 1500, iter = 6000,
                 save_pars = save_pars(all = TRUE))
saveRDS(bin_tanner, file = "./output/bin_tanner.rds")

bin_tanner = readRDS("./output/bin_tanner.rds")
summary(bin_tanner)
pp_check(bin_tanner, type = "dens_overlay", ndraws = 100)
check_hmc_diagnostics(bin_tanner$fit)
neff_lowest(bin_tanner$fit)
rhat_highest(bin_tanner$fit)
plot(bin_tanner, ask = FALSE)
bayes_R2(bin_tanner)

## Look at conditional effects
ce_tanner = brms::conditional_effects(bin_tanner, effect = "year")
plot(ce_tanner)
ce_tanner$year

## Generate conditional effects by hand for comparison
nd = data.frame(year = 2015:2017, n_total = c(1, 1, 1))
epred_tanner = posterior_epred(bin_tanner, newdata = nd)
apply(epred_tanner, 2, quantile, prob = c(0.025, 0.5, 0.975))


## Fit a binomial prevalence model using 0/1 PCR status
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date
  filter(species_name == "Chionoecetes bairdi",
         index_site %in% c(1, 2, 3),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  dplyr::select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth,
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
         depth = as.numeric(depth),
         fourth.root.cpue70 = tanner70under_cpue^0.25, #Transformed CPUE of tanner <70mm CW
         fouth.root.cpueimm = tannerimm_cpue^0.25) -> tanner.dat #Transformed CPUE of immature tanner (Jensen protocol cutline)

bin_tanner_pcr = brm(pcr ~ year,
                     data = tanner.dat,
                     family = bernoulli(link = "logit"),
                     cores = 4, chains = 4, 
                     warmup = 1500, iter = 6000)

saveRDS(bin_tanner_pcr, file = "./output/bin_tanner_pcr.rds")
bin_tanner_pcr = readRDS("./output/bin_tanner_pcr.rds")

pp_check(bin_tanner_pcr, type = "dens_overlay", ndraws = 100)
ce_tanner_pcr = brms::conditional_effects(bin_tanner_pcr, effect = "year")
ce_tanner_pcr$year  ## almost identical to previous model!!

#Given similarities in methods, let's present the binomial models using indv level PCR
  #data (essentially the same approach as fitting our drivers model with 
  #only a year effect) vrs site level proportions

## #Combine snow and tanner plots for Fig 5 in ms  -------------------------------------------
df_snow = ce_snow$year
df_snow$species = "Snow crab"
df_tanner = ce_tanner$year
df_tanner$species = "Tanner crab"
df_year = rbind(df_snow, df_tanner)

dodge = position_dodge(width = 0.5) ## to offset datapoints on plot
new_colors = c("#2171b5", "#238b45")

g = ggplot(df_year) +
    geom_point(aes(year, estimate__ * 100, color = species), size = 3, position = dodge) +
    geom_errorbar(aes(year, ymin = lower__ * 100, ymax = upper__ * 100, color = species),
                  width = 0.3, size = 0.5, position = dodge) +
    ylab("Prevalance (%)") + xlab("") +
    scale_colour_manual(values = new_colors) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank()) +
    theme(legend.title = element_blank())
print(g)
ggsave("./figs/Fig5_annual_binomial.png", width = 4, height = 3, dpi = 300)
