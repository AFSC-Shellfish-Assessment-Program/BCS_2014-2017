# notes ----
#Objective 2: Annual Hematodinium prevalence estimated via zero-inflated hurdle models 

#load
library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(marginaleffects)
library(emmeans)
library(MARSS)
library(corrplot)
library(tidybayes)  
library(broom)          
library(broom.mixed) 
source("./scripts/stan_utils.R")

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

##################################################
#Snow crab models 

#Data wrangling: station level prevalence 
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
  summarise(Prevalance = (sum(pcr)/n())*100,
            avg_size = mean(size)) -> prev.snow

#plot distribution of response
prev.snow %>%
  ggplot(aes(Prevalance)) +
  geom_histogram() #lots of zeros! 
#Hmmm..without the zeros, this distribution is not exponential, and we've got a good #
  #of stations with 100% prevalence. With no hurdle_gaussian() (or an even better choice)
  #family in brms though, we'll use a hurdled lognormal model. Not ideal, and may be 
  #worth exploring a custom brms family 

#Proportion of zeros in data
prev.snow %>%
  mutate(is_zero = Prevalance == 0) %>%
  ungroup() %>%
  count(is_zero) %>%
  mutate(prop = n/sum(n)) #30.1%

#Hurdle model with prevalence as response 
#Hurdle model incorporates information about both zeros and non-zeros by:
#1) Using a logistic regression model that predicts if prevalence is 0 or not (hurdle)
#2) use a lognormal model for outcomes that are not zero 
#i.e. this means that any predictions using the posterior distribution will reflect zero processes

#Model 1: Ignoring spatial variability, year effect only (Obj 2 analysis)
hurdle1_formula <- bf(
  #mu, mean part of formula
  Prevalance ~ year,
  #alpha, zero inflation part
  hu ~ year) 

hurdle1_snow <- brm(hurdle1_formula,
               data = prev.snow,
               family = hurdle_lognormal(),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(hurdle1_snow, file = "./output/hurdle1_snow.rds")
hurdle1_snow <- readRDS("./output/hurdle1_snow.rds")

tidy(hurdle1_snow)
pp_check(hurdle1_snow) 
#both zero and non zero processes incorporated into the posterior distribution

#MCMC convergence diagnostics 
check_hmc_diagnostics(hurdle1_snow$fit)
neff_lowest(hurdle1_snow$fit)
rhat_highest(hurdle1_snow$fit)
summary(hurdle1_snow)
bayes_R2(hurdle1_snow)

#Diagnostic Plots
plot(hurdle1_snow, ask = FALSE)
conditional_effects(hurdle1_snow) #default returns conditional means of both non-zero "mu" and zero "hu"
mcmc_plot(hurdle1_snow, prob = 0.95)
mcmc_plot(hurdle1_snow, transformations = "inv_logit_scaled")

#Conditional Effect for year 
conditional_effects(hurdle1_snow, effect = "year")

ce1s_1 <- conditional_effects(hurdle1_snow, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Snow crab") -> year_snow

#Average marginal effect of year 
years_ame <- hurdle1_snow %>% 
  emmeans(~ year,
          var = "year",
          epred = TRUE, re_formula = NA) %>% 
  gather_emmeans_draws()

years_ame %>%
  median_hdi()

ggplot(years_ame,aes(x = .value, fill=year)) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect",
       y = "Density") +
  theme_bw()

#Marginal effects: effect of year in the hurdling process 

#Combine hu_year term(s) with hurdle intercept and transform
hurdle_intercept <- tidy(hurdle1_snow) %>% 
  filter(term == "hu_(Intercept)") %>%  
  pull(estimate)

hurdle_lifeexp <- tidy(hurdle1_snow) %>%  
  filter(term == "hu_year2017") %>%  
  pull(estimate)

plogis(hurdle_intercept + hurdle_lifeexp) - plogis(hurdle_intercept)
#The probability of seeing 0% prevalence in 2017 decreased by 38% from 2015

##########################################
#Additional more complex snow crab models (just exploratory...)

#Model 2 using our opilio.final model structure for both zero and non-zero processes

hurdle2_formula <-  bf(
  #mu, mean part of formula
  Prevalance ~ s(avg_size, k = 4) + s(julian, k = 4) + + s(fourth.root.cpue70, k = 4) 
  + sex + s(depth, k = 4) + (1 | year/index),
  #alpha, zero inflation part
  hu ~ s(avg_size, k = 4) + s(julian, k = 4) + + s(fourth.root.cpue70, k = 4) 
  + sex + s(depth, k = 4) + (1 | year/index)) 


hurdle2 <- brm(hurdle2_formula,
               data = prev.dat,
               family = hurdle_lognormal(),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(hurdle2, file = "./output/hurdle2.rds")
hurdle2 <- readRDS("./output/hurdle2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(hurdle2$fit)
neff_lowest(hurdle2$fit)
rhat_highest(hurdle2$fit)
summary(hurdle2)
bayes_R2(hurdle2)

#Diagnostic Plots
plot(hurdle2, ask = FALSE)
plot(conditional_smooths(hurdle2), ask = FALSE)
conditional_effects(hurdle2)
mcmc_plot(hurdle2, type = "areas", prob = 0.95)
mcmc_rhat(rhat(hurdle2)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(hurdle2, pars = c("b_Intercept", "bs_ssize_1", "bs_sfourth.root.cpue70_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(hurdle2)) #Effective sample size: All ratios > 0.1

#Model #3 with year as a fixed effect to look at variation in prevalence across years 
hurdle3_formula <-  bf(
  #mu, mean part of formula
  Prevalance ~ s(avg_size, k = 4) + s(julian, k = 4) + + s(fourth.root.cpue70, k = 4) 
  + sex + s(depth, k = 4) + year + (1 | index),
  #alpha, zero inflation part
  hu ~ s(avg_size, k = 4) + s(julian, k = 4) + + s(fourth.root.cpue70, k = 4) 
  + sex + s(depth, k = 4) + year + (1 | index)) 


hurdle3 <- brm(hurdle3_formula,
               data = prev.dat,
               family = hurdle_lognormal(),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(hurdle3, file = "./output/hurdle3.rds")
hurdle3 <- readRDS("./output/hurdle3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(hurdle3$fit)
neff_lowest(hurdle3$fit)
rhat_highest(hurdle3$fit)
summary(hurdle3)
bayes_R2(hurdle3)

#Diagnostic Plots
plot(hurdle3, ask = FALSE)
plot(conditional_smooths(hurdle3), ask = FALSE)
conditional_effects(hurdle3)
mcmc_plot(hurdle3, prob = 0.95)
mcmc_plot(hurdle3, transformations = "inv_logit_scaled")

#Conditional Effect 
conditional_effects(hurdle3, effect = "year")

ce1s_1 <- conditional_effects(hurdle3, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Snow crab") -> year_snow2

#Average marginal effect of year 
years_ame <- hurdle3 %>% 
  emmeans(~ year,
          var = "year",
          epred = TRUE, re_formula = NA) %>% 
  gather_emmeans_draws()

years_ame %>%
  median_hdi()

ggplot(years_ame,aes(x = .value, fill=year)) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect",
       y = "Density") +
  theme_bw()

#Diagnostics and conditional effects plots look very similar to individual level models, as to be expected
#Note: this blog is great for script on extracting and interpreting coefficients from a hurdle model:
#https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/

#######################################################################
#Tanner Model

#data wrangling
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
  group_by(year, station, julian, depth, sex, index) %>%
  summarise(Prevalance = (sum(pcr)/n())*100,
            avg_size = mean(size)) -> prev.tanner

#plot distribution of response
prev.tanner %>%
  ggplot(aes(Prevalance)) +
  geom_histogram() #lots of zeros! 
#Without the zeros, this looks fairly normal. With no hurdle_gaussian()
  #family in brms though, we'll use a hurdled lognormal model again

#Proportion of zeros in data
prev.tanner %>%
  mutate(is_zero = Prevalance == 0) %>%
  ungroup() %>%
  count(is_zero) %>%
  mutate(prop = n/sum(n)) #47.2%

#Hurdle Model 1: Ignoring spatial variability, year effect only (Obj 2 analysis)
hurdle1_formula <- bf(
  #mu, mean part of formula
  Prevalance ~ year,
  #alpha, zero inflation part
  hu ~ year) 

hurdle1_tanner <- brm(hurdle1_formula,
                    data = prev.tanner,
                    family = hurdle_lognormal(),
                    cores = 4, chains = 4, iter = 2500,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(hurdle1_tanner, file = "./output/hurdle1_tanner.rds")
hurdle1_tanner <- readRDS("./output/hurdle1_tanner.rds")

tidy(hurdle1_tanner)
pp_check(hurdle1_tanner) 
#both zero and non zero processes incorporated into the posterior distribution

#MCMC convergence diagnostics 
check_hmc_diagnostics(hurdle1_tanner$fit)
neff_lowest(hurdle1_tanner$fit)
rhat_highest(hurdle1_tanner$fit)
summary(hurdle1_tanner)
bayes_R2(hurdle1_tanner)

#Diagnostic Plots
plot(hurdle1_tanner, ask = FALSE)
conditional_effects(hurdle1_tanner)
mcmc_plot(hurdle1_tanner, prob = 0.95)
mcmc_plot(hurdle1_tanner, transformations = "inv_logit_scaled")

#Conditional Effect for year 
conditional_effects(hurdle1_tanner, effect = "year")

ce1s_1 <- conditional_effects(hurdle1_tanner, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Tanner crab") -> year_tanner

######################################
#Combine snow and tanner plots for Fig 5 in ms 
dodge <- position_dodge(width=0.5) #to offset datapoints on plot 

new_colors <- c("#238b45","#2171b5")

year_tanner %>%
  full_join(year_snow) %>%
  #Combined conditional effect plot 
  ggplot() +
  geom_point(aes(year, estimate__, color=factor(species, 
                                                levels = c("Tanner crab", "Snow crab"))), size=3,
             position=dodge) +
  geom_errorbar(aes(year, ymin=lower__, ymax=upper__, color=factor(species, 
                                                                   levels = c("Tanner crab", "Snow crab"))), width=0.3, 
                size=0.5, position=dodge) +
  ylab("Prevalance (%)") + xlab("") +
  scale_colour_manual(values = new_colors) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.title= element_blank())
ggsave("./figs/annual_hurdle.png", dpi=300)
