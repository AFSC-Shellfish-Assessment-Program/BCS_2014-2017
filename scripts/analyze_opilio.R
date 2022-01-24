# analyze infection dynamics in C. opilio

library(tidyverse)
library(plyr)
library(rstan)
library(brms)
library(bayesplot)
library(MARSS)
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

# separate opilio data
opilio.dat <- dat %>%
  dplyr::filter(species_name == "Chionoecetes opilio",
                index_site %in% c(4, 5, 6),
                year %in% c(2015:2017),
                sex %in% c(1, 2),
                pcr_result %in% c(1, 0)) %>%
  dplyr::select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_latitude, bottom_depth, gear_temperature, snow70under_cpue) %>%
  dplyr::rename(pcr = pcr_result,
                station = gis_station,
                latitude = mid_latitude,
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

# appear to be some NAs - check?
# looks like one row of NA in 	year == 2015, station == I-23

check <- opilio.dat$year == 2015 & opilio.dat$station == "I-23"
opilio.dat[check,]

# fix

opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","julian"] <-
  opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","julian"][1]

opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","latitude"] <-
  opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","latitude"][1]

opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","depth"] <-
  opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","depth"][1]

opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","temperature"] <-
  opilio.dat[opilio.dat$year == 2015 & opilio.dat$station == "I-23","temperature"][1]

# and one NA size!

apply(is.na(opilio.dat), 2, which)  

opilio.dat <- opilio.dat[-apply(is.na(opilio.dat), 2, which)$size,]

# need dimension reduction for exogenous covariates (day, depth, longitude, temperature, cpue)

pca.dat <- opilio.dat %>%
  dplyr::group_by(station, year) %>%
  dplyr::summarise(julian = mean(julian),
                   depth = mean(depth),
                   latitude = mean(latitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70))

cor(pca.dat[,3:7]) # some lower correlations than for tanner


# plot - for the paper?
plot <- data.frame(Day_of_year = pca.dat$julian,
                   Depth_m = pca.dat$depth,
                   N_latitude = pca.dat$latitude,
                   Bottom_temperature_C = pca.dat$temperature,
                   Fourth_root_CPUE_70 = pca.dat$fourth.root.cpue70,
                   year = as.numeric(as.character(pca.dat$year))) %>%
  tidyr::pivot_longer(cols = c(-Day_of_year, -year))

ggplot(plot, aes(Day_of_year, value)) +
  geom_point(color = "grey100") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, color = "black", lwd = 0.3) +
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  geom_point(aes(color = as.factor(year))) +
  scale_color_manual(values = cb[c(2,4,6)]) +
  theme(legend.title = element_blank())

ggsave("./figs/opilio_julian_temp_depth_long_cpue.png", width = 4.5, height = 7.5, units = 'in')

# DFA is hard here b/c we want to include time as one of the time series, *and* we don't have continuous observations for DFA

# following Tanner analysis - fit a PCA
# use julian, latitude, temperature
PCA <- prcomp(pca.dat[,c(3,5,6)], scale = T, center = T)
PCA$rotation
PCA$x
pca.dat$pc1 <- PCA$x[,1]

# and join pc1 back in
pc1 <- pca.dat %>%
  dplyr::select(station, year, pc1)

opilio.dat <- left_join(opilio.dat, pc1)

## define model formula

opilio1_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + (1 | year/index/station)) # simple first model

## Show default priors

get_prior(opilio1_formula, opilio.dat, family = bernoulli(link = "logit"))

## fit binomial mode --------------------------------------
opilio1 <- brm(opilio1_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))


saveRDS(opilio1, file = "./output/opilio1.rds")

opilio1 <- readRDS("./output/opilio1.rds")
check_hmc_diagnostics(opilio1$fit)
neff_lowest(opilio1$fit)
rhat_highest(opilio1$fit)
summary(opilio1)
bayes_R2(opilio1)

# plot(opilio1$criteria$loo, "k")
plot(conditional_smooths(opilio1), ask = FALSE)
# posterior predictive test

y <- opilio.dat$pcr
yrep_opilio1  <- fitted(opilio1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_opilio1[sample(nrow(yrep_opilio1), 25), ]) +
  ggtitle("opilio1")


# pdf("./figs/trace_opilio1.pdf", width = 6, height = 4)
trace_plot(opilio1$fit)
# dev.off()

###################################################
# add depth as covariate
opilio2_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(depth, k = 4) + (1 | year/index/station))

opilio2 <- brm(opilio2_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(opilio2, file = "./output/opilio2.rds")

opilio2 <- readRDS("./output/opilio2.rds")
check_hmc_diagnostics(opilio2$fit)
neff_lowest(opilio2$fit)
rhat_highest(opilio2$fit)
summary(opilio2) 
bayes_R2(opilio2)

plot(conditional_smooths(opilio2), ask = FALSE)

# posterior predictive test

y <- opilio.dat$pcr
yrep_opilio2_hier_bernoulli  <- fitted(opilio2_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_opilio2_hier_bernoulli[sample(nrow(yrep_opilio2_hier_bernoulli), 25), ]) +
  ggtitle("opilio2_hier_bernoulli")

# model comparison
loo(opilio1, opilio2) # depth does not improve prediction



###############################################################################################
# examine cpue effect

opilio3_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(fourth.root.cpue70, k = 4) + (1 | year/index/station))                      

opilio3 <- brm(opilio3_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

# opilio3  <- add_criterion(opilio3, "loo",
#                                          moment_match = TRUE)

saveRDS(opilio3, file = "./output/opilio3.rds")

opilio3 <- readRDS("./output/opilio3.rds")

check_hmc_diagnostics(opilio3$fit)
neff_lowest(opilio3$fit)
rhat_highest(opilio3$fit)
summary(opilio3) 
bayes_R2(opilio3)
# plot(opilio3$criteria$loo, "k")

plot(conditional_smooths(opilio3), ask = FALSE)

# posterior predictive test

y <- opilio.dat$pcr
yrep_opilio3  <- fitted(opilio3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_opilio3[sample(nrow(yrep_opilio3), 25), ]) +
  ggtitle("opilio3")

# let's run the model comparison
loo(opilio1, opilio2, opilio3) # opilio3

###############################################################################################
# see if a sex effect improves model

opilio4_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(fourth.root.cpue70, k = 4) + sex + (1 | year/index/station))                      

opilio4 <- brm(opilio4_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

# opilio4  <- add_criterion(opilio4, "loo",
#                                          moment_match = TRUE)

saveRDS(opilio4, file = "./output/opilio4.rds")

opilio4 <- readRDS("./output/opilio4.rds")

check_hmc_diagnostics(opilio4$fit)
neff_lowest(opilio4$fit)
rhat_highest(opilio4$fit)
summary(opilio4) # increased indicdence for sex2 (female?)
bayes_R2(opilio4)
# plot(opilio4$criteria$loo, "k")

# posterior predictive test

y <- opilio.dat$pcr
yrep_opilio4  <- fitted(opilio4, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_opilio4[sample(nrow(yrep_opilio4), 25), ]) +
  ggtitle("opilio4")


# let's run the model comparison
loo(opilio1, opilio2, opilio3, opilio4)

######################################################
# finally, check for a year effect

opilio5_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + 
                         s(fourth.root.cpue70, k = 4) + sex + year + (1 | year/index/station))                      

opilio5 <- brm(opilio5_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 4000, # increasing iterations 
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

# opilio5  <- add_criterion(opilio5, "loo",
#                                          moment_match = TRUE)

saveRDS(opilio5, file = "./output/opilio5.rds")

opilio5 <- readRDS("./output/opilio5.rds")

check_hmc_diagnostics(opilio5$fit)
neff_lowest(opilio5$fit) # too low!
rhat_highest(opilio5$fit)
summary(opilio5) # no evidence of a year effect
bayes_R2(opilio5)
# plot(opilio5$criteria$loo, "k")

# posterior predictive test

y <- opilio.dat$pcr
yrep_opilio5  <- fitted(opilio5, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_opilio5[sample(nrow(yrep_opilio5), 25), ]) +
  ggtitle("opilio5")

# png("./figs/trace_opilio5.png", width = 6, height = 4, units = 'in', res = 300)
trace_plot(opilio5$fit)
# dev.off()

## additional suggested model - size, pc1, year
opilio6_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + year + (1 | year/index/station))                      

opilio6 <- brm(opilio6_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 4000, # increasing iterations 
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

# opilio6  <- add_criterion(opilio6, "loo",
#                                          moment_match = TRUE)

saveRDS(opilio6, file = "./output/opilio6.rds")

opilio6 <- readRDS("./output/opilio6.rds")

check_hmc_diagnostics(opilio6$fit)
neff_lowest(opilio6$fit) # too low!
rhat_highest(opilio6$fit)
summary(opilio6) # no evidence of a year effect
bayes_R2(opilio6)
# plot(opilio6$criteria$loo, "k")

# posterior predictive test

y <- opilio.dat$pcr
yrep_opilio6  <- fitted(opilio6, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_opilio6[sample(nrow(yrep_opilio6), 25), ]) +
  ggtitle("opilio6")

# png("./figs/trace_opilio6.png", width = 6, height = 4, units = 'in', res = 300)
trace_plot(opilio6$fit)
# dev.off()

#######

model.comp <- loo(opilio1, opilio2, opilio3, opilio4, opilio5, opilio6) # opilio5 marginally the best!

model.comp

# save model comparison for ms.
forms <- data.frame(formula=c(as.character(opilio5_formula)[1],
                              as.character(opilio4_formula)[1],
                              as.character(opilio3_formula)[1],
                              as.character(opilio6_formula)[1],
                              as.character(opilio1_formula)[1],
                              as.character(opilio2_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/opilio_model_comp.csv")

###########################
# plot predicted effects from best model (opilio5)
opilio5 <- readRDS("./output/opilio5.rds")

# first year
ce1s_1 <- conditional_effects(opilio5, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__)

ggplot(plot, aes(year, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Probability positive") +
  ggtitle("Opilio5 - posterior mean & 95% credible interval")

ggsave("./figs/opilio5_year_effect.png", width = 2, height = 3, units = 'in')

# then sex
ce1s_1 <- conditional_effects(opilio5, effect = "sex", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$sex %>%
  dplyr::select(sex, estimate__, lower__, upper__)

ggplot(plot, aes(sex, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Probability positive") +
  ggtitle("Opilio5 - posterior mean & 95% credible interval")

ggsave("./figs/opilio5_sex_effect.png", width = 2, height = 3, units = 'in')

# then size

## 95% CI
ce1s_1 <- conditional_effects(opilio5, effect = "size", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opilio5, effect = "size", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opilio5, effect = "size", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$size
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$size[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$size[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$size[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$size[["lower__"]]

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Carapace width (mm)", y = "Probability positive") +
  ggtitle("opilio5 - posterior mean & 80 / 90 / 95% credible intervals")

ggsave("./figs/opilio5_size_effect.png", width = 6, height = 4, units = 'in')


## cpue

## 95% CI
ce1s_1 <- conditional_effects(opilio5, effect = "fourth.root.cpue70", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opilio5, effect = "fourth.root.cpue70", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opilio5, effect = "fourth.root.cpue70", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue70
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue70[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue70[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue70[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue70[["lower__"]]

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "fourth.root.cpue70", y = "Probability positive") +
  ggtitle("opilio5 - posterior mean & 80 / 90 / 95% credible intervals")

ggsave("./figs/opilio5_fourth.root.cpue70_effect.png", width = 6, height = 4, units = 'in')


## predict overall prevalence

# first, create a combined year_index_station column to identify combinations that exist in the data
opilio.dat$yr_ind_st <- paste(opilio.dat$year, opilio.dat$index, opilio.dat$station, sep = "_")


new.dat <- data.frame(yr_ind_st = unique(opilio.dat$yr_ind_st),
                      size = 30,
                      pc1 = mean(unique(opilio.dat$pc1)),
                      fourth.root.cpue70 = mean(unique(opilio.dat$fourth.root.cpue70)),
                      sex = 2) 

new.dat$year <- map_chr(str_split(new.dat$yr_ind_st, "_"), 1)
new.dat$index <- map_chr(str_split(new.dat$yr_ind_st, "_"), 2)
new.dat$station <- map_chr(str_split(new.dat$yr_ind_st, "_"), 3)



posterior.predict <- posterior_epred(opilio5, newdata = new.dat)

opilio.estimate <- data.frame(species = "opilio",
                              estimate = mean(posterior.predict),
                              lower_95 = quantile(posterior.predict, probs = 0.025),
                              upper_95 = quantile(posterior.predict, probs = 0.975),
                              lower_90 = quantile(posterior.predict, probs = 0.05),
                              upper_90 = quantile(posterior.predict, probs = 0.95),
                              lower_80 = quantile(posterior.predict, probs = 0.1),
                              upper_80 = quantile(posterior.predict, probs = 0.9))

ggplot(opilio.estimate) +
  aes(x = species, y = estimate) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), color = "grey80") +
  geom_errorbar(aes(ymin = lower_90, ymax = upper_90), color = "grey60") +
  geom_errorbar(aes(ymin = lower_80, ymax = upper_80), color = "black") +
  geom_point(size = 3, color = "red3") +
  theme_classic()

ggsave("./figs/opilio estimate.png", width = 2, height = 2.5, units = 'in')
