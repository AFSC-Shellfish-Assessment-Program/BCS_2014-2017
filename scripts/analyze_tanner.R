# analyze infection dynamics in C. bairdi

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

# examine cpue distribution                                    
ggplot(dat, aes(tanner70under_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(dat, aes(tanner70under_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

# separate Tanner data
tanner.dat <- dat %>%
  dplyr::filter(species_name == "Chionoecetes bairdi",
                index_site %in% c(1, 2, 3),
                year %in% c(2015:2017),
                sex %in% c(1, 2),
                pcr_result %in% c(1, 0)) %>%
  dplyr::select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth, gear_temperature, tanner70under_cpue) %>%
  dplyr::rename(pcr = pcr_result,
                station = gis_station,
                longitude = mid_longitude,
                depth = bottom_depth,
                temperature = gear_temperature,
                index = index_site,
                fourth.root.cpue70 = tanner70under_cpue) %>%
  dplyr::mutate(year = as.factor(year),
                sex = as.factor(sex),
                index = as.factor(index),
                station = as.factor(station),
                fourth.root.cpue70 = fourth.root.cpue70^0.25) # transforming cpue here

nrow(tanner.dat) # 1285 samples!

# need dimension reduction for exogenous covariates (day, depth, longitude, temperature, cpue)

pca.dat <- tanner.dat %>%
  dplyr::group_by(station, year) %>%
  dplyr::summarise(julian = mean(julian),
                   depth = mean(depth),
                   longitude = -mean(longitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70))

cor(pca.dat[,3:7]) # temp no longer correlated with others! # cpue weakly collinear with depth

# plot - for the paper?
plot <- data.frame(Day_of_year = pca.dat$julian,
                   Depth_m = pca.dat$depth,
                   W_longitude = pca.dat$longitude,
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

ggsave("./figs/tanner_julian_temp_depth_long.png", width = 4.5, height = 7.5, units = 'in')

# DFA is hard here b/c we want to include time as one of the time series, *and* we don't have continuous observations for DFA

# could just fit a PCA!
PCA <- prcomp(pca.dat[,3:5], scale = T, center = T)
PCA$rotation
PCA$x
pca.dat$pc1 <- PCA$x[,1]

# and join pc1 back in
pc1 <- pca.dat %>%
  dplyr::select(station, year, pc1)

tanner.dat <- left_join(tanner.dat, pc1)

## define model formula
tanner1_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + (1 | year/index/station)) # simple first model

## Show default priors
get_prior(tanner1_formula, tanner.dat, family = bernoulli(link = "logit"))

## fit binomial mode --------------------------------------
tanner1 <- brm(tanner1_formula,
                         data = tanner.dat,
                         family = bernoulli(link = "logit"),
                         cores = 4, chains = 4, iter = 2500,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 14))


saveRDS(tanner1, file = "./output/tanner1.rds")

tanner1 <- readRDS("./output/tanner1.rds")
check_hmc_diagnostics(tanner1$fit)
neff_lowest(tanner1$fit)
rhat_highest(tanner1$fit)
summary(tanner1)
bayes_R2(tanner1)
plot(conditional_smooths(tanner1), ask = FALSE)

# plot(tanner1$criteria$loo, "k")
plot(conditional_smooths(tanner1), ask = FALSE)
# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner1  <- fitted(tanner1, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner1[sample(nrow(yrep_tanner1), 25), ]) +
  ggtitle("tanner1")

# this is not a good fit? I don't have experience with logistic models in brms

# pdf("./figs/trace_tanner1.pdf", width = 6, height = 4)
trace_plot(tanner1$fit)
# dev.off()

###################################################
# add temp as covariate
tanner2_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(temperature, k = 4) + (1 | year/index/station))

tanner2 <- brm(tanner2_formula,
                              data = tanner.dat,
                              family = bernoulli(link = "logit"),
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(tanner2, file = "./output/tanner2.rds")

tanner2 <- readRDS("./output/tanner2.rds")
check_hmc_diagnostics(tanner2$fit)
neff_lowest(tanner2$fit)
rhat_highest(tanner2$fit)
summary(tanner2) 
bayes_R2(tanner2)
# plot(tanner2_hier_bernoulli$criteria$loo, "k")

# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner2_hier_bernoulli  <- fitted(tanner2_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner2_hier_bernoulli[sample(nrow(yrep_tanner2_hier_bernoulli), 25), ]) +
  ggtitle("tanner2_hier_bernoulli")

# model comparison
loo(tanner1, tanner2) # temp does not improve prediction

# loo(tanner1, tanner2, moment_match = T) # moment matching crashes - need to try updating R / packages

###############################################################################################
# see if a sex effect improves model

tanner3_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + sex + (1 | year/index/station))                      

tanner3 <- brm(tanner3_formula,
                              data = tanner.dat,
                              family =bernoulli(link = "logit"),
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

# tanner3  <- add_criterion(tanner3, "loo",
#                                          moment_match = TRUE)

saveRDS(tanner3, file = "./output/tanner3.rds")

tanner3 <- readRDS("./output/tanner3.rds")

check_hmc_diagnostics(tanner3$fit)
neff_lowest(tanner3$fit)
rhat_highest(tanner3$fit)
summary(tanner3) # no evidence of a sex effect
bayes_R2(tanner3)
# plot(tanner3$criteria$loo, "k")

# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner3  <- fitted(tanner3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner3[sample(nrow(yrep_tanner3), 25), ]) +
  ggtitle("tanner3")

# still poor

# let's run the model comparison
loo(tanner1, tanner2, tanner3)

# finally, check for a year effect

tanner5_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + year + (1 | year/index/station))                      

tanner5 <- brm(tanner5_formula,
               data = tanner.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 5000, # increasing iterations - can reduce to 4000 for re-runs
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

# tanner5  <- add_criterion(tanner5, "loo",
#                                          moment_match = TRUE)

saveRDS(tanner5, file = "./output/tanner5.rds")

tanner5 <- readRDS("./output/tanner5.rds")

check_hmc_diagnostics(tanner5$fit)
neff_lowest(tanner5$fit) # too low!
rhat_highest(tanner5$fit)
summary(tanner5) # no evidence of a year effect
bayes_R2(tanner5)
# plot(tanner5$criteria$loo, "k")

# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner5  <- fitted(tanner5, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner5[sample(nrow(yrep_tanner5), 25), ]) +
  ggtitle("tanner5")

# png("./figs/trace_tanner5.png", width = 6, height = 4, units = 'in', res = 300)
trace_plot(tanner5$fit)
# dev.off()

loo(tanner1, tanner2, tanner3, tanner5) # tanner5 marginally the best!




tanner1_formula <-  bf(pcr ~ sex + maturity + index + year) # simple first model
tanner1.hier_formula <-  bf(pcr ~ sex + maturity + index + year + (1 | year/index/station)) # and hierarchical version of same                       
tanner2.hier_formula <-  bf(pcr ~ sex + maturity + index + (1 | year/index/station)) # and hierarchical version of same                       
tanner3.hier_formula <-  bf(pcr ~ sex + s(size, k = 3) + index + (1 | year/index/station)) # and hierarchical version of same                       
tanner4.hier_formula <-  bf(pcr ~ s(size, k = 3) + s(julian, k = 3) + (1 | year/index/station))
tanner5.hier_formula <-  bf(pcr ~ s(size, k = 3) + s(julian, k = 3) + sex + (1 | year/index/station))
tanner6.hier_formula <-  bf(pcr ~ s(size, k = 3) + s(julian, k = 3) + (1 | station_year))

############
# save model comparison for ms.
forms <- data.frame(formula=c(as.character(tanner4.hier_formula)[1],
                              as.character(tanner3.hier_formula)[1],
                              as.character(tanner6.hier_formula)[1],
                              as.character(tanner5.hier_formula)[1],
                              as.character(tanner1.hier_formula)[1],
                              as.character(tanner2.hier_formula)[1],
                              as.character(tanner1_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/growth_model_comp.csv")