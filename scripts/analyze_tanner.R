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


# load PCR data for each year
d1=read.csv("./data/bcs14.csv")
d2=read.csv("./data/bcs15.csv")
d3=read.csv("./data/bcs16.csv") 

names(d1); names(d2); names(d3)

name.check <- data.frame(names1 = names(d1),
                         names2 = names(d2),
                         names3 = names(d3))

name.check # so names are interchangeable, just some formatting differences

names(d2) <- names(d3) <- names(d1)
dat = rbind(d1, d2, d3) 

# add julian date

dat$date <- stringr::str_split_fixed(dat$START_TIME, " ", n = 2)[,1] # remove time, keep date

dat$julian <- lubridate::yday(chron::dates(dat$date))

# add index column

substr(dat$Specific_Location,10,10)

dat$index = as.numeric(substr(dat$Specific_Location,10,10))

# separate Tanner data
tanner.dat <- dat %>%
  dplyr::filter(Species_Name == "Chionoecetes bairdi",
                index %in% c(1, 2, 3),
                Year %in% c(2015, 2016),
                Sex %in% c(1, 2)) %>%
  dplyr::select(PCR_result, Size, Sex, index, Year, STATIONID, julian, START_LONGITUDE, BOTTOM_DEPTH, Bottom_Temp) %>%
  dplyr::rename(pcr = PCR_result,
                size = Size,
                sex = Sex,
                year = Year,
                station = STATIONID,
                longitude = START_LONGITUDE,
                depth = BOTTOM_DEPTH,
                temperature = Bottom_Temp) %>%
  dplyr::mutate(year = as.factor(year),
                sex = as.factor(sex),
                index = as.factor(index),
                station = as.factor(station)) 

nrow(tanner.dat) # 886 samples!

# need dimension reduction for exogenous covariates (day, depth, longitude, temperature)

pca.dat <- tanner.dat %>%
  dplyr::group_by(station, year) %>%
  dplyr::summarise(julian = mean(julian),
                   depth = mean(depth),
                   longitude = -mean(longitude),
                   temperature = mean(temperature))

cor(pca.dat[,3:6])

# plot - for the paper?
plot <- data.frame(Day_of_year = pca.dat$julian,
                   'Depth_m' = pca.dat$depth,
                   'W_longitude' = pca.dat$longitude,
                   'Bottom_temperature_C' = pca.dat$temperature,
                   year = as.numeric(as.character(pca.dat$year))) %>%
  tidyr::pivot_longer(cols = c(-Day_of_year, -year))
  
ggplot(plot, aes(Day_of_year, value)) +
  geom_point(color = "grey100") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, color = "black", lwd = 0.3) +
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  geom_point(aes(color = as.factor(year))) +
  scale_color_manual(values = cb[c(2,4)]) +
  theme(legend.title = element_blank())

ggsave("./figs/tanner_julian_temp_depth_long.png", width = 5, height = 6, units = 'in')

# DFA is hard here b/c we want to include time as one of the time series, *and* we don't have continuous observations for DFA

# could just fit a PCA!
PCA <- prcomp(dfa.dat[,3:6], scale = T, center = T)
PCA$rotation
PCA$x
pca.dat$pc1 <- PCA$x[,1]

# and join pc1 back in
tanner.dat <- left_join(tanner.dat, pca.dat)


## define model formula
tanner1_formula <-  bf(pcr ~ sex + maturity + index + year) # simple first model


tanner1.hier_formula <-  bf(pcr ~ sex + maturity + index + year + (1 | year/index/station)) # and hierarchical version of same                       

## Set model family
family <- bernoulli(link = "logit")


## Show default priors
get_prior(tanner1_formula, tanner.dat, family = family)

## fit binomial mode --------------------------------------
tanner1_bernoulli <- brm(tanner1_formula,
                         data = tanner.dat,
                         family = family,
                         cores = 4, chains = 4, iter = 2500,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 14))

# cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
# moment_match = TRUE)
saveRDS(tanner1_bernoulli, file = "./output/tanner1_bernoulli.rds")

tanner1_bernoulli <- readRDS("./output/tanner1_bernoulli.rds")
check_hmc_diagnostics(tanner1_bernoulli$fit)
neff_lowest(tanner1_bernoulli$fit)
rhat_highest(tanner1_bernoulli$fit)
summary(tanner1_bernoulli)
bayes_R2(tanner1_bernoulli)
# plot(tanner1_bernoulli$criteria$loo, "k")

# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner1_bernoulli  <- fitted(tanner1_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner1_bernoulli[sample(nrow(yrep_tanner1_bernoulli), 25), ]) +
  ggtitle("tanner1_bernoulli")

# this is not a good fit!

# pdf("./figs/trace_tanner1_bernoulli.pdf", width = 6, height = 4)
trace_plot(tanner1_bernoulli$fit)
# dev.off()

###################################################
# see if the hierarchical model is an improvement!

tanner1_hier_bernoulli <- brm(tanner1.hier_formula,
                              data = tanner.dat,
                              family = family,
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

# cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
# moment_match = TRUE)
saveRDS(tanner1_hier_bernoulli, file = "./output/tanner1_hier_bernoulli.rds")

tanner1_hier_bernoulli <- readRDS("./output/tanner1_hier_bernoulli.rds")
check_hmc_diagnostics(tanner1_hier_bernoulli$fit)
neff_lowest(tanner1_hier_bernoulli$fit)
rhat_highest(tanner1_hier_bernoulli$fit)
summary(tanner1_hier_bernoulli)
bayes_R2(tanner1_hier_bernoulli)
# plot(tanner1_hier_bernoulli$criteria$loo, "k")

# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner1_hier_bernoulli  <- fitted(tanner1_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner1_hier_bernoulli[sample(nrow(yrep_tanner1_hier_bernoulli), 25), ]) +
  ggtitle("tanner1_hier_bernoulli")

# the model still isn't great
# drop year as it seems unimportant

tanner2.hier_formula <-  bf(pcr ~ sex + maturity + index + (1 | year/index/station)) # and hierarchical version of same                       

tanner2_hier_bernoulli <- brm(tanner2.hier_formula,
                              data = tanner.dat,
                              family = family,
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

# cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
# moment_match = TRUE)
saveRDS(tanner2_hier_bernoulli, file = "./output/tanner2_hier_bernoulli.rds")

tanner2_hier_bernoulli <- readRDS("./output/tanner2_hier_bernoulli.rds")
check_hmc_diagnostics(tanner2_hier_bernoulli$fit)
neff_lowest(tanner2_hier_bernoulli$fit)
rhat_highest(tanner2_hier_bernoulli$fit)
summary(tanner2_hier_bernoulli)
bayes_R2(tanner2_hier_bernoulli)
# plot(tanner2_hier_bernoulli$criteria$loo, "k")

# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner2_hier_bernoulli  <- fitted(tanner2_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner2_hier_bernoulli[sample(nrow(yrep_tanner2_hier_bernoulli), 25), ]) +
  ggtitle("tanner2_hier_bernoulli")

# still poor 
# add size

tanner3.hier_formula <-  bf(pcr ~ sex + s(size, k = 3) + index + (1 | year/index/station)) # and hierarchical version of same                       

tanner3_hier_bernoulli <- brm(tanner3.hier_formula,
                              data = tanner.dat,
                              family = family,
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

tanner3_hier_bernoulli  <- add_criterion(tanner3_hier_bernoulli, "loo",
                                         moment_match = TRUE)

saveRDS(tanner3_hier_bernoulli, file = "./output/tanner3_hier_bernoulli.rds")

tanner3_hier_bernoulli <- readRDS("./output/tanner3_hier_bernoulli.rds")

check_hmc_diagnostics(tanner3_hier_bernoulli$fit)
neff_lowest(tanner3_hier_bernoulli$fit)
rhat_highest(tanner3_hier_bernoulli$fit)
summary(tanner3_hier_bernoulli)
bayes_R2(tanner3_hier_bernoulli)
# plot(tanner3_hier_bernoulli$criteria$loo, "k")

# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner3_hier_bernoulli  <- fitted(tanner3_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner3_hier_bernoulli[sample(nrow(yrep_tanner3_hier_bernoulli), 25), ]) +
  ggtitle("tanner3_hier_bernoulli")

# still poor

# let's run the model comparison
tanner1_bernoulli <- readRDS("./output/tanner1_bernoulli.rds")
tanner1_hier_bernoulli <- readRDS("./output/tanner1_hier_bernoulli.rds")
tanner2_hier_bernoulli <- readRDS("./output/tanner2_hier_bernoulli.rds")
tanner3_hier_bernoulli <- readRDS("./output/tanner3_hier_bernoulli.rds")

model.comp <- loo(tanner1_bernoulli, tanner1_hier_bernoulli, tanner2_hier_bernoulli, tanner3_hier_bernoulli)

model.comp # the three hierarchical models are very similar!

# try a model with size and season!

tanner4.hier_formula <-  bf(pcr ~ s(size, k = 3) + s(julian, k = 3) + (1 | year/index/station))

tanner4_hier_bernoulli <- brm(tanner4.hier_formula,
                              data = tanner.dat,
                              family = family,
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

# tanner4_hier_bernoulli  <- add_criterion(tanner4_hier_bernoulli, "loo",
# moment_match = TRUE)

saveRDS(tanner4_hier_bernoulli, file = "./output/tanner4_hier_bernoulli.rds")

tanner4_hier_bernoulli <- readRDS("./output/tanner4_hier_bernoulli.rds")

check_hmc_diagnostics(tanner4_hier_bernoulli$fit)
neff_lowest(tanner4_hier_bernoulli$fit)
rhat_highest(tanner4_hier_bernoulli$fit)
summary(tanner4_hier_bernoulli)
bayes_R2(tanner4_hier_bernoulli)

# plot(tanner4_hier_bernoulli$criteria$loo, "k")
plot(conditional_smooths(tanner4_hier_bernoulli), ask = FALSE)
# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner4_hier_bernoulli  <- fitted(tanner4_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner4_hier_bernoulli[sample(nrow(yrep_tanner4_hier_bernoulli), 25), ]) +
  ggtitle("tanner4_hier_bernoulli")

ggsave("./figs/tanner_4_bernoulli_ppc_overlay.png", width = 5, height = 4, units = 'in')

# and trace plot
png("./figs/trace_tanner4_hier_bernoulli.png", width = 6, height = 4, units = 'in', res = 300)
trace_plot(tanner4_hier_bernoulli$fit)
dev.off()




tanner1_bernoulli <- readRDS("./output/tanner1_bernoulli.rds")
tanner1_hier_bernoulli <- readRDS("./output/tanner1_hier_bernoulli.rds")
tanner2_hier_bernoulli <- readRDS("./output/tanner2_hier_bernoulli.rds")
tanner3_hier_bernoulli <- readRDS("./output/tanner3_hier_bernoulli.rds")
tanner4_hier_bernoulli <- readRDS("./output/tanner4_hier_bernoulli.rds")

model.comp <- loo(tanner1_bernoulli, tanner1_hier_bernoulli, tanner2_hier_bernoulli, 
                  tanner3_hier_bernoulli, tanner4_hier_bernoulli)

model.comp

# add sex - support for that covariate?
tanner5.hier_formula <-  bf(pcr ~ s(size, k = 3) + s(julian, k = 3) + sex + (1 | year/index/station))

tanner5_hier_bernoulli <- brm(tanner5.hier_formula,
                              data = tanner.dat,
                              family = family,
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

# tanner5_hier_bernoulli  <- add_criterion(tanner5_hier_bernoulli, "loo",
# moment_match = TRUE)

saveRDS(tanner5_hier_bernoulli, file = "./output/tanner5_hier_bernoulli.rds")

tanner5_hier_bernoulli <- readRDS("./output/tanner5_hier_bernoulli.rds")

check_hmc_diagnostics(tanner5_hier_bernoulli$fit)
neff_lowest(tanner5_hier_bernoulli$fit)
rhat_highest(tanner5_hier_bernoulli$fit)
summary(tanner5_hier_bernoulli) # no suggestion of a sex effect!
bayes_R2(tanner5_hier_bernoulli)

# plot(tanner5_hier_bernoulli$criteria$loo, "k")
plot(conditional_smooths(tanner5_hier_bernoulli), ask = FALSE)
# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner5_hier_bernoulli  <- fitted(tanner5_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner5_hier_bernoulli[sample(nrow(yrep_tanner5_hier_bernoulli), 25), ]) +
  ggtitle("tanner5_hier_bernoulli")

tanner1_bernoulli <- readRDS("./output/tanner1_bernoulli.rds")
tanner1_hier_bernoulli <- readRDS("./output/tanner1_hier_bernoulli.rds")
tanner2_hier_bernoulli <- readRDS("./output/tanner2_hier_bernoulli.rds")
tanner3_hier_bernoulli <- readRDS("./output/tanner3_hier_bernoulli.rds")
tanner4_hier_bernoulli <- readRDS("./output/tanner4_hier_bernoulli.rds")
tanner5_hier_bernoulli <- readRDS("./output/tanner5_hier_bernoulli.rds")

model.comp <- loo(tanner1_bernoulli, tanner1_hier_bernoulli, tanner2_hier_bernoulli, 
                  tanner3_hier_bernoulli, tanner4_hier_bernoulli, tanner5_hier_bernoulli)

model.comp # model 4 still the best

# try with a different random structure
tanner.dat$station_year <- as.factor(paste(tanner.dat$station, tanner.dat$year, sep = "_"))


tanner6.hier_formula <-  bf(pcr ~ s(size, k = 3) + s(julian, k = 3) + (1 | station_year))

tanner6_hier_bernoulli <- brm(tanner6.hier_formula,
                              data = tanner.dat,
                              family = bernoulli(link = "logit"),
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

# tanner6_hier_bernoulli  <- add_criterion(tanner6_hier_bernoulli, "loo",
# moment_match = TRUE)

saveRDS(tanner6_hier_bernoulli, file = "./output/tanner6_hier_bernoulli.rds")

tanner6_hier_bernoulli <- readRDS("./output/tanner6_hier_bernoulli.rds")

check_hmc_diagnostics(tanner6_hier_bernoulli$fit)
neff_lowest(tanner6_hier_bernoulli$fit)
rhat_highest(tanner6_hier_bernoulli$fit)
summary(tanner6_hier_bernoulli)
bayes_R2(tanner6_hier_bernoulli)

# plot(tanner6_hier_bernoulli$criteria$loo, "k")
plot(conditional_smooths(tanner6_hier_bernoulli), ask = FALSE)
# posterior predictive test

y <- tanner.dat$pcr
yrep_tanner6_hier_bernoulli  <- fitted(tanner6_hier_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner6_hier_bernoulli[sample(nrow(yrep_tanner6_hier_bernoulli), 25), ]) +
  ggtitle("tanner6_hier_bernoulli")

tanner1_bernoulli <- readRDS("./output/tanner1_bernoulli.rds")
tanner1_hier_bernoulli <- readRDS("./output/tanner1_hier_bernoulli.rds")
tanner2_hier_bernoulli <- readRDS("./output/tanner2_hier_bernoulli.rds")
tanner3_hier_bernoulli <- readRDS("./output/tanner3_hier_bernoulli.rds")
tanner4_hier_bernoulli <- readRDS("./output/tanner4_hier_bernoulli.rds")
tanner5_hier_bernoulli <- readRDS("./output/tanner5_hier_bernoulli.rds")
tanner6_hier_bernoulli <- readRDS("./output/tanner6_hier_bernoulli.rds")

model.comp <- loo(tanner1_bernoulli, tanner1_hier_bernoulli, tanner2_hier_bernoulli, 
                  tanner3_hier_bernoulli, tanner4_hier_bernoulli, tanner5_hier_bernoulli, tanner6_hier_bernoulli)

model.comp # model 4 still the best


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