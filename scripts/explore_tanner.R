# explore and analyze infection dynamics in C. bairdi

library(tidyverse)
library(plyr)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

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

str(dat)

# set up an index column
substr(dat$Specific_Location,10,10)
dat$index = as.numeric(substr(dat$Specific_Location,10,10))

#it appears there are some crabs from outside index areas?
unique(dat$index)
idxs = is.na(dat$index)
sum(idxs) # only one!
which(idxs == T)
dat = dat[!idxs,]

dat = dat[dat$PCR_result == 0 | dat$PCR_result == 1,] # limiting to PCR positive or negative - dropping code 3 (results uncertain)

unique(dat$Species_Name)

# removing species == NA
idxs = dat$Species_Name == ""
sum(idxs) # only one

dat=dat[!idxs,]


# check clutch codes 
check <- dat %>%
  dplyr::group_by(Species_Name, Clutch) %>%
  dplyr::summarise(count = dplyr::n())

check # that's a lot of data to throw away!


# select bairdi for analysis

bairdi <- dat %>%
  dplyr::filter(Species_Name == "Chionoecetes bairdi",
         index %in% c(1, 2, 3))

bairdi$index <- as.factor(bairdi$index)

# examine sample size by index / year

# now the same for index/year
plot <- bairdi %>%
  dplyr::select(Year, index, PCR_result) %>%
  dplyr::group_by(Year, index) %>%
  dplyr::summarise(PCR_0 = sum(PCR_result == 0),
                   PCR_1 = sum(PCR_result == 1)) %>%
  tidyr::pivot_longer(cols = c(-Year, -index))

ggplot(plot, aes(name, value)) +
  geom_bar(stat = "identity") +
  facet_grid(Year ~ index)

ggsave("./figs/bairdi_index_year_incidence.png", width = 6, height = 6, units = "in")

# plot surface temp - bottom temp relationship
ggplot(bairdi, aes(SURFACE_TEMP, Bottom_Temp)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F) +
  geom_point(aes(SURFACE_TEMP, Bottom_Temp, color = index))

ggsave("./figs/bairdi_surface_temp_bottom_temp.png", width = 6, height = 4, units = "in")

# and plot surface / bottom temps / bottom depth / dates by year / index

# first, calculate julian day
bairdi$date <- stringr::str_split_fixed(bairdi$START_TIME, " ", n = 2)[,1] # remove time, keep date

bairdi$julian <- lubridate::yday(chron::dates(bairdi$date))

plot <- bairdi %>%
  dplyr::group_by(Year, index, STATIONID) %>%
  dplyr::summarise(surface = mean(SURFACE_TEMP),
                   bottom = mean(Bottom_Temp),
                   depth = mean(BOTTOM_DEPTH),
                   julian = mean(julian)) 


ggplot(plot, aes(bottom)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(Year ~ index)

ggsave("./figs/bairdi_bottom_temp_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(depth)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(Year ~ index)

ggsave("./figs/bairdi_depth_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(Year ~ index) # big differences among index areas

ggsave("./figs/bairdi_date_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(julian, bottom)) +
  geom_point() # this is a big problem - julian day is collinear with temperature

ggsave("./figs/bairdi_date_vs_bottom_temp.png", width = 6, height = 4, units = "in")


# look at Julian day as a predictor of % positive
plot <- bairdi %>%
  dplyr::group_by(Year, STATIONID) %>%
  dplyr::summarise(julian = mean(julian), 
                   n = dplyr::n(),
                   proportion.positive = sum(PCR_result == 1) / dplyr::n())

ggplot(plot, aes(julian, proportion.positive)) +
  geom_point() + 
  facet_wrap(~Year) +
  geom_smooth(method = "gam")

ggsave("./figs/bairdi_date_vs_percent_positive.png", width = 6, height = 4, units = "in")


# and size as a predictor of % positive
plot <- bairdi %>%
  dplyr::group_by(Year, STATIONID) %>%
  dplyr::summarise(size = mean(Size), 
                   n = dplyr::n(),
                   proportion.positive = sum(PCR_result == 1) / dplyr::n())

ggplot(plot, aes(size, proportion.positive)) +
  geom_point() + 
  facet_wrap(~Year) +
  geom_smooth(method = "gam")

ggsave("./figs/bairdi_size_vs_percent_positive.png", width = 6, height = 4, units = "in")


# make a toy model with julian day, size, and year as predictors

toy.dat <- bairdi %>%
  dplyr::group_by(Year, STATIONID) %>%
  dplyr::summarise(size = mean(Size),
                   julian = mean(julian),
                   proportion.positive = sum(PCR_result == 1) / dplyr::n())


toy.dat$Year <- as.factor(toy.dat$Year)


mod <- mgcv::gam(proportion.positive ~ Year + s(julian, k = 4) + s(size, k = 4), dat = toy.dat)
summary(mod)

plot(mod)

## have a go at fitting simple models in brms -----------------------------------

## note that seasonal time of sampling is correlated with space (South to North)

## fit an initial model with index identity, year, size as covariates (and year/index/station as random term)

# clean up data again!

# and, based on conversation, we'll drop 2014, drop any data without sex and maturity data

# start by defining julian day
dat$date <- stringr::str_split_fixed(dat$START_TIME, " ", n = 2)[,1] # remove time, keep date

dat$julian <- lubridate::yday(chron::dates(dat$date))

tanner.dat <- dat %>%
  dplyr::filter(Species_Name == "Chionoecetes bairdi",
         index %in% c(1, 2, 3),
         Year %in% c(2015, 2016),
         Sex %in% c(1, 2),
         Maturity %in% c(0, 1)) %>%
  dplyr::select(PCR_result, Size, Sex, Maturity, index, Year, STATIONID, julian) %>%
  dplyr::rename(pcr = PCR_result,
                size = Size,
                sex = Sex,
                maturity = Maturity,
                year = Year,
                station = STATIONID) %>%
  dplyr::mutate(year = as.factor(year),
                sex = as.factor(sex),
                maturity = as.factor(maturity),
                index = as.factor(index),
                station = as.factor(station)) %>%
  na.omit()

nrow(tanner.dat) # 685 samples!



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

╠# cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
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

╠# cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
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

╠# cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
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

############
# save model comparison for ms.
forms <- data.frame(formula=c(as.character(growth2s_formula)[1],
                              as.character(growth2_formula)[1],
                              as.character(growth1_formula)[1],
                              as.character(growth0_formula)[1],
                              as.character(growth3_formula)[1],
                              as.character(growth4_formula)[1],
                              as.character(growth5_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/growth_model_comp.csv")