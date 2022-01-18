# explore and analyze infection dynamics in C. bairdi

library(tidyverse)

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
  group_by(Species_Name, Clutch) %>%
  summarise(count = n())

check # that's a lot of data to throw away!


# select bairdi for analysis

bairdi <- dat %>%
  filter(Species_Name == "Chionoecetes bairdi",
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
library(plyr)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

## note that seasonal time of sampling is correlated with space (South to North)

## fit an initial model with index identity, year, size as covariates (and year/index/station as random term)

# clean up data again!

tanner.dat <- dat %>%
  dplyr::select(PCR_result, Size, index, Year, STATIONID) %>%
  dplyr::rename(pcr = PCR_result, size = Size, year = Year, station = STATIONID) %>%
  dplyr::mutate(year = as.factor(year),
                index = as.factor(index),
                station = as.factor(station)) %>%
  na.omit()

## define model formula
tanner1_formula <-  bf(pcr ~ s(size, k = 3) + index + year + (1 | year/index/station))
                       
## Set model family
family <- bernoulli(link = "logit")


## Show default priors
get_prior(tanner1_formula, tanner.dat, family = family)

## fit binomial mode --------------------------------------
tanner1_bernoulli <- brm(tanner1_formula,
                    data = tanner.dat,
                    family = family,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))

# cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
                               # moment_match = TRUE)
saveRDS(tanner1_bernoulli, file = "./output/tanner1_bernoulli.rds")

tanner1_bernoulli <- readRDS("./output/tanner1_bernoulli.rds")
check_hmc_diagnostics(tanner1_bernoulli$fit)
neff_lowest(tanner1_bernoulli$fit)
rhat_highest(tanner1_bernoulli$fit)
summary(tanner1_bernoulli)
bayes_R2(tanner1_bernoulli)
plot(tanner1_bernoulli$criteria$loo, "k")
plot(conditional_smooths(tanner1_bernoulli), ask = FALSE)
y <- cod.data$cod
yrep_tanner1_bernoulli  <- fitted(tanner1_bernoulli, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_tanner1_bernoulli[sample(nrow(yrep_tanner1_bernoulli), 25), ]) +
  xlim(0, 500) +
  ggtitle("tanner1_bernoulli")
pdf("./figs/trace_tanner1_bernoulli.pdf", width = 6, height = 4)
trace_plot(tanner1_bernoulli$fit)
dev.off()

