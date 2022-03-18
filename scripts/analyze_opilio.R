# notes ----
#Analyze BCS infection dynamics in C.opilio using Bayesian multivariate models 

#load
library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(MARSS)
library(corrplot)
library(factoextra)
source("./scripts/stan_utils.R")

# set plot theme
theme_set(theme_bw())

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

########################################
#Data Manipulation
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(species_name == "Chionoecetes opilio",
         index_site %in% c(4, 5, 6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
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

nrow(opilio.dat) # 1511 samples!

#Check for missing data for PCA
opilio.dat %>%
  select(size, sex, year, index, station, julian, latitude, depth, temperature) %>%
  filter(!complete.cases(.)) 
 # Looks like just one missing size obsv. 

#Dimension reduction for exogenous covariates 
opilio.dat %>%
  group_by(year, station) %>%
  summarise(julian = mean(julian),
                   depth = mean(depth),
                   latitude = mean(latitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70)) -> pca.dat

cor(pca.dat[,3:7]) # some lower correlations than for tanner
corrplot(cor(pca.dat[,3:7])) 

# plot exogenous variables
pca.dat %>%
  rename(`Day of Year` = julian,
         `Depth (m)` = depth,
         `Latitude (N)` = latitude,
         `Bottom Temperature (C)` = temperature,
         `Fourth root CPUE` = fourth.root.cpue70) %>%
  pivot_longer(cols=`Depth (m)`:`Fourth root CPUE`) %>%
  ggplot(aes(`Day of Year`, value)) +
  geom_point(color = "grey100") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, color = "black", lwd = 0.3) +
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  geom_point(aes(color = as.factor(year))) +
  scale_color_manual(values = cb[c(2,4,6)]) +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank())
ggsave("./figs/opilio_julian_temp_depth_long_cpue.png", width = 4.5, height = 7.5, units = 'in')


#Fit a PCA with all exogenous covariates
PCA <- prcomp(pca.dat[,3:7], scale = T, center = T)
PCA$rotation #PC1 loads strongly/equally on all but depth 
get_eig(PCA)

#Fit PCA without Depth and extract PC1 for model runs 
PCA2 <- prcomp(pca.dat[,c(3,5:7)], scale = T, center = T)
PCA2$rotation #Variable loadings
get_eig(PCA2)
fviz_eig(PCA2) #Scree plot: PC1 explains ~72% of variance 
fviz_pca_var(PCA2, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)     
pca.dat$pc1 <- PCA2$x[,1]

#Join pc1 back in
pc1 <- pca.dat %>%
  dplyr::select(station, year, pc1)

#Final dataset for modeling 
opilio.dat <- left_join(opilio.dat, pc1)

####################################################
#Model 1: base model with size, pc1 and random year/index/station intercept 

## define model formula
opilio1_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + (1 | year/index/station)) 

## Show default priors
get_prior(opilio1_formula, opilio.dat, family = bernoulli(link = "logit"))

## fit binomial Model 1
opilio1 <- brm(opilio1_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio1, file = "./output/opilio1.rds")
opilio1 <- readRDS("./output/opilio1.rds")

#Model convergence diagnostics 
check_hmc_diagnostics(opilio1$fit)
neff_lowest(opilio1$fit)
rhat_highest(opilio1$fit)
summary(opilio1)
bayes_R2(opilio1)

#Plots
plot(conditional_smooths(opilio1), ask = FALSE)
plot(opilio1)
trace_plot(opilio1$fit)

#Posterior Predictive check 
pp_check(opilio1, nsamples = 100)

###################################################
#Model 2: Base model + depth 
  
opilio2_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(depth, k = 4) + (1 | year/index/station))

opilio2 <- brm(opilio2_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio2, file = "./output/opilio2.rds")
opilio2 <- readRDS("./output/opilio2.rds")

#Model convergence diagnostics 
check_hmc_diagnostics(opilio2$fit)
neff_lowest(opilio2$fit)
rhat_highest(opilio2$fit)
summary(opilio2)
bayes_R2(opilio2)

#Plot
plot(conditional_smooths(opilio2), ask = FALSE)
plot(opilio2)
trace_plot(opilio2$fit)

#Posterior Predictive check 
pp_check(opilio2, nsamples = 100)

# model comparison
loo(opilio1, opilio2) # oplio1 better, no support for adding depth

###############################################################################################
#Model 3: Base Model + Sex

opilio3_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + sex + (1 | year/index/station))                      

opilio3 <- brm(opilio3_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio3, file = "./output/opilio3.rds")
opilio3 <- readRDS("./output/opilio3.rds")

#Model convergence diagnostics 
check_hmc_diagnostics(opilio3$fit)
neff_lowest(opilio3$fit)
rhat_highest(opilio3$fit)
summary(opilio3)
bayes_R2(opilio3)

#Posterior Predictive check 
pp_check(opilio3, nsamples = 100)

#Plots
plot(conditional_smooths(opilio3), ask = FALSE)
plot(opilio3)

trace_plot(opilio3$fit) 

# model comparison
loo(opilio1, opilio2, opilio3) #Sex marginally improves model 

######################################################
# Model 4: Base model + Sex + Year

opilio4_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + sex +
                       year + (1 | year/index/station))                      

opilio4 <- brm(opilio4_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 4000, # increasing iterations 
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save Output
saveRDS(opilio4, file = "./output/opilio4.rds")
opilio4 <- readRDS("./output/opilio4.rds")

#Convergence Diagnostics 
check_hmc_diagnostics(opilio4$fit)
neff_lowest(opilio4$fit) # too low!
rhat_highest(opilio4$fit)
summary(opilio4) 
bayes_R2(opilio4)

#Plots
plot(conditional_smooths(opilio4), ask = FALSE)
plot(opilio4)
trace_plot(opilio4$fit)

#Model comparison
loo(opilio1, opilio2, opilio3, opilio4) #Including sex and year improves model 

#######################################################
#Final Model:  Run best opilio4 model with 10,000 iterations and set seed 

opiliofinal <- brm(opilio4_formula,
                   data = opilio.dat,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, iter = 10000,
                   save_pars = save_pars(all = TRUE), seed = 1, 
                   control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(opiliofinal, file = "./output/opiliofinal.rds")
opiliofinal <- readRDS("./output/opiliofinal.rds")

#Model convergence diagnostics
check_hmc_diagnostics(opiliofinal$fit)
neff_lowest(opiliofinal$fit)
rhat_highest(opiliofinal$fit)
summary(opiliofinal) 
bayes_R2(opiliofinal)

#Posterior Predictive check 
pp_check(opiliofinal, nsamples = 100) #PPC graphical check
pp_check(opiliofinal, type = "stat_grouped", stat = "mean", group = "year") #PPC for mean

#Trace plot 
trace_plot(opiliofinal$fit)
pdf("./figs/trace_opiliofinal.pdf", width = 6, height = 4)

#Area under the curve for final model 

draws_beta0 <- as.matrix(opiliofinal, variable = "b_Intercept")
logistic_beta0 <- plogis(draws_beta0)
# Summarize the posterior distribution
psych::describe(logistic_beta0)

mcmc_areas(opiliofinal, pars = "b_Intercept", 
           transformations = list("b_Intercept" = "plogis"), bw = "SJ")

# Compute AUC for predicting prevalence with the model
Prob <- predict(tanner1, type="response")
Prob <- Prob[,1]
Pred <- prediction(Prob, as.vector(pull(tanner.dat, pcr)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC 


############################################
#Save model comparison for ms.
forms <- data.frame(formula=c(as.character(opilio4_formula)[1],
                              as.character(opilio3_formula)[1],
                              as.character(opilio1_formula)[1],
                              as.character(opilio2_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/opilio_model_comp.csv")

###########################
# Plot predicted effects from best model 

opiliofinal  <- readRDS("./output/opiliofinal .rds")

# first year
ce1s_1 <- conditional_effects(opiliofinal, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__)

ggplot(plot, aes(year, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Probability positive") +
  ggtitle("Opilio Final - posterior mean & 95% credible interval")

ggsave("./figs/opiliofinal_year_effect.png", width = 2, height = 3, units = 'in')

# then sex
ce1s_1 <- conditional_effects(opiliofinal, effect = "sex", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$sex %>%
  dplyr::select(sex, estimate__, lower__, upper__)

ggplot(plot, aes(sex, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Probability positive") +
  ggtitle("Opilio Final - posterior mean & 95% credible interval")

ggsave("./figs/opiliofinal_sex_effect.png", width = 2, height = 3, units = 'in')

# then size

## 95% CI
ce1s_1 <- conditional_effects(opiliofinal, effect = "size", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opiliofinal, effect = "size", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opiliofinal, effect = "size", re_formula = NA,
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
  ggtitle("opilio Final - posterior mean & 80 / 90 / 95% credible intervals")

ggsave("./figs/opiliofinal_size_effect.png", width = 6, height = 4, units = 'in')


## cpue

## 95% CI
ce1s_1 <- conditional_effects(opiliofinal, effect = "fourth.root.cpue70", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opiliofinal, effect = "fourth.root.cpue70", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opiliofinal, effect = "fourth.root.cpue70", re_formula = NA,
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
  ggtitle("opilio Final - posterior mean & 80 / 90 / 95% credible intervals")

ggsave("./figs/opiliofinal_fourth.root.cpue70_effect.png", width = 6, height = 4, units = 'in')


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
