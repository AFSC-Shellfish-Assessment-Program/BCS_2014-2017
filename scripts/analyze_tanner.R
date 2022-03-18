# notes ----
#Analyze BCS infection dynamics in C. bairdi using Bayesian multivariate models 

# Authors: Mike Litzow & Erin Fedewa
# last updated: 2022/2/28

#load
library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(MARSS)
library(corrplot)
library(factoextra)
library(ROCR)
source("./scripts/stan_utils.R")

# set plot theme
theme_set(theme_bw())

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

########################################
#Data Manipulation

# data wrangling
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
         fourth.root.cpue70 = tanner70under_cpue^0.25,
         fouth.root.cpueimm = tannerimm_cpue^0.25) -> tanner.dat 

#Dimension reduction for exogenous covariates 
tanner.dat %>%
  group_by(station, year) %>%
  summarise(julian = mean(julian),
                   depth = mean(depth),
                   longitude = -mean(longitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70)) -> pca.dat

cor(pca.dat[,3:7]) 
corrplot(cor(pca.dat[,3:7])) #temp no longer correlated with others, cpue weakly collinear with depth

#Plot 
pca.dat %>%
  rename(`Day of Year` = julian,
         `Depth (m)` = depth,
         `Longitude (W)` = longitude,
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
ggsave("./figs/tanner_julian_temp_depth_long_cpue.png", width = 4.5, height = 7.5, units = 'in')

#Fit a PCA with all exogenous covariates
PCA <- prcomp(pca.dat[,3:7], scale = T, center = T)
PCA$rotation #Variable loadings

#PCA result plots 
fviz_eig(PCA) #Scree plot: PC1 explains ~60% of variance 
fviz_pca_var(PCA, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)     
get_eig(PCA)

#Refit PCA with just long/depth/julian day and extract pc1 only for models 
PCA2 <- prcomp(pca.dat[,3:5], scale = T, center = T)
PCA2$rotation #Variable loadings
fviz_eig(PCA2) #Scree plot
PCA2$x
pca.dat$pc1 <- PCA2$x[,1]

#Join pc1 back in
pc1 <- pca.dat %>%
  dplyr::select(station, year, pc1)

#Final dataset for modeling 
tanner.dat <- left_join(tanner.dat, pc1)

#####################################################
#Model 1: base model with size, pc1 and random year/index/station intercept 

## define model formula
tanner1_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4)  + (1 | year/index/station))

## Show default priors
get_prior(tanner1_formula, tanner.dat, family = bernoulli(link = "logit"))

## fit binomial model 1
tanner1 <- brm(tanner1_formula,
                         data = tanner.dat,
                         family = bernoulli(link = "logit"),
                         cores = 4, chains = 4, iter = 2500,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output
saveRDS(tanner1, file = "./output/tanner1.rds")
tanner1 <- readRDS("./output/tanner1.rds")

#Model convergence diagnostics 
check_hmc_diagnostics(tanner1$fit)
neff_lowest(tanner1$fit)
rhat_highest(tanner1$fit)
summary(tanner1)
bayes_R2(tanner1)

#Plots
plot(conditional_smooths(tanner1), ask = FALSE)
plot(tanner1)
trace_plot(tanner1$fit) 
mcmc_plot(tanner1, type = "areas", prob = 0.95)


#Compute posterior predictive probabilities to assess predictive performance 


#Posterior Predictive check 
pp_check(tanner1, nsamples = 100)

###################################################
# Model 2: Add temperature to base model 1 

tanner2_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(temperature, k = 4) + (1 | year/index/station))

tanner2 <- brm(tanner2_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner2, file = "./output/tanner2.rds")
tanner2 <- readRDS("./output/tanner2.rds")

#Convergence diagnostics
check_hmc_diagnostics(tanner2$fit)
neff_lowest(tanner2$fit)
rhat_highest(tanner2$fit)
summary(tanner2) 
bayes_R2(tanner2)

#Plots
plot(conditional_smooths(tanner2), ask = FALSE)
plot(tanner2)
trace_plot(tanner2$fit)

# model comparison
loo(tanner1, tanner2, moment_match = TRUE) # Temp does NOT improve prediction

###############################################################################################
# Model 3: Add CPUE to base model 

tanner3_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(fourth.root.cpue70, k = 4) + (1 | year/index/station))                      

tanner3 <- brm(tanner3_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output 
saveRDS(tanner3, file = "./output/tanner3.rds")
tanner3 <- readRDS("./output/tanner3.rds")

#Convergence Diagnostics 
check_hmc_diagnostics(tanner3$fit)
neff_lowest(tanner3$fit)
rhat_highest(tanner3$fit)
summary(tanner3) # no evidence of a CPUE effect
bayes_R2(tanner3)
plot(conditional_smooths(tanner3), ask = FALSE)

#Plots
plot(conditional_smooths(tanner3), ask = FALSE)
plot(tanner3)
trace_plot(tanner3$fit)

#Model comparison
loo(tanner1, tanner2, tanner3, moment_match = TRUE) #Same elpd b/w tanner3 and base model 
#No difference between predictive power of the two models so let's go with most 
  #parsimonious base model 

###############################################################################################
#Model 4: Add sex to base model 

tanner4_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + sex + (1 | year/index/station))                      

tanner4 <- brm(tanner4_formula,
                              data = tanner.dat,
                              family =bernoulli(link = "logit"),
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner4, file = "./output/tanner4.rds")
tanner4 <- readRDS("./output/tanner4.rds")

#Convergence diagnostics 
check_hmc_diagnostics(tanner4$fit)
neff_lowest(tanner4$fit)
rhat_highest(tanner4$fit)
summary(tanner4) # no evidence of a sex effect
bayes_R2(tanner4)

#Plots
plot(conditional_smooths(tanner4), ask = FALSE)
plot(tanner4)
trace_plot(tanner4$fit)

#Model comparison
loo(tanner1, tanner2, tanner3, tanner4, moment_match = TRUE) #no improvement by adding sex 

######################################################
# Model 5: Add year effect to base model 

tanner5_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + year + (1 | year/index/station))                      

tanner5 <- brm(tanner5_formula,
               data = tanner.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500, 
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner5, file = "./output/tanner5.rds")
tanner5 <- readRDS("./output/tanner5.rds")

#Convergence Diagnostics 
check_hmc_diagnostics(tanner5$fit)
neff_lowest(tanner5$fit) # too low!
rhat_highest(tanner5$fit)
summary(tanner5) # no evidence of a year effect
bayes_R2(tanner5)

#Plots
plot(conditional_smooths(tanner5), ask = FALSE)
plot(tanner5)
trace_plot(tanner5$fit)

#Model comparison
loo(tanner1, tanner2, tanner3, tanner4, tanner5, moment_match = TRUE) # tanner5 marginally the best

##########################################
#Model 6:  Compare best population-level model (tanner5) with index/station group-level term 

tanner6_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + year + (1 | index/station)) 

tanner6 <- brm(tanner6_formula,
                 data = tanner.dat,
                 family = bernoulli(link = "logit"),
                 cores = 4, chains = 4, iter = 2500,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(tanner6, file = "./output/tanner6.rds")
tanner6 <- readRDS("./output/tanner6.rds")

#Model convergence diagnostics
check_hmc_diagnostics(tanner6$fit)
neff_lowest(tanner6$fit)
rhat_highest(tanner6$fit)
summary(tanner6) 
bayes_R2(tanner6)

#Model comparison between random effects structures
loo(tanner5, tanner6, moment_match = TRUE) #year/site/station group level term much better 

#########################################
#Full model comparison

#All models 
  #Not including model 6 here
model.comp <- loo(tanner1, tanner2, tanner3, tanner4, tanner5, moment_match=TRUE)
  model.comp
pbma_wts <- loo_model_weights(tanner1, tanner2, tanner3, tanner4, tanner5, method = "pseudobma")
  pbma_wts #tanner 5 gets highest weight 

#Save model output 
sjPlot::tab_model(tanner1, tanner2, tanner3, tanner4, tanner5)

forms <- data.frame(formula=c(as.character(tanner5_formula)[1],
                              as.character(tanner3_formula)[1],
                              as.character(tanner1_formula)[1],
                              as.character(tanner2_formula)[1],
                              as.character(tanner4_formula)[1],
                              as.character(tanner6_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/tanner_model_comp.csv")

##########################################
#Exploring model averaging/stacking to incorporate all models for prediction/propagating
  #uncertainty b/c all models are so similar 

# Prediction from model stacking by Yao et al. (2018)
pred_stacking <- pp_average(tanner1, tanner2, tanner3, tanner4, tanner5, method = "predict")

# Prediction from pseudo BMA
# 1. Obtain predictions from each model
pred_m12345 <- map(list(tanner1, tanner2, tanner3, tanner4, tanner5), posterior_predict)
# 2. Obtain model weights based on Pseudo-BMA (with Bayesian bootstrap)
pbma_wts <- loo_model_weights(tanner1, tanner2, tanner3, tanner4, tanner5, method = "pseudobma") #can change to method="waic"
# 3. Obtain weighted predictions
pred_pbma <- map2(pred_m12345, pbma_wts, `*`) %>% 
  reduce(`+`) %>% 
  posterior_summary()
# Compare the weights
ggplot(tibble(stacking = pred_stacking[ , "Estimate"], 
              pbma = pred_pbma[ , "Estimate"]), aes(x = pbma, y = stacking)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) #two methods give very similar predictions

###########################
#Final Model:  Run tanner1 model with 10,000 iterations and set seed 

tannerfinal <- brm(tanner1_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 10000,
               save_pars = save_pars(all = TRUE), seed = 1, 
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(tannerfinal, file = "./output/tannerfinal.rds")
tannerfinal <- readRDS("./output/tannerfinal.rds")

#Model convergence diagnostics
check_hmc_diagnostics(tannerfinal$fit)
neff_lowest(tannerfinal$fit)
rhat_highest(tannerfinal$fit)
summary(tannerfinal) 
bayes_R2(tannerfinal)

#Posterior Predictive check 
pp_check(tannerfinal, nsamples = 100) #PPC graphical check
pp_check(tannerfinal, type = "stat_grouped", stat = "mean", group = "year") #PPC for mean

#Trace plot 
trace_plot(tannerfinal$fit)
pdf("./figs/trace_tannerfinal.pdf", width = 6, height = 4)

#Area under the curve for final model 

draws_beta0 <- as.matrix(tannerfinal, variable = "b_Intercept")
logistic_beta0 <- plogis(draws_beta0)
# Summarize the posterior distribution
psych::describe(logistic_beta0)

mcmc_areas(tannerfinal, pars = "b_Intercept", 
           transformations = list("b_Intercept" = "plogis"), bw = "SJ")

# Compute AUC for predicting prevalence with the model
Prob <- predict(tanner1, type="response")
Prob <- Prob[,1]
Pred <- prediction(Prob, as.vector(pull(tanner.dat, pcr)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC #Looks good!

################################
#Plot predicted effects from best model 
tannerfinal <- readRDS("./output/tannerfinal.rds")

#Size
## 95% CI
ce1s_1 <- conditional_effects(tanner1 , effect = "size", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tanner1 , effect = "size", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tanner1 , effect = "size", re_formula = NA,
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
  ggtitle("Tanner Final - posterior mean & 80 / 90 / 95% credible intervals")
ggsave("./figs/tannerfinal_size_effect.png", width = 6, height = 4, units = 'in')

##PC1 (day of year, depth, longitude)
## 95% CI
ce1s_1 <- conditional_effects(tannerfinal , effect = "pc1", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tannerfinal , effect = "pc1", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tannerfinal , effect = "pc1", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$pc1
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$pc1[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$pc1[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$pc1[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$pc1[["lower__"]]

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "PC1 (day of year, depth, longitude)", y = "Probability positive") +
  ggtitle("Tanner Final - posterior mean & 80 / 90 / 95% credible intervals")
ggsave("./figs/tannerfinal_pc1_effect.png", width = 6, height = 4, units = 'in')

##Year
ce1s_1 <- conditional_effects(tannerfinal, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__)

ggplot(plot, aes(year, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Probability positive") +
  ggtitle("Tanner Final - posterior mean & 95% credible interval")
ggsave("./figs/tannerfinal_year_effect.png", width = 6, height = 4, units = 'in')

#####################################################
#Predictions for best model 

#Predict for size < 60mm
new.dat <- tanner.dat %>%
  dplyr::filter(size < 60)

##Predict overall prevalence

# first, create a combined year_index_station column to identify combinations that exist in the data
tanner.dat$yr_ind_st <- paste(tanner.dat$year, tanner.dat$index, tanner.dat$station, sep = "_")


new.dat <- data.frame(yr_ind_st = unique(tanner.dat$yr_ind_st),
                      size = 30,
                      pc1 = mean(unique(tanner.dat$pc1))
                      )

new.dat$year <- map_chr(str_split(new.dat$yr_ind_st, "_"), 1)
new.dat$index <- map_chr(str_split(new.dat$yr_ind_st, "_"), 2)
new.dat$station <- map_chr(str_split(new.dat$yr_ind_st, "_"), 3)

posterior.predict <- posterior_epred(tannerfinal, newdata = new.dat)

tanner.estimate <- data.frame(species = "tanner",
                              estimate = mean(posterior.predict),
                              lower_95 = quantile(posterior.predict, probs = 0.025),
                              upper_95 = quantile(posterior.predict, probs = 0.975),
                              lower_90 = quantile(posterior.predict, probs = 0.05),
                              upper_90 = quantile(posterior.predict, probs = 0.95),
                              lower_80 = quantile(posterior.predict, probs = 0.1),
                              upper_80 = quantile(posterior.predict, probs = 0.9))

ggplot(tanner.estimate) +
  aes(x = species, y = estimate) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), color = "grey80") +
  geom_errorbar(aes(ymin = lower_90, ymax = upper_90), color = "grey60") +
  geom_errorbar(aes(ymin = lower_80, ymax = upper_80), color = "black") +
  geom_point(size = 3, color = "red3") +
  theme_classic()

ggsave("./figs/tanner estimate.png", width = 2, height = 2.5, units = 'in')

############################################################
#Additional models tested:

#Tanner7- Base model plus tanner imm_cpue (Pam's cutoff) vrs <70mm CPUE 
tanner7_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(fouth.root.cpueimm, k = 4) + (1 | year/index/station))

tanner7 <- brm(tanner7_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(tanner7, file = "./output/tanner7.rds")
tanner7 <- readRDS("./output/tanner7.rds")

#Model convergence diagnostics
check_hmc_diagnostics(tanner7$fit)
neff_lowest(tanner7$fit)
rhat_highest(tanner7$fit)
summary(tanner7) 
bayes_R2(tanner7)

#Model comparison
loo(tanner1, tanner3, tanner5, tanner7, moment_match=TRUE)

############
#Tanner 8: Base model plus tanner CPUE 70mm UNTRANSFORMED

tanner8_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + s(tanner70under_cpue , k = 4) + (1 | year/index/station))

tanner8 <- brm(tanner8_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Divergent transitions after multiple model run attempts 

