# notes ----
#Analyze BCS infection dynamics in C. bairdi using Bayesian multivariate models
#Investigate factors that may be important in driving disease prevalence 
  #(host size/sex, depth, temperature, lat/long, immature crab density, date of sampling)

#load
library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(MARSS)
library(corrplot)
library(factoextra)
library(patchwork)
library(broom.mixed)
library(pROC)
library(tidybayes)
library(knitr)
library(loo)
library(sjPlot)
source("./scripts/stan_utils.R")

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################################
#Functions 

# Fisher Pearson skew 
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}

#Plot PPC test statistic plots
posterior_stat_plot <- function(obs, model, samples=1000, statistic="mean"){
  fig <- ppc_stat(y = obs, 
                  yrep = posterior_predict(model, ndraws = samples),  #matrix of draws from posterior distribution
                  stat = statistic)
  
  return(fig)
}

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
         fourth.root.cpue70 = tanner70under_cpue^0.25, #Transformed CPUE of tanner <70mm CW
         fouth.root.cpueimm = tannerimm_cpue^0.25) -> tanner.dat #Transformed CPUE of immature tanner (Jensen protocol cutline)

#Dimension reduction for highly correlated exogenous covariates 
tanner.dat %>%
  group_by(station, year) %>%
  summarise(julian = mean(julian),
                   depth = mean(depth),
                   longitude = -mean(longitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70)) -> pca.dat

cor(pca.dat[,3:7]) 
corrplot(cor(pca.dat[,3:7])) 

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
#Model 1: base model with crab size, pc1 and random year/index intercept 

## define model formula
tanner1_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4)  + (1 | year/index))

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

#MCMC convergence diagnostics 
check_hmc_diagnostics(tanner1$fit)
neff_lowest(tanner1$fit)
rhat_highest(tanner1$fit)
summary(tanner1)
bayes_R2(tanner1)

#Diagnostic Plots
plot(tanner1)
plot(conditional_smooths(tanner1), ask = FALSE)
mcmc_plot(tanner1, type = "areas", prob = 0.95)
mcmc_rhat(rhat(tanner1)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tanner1, pars = c("b_Intercept", "bs_ssize_1", "bs_spc1_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tanner1)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance 
pp_check(tanner1, ndraws = 100) #Not meaningful b/c pcr can only be 0 or 1

#PPC: Mean and skewness summary statistics - run functions above 
y_obs <- tanner.dat$pcr #Observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, tanner1) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,tanner1, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1 

#Histograms show distribution of the fraction of pcr outcomes given all values explored 
  #in the posterior draws. Observed results (black/red line) are near the most likely results 
  #that the model predicts- look good!

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(tanner1)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good 

# Compute AUC for predicting prevalence with the model
y_obs <- tanner.dat$pcr #observed values
preds <- posterior_epred(tanner1) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Looks like our model discriminates fairly well 

#################################################
# Model 2: Base model with crab size, julian day and random year/index intercept

tanner2_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + (1 | year/index))

tanner2 <- brm(tanner2_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner2, file = "./output/tanner2.rds")
tanner2 <- readRDS("./output/tanner2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tanner2$fit)
neff_lowest(tanner2$fit)
rhat_highest(tanner2$fit)
summary(tanner2)
bayes_R2(tanner2)

#Diagnostic Plots
plot(tanner2, ask = FALSE)
plot(conditional_smooths(tanner2), ask = FALSE)
mcmc_plot(tanner2, type = "areas", prob = 0.95)
mcmc_rhat(rhat(tanner2)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tanner2, pars = c("b_Intercept", "bs_ssize_1", "bs_sjulian_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tanner2)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check: Mean and skewness summary statistics 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, tanner2) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,tanner2, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(tanner2)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good 

# Compute AUC for predicting prevalence with the model
y_obs <- tanner.dat$pcr #observed values
preds <- posterior_epred(tanner2) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well 

# model comparison
loo(tanner1, tanner2, moment_match = TRUE) #Although tanner1 has better predictive capacity, let's move
  #forward with julian day to capture seasonal progression for ease of interpretation 

#################################################
# Model 3: Add temperature to base model 2 

tanner3_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) +  s(temperature, k = 4) + (1 | year/index))

tanner3 <- brm(tanner3_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner3, file = "./output/tanner3.rds")
tanner3 <- readRDS("./output/tanner3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tanner3$fit)
neff_lowest(tanner3$fit)
rhat_highest(tanner3$fit)
summary(tanner3)
bayes_R2(tanner3)

#Diagnostic Plots
plot(tanner3, ask = FALSE)
plot(conditional_smooths(tanner3), ask = FALSE)
mcmc_plot(tanner3, type = "areas", prob = 0.95)
mcmc_rhat(rhat(tanner3)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tanner3, pars = c("b_Intercept", "bs_ssize_1", "bs_stemperature_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tanner3)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check: Mean and skewness summary statistics 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, tanner3) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,tanner3, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(tanner3)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good 

# Compute AUC for predicting prevalence with the model
y_obs <- tanner.dat$pcr #observed values
preds <- posterior_epred(tanner3) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well 

# model comparison
loo(tanner2, tanner3, moment_match = TRUE) # Temp improves prediction, though predictive performance is very similar

###############################################################################################
# Model 4: Add CPUE (fourth-root transformed) to model 3 

tanner4_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) +  s(temperature, k = 4) +  s(fourth.root.cpue70, k = 4) + (1 | year/index))                    

tanner4 <- brm(tanner4_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output 
saveRDS(tanner4, file = "./output/tanner4.rds")
tanner4 <- readRDS("./output/tanner4.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tanner4$fit)
neff_lowest(tanner4$fit)
rhat_highest(tanner4$fit)
summary(tanner4)
bayes_R2(tanner4)

#Diagnostic Plots
plot(tanner4, ask = FALSE)
plot(conditional_smooths(tanner4), ask = FALSE)
mcmc_plot(tanner4, type = "areas", prob = 0.95)
mcmc_rhat(rhat(tanner4)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tanner4, pars = c("b_Intercept", "bs_ssize_1", "bs_sfourth.root.cpue70_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tanner4)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check: Mean and skewness summary statistics 
y_obs <- tanner.dat$pcr #observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, tanner4) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,tanner4, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(tanner4)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(tanner4) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well 

#Model comparison
loo(tanner3, tanner4, moment_match = TRUE) #CPUE doesn't improve model predictions, though very similar  
###############################################################################################
#Model 5: Add sex to model 3 

tanner5_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) +  s(temperature, k = 4) +  sex + (1 | year/index))                     

tanner5 <- brm(tanner5_formula,
                              data = tanner.dat,
                              family =bernoulli(link = "logit"),
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner5, file = "./output/tanner5.rds")
tanner5 <- readRDS("./output/tanner5.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tanner5$fit)
neff_lowest(tanner5$fit)
rhat_highest(tanner5$fit)
summary(tanner5)
bayes_R2(tanner5)

#Diagnostic Plots
plot(tanner5, ask = FALSE)
plot(conditional_smooths(tanner5), ask = FALSE)
mcmc_plot(tanner5, type = "areas", prob = 0.95)
mcmc_rhat(rhat(tanner5)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tanner5, pars = c("b_Intercept", "bs_ssize_1", "b_sex2"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tanner5)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check: Mean and skewness summary statistics 
y_obs <- tanner.dat$pcr #observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, tanner5) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,tanner5, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(tanner5)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(tanner5) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well

#Model comparison
loo(tanner3, tanner5, moment_match = TRUE) #no improvement by adding sex, but models very similar 

#########################################
#Full Model Comparison

#LOO-CV
model.comp <- loo(tanner2, tanner3, tanner4, tanner5, moment_match=TRUE)
  model.comp
#All models very similar: Cross-validation is uncertain whether estimating additional 
  #parameters improves predictive performance 
  
#Table of Rsq Values 
rbind(bayes_R2(tanner2), 
        bayes_R2(tanner3), 
        bayes_R2(tanner4), 
        bayes_R2(tanner5)) %>%
  as_tibble() %>%
  mutate(model = c("tanner2", "tanner3", "tanner4", "tanner5"),
           r_square_posterior_mean = round(Estimate, digits = 2)) %>%
  select(model, r_square_posterior_mean) 

#Model weights 
  #PSIS-LOO
loo1 <- loo(tanner2)
loo2 <- loo(tanner3)
loo3 <- loo(tanner4)
loo4 <- loo(tanner5)

loo_list <- list(loo1, loo2, loo3, loo4)

#Compute and compare Pseudo-BMA weights without Bayesian bootstrap, 
  #Pseudo-BMA+ weights with Bayesian bootstrap, and Bayesian stacking weights
stacking_wts <- loo_model_weights(loo_list, method="stacking")
pbma_BB_wts <- loo_model_weights(loo_list, method = "pseudobma")
pbma_wts <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
round(cbind(stacking_wts, pbma_wts, pbma_BB_wts),2)
#Tanner3 consistently highest weighted model...enough rationale to select?
  
#Save model output 
tab_model(tanner2, tanner3, tanner4, tanner5)

forms <- data.frame(formula=c(as.character(tanner3_formula)[1],
                              as.character(tanner4_formula)[1],
                              as.character(tanner5_formula)[1],
                              as.character(tanner2_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/tanner_model_comp.csv")

##########################################
#Exploring model averaging/stacking to incorporate all models for prediction/propagating uncertainty
  #b/c all models are so similar and LOO-CV is less accurate for separating between small effect sizes 

# Predictions from model stacking by Yao et al. (2018)
pred_stacking <- pp_average(tanner2, tanner3, tanner4, tanner5, method = "posterior_predict")
  #Stacking is the default method here... 

# Prediction from pseudo-Bayesian model averaging (BMA)
# 1. Obtain predictions from each model
pred_m2345 <- map(list(tanner2, tanner3, tanner4, tanner5), posterior_predict)
# 2. Obtain model weights based on Pseudo-BMA (with Bayesian bootstrap)
pbma_wts <- loo_model_weights(tanner2, tanner3, tanner4, tanner5, method = "pseudobma") 
# 3. Obtain weighted predictions
pred_pbma <- map2(pred_m2345, pbma_wts, `*`) %>% 
  reduce(`+`) %>% 
  posterior_summary()
# Compare the weights
ggplot(tibble(stacking = pred_stacking[ , "Estimate"], 
              pbma = pred_pbma[ , "Estimate"]), aes(x = pbma, y = stacking)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x="Pseudo-Bayesian Model Averaging", y="Model Stacking") #two methods give very similar predictions

###############################
#Explore tanner3 vrs tanner4....top two models (base + temp, base + temp + CPUE) 

#1) Compare LOO ELPC pointwise comparisons

#Compute LOO-CV for each model 
loo1 <- loo(tanner3, save_psis = TRUE, cores = 2)
loo2 <- loo(tanner4, save_psis = TRUE, cores = 2)

# Obtain pointwise ELPD values
elpd1 <- loo1$pointwise[,"elpd_loo"]
elpd2 <- loo2$pointwise[,"elpd_loo"]

# Build differences dataframe
elpd_df <- tibble(diff12 = elpd2 - elpd1) %>% 
  mutate(idx = 1:n())

# Plot LOO ELPD Point-wise Model Comparisons 
ggplot(elpd_df, aes(x = idx, y = diff12)) +
  geom_point(alpha=0.7) +
  geom_hline(aes(yintercept=0)) +
  theme_bw() + 
  labs(x = "Index", y = "ELPD Difference",
       title = "Tanner3 - Tanner4 ") +
  theme(title = element_text(size=8))
#Positive ELPD differences would indicate tanner3 is a better performing model...
  #but no apparent trend in differences b/w models across the range of observations 

#2) Compare Pareto k diagnostics

# Get khat values
k_tanner3 <- loo1$psis_object$diagnostics$pareto_k
k_tanner4 <- loo2$psis_object$diagnostics$pareto_k

# Plot values
tibble(idx = seq_along(k_tanner3), 
       Tanner3 = k_tanner3,
       Tanner4 = k_tanner4) %>% 
  pivot_longer(cols = -idx) %>% 
  ggplot(., aes(x = idx, y = value, color = name)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0)) +
  theme_bw() +
  theme(legend.title= element_blank()) +
  ylim(-1,1) +
  facet_wrap(~ name) +
  labs(y = expression(hat(k)), x ="Observation Index") +
  ggtitle("Pareto-k diagnostic (PSIS diagnostic)")

#No influence points (all khat <0.7). One khat value > 0.5 in tanner3

###########################
#Below is a preliminary final model using tanner3, though all models tested are very similar
  #Should we be using predictions from model averaging, or is there enough rationale for tanner3? 

#Final Model:  Run tanner3 model with 10,000 iterations and set seed for reproducibility 
tannerfinal <- brm(tanner3_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 10000,
               save_pars = save_pars(all = TRUE), seed = 3, 
               control = list(adapt_delta = 0.9999, max_treedepth = 14))

#Save model output 
saveRDS(tannerfinal, file = "./output/tannerfinal.rds")
tannerfinal <- readRDS("./output/tannerfinal.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tannerfinal$fit)
neff_lowest(tannerfinal$fit)
rhat_highest(tannerfinal$fit)
summary(tannerfinal)
bayes_R2(tannerfinal)

#Diagnostic Plots
plot(tannerfinal, ask = FALSE)
plot(conditional_smooths(tannerfinal), ask = FALSE)
mcmc_plot(tannerfinal, type = "areas", prob = 0.95)
mcmc_rhat(rhat(tannerfinal)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tannerfinal, pars = c("b_Intercept", "bs_ssize_1", "bs_stemperature_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tannerfinal)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check: Mean and skewness summary statistics 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, tannerfinal) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,tannerfinal, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(tannerfinal)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good

# Compute AUC for predicting prevalence with the model
y_obs <- tanner.dat$pcr
preds <- posterior_epred(tannerfinal)
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Looks like our model discriminates fairly well 

#Get posterior estimates with tidybayes 
get_variables(tannerfinal)

post_beta_means <- tannerfinal %>% 
  gather_draws(b_Intercept, bs_ssize_1, bs_spc1_1,
               sd_year__Intercept, `sd_year:index__Intercept`, 
               `sd_year:index:station_Intercept`, sds_ssize_1, sds_spc1_1) %>% 
  median_hdi() #Huh, stumped on how to use spread/gather_draws with multilevel models, 
              #Might have to come back to this one 

################################
#Plot predicted effects from best model...should be much easier with tidybayes

#Size
## 95% CI
ce1s_1 <- conditional_effects(tannerfinal , effect = "size", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tannerfinal , effect = "size", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tannerfinal , effect = "size", re_formula = NA,
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

##Julian Day
## 95% CI
ce1s_1 <- conditional_effects(tannerfinal , effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tannerfinal , effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tannerfinal , effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Julian Day", y = "Probability positive") +
  ggtitle("Tanner Final - posterior mean & 80 / 90 / 95% credible intervals")
ggsave("./figs/tannerfinal_julian_effect.png", width = 6, height = 4, units = 'in')

##Temperature 
## 95% CI
ce1s_1 <- conditional_effects(tannerfinal , effect = "temperature", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tannerfinal , effect = "temperature", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tannerfinal , effect = "temperature", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$temperature
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temperature[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temperature[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temperature[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temperature[["lower__"]]

ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Temperature (C)", y = "Probability positive") +
  ggtitle("Tanner Final - posterior mean & 80 / 90 / 95% credible intervals")
ggsave("./figs/tannerfinal_temp_effect.png", width = 6, height = 4, units = 'in')

#####################################################
#Predictions for final model ...need to follow up on this 

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

