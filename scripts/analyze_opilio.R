# notes ----
#Objective 3: 
#Investigate drivers of Hematodinium infection in snow crab using Bayesian multivariate models  
  #(i.e. host size/sex, depth, temperature, lat/long, immature crab density, date of sampling)

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
library(factoextra)
library(patchwork)
library(modelr)
library(broom.mixed)
library(pROC)
library(ggthemes)
library(tidybayes)
library(RColorBrewer)
library(knitr)
library(loo)
library(sjPlot)
source("./scripts/stan_utils.R")


# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

# color palettes
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
my_colors <- RColorBrewer::brewer.pal(7, "GnBu")[c(3,5,7)]

########################################
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

nrow(opilio.dat) #1510 samples

#Check for missing data for PCA
opilio.dat %>%
  select(size, sex, year, index, station, julian, latitude, depth, temperature) %>%
  filter(!complete.cases(.)) #Looks good  

#Assess collinearity b/w covariates 
opilio.dat %>%
  group_by(year, station) %>%
  summarise(julian = mean(julian),
                   depth = mean(depth),
                   latitude = mean(latitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70)) -> pca.dat.opilio

cor(pca.dat.opilio[,3:7]) #Temperature, latitude, CPUE and julian day r > 0.6
corrplot(cor(pca.dat.opilio[,3:7]), method = 'number') 

#Assess the Variance Inflation Factor of each predictor 
model <- lm(pcr ~ julian + latitude + temperature + snow70under_cpue, data=opilio.dat ) 
summary(model) 
vif(model) #CPUE VIF fairly low (<2.5)

#Given that latitude, temp and Julian day are all intracorrelated and 
#CPUE/temp correlation may just be spurious, let's drop CPUE from PCA and test as
#a covariate in the Bayesian regression models 

# plot exogenous variables
pca.dat.opilio %>%
  rename("Depth (m)" = depth,
         "Latitude (N)" = latitude,
         "Bottom Temperature (C)" = temperature,
         "Snow Crab Density (CPUE)" = fourth.root.cpue70) %>%
  pivot_longer(4:7, names_to = "variable", values_to = "data") %>% 
  ggplot(aes(julian, data)) +
  geom_point(aes(color = as.factor(year))) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = T, alpha = 0.2, 
              color = "black", lwd = 0.3) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(x= "Day of Year") -> snow_plot
#See analyze_tanner.R script to combine plots for Fig 2 of ms

#Dimension reduction for temp/lat/day using PCA
pca.dat.opilio %>%
  ungroup() %>%
  select(julian, latitude, temperature) %>%
  prcomp(scale = T, center = T) -> PCA

get_eig(PCA)
fviz_eig(PCA) #Scree plot: PC1 explains ~82% of variance 
fviz_pca_var(PCA, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE) + labs(title="") -> plot1

#Extract rotation matrix and plot loadings of pc1
PCA %>%
  tidy(matrix = "rotation") %>%
  filter(PC == 1) %>%
  mutate(covariate = case_when(column == "julian" ~ "Julian day",
                               column == "latitude" ~ "Latitude",
                               column == "temperature" ~ "Temperature")) %>%
  select(-PC, -column) %>%
  ggplot(aes(covariate, value)) +
    geom_bar(stat='identity') +
  ylab("Loading") + xlab("") +
  theme_bw() -> plot2

#Figure S1b for MS
plot1 + plot2

#Extract pc1 for model runs and join to opilio dataset
pca.dat.opilio$pc1 <- PCA$x[,1] 

pc1 <- pca.dat.opilio %>%
  select(station, year, pc1)

opilio.dat <- left_join(opilio.dat, pc1) #Final dataset for modeling 

####################################################
#Model 1: base model with size, pc1 and random year/index intercept 

## define model formula
opilio1_formula <-  bf(pcr ~ s(size, k = 4) + s(pc1, k = 4) + (1 | year/index)) 

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

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio1$fit)
neff_lowest(opilio1$fit)
rhat_highest(opilio1$fit)
summary(opilio1)
bayes_R2(opilio1)

#Diagnostic Plots
plot(opilio1, ask = FALSE)
plot(conditional_smooths(opilio1), ask = FALSE)
mcmc_plot(opilio1, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio1)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio1, pars = c("b_Intercept", "bs_ssize_1", "bs_spc1_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio1)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
y_obs <- opilio.dat$pcr #Observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio1) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio1, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1 

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio1)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio1) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well (though should preferably be >0.80)

#Final diagnostic: Correct classification rate
#calculate the predicted probabilities of 0/1 in the original data from the fitted model
  Pred <- predict(opilio1, type = "response")
  Pred <- if_else(Pred[,1] > 0.5, 1, 0)
ConfusionMatrix <- table(Pred, pull(opilio.dat, pcr)) #`pull` results in a vector
#correct classification rate
  sum(diag(ConfusionMatrix))/sum(ConfusionMatrix) #Model correctly classifies 77% of observations
  ConfusionMatrix

###################################################
#Model 2: Base model with more parsimonious Julian day vrs pc1 

opilio2_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + (1 | year/index)) 

opilio2 <- brm(opilio2_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio2, file = "./output/opilio2.rds")
opilio2 <- readRDS("./output/opilio2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio2$fit)
neff_lowest(opilio2$fit)
rhat_highest(opilio2$fit)
summary(opilio2)
bayes_R2(opilio2)

#Diagnostic Plots
plot(opilio2, ask = FALSE)
plot(conditional_smooths(opilio2), ask = FALSE)
mcmc_plot(opilio2, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio2)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio2, pars = c("b_Intercept", "bs_ssize_1", "bs_sjulian_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio2)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
y_obs <- opilio.dat$pcr #Observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio2) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio2, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1 

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio2)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio2) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well 

#Final diagnostic: Correct classification rate (similar to above)
#calculate the predicted probabilities of 0/1 in the original data from the fitted model
Pred <- predict(opilio2, type = "response")
Pred <- if_else(Pred[,1] > 0.5, 1, 0)
ConfusionMatrix <- table(Pred, pull(opilio.dat, pcr)) #`pull` results in a vector
#correct classification rate
sum(diag(ConfusionMatrix))/sum(ConfusionMatrix) #Model correctly classifies 78% of observations
ConfusionMatrix

# model comparison
loo(opilio1, opilio2) #More parsimonious model with Julian day improves predictive accuracy

#Smooths of pc1 vrs Julian
conditional_smooths(opilio1, effects = "pc1")
conditional_smooths(opilio2, effects = "julian")
#Smooths look very similar, i.e. Julian day is capturing the pc1 pattern, and we're
  #not losing much by dropping temperature and latitude. Let's use julian day
  #as base model and add covariates from there 

###############################################################################################
#Model 3: Add CPUE to base model 2
#Not testing an opilio model with temperature b/c highly correlated with julian day
 
opilio3_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(fourth.root.cpue70, k = 4) + (1 | year/index))                       

opilio3 <- brm(opilio3_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio3, file = "./output/opilio3.rds")
opilio3 <- readRDS("./output/opilio3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio3$fit)
neff_lowest(opilio3$fit)
rhat_highest(opilio3$fit)
summary(opilio3)
bayes_R2(opilio3)

#Diagnostic Plots
plot(opilio3, ask = FALSE)
plot(conditional_smooths(opilio3), ask = FALSE)
mcmc_plot(opilio3, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio3)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio3, pars = c("b_Intercept", "bs_ssize_1", "bs_sfourth.root.cpue70_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio3)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio3) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio3, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1  

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio3)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio3) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well 

# model comparison
loo(opilio2, opilio3) #CPUE improves predictive capacity of base model, though fairly
#small improvement in SE suggests models are very similar 

######################################################
# Model 4: Add sex to model 2 

opilio4_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(fourth.root.cpue70, k = 4) + sex + (1 | year/index))                      

opilio4 <- brm(opilio4_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500, 
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save Output
saveRDS(opilio4, file = "./output/opilio4.rds")
opilio4 <- readRDS("./output/opilio4.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio4$fit)
neff_lowest(opilio4$fit)
rhat_highest(opilio4$fit)
summary(opilio4)
bayes_R2(opilio4)

#Diagnostic Plots
plot(opilio4, ask = FALSE)
plot(conditional_smooths(opilio4), ask = FALSE)
mcmc_plot(opilio4, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio4)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio4, pars = c("b_Intercept", "bs_ssize_1", "b_sex2"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio4)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio4) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio4, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1  

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio4)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio4) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well

#Model comparison
loo(opilio3, opilio4)  #Opilio4 marginally better, although again- models very similar 

######################################################
# Model 5: Add depth to model 4

opilio5_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(fourth.root.cpue70, k = 4) 
                       + sex + s(depth, k = 4) + (1 | year/index))                      

opilio5 <- brm(opilio5_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500, 
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save Output
saveRDS(opilio5, file = "./output/opilio5.rds")
opilio5 <- readRDS("./output/opilio5.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio5$fit)
neff_lowest(opilio5$fit)
rhat_highest(opilio5$fit)
summary(opilio5)
bayes_R2(opilio5)

#Diagnostic Plots
plot(opilio5, ask = FALSE)
plot(conditional_smooths(opilio5), ask = FALSE)
mcmc_plot(opilio5, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio5)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio5, pars = c("b_Intercept", "bs_ssize_1", "b_sex2"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio5)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio5) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio5, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1  

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio5)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio5) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well

#Model comparison
loo(opilio4, opilio5)  #Opilio5 marginally better 

############################################
#Full Model Comparison

#LOO-CV
model.comp <- loo(opilio2, opilio3, opilio4, opilio5, moment_match=TRUE)
model.comp #Model 5 marginally best 

#Table of Rsq Values 
rbind(bayes_R2(opilio2), 
      bayes_R2(opilio3), 
      bayes_R2(opilio4),
      bayes_R2(opilio5)) %>%
  as_tibble() %>%
  mutate(model = c("opilio2", "opilio3", "opilio4", "opilio5"),
         r_square_posterior_mean = round(Estimate, digits = 2)) %>%
  select(model, r_square_posterior_mean) 

#Model weights 
#PSIS-LOO
loo2 <- loo(opilio2)
loo3 <- loo(opilio3)
loo4 <- loo(opilio4)
loo5 <- loo(opilio5)

loo_list <- list(loo2, loo3, loo4, loo5)

#Compute and compare Pseudo-BMA weights without Bayesian bootstrap, 
#Pseudo-BMA+ weights with Bayesian bootstrap, and Bayesian stacking weights
stacking_wts <- loo_model_weights(loo_list, method="stacking")
pbma_BB_wts <- loo_model_weights(loo_list, method = "pseudobma")
pbma_wts <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
round(cbind(stacking_wts, pbma_wts, pbma_BB_wts),2)
#Opilio5 has highest weight across all 3 methods 

#Save model output  
tab_model(opilio2, opilio3, opilio4, opilio5)

forms <- data.frame(formula=c(as.character(opilio5_formula)[1],
                              as.character(opilio4_formula)[1],
                              as.character(opilio3_formula)[1],
                              as.character(opilio2_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/snow_model_comp.csv")

#######################################################
#Final Model:  Run best opilio5 model with 10,000 iterations and set seed 

opiliofinal <- brm(opilio5_formula,
                   data = opilio.dat,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, iter = 10000,
                   save_pars = save_pars(all = TRUE), seed = 1, 
                   control = list(adapt_delta = 0.9999, max_treedepth = 14))

#Save model output 
saveRDS(opiliofinal, file = "./output/opiliofinal.rds")
opiliofinal <- readRDS("./output/opiliofinal.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opiliofinal$fit)
neff_lowest(opiliofinal$fit)
rhat_highest(opiliofinal$fit)
summary(opiliofinal)
bayes_R2(opiliofinal)

#Diagnostic Plots
plot(opiliofinal, ask = FALSE)
plot(conditional_smooths(opiliofinal), ask = FALSE)
mcmc_plot(opiliofinal, type = "areas", prob = 0.95) #all credible intervals but sex contain zero 
mcmc_rhat(rhat(opiliofinal)) #Potential scale reduction: All rhats < 1.1
mcmc_neff(neff_ratio(opiliofinal)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check: Mean and skewness summary statistics 
y_obs <- opilio.dat$pcr #observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opiliofinal) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opiliofinal, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opiliofinal)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opiliofinal) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well 

###########################
#Extract and plot conditional effects of each predictor from best model
#conditioning on the mean for all other predictors, yr/site effects ignored 

#Sex effect plot 
#Need to save settings from conditional effects as an object to plot in ggplot
ce1s_1 <- conditional_effects(opiliofinal, effect = "sex", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$sex %>%
  dplyr::select(sex, estimate__, lower__, upper__) %>%
  mutate(sex = case_when(sex == 1 ~ "Male",
                         sex == 2 ~ "Female")) %>%
ggplot(aes(factor(sex, levels = c("Male", "Female")), estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  labs(y="Probability of infection", x="") +
  theme_bw() -> sexplot

#Size
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

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = size, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Carapace width (mm)", y = "") +
  theme_bw() -> sizeplot

#Julian Day
## 95% CI
ce1s_1 <- conditional_effects(opiliofinal, effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opiliofinal, effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opiliofinal, effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = julian, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Julian Day", y = "Probability of infection") +
  theme_bw() -> dayplot

#Depth
## 95% CI
ce1s_1 <- conditional_effects(opiliofinal, effect = "depth", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opiliofinal, effect = "depth", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opiliofinal, effect = "depth", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$depth
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$depth[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$depth[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$depth[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$depth[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = depth, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Depth (m)", y = "Probability of infection") +
  theme_bw() -> depthplot

#CPUE
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

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = fourth.root.cpue70, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Snow crab CPUE", y = "Probability of infection") +
  theme_bw() -> cpueplot

#Combine plots for Fig 7 of MS
(sexplot + sizeplot) / (dayplot + depthplot) / (cpueplot + plot_spacer()) +
  plot_annotation(tag_levels = 'a')
ggsave("./figs/opilioFig7.png", height=9)

###################################
#Extract Marginal Effects: instantaneous slope of one explanatory value with all 
  #other values held constant 

#Marginal effect at the mean: julian day 
opiliofinal %>%
  emtrends(~ julian, 
           var = "julian", 
           regrid = "response", re_formula = NA)
#on average, a one-day increase in Julian day is associated with a 3.5% increase in 
#the probability of infection

#Marginal effect at various levels of julian day 
opiliofinal %>% 
  emtrends(~ julian, var = "julian",
           at = list(julian = 
                       seq(min(opilio.dat$julian), 
                           max(opilio.dat$julian), 1)),
           re_formula = NA) %>%
  as_tibble() %>%
#and plot 
ggplot(aes(x = julian, y = julian.trend )) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.1) +
  geom_line(size = 1) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Julian Day", y = "Marginal effect of julian day on probability of infection") +
  theme_bw() 

#Marginal effect at the mean: sex
opiliofinal %>%
  emmeans(~ sex, 
           var = "sex", 
           regrid = "response", re_formula = NA)
#Females have a 9.3% increase in the probability of infection than males 

# Posterior predictions by sex
grand_mean <- opiliofinal %>% 
  epred_draws(newdata = expand_grid(sex = c(1, 2),
                                    julian = mean(opilio.dat$julian),
                                    depth = mean(opilio.dat$depth),
                                    size = mean(opilio.dat$size), 
                                    fourth.root.cpue70 = mean(opilio.dat$fourth.root.cpue70)), 
              re_formula = NA)

#plot
ggplot(grand_mean, aes(x = .epred, y = "Grand mean", fill = sex)) +
  stat_halfeye() +
  labs(x = "Predicted probability of infection", y = NULL,
       fill = "Sex",
       subtitle = "Posterior predictions") +
  theme_bw() 

######################################################
#Generating posterior predictions for final model 

#global julian mean-ignoring year/site specific deviations 
grand_mean <- opiliofinal %>% 
  #create dataset across a range of observed julian days sampled
  epred_draws(newdata = expand_grid(julian = range(opilio.dat$julian),
                                    depth = mean(opilio.dat$depth),
                                    size = mean(opilio.dat$size), 
                                    fourth.root.cpue70 = mean(opilio.dat$fourth.root.cpue70),
                                    sex = levels(opilio.dat$sex)),  
              re_formula = NA) #ignoring random effects 
#plot
ggplot(grand_mean, aes(x = julian, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Julian day", y = "Probability of infection",
       fill = "Credible interval") +
  theme_bw() +
  theme(legend.position = "bottom")

##########
#Year-specific posterior predictions across day 
all_years <- opiliofinal %>% 
  epred_draws(newdata = expand_grid(julian = range(opilio.dat$julian),
                                    depth = mean(opilio.dat$depth),
                                    size = mean(opilio.dat$size), 
                                    fourth.root.cpue70 = mean(opilio.dat$fourth.root.cpue70),
                                    sex = levels(opilio.dat$sex), 
                                    year = levels(opilio.dat$year)), 
              re_formula = ~ (1 | year)) #only predict using yr effects, not site too 

ggplot(all_years, aes(x = julian, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Day of year", y = "Probability of Infection",
       fill = "Credible interval") +
  facet_wrap(vars(year)) +
  theme_bw() +
  theme(legend.position = "bottom")

#average marginal effect by year
all_years_ame <- opiliofinal %>% 
  emtrends(~ julian + year,
           var = "julian",
           at = list(year = levels(opilio.dat$year)),
           epred = TRUE, re_formula = ~ (1 | year)) %>% 
  gather_emmeans_draws()

ggplot(all_years_ame,aes(x = .value)) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a 1 day increase in sampling",
       y = "Density") +
  facet_wrap(~year) +
  theme_bw()

#post and interval summaries of draws from size effect 
all_years_ame %>% median_hdi()
#Probability of infection varies by year 

###########################################################
#Lastly, to compare annual probability of infection, lets use best model
#with year as a fixed effect 

## fit snow model
snow_year_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(fourth.root.cpue70, k = 4) 
                         + sex + s(depth, k = 4) + year + (1 | index))


snow_year <- brm(snow_year_formula,
                   data = opilio.dat,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, iter = 2500,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(snow_year, file = "./output/snow_year.rds")
opilio_year <- readRDS("./output/snow_year.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio_year$fit)
neff_lowest(opilio_year$fit)
rhat_highest(opilio_year$fit)
summary(opilio_year)
bayes_R2(opilio_year)

#Diagnostic Plots
plot(opilio_year, ask = FALSE)
plot(conditional_smooths(opilio_year), ask = FALSE)
mcmc_plot(opilio_year, prob = 0.95)
mcmc_plot(opilio_year, transformations = "inv_logit_scaled")

#Conditional Effect 
conditional_effects(opilio_year, effect = "year")

ce1s_1 <- conditional_effects(opilio_year, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Snow crab") -> year_snow

#Average marginal effect of year 
years_ame <- opilio_year %>% 
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

#######################################################################
#Additional exploratory model runs for determine base model structure
#Model 2.5: Base model with Julian day plus temperature vrs pc1 

opilio2.5_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(temperature, k = 4) + (1 | year/index)) 

opilio2.5 <- brm(opilio2.5_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio2.5, file = "./output/opilio2.5.rds")
opilio2.5 <- readRDS("./output/opilio2.5.rds")

#MCMC convergence diagnostics 
summary(opilio2.5)
bayes_R2(opilio2.5)

#Diagnostic Plots
plot(opilio2.5, ask = FALSE)
plot(conditional_smooths(opilio2.5), ask = FALSE)

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio2.5)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio2.5) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well 

# model comparison
loo(opilio2, opilio2.5) #Not much difference, though looks like temp slightly improves predictive accuracy

#Smooths of Julian vrs Julian with temp
plot(conditional_smooths(opilio2, effects = "pc1"), ask=FALSE)
plot(conditional_smooths(opilio2.5, effects = "julian"), ask=FALSE)
#Julian day smooths look very similar, i.e. We're not losing much by dropping temperature 

