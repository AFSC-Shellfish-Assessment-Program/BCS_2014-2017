# notes ----
#Testing preliminary GLMM models with tanner crab infection status/covariates

# Author: Erin Fedewa
# last updated: 2022/1/20

# load ----
library(tidyverse)
library(lubridate)
library(corrplot)
library(MASS)
library(lme4)
library(car)
library(visreg)
library(splines)
library(sjPlot)

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

#############################

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")
  head(dat)
  
# data wrangling
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(species_name == "Chionoecetes bairdi",
                index_site %in% c(1, 2, 3),
                year %in% c(2015:2017),
                sex %in% c(1, 2),
                pcr_result %in% c(1, 0)) %>%
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth, gear_temperature, tanner70under_cpue) %>%
  rename(pcr = pcr_result,
                station = gis_station,
                longitude = mid_longitude,
                depth = bottom_depth,
                temperature = gear_temperature,
                index = index_site,
                fourth.root.cpue70 = tanner70under_cpue) %>%
  mutate(year = as.factor(year),
                sex = as.factor(sex),
                index = as.factor(index),
                station = as.factor(station),
                depth = as.numeric(depth),
                fourth.root.cpue70 = fourth.root.cpue70^0.25) -> tanner.dat # transforming cpue here

nrow(tanner.dat) # 1285 samples!

#Data exploration
tanner.dat %>%
  group_by(year, station) %>%
  summarise(size = mean(size), 
            temperature = mean(temperature),
            CPUE = mean(fourth.root.cpue70),
            proportion.positive = sum(pcr == 1) / n()) -> plot

#Size plot
ggplot(plot, aes(size, proportion.positive)) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

#Temperature plot
ggplot(plot, aes(temperature, proportion.positive)) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

#CPUE plot 
ggplot(plot, aes(CPUE, proportion.positive)) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

#Sample sizes
tanner.dat %>%
  group_by(year, index, station) %>%
  count() #Sometimes only one crab sampled at a station

################################################

#Check for correlation between continuous covariates 
tanner.dat %>%
  select(size, julian, longitude, depth, temperature, fourth.root.cpue70) %>%
  as.data.frame() -> cov
  corrplot(cor(cov)) #depth, longitude and julian day are problematic 

# need dimension reduction for exogenous covariates (day, depth, longitude, temperature, cpue)
tanner.dat %>%
  group_by(station, year) %>%
  summarise(julian = mean(julian),
                   depth = mean(depth),
                   longitude = -mean(longitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70)) -> pca.dat

cor(pca.dat[,3:7]) # temp no longer correlated with others, cpue weakly collinear with depth

# plot 
plot <- data.frame(Day_of_year = pca.dat$julian,
                   Depth_m = pca.dat$depth,
                   W_longitude = pca.dat$longitude,
                   Bottom_temperature_C = pca.dat$temperature,
                   Fourth_root_CPUE_70 = pca.dat$fourth.root.cpue70,
                   year = as.numeric(as.character(pca.dat$year))) %>%
  pivot_longer(cols = c(-Day_of_year, -year))

ggplot(plot, aes(Day_of_year, value)) +
  geom_point(color = "grey100") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, color = "black", lwd = 0.3) +
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  geom_point(aes(color = as.factor(year))) +
  scale_color_manual(values = cb[c(2,4,6)]) +
  theme(legend.title = element_blank())

# DFA is hard here b/c we want to include time as one of the time series, *and* 
  #we don't have continuous observations for DFA

# Fit a PCA
PCA <- prcomp(pca.dat[,3:5], scale = T, center = T)
summary(PCA) #PC1 explains most of variation
PCA$rotation #Variable loadings 
PCA$x #Value of rotated data 
pca.dat$pc1 <- PCA$x[,1]

# and join pc1 back in
pca.dat %>%
  dplyr::select(station, year, pc1) -> pc1

tanner.dat <- left_join(tanner.dat, pc1)

####################################################
#Approach 1: use GLM's to assess whether occurrence varies as a function of size/pc1/year/site/cpue 
  #No random effects here- not accounting for variation across index sites/nested design

#Full model for spatiotemporal effects: across year and across season (depth/day/long)
glm.full <-  glm(pcr ~ size + year + pc1 + temperature + fourth.root.cpue70, data=tanner.dat,family=binomial(link="logit"))
summary(glm.full)

#Drop CPUE
glm.1 <- update(glm.full, .~. -fourth.root.cpue70)
summary(glm.1)

#Drop Temperature
glm.2 <- update(glm.1, .~. -temperature)
summary(glm.2)

#Final model using stepwise AIC
mod.glm.final <- stepAIC(glm.full, direction="both")
formula(mod.glm.final)
summary(mod.glm.final)

#Visualize marginal effects in final model
visreg(mod.glm.final, "pc1", scale="response") #inverse-logit transformed probability
visreg(mod.glm.final, "year", scale="response")
visreg(mod.glm.final, "size", scale="response")

#Full model w/ cubic splines for nonlinear effects
glm.smooth <-  glm(pcr ~ ns(size,3) + year + ns(pc1,3) + ns(temperature,3) + ns(fourth.root.cpue70,3), 
              data=tanner.dat,family=binomial(link="logit"))
summary(glm.smooth)

#Final model using stepwise AIC
mod.glm.smooth <- stepAIC(glm.smooth, direction="both")
formula(mod.glm.smooth)
summary(mod.glm.smooth)

#Compare linear vrs nonlinear fixed effects structure in full model 
AIC(glm.full,glm.smooth) 

###################################################
#Approach 2: Use GLMM's to specify 1) a random structure for 
  #year/site/station nested design and 2) find optimal fixed structure 

### Random effects structure:
#Model 1:  Size and PC1 as linear fixed effect, nested year-index as random intercept 
glmm.1 <- glmer(pcr ~ size + pc1 + (1 | year/index), family=binomial(link = "logit"), data=tanner.dat)
  summary(glmm.1) #Warning: random effects are very small 

###Fixed effects structure:
#Full model with all fixed effects
    #NOTE: Should be using ML to estimate fixed effects but you can't specify method= in glmer?? 
glmm.full <- glmer(pcr ~ size + year + pc1 + temperature + fourth.root.cpue70 + (1 | year),
                   family=binomial(link = "logit"), data=tanner.dat)
summary(glmm.full) #Convergence issues...running out of df? 

#Full model with nonlinear fixed effects
glmm.full.smooth <- glmer(pcr ~ ns(size,3) + year + ns(pc1,3) + ns(temperature,3) + ns(fourth.root.cpue70,3) + (1 | year/index),
                   family=binomial(link = "logit"), data=tanner.dat) 
summary(glmm.full.smooth) #Convergence issues...running out of df?

#Compare fixed effects structure
anova(glmm.full, glmm.full.smooth) #Linear fixed effects favored 

#Drop temperature and apply likelihood ratio test 
glmm.3 <- update(glmm.full, .~. -temperature)
anova(glmm.full, glmm.3)
summary(glmm.3)

#Drop CPUE and apply likelihood ratio test 
glmm.4 <- update(glmm.3, .~. -fourth.root.cpue70)
anova(glmm.3, glmm.4) #Very little difference in two models 
summary(glmm.4)

#Drop year and apply likelihood ratio test 
glmm.5 <- update(glmm.4, .~. -year)
anova(glmm.4, glmm.5)
summary(glmm.5)

#Refit best model with REML....but no method argument
glmm.final <- glmer(pcr ~ size + pc1 + (1 | year/index), family=binomial(link = "logit"), data=tanner.dat)
  
# Explore residuals
plot(glmm.final) #Yikes....several high leverage observations 
visreg(glmm.final)

#Visualize marginal effects in final model
visreg(glmm.final, "pc1", scale="response") #inverse-logit transformed probability
visreg(glmm.final, "size", scale="response")

#Effect size 
sjPlot::plot_model(glmm.final, show.values = TRUE, value.offset = .3)

#Moving to Bayesian hierarchical modeling approach 



