# notes ----
#Testing preliminary GLM models with tanner crab data

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
  dplyr::select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth, gear_temperature, tanner70under_cpue) %>%
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
  dplyr::select(size, julian, longitude, depth, temperature, fourth.root.cpue70) %>%
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

# DFA is hard here b/c we want to include time as one of the time series, *and* we don't have continuous observations for DFA

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
#Modeling 

#plot logit transformation
#Size plot
ggplot(tanner.dat, aes(logit(size), logit(pcr))) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

#Temperature plot
ggplot(tanner.dat, aes(logit(temperature), logit(pcr))) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

#CPUE plot 
ggplot(tanner.dat, aes(logit(fourth.root.cpue70), logit(pcr))) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

#####################################################

#Approach 1: use GLM's to assess whether prevalence varies as a function of size/pc1/year/site/cpue 
  #No random effects here- not accounting for variation across index sites/nested design

#Full model for spatiotemporal effects: across year and across season (depth/day/long)
glm.1 <-  glm(pcr ~ size + year + pc1, data=tanner.dat,family=binomial(link="logit"))
summary(glm.1)

#Full model for spatiotemporal effects w/ cubic splines
glm.2 <-  glm(pcr ~ ns(size,3) + year + ns(pc1,3), data=tanner.dat,family=binomial(link="logit"))
summary(glm.2)
AIC(glm.1,glm.2) #Let's keep fixed effects linear moving forward 

#Add temperature
glm.3 <-  glm(pcr ~ size + year + pc1 +temperature, data=tanner.dat,family=binomial(link="logit"))
summary(glm.3)
AIC(glm.1,glm.3) #temperature doesn't improve model 

#Add CPUE
glm.4 <-  glm(pcr ~ size + year + pc1 + fourth.root.cpue70, data=tanner.dat,family=binomial(link="logit"))
summary(glm.4)
AIC(glm.1,glm.4) #CPUE doesn't improve model 

#Visualize final model
visreg(glm.1, "pc1", scale="response")
visreg(glm.1, "year", scale="response")
visreg(glm.1, "size", scale="response")

###################################################
#Approach 2: Use GLMM's to specify a random effect for year/site/station nested design

#Model 1:  Size as linear fixed effect, nested year-index-station as random intercept 
glmm.1 <- glmer(pcr ~ size + (1 | index/station), family=binomial(link = "logit"), data=tanner.dat)
summary(glmm.1)
#Model does not converge, not enough df due to low per station sample size

#Model 1:  Size as linear fixed effect, nested year-index-station as random intercept 
glmm.1 <- glmer(pcr ~ size + (1 | year/index/station), family=binomial(link = "logit"), data=tanner.dat)
  #Model does not converge, not enough df due to low per station sample size 

#Model 2: Site as random effect only 
glmm.2 <- glmer(pcr ~ size + (1 | index), family=binomial(link = "logit"), data=tanner.dat)
summary(glmm.2)

#Model 3: Add year
glmm.3 <- glmer(pcr ~ size + year + (1 | index), family=binomial(link = "logit"), data=tanner.dat)
summary(glmm.3)

#Add pc1
glmm.4 <- glmer(pcr ~ size + year + pc1 + (1 | index), family=binomial(link = "logit"), data=tanner.dat)
summary(glmm.4)


#Does random effect improve model?
m.glm <- glm(pcr ~ size + pc1 + year, data=tanner.dat,family=binomial(link="logit"))
summary(m.glm)
m.glmm <- glmer(pcr ~ size + pc1 + year + (1 | index), family=binomial(link = "logit"), data=tanner.dat)
summary(m.glmm)
AIC(m.glm, m.glmm)

####Do the same for year!!!


#don't assume PC1 represents seasonal - 1 of factors is driving but don't know
  #vrs using index site instead of pc1 variation: could be due to a host of different things ____
####Try models with site instead of pc1 
#If using PC1, frame as site differences due to depth/day sampled/log...vrs
  #differences in index site might be due to temp, depth, 

#Is seasonality due to increased detection as disease progresses vrs transmission (more crab getting it)
    #Email Hamish 
#include 2014 but look at patterns vrs true 

#Cite morado paper (very low prevalance at lg sizes) so used non-linear effect 


#spline terms
#random slope?--convergence issues.... no random effect? 
#random effects-partioning variance....but interested in year as fixed effect and don't really care
  #about station 

#Future: assess visual vrs PCR from historic data---if follows same pattern then have entire dataset to run covariates 
  #Lab exp/temp effects paper
  #Long timeseries BCS drivers (visual, PCR)
  #recent declines fit to BCS
  #NBS fit into this: 2017-2019


m1 <- glmer(pcr ~ size + year + (1 | index), family=binomial(link = "logit"), data=tanner.dat)
summary(m1)
m1 <- glmer(pcr ~ size + pc1 + year + (1 | index), family=binomial(link = "logit"), data=tanner.dat)
summary(m1)
#Does the slope of index site vary based on size of crab sampled there?
m3 <- glmer(pcr ~ size + pc1 + year + (1 + size | index), family=binomial(link = "logit"), data=tanner.dat)
#model doesn't converge 

#Random rationale: don't care about variation across index site becasue that's captured by pc1
# and interest lies in underlying pop prevalance, not variation in prev across sites 
  #year is of interest though
#Treating factors with small numbers of levels as random will in the best case lead to very small and/or imprecise estimates of random effects; 
#in the worst case it will lead to various numerical difficulties such as lack of convergence, zero variance estimates


#Model 1:  Size as linear fixed effect, nested year-index-station as random intercept 
#Model 1:  Size as linear fixed effect, nested year-index-station as random intercept 
m1 <- glmer(pcr ~ size + (1 | year/index/station), family=binomial, data=tanner.dat)



#Model 1: size, PC1 and random effect for year-index-station
model1.PQL <- glmmPQL(pcr ~ size + pc1, random= ~1 | year/index/station, data=tanner.dat, family=binomial)
  summary(model1.PQL)

model1.lme4 <- lmer(pcr ~ size + pc1 + (1 | year/index/station), data=tanner.dat, family=binomial)
summary(model1.lme4)

spline and random effect 

yearXstation
size spline, station fixed effect, random effect 
-Need to look at prev vrs size rxn-does this warrant nonlinear rxn
    keep size as fixed/linear then non linear 
Always going to get better fits with more parameters-i.e. spline

















