notes ----
  #Testing preliminary GLM models with opilio crab data
  
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
  filter(species_name == "Chionoecetes opilio",
         index_site %in% c(4, 5, 6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_latitude, bottom_depth, gear_temperature, opilio70under_cpue) %>%
  rename(pcr = pcr_result,
                station = gis_station,
                latitude = mid_latitude,
                depth = bottom_depth,
                temperature = gear_temperature,
                index = index_site,
                fourth.root.cpue70 = opilio70under_cpue) %>%
  mutate(year = as.factor(year),
                sex = as.factor(sex),
                index = as.factor(index),
                station = as.factor(station),
                fourth.root.cpue70 = fourth.root.cpue70^0.25) -> opilio.dat # transforming cpue here

nrow(opilio.dat) # 1511 samples!

#Check for missing data 
opilio.dat %>%
  select(size, sex, index, julian, latitude, depth, temperature) %>%
  filter(!complete.cases(.)) 
  summarise(num_na = sum(is.na(.)))
# Looks like one row of NA in year == 2015, station == I-23

  #Fix and overwrite 
  opilio.dat %>%
    mutate(julian = replace_na(189),
           latitude = replace_na(57.67492),
           depth = replace_na(100),
           temperature = replace_na(2.3)) %>%
    filter(!is.na(size)) -> opilio.dat

#Data exploration
opilio.dat %>%
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
opilio.dat %>%
  group_by(year, index, station) %>%
  count() %>%
  print(n=100) #Sometimes only 1 obsv within a level

################################################

#Check for correlation between continuous covariates 
opilio.dat %>%
  select(size, julian, latitude, depth, temperature, fourth.root.cpue70) %>%
  as.data.frame() -> cov
corrplot(cor(cov)) #temperature, lat and julian day are problematic 

# need dimension reduction for exogenous covariates (day, depth, lat, temperature, cpue)
pca.dat <- opilio.dat %>%
  group_by(station, year) %>%
  summarise(julian = mean(julian),
                   depth = mean(depth),
                   latitude = mean(latitude),
                   temperature = mean(temperature),
                   fourth.root.cpue70 = mean(fourth.root.cpue70))

cor(pca.dat[,3:7]) # some lower correlations than for tanner

# plot 
plot <- data.frame(Day_of_year = pca.dat$julian,
                   Depth_m = pca.dat$depth,
                   N_latitude = pca.dat$latitude,
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

ggsave("./figs/opilio_julian_temp_depth_long_cpue.png", width = 4.5, height = 7.5, units = 'in')

# DFA is hard here b/c we want to include time as one of the time series, *and* we don't have continuous observations for DFA

# following Tanner analysis - fit a PCA
# use julian, latitude, temperature
PCA <- prcomp(pca.dat[,c(3,5,6)], scale = T, center = T)
PCA$rotation
PCA$x
pca.dat$pc1 <- PCA$x[,1]

# and join pc1 back in
pc1 <- pca.dat %>%
  dplyr::select(station, year, pc1)

#Final dataset for modeling
opilio.dat <- left_join(opilio.dat, pc1)

####################################################
#Approach 1: use GLM's to assess whether prevalence varies as a function of size/pc1/year/depth/cpue 
#No random effects here- not accounting for variation across index sites/nested design

#Full model with size/pc1/year/depth/cpue
glm.full <-  glm(pcr ~ size + year + pc1 + depth + fourth.root.cpue70, data=opilio.dat,family=binomial(link="logit"))
summary(glm.full)

#Full model w/ cubic splines for nonlinear effects
glm.smooth <-  glm(pcr ~ ns(size,3) + year + ns(pc1,3) + ns(depth,3) + ns(fourth.root.cpue70,3), 
                   data=opilio.dat,family=binomial(link="logit"))
summary(glm.smooth)

#Compare linear vrs nonlinear fixed effects structure in full model 
AIC(glm.full,glm.smooth) #nonlinear fixed effects model considerably better 

#Final model using stepwise AIC
mod.glm.final <- stepAIC(glm.smooth, direction="both")
formula(mod.glm.final) #full model is best fit 
summary(mod.glm.final)

#Visualize marginal effects in final model
visreg(mod.glm.final, "pc1", scale="response") #inverse-logit transformed probability
visreg(mod.glm.final, "year", scale="response")
visreg(mod.glm.final, "size", scale="response")
visreg(mod.glm.final, "depth", scale="response")
visreg(mod.glm.final, "fourth.root.cpue70", scale="response") 

###################################################
#Approach 2: Use GLMM's to specify 1) a random structure for 
#year/site/station nested design and 2) find optimal fixed structure 

### Random effects structure:
#Model 1:  Size and PC1 as linear fixed effect, nested year-index-station as random intercept 
glmm.1 <- glmer(pcr ~ size + pc1 + (1 | year/index/station), family=binomial(link = "logit"), data=opilio.dat)
summary(glmm.1)  

#Model 2: Size and PC1 as linear fixed effect, nested index-station as random intercept 
glmm.2 <- glmer(pcr ~ size + pc1 + (1 | index/station), family=binomial(link = "logit"), data=opilio.dat)
summary(glmm.2) #Warning: random effects are very small 

anova(glmm.1,glmm.2) #year/index/site structure preferred 

###Fixed effects structure:
#Full model with all fixed effects
#NOTE: Should be using ML to estimate fixed effects but you can't specify method= in glmer?? 
glmm.full <- glmer(pcr ~ size + year + pc1 + depth + fourth.root.cpue70 + (1 | year/index/station),
                   family=binomial(link = "logit"), data=opilio.dat)
summary(glmm.full) #Convergence issues...running out of df? 

#Full model with nonlinear fixed effects
glmm.full.smooth <- glmer(pcr ~ ns(size,3) + year + ns(pc1,3) + ns(depth,3) + ns(fourth.root.cpue70,3) + (1 | year/index/station),
                          family=binomial(link = "logit"), data=opilio.dat) 
summary(glmm.full.smooth) #Convergence issues...running out of df?

#Compare fixed effects structure
anova(glmm.full, glmm.full.smooth) #Nonlinear effects favored 

#Drop depth and apply likelihood ratio test 
glmm.3 <- update(glmm.full.smooth, .~. -ns(depth,3))
anova(glmm.full.smooth, glmm.3)
summary(glmm.3)

#Drop CPUE and apply likelihood ratio test 
glmm.4 <- update(glmm.3, .~. -ns(fourth.root.cpue70,3))
anova(glmm.3, glmm.4) #Very little difference in two models 
summary(glmm.4)

#Drop size and apply likelihood ratio test 
glmm.5 <- update(glmm.4, .~. -ns(size,3))
anova(glmm.4, glmm.5) #Need to keep size in 

#Refit best model with REML....but no method argument
glmm.final <- glmer(pcr ~ ns(size, 3) + year + ns(pc1, 3) + (1 | year/index/station), family=binomial(link = "logit"), data=opilio.dat)

# Explore residuals
plot(glmm.final) #Yikes....several high leverage observations 
visreg(glmm.final)

#Visualize marginal effects in final model
visreg(glmm.final, "pc1", scale="response") #inverse-logit transformed probability
visreg(glmm.final, "size", scale="response")
visreg(glmm.final, "year", scale="response")

#Effect size 
sjPlot::plot_model(glmm.final, show.values = TRUE, value.offset = .3)
