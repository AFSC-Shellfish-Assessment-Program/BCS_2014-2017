# notes ----
#Analyze BCS infection dynamics in C. bairdi using Bayesian multivariate models
#Investigate factors that may be important in driving disease occurrence 
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

#Check for missing data for PCA
tanner.dat %>%
  select(size, sex, year, index, station, julian, longitude, depth, temperature) %>%
  filter(!complete.cases(.)) #Looks good 

#Assess collinearity b/w covariates 
tanner.dat %>%
  group_by(station, year) %>%
  summarise(julian = mean(julian),
            depth = mean(depth),
            longitude = -mean(longitude),
            temperature = mean(temperature),
            fourth.root.cpue70 = mean(fourth.root.cpue70)) -> pca.dat

cor(pca.dat[,3:7]) # depth, longitude and julian day highly correlated
corrplot(cor(pca.dat[,3:7]), method = 'number') 

#Plot 
pca.dat %>%
  rename("Depth (m)" = depth,
         "Longitude (W)" = longitude,
         "Bottom Temperature (C)" = temperature,
         "Fourth root CPUE" = fourth.root.cpue70) %>%
  pivot_longer(4:7, names_to = "variable", values_to = "data") %>% 
  ggplot(aes(julian, data)) +
  geom_point(aes(color = as.factor(year))) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = T, alpha = 0.2, 
              color = "black", lwd = 0.3) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position="none") +
  labs(x= "Day of Year") -> tanner_plot

#Combine plots for Fig 2 of ms (run "analyze_opilio.R" lines 1-124 first)
tanner_plot + snow_plot + plot_annotation(tag_levels = 'a')
ggsave("./figs/exog_variables.png")

#Dimension reduction for depth/long/day using PCA
pca.dat %>%
  ungroup() %>%
  select(julian, longitude, depth) %>%
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
                               column == "longitude" ~ "Longitude",
                               column == "depth" ~ "Depth")) %>%
  select(-PC, -column) %>%
  ggplot(aes(covariate, value)) +
  geom_bar(stat='identity') +
  ylab("Loading") + xlab("") +
  theme_bw() -> plot2

#Figure S1b for MS
plot1 + plot2

#Extract pc1 for model runs and join to opilio dataset
pca.dat$pc1 <- PCA$x[,1] 

pc1 <- pca.dat %>%
  select(station, year, pc1)

tanner.dat <- left_join(tanner.dat, pc1) #Final dataset for modeling 

## fit Tanner model
dat$year_fac <- as.factor(dat$year)
tanner_year_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(temperature, k = 4) + year + (1 | index))

tanner_year <- brm(tanner_year_formula,
               data = tanner.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner_year, file = "./output/tanner_year.rds")
tanner_year <- readRDS("./output/tanner_year.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tanner_year$fit)
neff_lowest(tanner_year$fit)
rhat_highest(tanner_year$fit)
summary(tanner_year)
bayes_R2(tanner_year)

plot_tanner <- conditional_effects(tanner_year, effect = "year")
