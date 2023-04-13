# notes ----
#Objective 2: Baseline estimates of Hematodinium infection prevalence in EBS monitoring sites
  
# Author: Erin Fedewa
# last updated: 2022/9/11

# load ----
library(tidyverse)
library(patchwork)
library(ggthemes)
library(lubridate)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

###########################################
#Data wrangling 
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  select(species_name, pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth, 
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
         depth = as.numeric(depth)) -> dat2

#Size range of crab sampled by species
dat2 %>%
  group_by(species_name, year) %>%
  summarise(Min_size = min(size, na.rm=T),
            Max_size = max(size, na.rm=T),
            Avg_size = mean(size, na.rm=T))

#Calculate prevalence by year
dat2 %>%
  group_by(species_name, year) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) 

#Calculate prevalence and sample sizes within index sites
dat2 %>%
  group_by(species_name, index, year) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) -> prev

#Calculate sample size and binomial CI's
  a <- .05 #significance level
dat2 %>%
  group_by(species_name, index, year) %>%
  summarise(n = n()) %>%
  full_join(prev) %>%
  mutate(CI_upper = Prevalance + (qnorm(1-a/2))*sqrt((1/n)*Prevalance*(100-Prevalance)),
         CI_lower = Prevalance + (-qnorm(1-a/2))*sqrt((1/n)*Prevalance*(100-Prevalance))) -> prev_n
  

#We've got snow crab in the dataset for index site 2 so we'll need 
  #to filter out for plots (confirmed correct species with raw data sheets)

#Prevalence Plot by species/site/year
my_colors <- RColorBrewer::brewer.pal(7, "GnBu")[2:7]

prev_n %>%
  filter(species_name=="Chionoecetes opilio" &
           index %in% c(4:6)|
           species_name=="Chionoecetes bairdi" &
           index %in% c(1:3)) %>%
  mutate(species = case_when(species_name == "Chionoecetes opilio" ~ "Snow crab",
                             species_name == "Chionoecetes bairdi" ~ "Tanner crab")) %>%
  ggplot(aes(x=index, y=Prevalance)) +
  geom_bar(aes(fill=index, group=index),stat='identity', show.legend = FALSE) +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.2, color="grey") + 
  geom_text(aes(label=paste0("n=",n), group=index), vjust = -2.2, size=3) +
  scale_y_continuous(limits = c(0,70)) +
  labs(x="Monitoring Site", y="Prevalance (%)") +
  facet_grid(year ~ factor(species, levels = c("Tanner crab", "Snow crab")), scales = "free_x") +
  theme_bw() 
ggsave("./figs/opilio_tanner_prev.png")

#Prevelance plot by species/site/size version 2 
#tanner
tanner_colors <- RColorBrewer::brewer.pal(8, "Greens")[c(3,5,8)]

prev_n %>%
  filter(species_name=="Chionoecetes bairdi" &
           index %in% c(1:3)) %>%
  mutate(species = case_when(species_name == "Chionoecetes bairdi" ~ "Tanner crab")) %>%
  ggplot(aes(x=year, y=Prevalance, fill=index)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75) + 
  scale_fill_manual(values = tanner_colors) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.2, color="dark grey", position = position_dodge(width = .9)) + 
  geom_text(aes(year, Prevalance, label=paste0("n=",n)), size=2.4,
            position = position_dodge(width = .9), vjust=-5) +
  scale_y_continuous(limits = c(0,70)) +
  labs(x="", y="Prevalance (%)", fill="Monitoring Site", title="Tanner Crab") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(panel.grid.major = element_line()) -> tanner

#snow
snow_colors <- RColorBrewer::brewer.pal(8, "Blues")[c(3,5,8)]

prev_n %>%
  filter(species_name=="Chionoecetes opilio" &
           index %in% c(4:6)) %>%
  mutate(species = case_when(species_name == "Chionoecetes opilio" ~ "Snow crab")) %>%
  ggplot(aes(x=year, y=Prevalance, fill=index)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75) + 
  scale_fill_manual(values = snow_colors) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.2, color="dark grey", position = position_dodge(width = .9)) + 
  geom_text(aes(year, Prevalance, label=paste0("n=",n)), size=2.4,
            position = position_dodge(width = .9), vjust=-5) +
  scale_y_continuous(limits = c(0,70)) +
  labs(x="", y="Prevalance (%)", fill="Monitoring Site", title="Snow Crab") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(panel.grid.major = element_line()) -> snow

#Combine and save figure
tanner / snow
ggsave("./figs/opilio_tanner_prev2.png",dpi=300)






 
  
  











