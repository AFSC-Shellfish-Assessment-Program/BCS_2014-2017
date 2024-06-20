# notes ----
#Data exploration: BCS infection dynamics in C. bairdi 

# Author: Mike Litzow & Erin Fedewa
# last updated: 2022/11/7

# load ----
library(tidyverse)
library(lubridate)

#load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

#colors
my_colors <- RColorBrewer::brewer.pal(7, "GnBu")[c(3,7)]

###########################################################
# data wrangling
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(species_name == "Chionoecetes bairdi",
         index_site %in% c(1, 2, 3),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_longitude, bottom_depth, 
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
         depth = as.numeric(depth)) -> tanner.dat 

#examine <70mm Tanner cpue distribution to decide on transformation                                    
ggplot(tanner.dat, aes(tanner70under_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

  #fourth root
ggplot(tanner.dat, aes(tanner70under_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

  #log
ggplot(tanner.dat, aes(log(tanner70under_cpue))) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#examine protocol target size Tanner cpue distribution to decide on transf.
ggplot(tanner.dat, aes(tannerimm_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(tanner.dat, aes(tannerimm_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

ggplot(tanner.dat, aes(log(tannerimm_cpue))) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#Transform both CPUE indices with 4th root
tanner.dat %>%
  mutate(fourth.root.cpue70 = tanner70under_cpue^0.25,
         fouth.root.cpueimm = tannerimm_cpue^0.25) -> tanner.dat 

###############################################
#Data exploration

nrow(tanner.dat) # 1285 tanner samples

#Sample sizes
tanner.dat %>%
  group_by(year, index, station) %>%
  count() %>%
  print(n=100) #Only one crab sampled at  some stations so attributing
#variance to station level will be difficult 

#total sample size of uncertain diagnosis
dat %>%
 filter(pcr_result == 3,
        year %in% c(2015:2017)) %>%
  count() #220 out of 2814 samples 

# BCS+/- occurrence by index/year
tanner.dat %>%
  select(year, index, pcr) %>%
  group_by(year, index) %>%
  summarise(PCR_0 = sum(pcr == 0),
            PCR_1 = sum(pcr == 1)) %>%
  pivot_longer(cols = c(-year, -index)) -> plot2

#% +/- stacked barplot
ggplot(plot2, aes(fill=name, y=value, x=year)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ index)

#Sample sizes by maturity
dat %>%
  #filter(maturity != "NA") %>%
  group_by(year,maturity) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(maturity), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Maturity Status", y = "Sample size") #lots of missing maturity info in 14/15

#Sample sizes by maturity/sex
dat %>%
  filter(sex != "NA") %>%
  group_by(year,maturity, sex) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(maturity), y= n), stat='identity') +
  facet_grid(sex~year) +
  theme_bw() +
  labs(x= "Maturity Status", y = "Sample size") 

#Sample sizes by shell condition/sex
dat %>%
  filter(sex != "NA") %>%
  filter(shell_cond != "NA") %>%
  group_by(year,shell_cond, sex) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(shell_cond), y= n), stat='identity') +
  facet_grid(sex~year) +
  theme_bw() +
  labs(x= "Shell Condition", y = "Sample size")
#We're missing maturity info, but prevalence of old shell 3 & 4 crab suggests that 
  #collections were not targeting immature crab only 

#Size range sampled across years
tanner.dat %>% 
  summarize(avg_cw = mean(size, na.rm=T), 
            max_cw = max(size, na.rm=T), 
            min_cw = min(size, na.rm=T))

#Size composition sampled by index site/yr
tanner.dat %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  group_by(year, index) %>%
  ggplot() +
  geom_histogram(aes(x=size, fill=Sex), position = "stack", bins=50) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Tanner crab carapace width (mm)", y = "Count")
ggsave("./figs/tanner_size.png", width=6.75)

#Combined tanner/snow figure for Fig. S1 
dat %>%
  filter(index_site %in% c(1, 2, 3, 4, 5, 6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F"), 
         species = recode_factor(species_name, "Chionoecetes bairdi" = "Tanner crab",
                                      "Chionoecetes opilio" = "Snow crab")) %>%
  group_by(species, year, index_site) %>%
  ggplot() +
  geom_density(aes(x=size, fill=Sex), position = "stack", alpha=.6) +
  scale_fill_manual(values=my_colors) +
  facet_grid(species~year) +
  theme_bw() +
  xlim(0, NA) +
  labs(x= "Carapace width (mm)", y = "Density")
ggsave("./figs/Fig S1.png")

#Size-frequency distribution of uninfected vrs infected
tanner.dat %>%
  ggplot(aes(size, fill=as.factor(pcr), color=as.factor(pcr))) +  
  geom_histogram(position="identity",alpha=0.5) +
  theme_bw()

#Size-frequency distribution of maturity status- only those w/ chela measurements! 
dat %>%
  filter(species_name == "Chionoecetes bairdi",
         index_site %in% c(1, 2, 3),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         maturity %in% c(0,1),
         pcr_result %in% c(1, 0)) %>%
  ggplot(aes(size, fill=as.factor(maturity), color=as.factor(maturity))) +  
  geom_histogram(position="identity",alpha=0.5) +
  theme_bw()

#Percent prevalence by size bin
tanner.dat %>% 
  mutate(size_bin = cut(size, breaks=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140))) %>%
  group_by(size_bin) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) %>%
  filter(size_bin != "NA") -> size
ggplot(size, aes(as.factor(size_bin), Prevalance)) +
  geom_col()

#This is something to keep in mind when interpreting changes in prev. across
#site and year- in 2016 and 2017 prevalence was very low, but is likely due to 
#most samples being taken from mature males (despite protocol specifying imm).
#Without pulling them, it's difficult to interpret changes in prevalence as 
#a true measure of the disease prev, vrs. just a factor of the subset of population 
#being sampled. 

###########################################
#Covariate data exploration

#Plot range of observed data by year/index site 
tanner.dat %>%
  group_by(year, index, station) %>%
  summarise(temperature = mean(temperature),
            depth = mean(depth),
            julian = mean(julian),
            cpue = mean(tanner70under_cpue)) -> plot 

#Temperature
ggplot(plot, aes(temperature)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

#Depth
ggplot(plot, aes(depth)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

#CPUE
ggplot(plot, aes(cpue)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)


#Julian Day 
ggplot(plot, aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index) #differences among index areas

#Julian date vrs temperature 
ggplot(plot, aes(julian, temperature)) +
  geom_point() 

#Plot explanatory variables as predictors of proportion BCS+ by year/station 
tanner.dat %>%
  group_by(year, station) %>%
  summarise(size = mean(size), 
            julian = mean(julian),
            temperature = mean(temperature),
            CPUE70 = mean(fourth.root.cpue70),
            CPUEimm = mean(fouth.root.cpueimm),
            proportion.positive = sum(pcr == 1) / n()) -> plot3

# Julian day vrs  %positive plot
ggplot(plot3, aes(julian, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") #Large seasonal effect!

#Mean size-at-station vrs %positive plot 
ggplot(plot3, aes(size, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") #Large size effect!

#Temp-at-station vrs %positive plot 
ggplot(plot3, aes(temperature, proportion.positive)) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

# <70mm CPUE-at-station vrs %positive plot 
ggplot(plot3, aes(CPUE70, proportion.positive)) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

# Immature CPUE-at-station vrs %positive plot 
ggplot(plot3, aes(CPUEimm, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") 
#Very similar to above plot- CPUE metrics likely very similar 



















