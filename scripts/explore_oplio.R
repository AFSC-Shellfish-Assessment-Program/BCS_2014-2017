# explore and analyze infection dynamics in C. opilio

library(tidyverse)

# load PCR data for each year
d1=read.csv("./data/bcs14.csv")
d2=read.csv("./data/bcs15.csv")
d3=read.csv("./data/bcs16.csv") 

names(d1); names(d2); names(d3)

name.check <- data.frame(names1 = names(d1),
                         names2 = names(d2),
                         names3 = names(d3))

name.check # so names are interchangeable, just some formatting differences

names(d2) <- names(d3) <- names(d1)
dat = rbind(d1, d2, d3) 

str(dat)

# set up an index column
substr(dat$Specific_Location,10,10)
dat$index = as.numeric(substr(dat$Specific_Location,10,10))

#it appears there are some crabs from outside index areas?
unique(dat$index)
idxs = is.na(dat$index)
sum(idxs) # only one!
which(idxs == T)
dat = dat[!idxs,]

dat = dat[dat$PCR_result == 0 | dat$PCR_result == 1,] # limiting to PCR positive or negative - dropping code 3 (results uncertain)

unique(dat$Species_Name)

# removing species == NA
idxs = dat$Species_Name == ""
sum(idxs) # only one

dat=dat[!idxs,]


# check clutch codes 
check <- dat %>%
  group_by(Species_Name, Clutch) %>%
  summarise(count = n())

check # that's a lot of data to throw away!


# select opilio for analysis

opilio <- dat %>%
  filter(Species_Name == "Chionoecetes opilio",
         index %in% c(4, 5, 6))

opilio$index <- as.factor(opilio$index)

# examine sample size by index / year

# now the same for index/year
plot <- opilio %>%
  dplyr::select(Year, index, PCR_result) %>%
  dplyr::group_by(Year, index) %>%
  dplyr::summarise(PCR_0 = sum(PCR_result == 0),
                   PCR_1 = sum(PCR_result == 1)) %>%
  tidyr::pivot_longer(cols = c(-Year, -index))

ggplot(plot, aes(name, value)) +
  geom_bar(stat = "identity") +
  facet_grid(Year ~ index)

ggsave("./figs/opilio_index_year_incidence.png", width = 6, height = 6, units = "in")

# plot surface temp - bottom temp relationship
ggplot(opilio, aes(SURFACE_TEMP, Bottom_Temp)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F) +
  geom_point(aes(SURFACE_TEMP, Bottom_Temp, color = index))

# very different from Tanner!
ggsave("./figs/opilio_surface_temp_bottom_temp.png", width = 6, height = 4, units = "in")

# and plot surface / bottom temps / bottom depth / dates by year / index

# first, calculate julian day
opilio$date <- stringr::str_split_fixed(opilio$START_TIME, " ", n = 2)[,1] # remove time, keep date

opilio$julian <- lubridate::yday(chron::dates(opilio$date))

plot <- opilio %>%
  dplyr::group_by(Year, index, STATIONID) %>%
  dplyr::summarise(surface = mean(SURFACE_TEMP),
                   bottom = mean(Bottom_Temp),
                   depth = mean(BOTTOM_DEPTH),
                   julian = mean(julian)) 


ggplot(plot, aes(bottom)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(Year ~ index)

ggsave("./figs/opilio_bottom_temp_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(depth)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(Year ~ index)

ggsave("./figs/opilio_depth_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(Year ~ index) # big differences among index areas

ggsave("./figs/opilio_date_by_index.png", width = 5, height = 4, units = "in")

ggplot(plot, aes(julian, bottom)) +
  geom_point() # this is a big problem - julian day is collinear with temperature

ggsave("./figs/opilio_date_vs_bottom_temp.png", width = 6, height = 4, units = "in")


# look at Julian day as a predictor of % positive
plot <- opilio %>%
  dplyr::group_by(Year, STATIONID) %>%
  dplyr::summarise(julian = mean(julian), 
                   n = dplyr::n(),
                   proportion.positive = sum(PCR_result == 1) / dplyr::n())

ggplot(plot, aes(julian, proportion.positive)) +
  geom_point() + 
  facet_wrap(~Year) +
  geom_smooth(method = "gam")

ggsave("./figs/opilio_date_vs_percent_positive.png", width = 6, height = 4, units = "in")


# and size as a predictor of % positive
plot <- opilio %>%
  dplyr::group_by(Year, STATIONID) %>%
  dplyr::summarise(size = mean(Size), 
                   n = dplyr::n(),
                   proportion.positive = sum(PCR_result == 1) / dplyr::n())

ggplot(plot, aes(size, proportion.positive)) +
  geom_point() + 
  facet_wrap(~Year) +
  geom_smooth(method = "gam")

ggsave("./figs/opilio_size_vs_percent_positive.png", width = 6, height = 4, units = "in")


# make a toy model with julian day, size, and year as predictors

toy.dat <- opilio %>%
  dplyr::group_by(Year, STATIONID) %>%
  dplyr::summarise(size = mean(Size),
                   julian = mean(julian),
                   proportion.positive = sum(PCR_result == 1) / dplyr::n())


toy.dat$Year <- as.factor(toy.dat$Year)


mod <- mgcv::gam(proportion.positive ~ Year + s(julian, k = 4) + s(size, k = 4), dat = toy.dat)
summary(mod)

plot(mod)
