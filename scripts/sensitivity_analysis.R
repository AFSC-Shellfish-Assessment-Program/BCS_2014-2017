# notes ----
#Sensitivity and specificity analysis between PCR and visual diagnosis

# Author: Erin Fedewa
# last updated: 2022/9/11

#Follow ups: add CI's for sensitivity and specificity 

# load ----
library(tidyverse)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

############################################
dat %>%
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  mutate(sen_spec = case_when(pcr_result==1 & visual_positive==1 ~ "TP", #true positive
                              pcr_result==0 & visual_positive==1 ~ "FP", #false positive
                              pcr_result==1 & visual_positive==0 ~ "FN", #false negative
                              pcr_result==0 & visual_positive==0 ~ "TN")) %>% #true negative
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#Wow...7% probability of diagnosing visually + when hematodinium is present 
  #99% probability of diagnosing visually neg when hemato is absent 

###############################################
#Does specificity/sensitivity differ by crab size/seasonality? May need to follow up
#on this approach, but for now, let's assess variation in sensitivity in subgroups of 
#the population (i.e. just using the median julian date/size to split data)

dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  mutate(sen_spec = case_when(pcr_result==1 & visual_positive==1 ~ "TP", #true positive
                              pcr_result==0 & visual_positive==1 ~ "FP", #false positive
                              pcr_result==1 & visual_positive==0 ~ "FN", #false negative
                              pcr_result==0 & visual_positive==0 ~ "TN")) -> dat2 #true negative
range(dat2$julian)
median(dat2$julian) #198
range(dat2$size, na.rm=T)
median(dat2$size, na.rm=T) #55mm

######Effect of Seasonality on sensitivity: 

#Crab sampled on day 198 or earlier 
dat2 %>%
  filter(julian <= 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
          FP_n = nrow(filter(.,sen_spec=="FP")),
          FN_n = nrow(filter(.,sen_spec=="FN")),
          TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#Crab sampled post day 198 
dat2 %>%
  filter(julian > 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

######Effect of crab size on sensitivity: 

#Crab <=55mm carapace width
dat2 %>%
  filter(size <= 55) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#Crab >55mm carapace width 
dat2 %>%
  filter(size > 55) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#Sample size of visually positive crab is so small that we really don't see much 
  #difference in sensitivity by date sampled or crab size 
