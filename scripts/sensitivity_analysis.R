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

#Calculate sample size and binomial CI's
a <- .05 #significance level
dat2 %>%
  group_by(species_name, index, year) %>%
  summarise(n = n()) %>%
  full_join(prev) %>%
  mutate(CI_upper = Prevalance + (qnorm(1-a/2))*sqrt((1/n)*Prevalance*(100-Prevalance)),
         CI_lower = Prevalance + (-qnorm(1-a/2))*sqrt((1/n)*Prevalance*(100-Prevalance))) -> prev_n

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

######################################################
#Exploration of a bayesian latent-class model for sensitivity & specificity
  #in the absence of a gold standard 

#Approach 1: JAGS model
library(readxl)
library(runjags)
library(rjags)
testjags()

Sys.setenv(JAGS_HOME = "C:/Users/erin.fedewa/AppData/Local/Programs/JAGS/JAGS-4.3.1/x64/bin/jags-terminal.exe")
 
#Reference: https://github.com/paoloeusebi/BLCM-Covid19/blob/master/covid_r1_ind_I.R
#https://academic.oup.com/aje/article/190/8/1689/6206818

#Model structure for run.jags() is a little over my head....

#Approach 2: R Shiny App for latent class models
#https://www.nandinidendukuri.com/how-to-guide-on-a-r-shiny-app-to-do-bayesian-inference-for-diagnostic-meta-analysis/
dat2 %>%
  filter(julian <= 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN")))

dat2 %>%
  filter(julian > 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN")))

tp <- c(20,36) 
fp <- c(1,1) 
fn <- c(242,477) 
tn <- c(1316,715) 
cell <- cbind(tp, fp, fn, tn)
n <- length(tp)  #  Number of studies      
write("tp fp fn tn","output/Sens_Data.txt")
for (j in 1:n) write(cell[j,],"output/Sens_Data.txt",append=T)

#Output used to run model on R shiny app 
