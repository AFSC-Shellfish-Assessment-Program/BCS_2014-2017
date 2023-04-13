# notes ----
#Manuscript Objective 1:compare testing accuracy of diagnostic methods
#Sensitivity and specificity analysis between PCR and visual diagnosis

# Author: Erin Fedewa
# last updated: 2023/1/3

# load ----
library(tidyverse)
library(lubridate)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

#Load color palette
my_colors <- RColorBrewer::brewer.pal(7, "GnBu")[4:7]

############################################
#Sample sizes
dat %>%
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  summarize(total = n(),
            pcr_pos = sum(pcr_result==1),
            vis_pos = sum(visual_positive==1))

#Calculate sensitivity and specificity
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

#7% probability of diagnosing visually + when PCR is +
  #99% probability of diagnosing visually neg when hemato is absent 

###############################################
#Does specificity/sensitivity differ by crab size/seasonality? 
#let's assess variation in sensitivity in subgroups of 
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
#Plot Classification Accuracy  

#All data
dat %>%
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  mutate(sen_spec = case_when(pcr_result==1 & visual_positive==1 ~ "TP", #true positive
                              pcr_result==0 & visual_positive==1 ~ "FP", #false positive
                              pcr_result==1 & visual_positive==0 ~ "FN", #false negative
                              pcr_result==0 & visual_positive==0 ~ "TN")) %>% #true negative
  summarise(TP_n = sum(sen_spec=="TP")/n()*100,
            FP_n = sum(sen_spec=="FP")/n()*100,
            FN_n = sum(sen_spec=="FN")/n()*100,
            TN_n = sum(sen_spec=="TN")/n()*100) %>%
  mutate(Test = print("All Samples")) -> total

#Crab sampled on day 198 or earlier 
dat2 %>%
  filter(julian <= 198) %>%
  summarise(TP_n = sum(sen_spec=="TP")/n()*100,
            FP_n = sum(sen_spec=="FP")/n()*100,
            FN_n = sum(sen_spec=="FN")/n()*100,
            TN_n = sum(sen_spec=="TN")/n()*100) %>%
  mutate(Test = print("June-July Samples")) -> pre

#Crab sampled post day 198 
dat2 %>%
  filter(julian > 198) %>%
  summarise(TP_n = sum(sen_spec=="TP")/n()*100,
            FP_n = sum(sen_spec=="FP")/n()*100,
            FN_n = sum(sen_spec=="FN")/n()*100,
            TN_n = sum(sen_spec=="TN")/n()*100) %>%
  mutate(Test = print("July-Aug Samples")) -> post 
  
#Combine datasets and plot 
total %>%
  rbind(post) %>%
  rbind(pre) %>%
  rename("True Positive" = "TP_n", "False Positive" = "FP_n",
         "True Negative" = "TN_n", "False Negative" = "FN_n") %>%
  pivot_longer(c(1:4), names_to="Classification", values_to = "Percent") %>%
  mutate(Diagnosis = case_when(Classification=="False Negative" | Classification=="True Positive" ~ "Infected", 
                               Classification=="False Positive" | Classification=="True Negative" ~ "Not Infected")) %>%
  ggplot(aes(Diagnosis, Percent, fill=Classification)) +
    geom_bar(position = "stack", stat="identity") +
    theme_bw() +
    scale_fill_manual(values = my_colors) +
    facet_wrap(~factor(Test, levels=c("All Samples", "June-July Samples", "July-Aug Samples")))
ggsave("./figs/sens_spec.png", width=7)


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
#https://cran.r-project.org/web/packages/runjags/vignettes/quickjags.html

#Model structure for run.jags() is a little over my head....side burner for now

#Approach 2: R Shiny App for latent class models
#https://www.nandinidendukuri.com/how-to-guide-on-a-r-shiny-app-to-do-bayesian-inference-for-diagnostic-meta-analysis/
dat2 %>%
  filter(julian <= 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) -> early

dat2 %>%
  filter(julian > 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) -> late

#Set up data structure to feed to R shiny model 
tp <- c(early$TP_n, late$TP_n) 
fp <- c(early$FP_n, late$FP_n) 
fn <- c(early$FN_n, late$FN_n) 
tn <- c(early$TN_n, late$TN_n) 
cell <- cbind(tp, fp, fn, tn)
n <- length(tp)  #  Number of studies      
write("tp fp fn tn","output/Sens_Data.txt")
for (j in 1:n) write(cell[j,],"output/Sens_Data.txt",append=T)

#Output used to run model on R shiny app 
