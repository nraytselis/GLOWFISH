setwd("~/Desktop/Whitefish/LakeGeneva")
library(readr)
library(dplyr)
library(tidyr)
library(lubridate) #add rows for missing dates
library(zoo) #interpolate values
library(mgcv)

####Bring in the data####

####Recruitment - 2002-2018
#measures of CPUE for year 3 recruits
  Recruitment = read_csv("clean_fish_recruit_Lake_Geneva.csv") 
  
  Recruitment = Recruitment %>% mutate(Age = year_yyyy - year_class) #confirm that all data is for age 3 recruits - it is!

####Growth - 1974-2024
#I generated growth rate data based on theoretical fish under mechanistic whitefish larval bioenergetics model
#See ~/Desktop/Whitefish/LakeGeneva/Mechanistic_Whitefish_Larval_Bioenergetics_Model_LakeGeneva_8-3-26
  Growth = read_csv("Lake_Geneva_whitefish_yearly_growth.csv")

####Water Chem Parameters 
  WaterChem = read_csv("clean_water_chemistry.csv")

####Total Phosphorous - 1974-2024
  TP = WaterChem %>% filter(W_parameter == "TP") 
  
  TP_average = TP %>% group_by(year_yyyy) %>% summarise(mean(W_value))
  colnames(TP_average) = c(c("year_yyyy", "TP"))
  
####Temperature - 1974-2024
  Temperature = WaterChem %>% filter(W_parameter == "temperature") #temperature never drops below freezing (ignore winter severity?)
  
  Temp_average = Temperature %>% group_by(year_yyyy) %>% summarise(mean(W_value))
  colnames(Temp_average) = c(c("year_yyyy", "Temp_C"))
####Mussels 
  #Quagga didn't appear until 2015, recruitment data doesn't go beyond 2012 
  

####Combine datasets####
  RecruitSelect = Recruitment %>% select(c(4,18))
  colnames(RecruitSelect) = c("year_yyyy", "CPUE")

  GrowthSelect = Growth %>% select(c(1,3))
  colnames(GrowthSelect) = c("year_yyyy", "growth_mm_day")
  
  Geneva = list(RecruitSelect, GrowthSelect, TP_average, Temp_average) %>% 
    reduce(left_join, by = "year_yyyy")
  
  GenevaYCS = gam(CPUE ~ s(growth_mm_day) + s(TP), 
      data = Geneva, family = nb())  
  
  summary(GenevaYCS)
  
  plot(GenevaYCS)
  