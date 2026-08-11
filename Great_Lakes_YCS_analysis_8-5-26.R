library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("~/Desktop/Whitefish/GreatLakes")

recruitmentData = read_csv("GLOWFISH_SuperiorHuronMichiganErieOntarioSimcoe_Brown.csv")
TP = read_csv("TP_GLNPO_Clean.csv") #total phosphorous
WS = read_csv("Percent_Ice_Cover_Great_Lakes.csv") #winter severity - % ice cover data 

TP$Season = as.factor(TP$Season)
####Separate water bodies into new dfs

###YCS - recruitment data for each of the Great Lakes + Lake Simcoe 
Lake_Ontario = recruitmentData %>% filter(waterbody_name == "Ontario_Brown")

Lake_Simcoe = recruitmentData %>% filter(waterbody_name == "Simcoe_Brown") #no TP or WS data currently 

Lake_Erie = recruitmentData %>% filter(waterbody_name == "Erie_Brown")

Lake_Huron = recruitmentData %>% filter(waterbody_name == "Huron_Brown")

Lake_Michigan = recruitmentData %>% filter(waterbody_name == "Michigan_Brown")

Lake_Superior = recruitmentData %>% filter(waterbody_name == "Superior_Brown")

###Total Phosphorous (data are for Spring and Summer, taken across a few days, and cover multiple stations)
TP_Ontario = TP %>% filter(Lake == "Ontario")
  
TP_Erie = TP %>% filter(Lake == "Erie")
  
TP_Huron = TP %>% filter(Lake == "Huron")

TP_Michigan = TP %>% filter(Lake == "Michigan")

TP_Superior = TP %>% filter(Lake == "Superior")


###Winter Severity (% Ice cover)

WS_Ontario = WS[, c(1,10,11)]

WS_Erie = WS[, c(1,8,9)]

WS_Huron = WS[, c(1,6,7)]
  
WS_Michigan = WS[, c(1,4,5)]
  
WS_Superior = WS[, c(1:3)]


####YCS - over time####

recruitmentData <- recruitmentData[!is.na(recruitmentData$data_name), ]

ggplot(recruitmentData, aes(x = year_class, y = fish_abundance_value, group = data_name, color = data_name)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ data_name) + theme_minimal() + labs(x = "Year Class", y = "YCS")


####Look at how much variation in Total Phosphorous there is across sampling sites 
#to see if we can take an average across all sites for a given lake/year

####Lake Ontario####
length(unique(TP_Ontario$Station)) #24 different stations 

TP_Ontario_Summary = TP_Ontario %>% group_by(Year,Season) %>% summarise(AcrossLakeMeanTP = mean(Value), #mean across multiple stations
                                                                 SD = sd(Value),
                                                                 n = n(),
                                                                 SE_TP = sd(Value)/sqrt(n))
ggplot(TP_Ontario_Summary,
       aes(x = Year,
           y = AcrossLakeMeanTP, group = Season, color = Season)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = AcrossLakeMeanTP - SE_TP,
    ymax = AcrossLakeMeanTP + SE_TP)) + theme_minimal() + labs(y = "Mean TP across all stations", title = "Lake Ontario")


ggplot(WS_Ontario,aes(x = Year, y = Ontario_Cover)) + geom_point() + geom_line() + theme_minimal() +
  labs(y = "% Ice Cover", title = "Lake Ontario")


####Key takeaways for Lake Ontario TP####
      ####total phosphorous has decreased over time across all stations for both spring and summer 
          #for now 

####Combine dataframes to do analysis####

Lake_Ontario = Lake_Ontario %>% mutate(Year = year_class) #copy over year class column into new year column 

Lake_Ontario_EnvioData = merge(TP_Ontario_Summary,WS_Ontario,by = "Year")

Lake_Ontario_alldata = merge(Lake_Ontario_EnvioData,Lake_Ontario, by = "Year")

Lake_Ontario_alldata_select = Lake_Ontario_alldata %>% select(c(1,2,3,6,7,23,24,25))


#GAM: YCS ~ Larval Growth + Mussels + TP + Winter Severity
library(mgcv)


Ontario_gam <- gam(fish_abundance_value ~ s(AcrossLakeMeanTP) + s(Ontario_Cover) + s(Season, bs = "re"), 
             data = Lake_Ontario_alldata_select)

summary(Ontario_gam)

plot(Ontario_gam)

####Key Takeaways - total phosphorous, but not winter severity, is a strong predictor of YCS in Lake Ontario####


####Lake Huron####
TP_Huron_Summary = TP_Huron %>% group_by(Year,Season) %>% summarise(AcrossLakeMeanTP = mean(Value), #mean across multiple stations
                                                                        SD = sd(Value),
                                                                        n = n(),
                                                                        SE_TP = sd(Value)/sqrt(n))
#Unlabled rows are Winter sampling
levels(TP_Huron_Summary$Season) <- c(levels(TP_Huron_Summary$Season), "Winter")
TP_Huron_Summary$Season[is.na(TP_Huron_Summary$Season)] <- "Winter"

ggplot(TP_Huron_Summary,
       aes(x = Year,
           y = AcrossLakeMeanTP, group = Season, color = Season)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = AcrossLakeMeanTP - SE_TP,
    ymax = AcrossLakeMeanTP + SE_TP)) + theme_minimal() +
  labs(y = "Mean TP cross all stations", title = "Lake Huron")


ggplot(WS_Huron,aes(x = Year, y = Huron_Cover)) + geom_point() + geom_line() + theme_minimal() +
  labs(y = "% Ice Cover", title = "Lake Huron")



Lake_Huron = Lake_Huron %>% mutate(Year = year_class) #copy over year class column into new year column 

Lake_Huron_EnvioData = merge(TP_Huron_Summary,WS_Huron,by = "Year")

Lake_Huron_alldata = merge(Lake_Huron_EnvioData,Lake_Huron, by = "Year")

Lake_Huron_alldata_select = Lake_Huron_alldata %>% select(c(1,2,3,6,7,23,24,25))

Huron_gam <- gam(fish_abundance_value ~ s(AcrossLakeMeanTP) + s(Huron_Cover) + s(Season, bs = "re"), 
                   data = Lake_Huron_alldata_select)

summary(Huron_gam) #both ice cover and TP significant
plot(Huron_gam)


#####Lake Michigan####
TP_Michigan_Summary = TP_Michigan %>% group_by(Year,Season) %>% summarise(AcrossLakeMeanTP = mean(Value), #mean across multiple stations
                                                                        SD = sd(Value),
                                                                        n = n(),
                                                                        SE_TP = sd(Value)/sqrt(n))
levels(TP_Michigan_Summary$Season) <- c(levels(TP_Michigan_Summary$Season), "Winter")
TP_Michigan_Summary$Season[is.na(TP_Michigan_Summary$Season)] <- "Winter"

ggplot(TP_Michigan_Summary,
       aes(x = Year,
           y = AcrossLakeMeanTP, group = Season, color = Season)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = AcrossLakeMeanTP - SE_TP,
    ymax = AcrossLakeMeanTP + SE_TP)) + theme_minimal() + labs(y = "Mean TP across all stations", title = "Lake Michigan")


ggplot(WS_Michigan,aes(x = Year, y = Michigan_Cover)) + geom_point() + geom_line() + theme_minimal() +
  labs(y = "% Ice Cover", title = "Lake Michigan")


ggplot(TP_Michigan_Summary,
       aes(x = Year,
           y = AcrossLakeMeanTP, group = Season, color = Season)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = AcrossLakeMeanTP - SE_TP,
    ymax = AcrossLakeMeanTP + SE_TP)) + theme_minimal() +
  labs(y = "Mean TP cross all stations", title = "Lake Michigan")


ggplot(WS_Michigan,aes(x = Year, y = Michigan_Cover)) + geom_point() + geom_line() + theme_minimal() +
  labs(y = "% Ice Cover", title = "Lake Michigan")



Lake_Michigan = Lake_Michigan %>% mutate(Year = year_class) #copy over year class column into new year column 

Lake_Michigan_EnvioData = merge(TP_Michigan_Summary,WS_Michigan,by = "Year")

Lake_Michigan_alldata = merge(Lake_Michigan_EnvioData,Lake_Michigan, by = "Year")

Lake_Michigan_alldata_select = Lake_Michigan_alldata %>% select(c(1,2,3,6,7,23,24,25))

Michigan_gam <- gam(fish_abundance_value ~ s(AcrossLakeMeanTP) + s(Michigan_Cover) + s(Season, bs = "re"), 
                 data = Lake_Michigan_alldata_select)

summary(Michigan_gam) #TP and Season significant 
plot(Michigan_gam)

#Spring and Summer follow a similar pattern 
TP_Michigan_Summary2 = TP_Michigan %>% group_by(Year) %>% summarise(AcrossLakeSeasonMeanTP = mean(Value), #mean across multiple stations
                                                                          SD = sd(Value),
                                                                          n = n(),
                                                                          SE_TP = sd(Value)/sqrt(n))
TP_Michigan_Summary2$Year = as.factor(TP_Michigan_Summary2$Year)

Lake_Michigan_EnvioData2 = merge(TP_Michigan_Summary2,WS_Michigan,by = "Year")

Lake_Michigan_alldata2 = merge(Lake_Michigan_EnvioData2,Lake_Michigan, by = "Year")

Lake_Michigan_alldata_select2 = Lake_Michigan_alldata2 %>% select(c(1,2,3,6,7,23,24,25))


####Lake Superior####
TP_Superior_Summary = TP_Superior %>% group_by(Year,Season) %>% summarise(AcrossLakeMeanTP = mean(Value), #mean across multiple stations
                                                                          SD = sd(Value),
                                                                          n = n(),
                                                                          SE_TP = sd(Value)/sqrt(n))
ggplot(TP_Superior_Summary,
       aes(x = Year,
           y = AcrossLakeMeanTP, group = Season, color = Season)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = AcrossLakeMeanTP - SE_TP,
    ymax = AcrossLakeMeanTP + SE_TP)) + theme_minimal() + labs(y = "Mean TP across all stations", title = "Lake Superior")


ggplot(WS_Superior,aes(x = Year, y = Superior_Cover)) + geom_point() + geom_line() + theme_minimal() +
  labs(y = "% Ice Cover", title = "Lake Superior")


Lake_Superior = Lake_Superior %>% mutate(Year = year_class) #copy over year class column into new year column 

Lake_Superior_EnvioData = merge(TP_Superior_Summary,WS_Superior,by = "Year")

Lake_Superior_alldata = merge(Lake_Superior_EnvioData,Lake_Superior, by = "Year")

Lake_Superior_alldata_select = Lake_Superior_alldata %>% select(c(1,2,3,6,7,23,24,25))

Superior_gam <- gam(fish_abundance_value ~ s(AcrossLakeMeanTP) + s(Superior_Cover) + s(Season, bs = "re"), 
                    data = Lake_Superior_alldata_select)

summary(Superior_gam) #TP and Season significant 
plot(Superior_gam)


####Lake Erie####
TP_Erie_Summary = TP_Erie %>% group_by(Year,Season) %>% summarise(AcrossLakeMeanTP = mean(Value), #mean across multiple stations
                                                                          SD = sd(Value),
                                                                          n = n(),
                                                                          SE_TP = sd(Value)/sqrt(n))
levels(TP_Erie_Summary$Season) <- c(levels(TP_Erie_Summary$Season), "Winter")
TP_Erie_Summary$Season[is.na(TP_Erie_Summary$Season)] <- "Winter"


ggplot(TP_Erie_Summary,
       aes(x = Year,
           y = AcrossLakeMeanTP, group = Season, color = Season)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = AcrossLakeMeanTP - SE_TP,
    ymax = AcrossLakeMeanTP + SE_TP)) + theme_minimal() + labs(y = "Mean TP across all stations", title = "Lake Erie")


ggplot(WS_Erie,aes(x = Year, y = Erie_Cover)) + geom_point() + geom_line() + theme_minimal() +
  labs(y = "% Ice Cover", title = "Lake Erie")


Lake_Erie = Lake_Erie %>% mutate(Year = year_class) #copy over year class column into new year column 

Lake_Erie_EnvioData = merge(TP_Erie_Summary,WS_Erie,by = "Year")

Lake_Erie_alldata = merge(Lake_Erie_EnvioData,Lake_Erie, by = "Year")

Lake_Erie_alldata_select = Lake_Erie_alldata %>% select(c(1,2,3,6,7,23,24,25))

Erie_gam <- gam(fish_abundance_value ~ s(AcrossLakeMeanTP) + s(Erie_Cover) + s(Season, bs = "re"), 
                    data = Lake_Erie_alldata_select)

summary(Erie_gam) #TP and Season significant 
plot(Erie_gam)




####Mussel Data - Michigan####

MichiganMussels = read_csv("Mussel_Southern_Basin_Lake_Michigan.csv")

MichiganMussels = MichiganMussels[!is.na(MichiganMussels$Biomass_Density_g_per_m2), ]

#sum across depths

MichiganMusselsSummarized = MichiganMussels %>% group_by(Year,Species) %>% summarise(MusselDensity = sum(Biomass_Density_g_per_m2))

ggplot(MichiganMusselsSummarized, aes(x = as.factor(Year), y = MusselDensity, group = Species, color = Species)) + 
  geom_point() + geom_line() + theme_minimal() +
  labs(x = "Year", y = "Mussel Density (g/m2)", title = "Lake Michigan Southern Basin")


MichiganMusselsSummarizedwide = MichiganMusselsSummarized %>% pivot_wider(
  names_from = Species, 
  values_from = MusselDensity) 


MichiganALL = merge(Lake_Michigan_alldata_select2,MichiganMusselsSummarizedwide, by = "Year") 


MichiganMussel_gam <- gam(fish_abundance_value ~ s(AcrossLakeSeasonMeanTP) + s(Michigan_Cover) + s(Quagga, k = 3) + s(Zebra, k = 3), 
                          data = MichiganALL)
summary(MichiganMussel_gam) #I don't really trust this result because of reduce k (degress of freedom)

plot(MichiganMussel_gam)


MichiganALLSeason = merge(Lake_Michigan_alldata_select, MichiganMusselsSummarizedwide, by = "Year")

MichiganMussel_gam2 <- gam(fish_abundance_value ~ s(AcrossLakeMeanTP) + s(Michigan_Cover) + s(Quagga, k = 3) + s(Zebra, k = 3) + s(Season, bs = "re"), 
                          data = MichiganALLSeason)

summary(MichiganMussel_gam2)


####ordinal scale - correct way to go about this ####

# Zebra Mussels (Max density around 10)
MichiganMusselsSummarizedwide$Zebra_ordinal <- cut(MichiganMusselsSummarizedwide$Zebra, 
                         breaks = c(0, 3, 7, 10), 
                         labels = c("Low", "Medium", "High"), 
                         include.lowest = TRUE)

# Quagga Mussels (Max density around 65)
MichiganMusselsSummarizedwide$Quagga_ordinal <- cut(MichiganMusselsSummarizedwide$Quagga, 
                         breaks = c(0, 20, 40, 65), 
                         labels = c("Low", "Medium", "High"), 
                         include.lowest = TRUE)


#Is it more accurate to take the mean TP across both seasons since the YCS is based on the whole year and the other predictors are not season specific???

MichiganALLoridinal2 = merge(Lake_Michigan_alldata_select2,MichiganMusselsSummarizedwide, by = "Year") 

# Convert to ordered factors first
MichiganALLoridinal2$Quagga_ordered <- ordered(MichiganALLoridinal2$Quagga_ordinal)
MichiganALLoridinal2$Zebra_ordered  <- ordered(MichiganALLoridinal2$Zebra_ordinal)

#Quagga and Zebra cannot be include as smooths because they are ordinal


MichiganMussel_gam2 <- gam(fish_abundance_value ~ s(AcrossLakeSeasonMeanTP) + s(Michigan_Cover) + Quagga_ordered + Zebra_ordered, 
                          data = MichiganALLoridinal2)

summary(MichiganMussel_gam2)

plot(MichiganMussel_gam2)

#How to interpet parametrix coefficients
library(marginaleffects)
plot_predictions(MichiganMussel_gam2, condition = "Quagga_ordered")
plot_predictions(MichiganMussel_gam2, condition = "Zebra_ordered")



plot(MichiganMussel_gam2)
####Erie Mussel Data####

ErieMussels = read_csv("Erie_mussel_data_8-10-26.csv")

#Interpolate missing data for years in between sampling dates 
library(zoo)

ErieMussels_interp <- ErieMussels %>%
  arrange(Year) %>%
  mutate(
    interpolated_density = na.approx(
      Dress_spp_density,
      x = Year,
      rule = 2, #extends the nearest observed value (key for the beginning)
      na.rm = FALSE #instead of removing NAs it keeps them in the output
    )
  )


#creating the df for gamm

TP_Erie_Summary2 = TP_Erie %>% group_by(Year) %>% summarise(AcrossLakeSeasonMeanTP = mean(Value), #mean across multiple stations
                                                                    SD = sd(Value),
                                                                    n = n(),
                                                                    SE_TP = sd(Value)/sqrt(n))

TP_Erie_Summary2$Year = as.factor(TP_Erie_Summary2$Year)

Lake_Erie_EnvioData2 = merge(TP_Erie_Summary2,WS_Erie,by = "Year")

Lake_Erie_alldata2 = merge(Lake_Erie_EnvioData2,Lake_Erie, by = "Year")

Lake_Erie_alldata_select2 = Lake_Erie_alldata2 %>% select(c(1,2,3,6,7,23,24,25))

ErieALL = merge(Lake_Erie_alldata_select2, ErieMussels_interp, by = "Year")

#convert to ordinal

ErieALL$interpolated_Dress_spp_ordinal <- cut(ErieALL$interpolated_density, 
                                                    breaks = c(350, 1500, 2650, 3800), 
                                                    labels = c("Low", "Medium", "High"), 
                                                    include.lowest = TRUE)

ErieALL$interpolated_Dress_spp_ordered <- ordered(ErieALL$interpolated_Dress_spp_ordinal)


ErieMussel_gam <- gam(fish_abundance_value ~ s(AcrossLakeSeasonMeanTP) + s(Erie_Cover) + interpolated_Dress_spp_ordered, 
                           data = ErieALL)

summary(ErieMussel_gam)

plot_predictions(ErieMussel_gam, condition = "interpolated_Dress_spp_ordered")

plot(ErieMussel_gam)
