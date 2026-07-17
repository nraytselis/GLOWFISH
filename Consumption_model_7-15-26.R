#Play with Lake Geneva Data
setwd("~/Desktop/Whitefish_Test/LakeGeneva")
library(readr)
library(dplyr)
GenevaZoops = read_csv("clean_zooplankton_abundance_LakeGeneva.csv")
GenevaLarvFish = read_csv("clean_fish_larval_LakeGeneva.csv")
GenevaZoopLength = read_csv("clean_zooplankton_length_LakeGeneva.csv")
GenevaLarvLength = read_csv("clean_larval_length_LakeGeneva.csv")
StationData = read_csv("clean_stationid_LakeGeneva.csv")


#ER = SV * Prey Density * Light

####Surface Volume####
#SV = swim speed * RA (reactive area)
#swim speed is standardized (~1 body length/s)
#RA = RD^2 * pi * Prop where Prop = 0.5
#RD = Prey Length /2*tan(αl/2)
##αl = 0.0167 ⋅ e^(9.14 − 2.4 ⋅ ln(l) + 0.229 ⋅ (ln(l))) 
#l is the total length of the fish larvae (mm)

#does prey length (zoops) or pred length (fish larvae) vary with...
  #DOY
  #temp 
  #species (zoops)

####Prey density####
#does this vary with...
#day of year?
#prey species?


#####Light#### 
#varies with DOY (day of year)

library(mgcv) #gams


####Pt. 1 - Surface Volume####

#Two parameters that depend on data in calculating SV are zooplankton length and fish larvae length

#Do zooplankton species differ in length? - YES
plot_zoop_length <- GenevaZoopLength %>%
  group_by(taxa_name_orig) %>%
  summarize(
    n = n(),
    mean_val = mean(length_mm),
    se_val = sd(length_mm) / sqrt(n)
  )

ggplot(plot_zoop_length, aes(x = taxa_name_orig, y = mean_val)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.1) +
  labs(x = "Zooplankton Taxa", y = "Density per L")

#Do individual species differ in size(length) over time? - Unlikely to be major differences based on zoop reproductive times
Daphnia = GenevaZoopLength %>% filter(taxa_name_orig == "Daphnia_sp")
ggplot(Daphnia, aes(x=month_mm,y=length_mm)) + geom_point() #visually doesn't appear to differ with month
Leptodora_kindtii = GenevaZoopLength %>% filter(taxa_name_orig == "Leptodora_kindtii")
ggplot(Leptodora_kindtii, aes(x=month_mm,y=length_mm)) + geom_point() #visually doesn't appear to differ with month


#Does larval fish size differ with time of year? 
#Probably but all fish length data is from April - cannot answer this question. Estimate mean length due to limited data. 
PredLength = mean(GenevaLarvLength$length_mm) 

####total encounter rate based on the proportion of each zooplankton species at a given time####
#Estimate individual encounter rates based on prey densities at specific times and prey size (which is relatively constant)
#RD equation is the one impacted by this
#RD (mm) = Prey Length /2*tan(αl/2)
#αl = 0.0167 ⋅ e^(9.14 − 2.4 ⋅ ln(PredLength) + 0.229 ⋅ (ln(PredLength))) 
acuityaverage = 0.0167 * exp(9.14 - 2.4 * log(PredLength) + 0.229 * (log(PredLength))) 

#Prey length differs by species 
GenevaZoopRDValues <- GenevaZoopLength %>% 
  group_by(taxa_name_orig) %>% 
  summarise(
    RD = mean(length_mm) / 2 * tan(mean(acuityaverage) / 2)
  )

#RA values, by species
#RA = RD^2 * pi * 0.5 
GenevaZoopRDValues$RA = (GenevaZoopRDValues$RD)^2 * pi * 0.5

#rename because now multiple counter rate parameters
ERparms = GenevaZoopRDValues

#SV (m^3/s) = swimming speed * reactive area
#Because swimming speeds can vary among species and we are not aware of length-based estimates of swimming speed for larval whitefish species, 
#we will use a generalized relationship from Miller et al. (1988), which is ~1 body length/s
ERparms$SV = ERparms$RA #SV values for each zoop species

####Pt. 2 - Prey density####

#convert dates to day of the year in Geneva dataframe
library(lubridate)
library(tidyverse)

GenevaZoops <- GenevaZoops %>% 
  rename(
    year = year_yyyy,
    month = month_mm,
    day = day_of_month_dd
  )


GenevaZoops <- GenevaZoops %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-")), day_of_year = yday(date))

#why are there negative values for densities for bosmina?
#for now, take absolute value of zoop densities
plot(gam(abs(zoop_value) ~ s(day_of_year),
                 data = GenevaZoops,
                 family = nb()))

summary(gam(abs(zoop_value) ~ s(day_of_year),
            data = GenevaZoops,
            family = nb()))

#highly significant relationship between zooplankton abundance and day of the year 

#how many different species?
unique(GenevaZoops$taxa_name_orig)

#Do species differ in density? - YES
plot_zoop_value <- GenevaZoops %>%
  group_by(taxa_name_orig) %>%
  summarize(
    n = n(),
    mean_val = mean(zoop_value),
    se_val = sd(zoop_value) / sqrt(n)
  )

ggplot(plot_zoop_value, aes(x = taxa_name_orig, y = mean_val)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.1) +
  labs(x = "Zooplankton Taxa", y = "Density per L")


#Do dominant zoop taxa differ in density and with time of year? - somewhat, seems like overall density varies more
species_summary_month <- GenevaZoops %>% 
  group_by(month, taxa_name_orig) %>% 
  summarise(mean_density = mean(zoop_value)) 

ggplot(data = species_summary_month, aes(x = month, y = mean_density, fill = taxa_name_orig)) + 
  geom_col()


#SV (prey species) * prey density (@ each day of year)
unique(GenevaZoops$taxa_name_orig) #different taxa names - how to deal with this?
unique(ERparms$taxa_name_orig)
#different taxa names - how to deal with this? Overlap of most dominant species: Eudiaptomus gracilis,Daphnia, Cyclops_prealpinus
#filter by these taxa for now
GenevaZoopsFiltered = GenevaZoops %>% filter(taxa_name_orig == "Eudiaptomus gracilis" | 
                                               taxa_name_orig == "Daphnia" |
                                               taxa_name_orig == "Cyclops")

ERparmsFiltered = ERparms %>% filter(taxa_name_orig == "Eudiaptomus_gracilis" | 
                                               taxa_name_orig == "Daphnia_sp" |
                                               taxa_name_orig == "Cyclops_prealpinus")

#give common names, note that I am unsure if Cyclops_prealpinus is the cyclops in the dataframe with zoop densities
ERparmsFilteredRename = ERparmsFiltered %>%
  mutate(taxa_name_orig = case_when(
    taxa_name_orig == "Daphnia_sp" ~ "Daphnia",
    taxa_name_orig == "Cyclops_prealpinus" ~ "Cyclops",
    taxa_name_orig == "Eudiaptomus_gracilis" ~ "Eudiaptomus"))  
  
GenevaZoopsFilteredRename = GenevaZoopsFiltered %>%
  mutate(taxa_name_orig = case_when(
    taxa_name_orig == "Daphnia" ~ "Daphnia",
    taxa_name_orig == "Cyclops" ~ "Cyclops",
    taxa_name_orig == "Eudiaptomus gracilis" ~ "Eudiaptomus"))  

zoopsdf = merge(ERparmsFilteredRename,GenevaZoopsFilteredRename, by = "taxa_name_orig")

#convert dates to day of the year in Geneva dataframe
library(lubridate)
library(tidyverse)

zoopsdf <- zoopsdf %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-")), day_of_year = yday(date))


####Pt. 3 - Light####
#Station Data has lat long coords. Focus on SHL2 because that's what the 2016 data is from
Lats = StationData$stationid_lat_decdeg
lat = Lats[1]
#DayLength = 2/15*arcos*(-tan(Lat)*tan(solardeclination))
#solardeclination based on day of the year
#arcos provides angle in degrees which is multiplied by 2/15 to get hours, since earth rotates 15 degrees per hour
#solardeclination = 23.45*sin(360/365*(284+d))
#full formula:    dayLength = 2/15*arcos*(-tan(Lat)*tan(23.45*sin(360/365*(284+day)))) 

deg2rad <- function(deg) deg * pi / 180
rad2deg <- function(rad) rad * 180 / pi

solardec <- 23.45 * sin(deg2rad(360/365 * (284 + 1))) #get january 1st 

dayLength <- acos(-tan(deg2rad(lat)) * tan(deg2rad(solardec)))
dayLengthdeg <- rad2deg(dayLength)

day_length <- (2/15) * dayLengthdeg #8 hrs 46 min daylight

#loop through rest of the year
day_of_year <- 1:365

resultsLight <- data.frame(
  day_of_year = numeric(),
  hoursdaylight = numeric()
)

for (j in seq_along(day_of_year)) {
  
  solardec <- 23.45 * sin(deg2rad(360/365 * (284 + day_of_year[j])))
  
  dayLength <- acos(-tan(deg2rad(lat)) * tan(deg2rad(solardec)))
  
  dayLengthdeg <- rad2deg(dayLength)
  
  hoursdaylight <- (2/15) * dayLengthdeg
  
  resultsLight <- rbind(
    resultsLight,
    data.frame(
      day_of_year = day_of_year[j],
      hoursdaylight = hoursdaylight
    )
  )
}

ggplot(resultsLight, aes(x=day_of_year,y=hoursdaylight)) + geom_point()

zoopsdf = merge(resultsLight,zoopsdf,by = "day_of_year")


####Connecting all the pieces####

ERdf = merge(zoopsdf,resultsLight, by = "day_of_year") 

ERdf$ER = ERdf$SV * ERdf$zoop_value * ERdf$hoursdaylight.x

####plot of encounter rates based on day of the year for three different zooplankton taxa###
ggplot(ERdf, aes(x=day_of_year,y=ER,group=taxa_name_orig,color=taxa_name_orig)) + geom_line()

#take encounter rate averages for each species across the whole year 
EncounterRateSummary = ERdf %>%
  group_by(taxa_name_orig) %>%
  summarize(
    n = n(),
    mean_ER = mean(ER),
    se_ER = sd(ER) / sqrt(n)
  )

ggplot(EncounterRateSummary, aes(x = taxa_name_orig, y = mean_ER)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_ER - se_ER, ymax = mean_ER + se_ER), width = 0.1) +
  labs(x = "Zooplankton Taxa", y = "Encounter Rate") + theme_classic()

#encounter rate units are individuals/second


#total encounter rate per day 
ERdf$encounter = ERdf$ER*60*60*24 

#Ignore handling time
#Letcher et al. (1996) also included handling time to inform “decisions” for individuals across encountered prey types, to inform a stochastic process through an individual based 
#model, but this level of detail is beyond the scope of our objective to estimate daily consumption of key prey under varying densities.


#Because foraging theory posits that larval fish do not attack each prey encountered, 
#capture success will be modeled to increase with larval fish length and decrease with prey size.
#Numbers of zooplankton consumed will be converted to biomass and then to energy density 
#for coupling with metabolic costs

#Capture Success = CSNum*fishLength/CSDen+fishLength 
#CSNum/CSDen = Prey trickiness, different zoops have different escape behaviors 
#We will assume that bigger zoops are trickier to catch as done in Letcher: 

#Capture successes from Letcher for: rotifer, nauplii, copepodite, copepod
#95% capture success for fish feeding on rotifers
#Capture success numerator 0.95; 0.90; 0.70; 0.90 
#Capture success denominator 10.0; 750.0; 5×E7; 5×E8

#Use Letcher's estimate for copepods

#Calculate capture success 
CS = CSNum*fishLength/CSDen+fishLength 

ERdf$captured      <- ERdf$dailyEncounter * ERdf$CS
ERdf$consumption   <- ERdf$captured * ERdf$Weight 

TotalConsumption <- sum(ERdf$consumption)


#Whitefish bioenergetic parameters (Huuskonen et al. 1998)

#Metabolic Rate (R), egestion (F), excretion (U) 

Ractive = 0.00584
Rbasal = 0.05341
egestion = 0.19
excretion = 0.07

#daily growth =  G = C - (R + F + U) 
#consumption estimated using foraging model 


#Wisconsin Style Bioenergetics Model
#larval whitefish have an energy density of 4,500 j/g wet and their zooplankton prey have an energy density of 2,000 j/g wet
#https://www.sciencedirect.com/science/article/abs/pii/S0304380023003058

#Growth = (Consumption*Energy Density of Prey - Total Metabolic Costs)/Energy Density of Fish






