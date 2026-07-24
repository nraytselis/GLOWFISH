setwd("~/Desktop/Whitefish/LakeGeneva")
library(readr)
library(dplyr)
library(tidyr)
library(lubridate) #add rows for missing dates
library(zoo) #interpolate values 



#bring in the data 
GenevaZoops = read_csv("clean_zooplankton_abundance_LakeGeneva.csv")
GenevaLarvFish = read_csv("clean_fish_larval_LakeGeneva.csv")
GenevaZoopLength = read_csv("clean_zooplankton_length_LakeGeneva.csv")
GenevaLarvLength = read_csv("clean_larval_length_LakeGeneva.csv")
StationData = read_csv("clean_stationid_LakeGeneva.csv")


####Light####
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

##interpolate zooplankton densities 

#add rows for missing dates
#add day of year column
GenevaZoops <- GenevaZoops %>% 
  rename(
    year = year_yyyy,
    month = month_mm,
    day = day_of_month_dd
  )


GenevaZoops <- GenevaZoops %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-")), day_of_year = yday(date))

#filter for daphnia only
GenevaZoops = GenevaZoops %>% filter(taxa_name_orig == "Daphnia")

GenevaZoops_interp <- GenevaZoops %>%
  mutate(date = as.Date(date),
         year = year(date)) %>%
  group_by(year) %>%
  complete(
    date = seq.Date(
      floor_date(first(date), "year"),
      ceiling_date(first(date), "year") - days(1),
      by = "day"
    )
  ) %>%
  ungroup()

#interpolate densities 
GenevaZoops_interp <- GenevaZoops_interp %>%
  arrange(date) %>%
  mutate(
    interpolated_zoopValue = na.approx(
      zoop_value,
      x = date,
      rule = 2, #extends the nearest observed value (key for the beginning)
      na.rm = FALSE #instead of removing NAs it keeps them in the output
    ),
    day_of_year = yday(date) #populates interpolate rows with days of the year 
  )

#convert to density in m^3
GenevaZoops_interp <- GenevaZoops_interp %>%
  mutate(
    density_m3 = interpolated_zoopValue * 1000
  )

#combine daylight and interpolated dateframes 
GenevaZoopsFull <- GenevaZoops_interp %>%
  left_join(resultsLight, by = "day_of_year")

################# Daily fish growth simulation

## Daily daylight hours (changes every day of the year) 
daylight_hours <- GenevaZoopsFull$hoursdaylight
n_days <- nrow(GenevaZoopsFull)

## Initial fish length (mm)
Pred_Length <- numeric(n_days + 1)
Pred_Length[1] <- 10 #assume fish starts off at 10mm

## Average Daphnia length (mm)
zoop_length_mm_df <- GenevaZoopLength %>%
  group_by(taxa_name_orig) %>%
  summarise(mean_length = mean(length_mm))

zoop_length_mm <- zoop_length_mm_df$mean_length[4] #select specific the daphnia 

## Convert prey length to weight (µg)

#regression terms from Watkins paper
#Ln(W) = Ln(alpha) + Beta*Ln(L), W in ug and L in mm
#coefficients for copepods and cladocern:
#Copepods...Ln(alpha) = 1.953
#Copepods...Beta = 2.399
#Daphnia (Cladoceran)...Ln(alpha) = 1.468 
#Daphnia...Beta = 2.829

Zoop_Weight_ug <- exp(1.468 + 2.829 * log(zoop_length_mm))

## Convert zooplankton density from density per L to density per m^3
GenevaZoops_interp$density_m3 <-GenevaZoops_interp$interpolated_zoopValue * 1000 

## Activity coefficient
A <- 1


## store results from simulation

Weight_g             <- numeric(n_days + 1)
Fish_Speed           <- numeric(n_days)
acuityaverage        <- numeric(n_days)
RD_mm                <- numeric(n_days)
RD_m                 <- numeric(n_days)
SV                   <- numeric(n_days)
ER_s                 <- numeric(n_days)
ERdaylight           <- numeric(n_days)
CS                   <- numeric(n_days)
captured             <- numeric(n_days)
consumption_ugPerday <- numeric(n_days)
F_egest              <- numeric(n_days)
U_excret             <- numeric(n_days)
SDA                  <- numeric(n_days)
R                    <- numeric(n_days)
growth_ugPerday      <- numeric(n_days)
growth_g             <- numeric(n_days)

## Initial fish weight (g)
Weight_g[1] <- 0.00000188 * Pred_Length[1]^3.296


## Daily growth model - currently runs for 1 year 
for(i in seq_len(n_days)) {
  
  ## Environmental conditions for this day
  daylight_today <- GenevaZoopsFull$hoursdaylight[i]
  density_today  <- GenevaZoopsFull$density_m3[i]
  
  ## Swimming speed (m/s)
  Fish_Speed[i] <- Pred_Length[i] / 1000
  
  ## Visual acuity
  acuityaverage[i] <-
    0.0167 *
    exp(
      9.14 -
        2.4 * log(Pred_Length[i]) +
        0.229 * log(Pred_Length[i])^2
    )
  
  ## Reactive distance (mm and m)
  RD_mm[i] <- (zoop_length_mm / 2) *
    tan(acuityaverage[i] / 2)
  
  RD_m[i] <- RD_mm[i] / 1000
  
  ## Search volume (m^3/s)
  SV[i] <- Fish_Speed[i] * pi * RD_m[i]^2 * 0.5
  
  ## Encounter rate (prey/s)
  ER_s[i] <- SV[i] * density_today
  
  ## Encounter rate during daylight
  ERdaylight[i] <- ER_s[i] * daylight_today * 3600
  
  ## Capture success - Anneville 2010 https://link.springer.com/article/10.1007/s10641-010-9755-1
  #LDS<3, 12.5-13.5mm
  #LDS>3, 15-22 mm
  #capture success of daphnia by LDS<3 = 0.15382444902005643
  #capture success of daphnia by LDS>3 = 0.5467894455883474
  
  if(Pred_Length[i] < 15){
    CS[i] <- 0.15382444902005643
  } else {
    CS[i] <- 0.5467894455883474
  }
  
  
  ## Number of prey captured
  captured[i] <- ERdaylight[i] * CS[i]
  
  ## Consumption (µg/day)
  consumption_ugPerday[i] <- captured[i] * Zoop_Weight_ug
  
  ## Bioenergetics losses
  F_egest[i]  <- 0.19 * consumption_ugPerday[i]
  U_excret[i] <- 0.07 * consumption_ugPerday[i]
  SDA[i]      <- 0.17 * (consumption_ugPerday[i] - F_egest[i])
  
  ## Respiration
  R[i] <-
    0.00584 *
    Weight_g[i]^(-0.05341) *
    exp(0.5060)
  
  ## Daily growth (µg/day) 
  growth_ugPerday[i] <-
    consumption_ugPerday[i] -
    F_egest[i] -
    U_excret[i] -
    SDA[i] -
    (R[i] * A)
  
  ## Convert growth to grams
  growth_g[i] <- growth_ugPerday[i] /1000000 #growth shouldn't be negative 
  
  ## Update fish weight
  Weight_g[i + 1] <- Weight_g[i] + growth_g[i]
  
  ## Convert weight back to length (mm)
  #calculate weight from zooplankton lengths using regression terms from Watkins 
  #Ln(W) = Ln(alpha) + Beta*Ln(L), W in ug and L in mm
  #coefficients for copepods and cladocern:
  #Copepods...Ln(alpha) = 1.953
  #Copepods...Beta = 2.399
  #Daphnia (Cladoceran)...Ln(alpha) = 1.468 
  #Daphnia...Beta = 2.829
  #Ln(W) = Ln(1.468) + 2.829*Ln(Pred_Length[i])
  #Pred_length[i+1] = (W/1.468)^(1/2.892)
  
  Pred_Length[i + 1] <-(Weight_g[i + 1]/1.468)^(1/2.892)
}


## Results
results <- data.frame(
  Date = GenevaZoopsFull$date,
  Year = GenevaZoopsFull$year,
  Day_of_Year = GenevaZoopsFull$day_of_year,
  Daylight_hours = GenevaZoopsFull$hoursdaylight,
  Zooplankton_density_m3 = GenevaZoopsFull$density_m3,
  Fish_length_mm = Pred_Length[-(n_days + 1)],
  Fish_weight_g = Weight_g[-(n_days + 1)],
  Consumption_ug = consumption_ugPerday,
  Growth_ug = growth_ugPerday,
  Growth_g = growth_g
)

head(results)