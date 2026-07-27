setwd("~/Desktop/Whitefish/LakeGeneva")
library(readr)
library(dplyr)
library(tidyr)
library(lubridate) #add rows for missing dates
library(zoo) #interpolate values 


#For each model day, we will use larval size and water temperature to growth and respiration (losses)
#consumption (Fig. 4). For each lake and year, we will run the model until the larval fish undergoes metamorphosis to the 
#juvenile stage (i.e., reaches 25 mm). The daily growth estimates during this larval stage will be used to estimate the mean 
#daily growth rate from hatch to metamorphosis for each year (cohort) in each lake.

#bring in the data 
GenevaZoops = read_csv("clean_zooplankton_abundance_LakeGeneva.csv")
GenevaLarvFish = read_csv("clean_fish_larval_LakeGeneva.csv")
GenevaZoopLength = read_csv("clean_zooplankton_length_LakeGeneva.csv")
GenevaLarvLength = read_csv("clean_larval_length_LakeGeneva.csv")
StationData = read_csv("clean_stationid_LakeGeneva.csv")
WaterChem = read_csv("clean_water_chemistry.csv")

###Filter Water Chemistry Data for just temperatures rows.
#multiple parameters exist in the W_value column because of drop_downs in excel sheets
Temp = WaterChem %>% filter(W_units == "celsius")

#add date column 
Temp <- Temp %>% 
  rename(
    year = year_yyyy,
    month = month_mm,
    day = day_of_month_dd
  )

Temp <- Temp %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-")), day_of_year = yday(date))

#add rows for additional dates
Temp_interp <- Temp_interp %>%
  group_by(year) %>%
  arrange(date) %>%
  mutate(
    interpolated_temp = na.approx(
      W_value,
      x = date,
      rule = 2,
      na.rm = FALSE
    ),
    day_of_year = yday(date)
  ) %>%
  ungroup()

#interpolate temperature 
Temp_interp <- Temp_interp %>%
  arrange(date) %>%
  mutate(
    interpolated_temp = na.approx(
      W_value,
      x = date,
      rule = 2, #extends the nearest observed value (key for the beginning)
      na.rm = FALSE #instead of removing NAs it keeps them in the output
    ),
    day_of_year = yday(date) #populates interpolate rows with days of the year 
  )

#select columns for specifically date, temp, and temp units
Temp_interp_select <- Temp_interp %>%
  select(date, interpolated_temp)

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
  group_by(year) %>%
  arrange(date) %>%
  mutate(
    interpolated_zoopValue = na.approx(
      zoop_value,
      x = date,
      rule = 2,
      na.rm = FALSE
    ),
    day_of_year = yday(date)
  ) %>%
  ungroup()

#convert to density in m^3
GenevaZoops_interp <- GenevaZoops_interp %>%
  mutate(
    density_m3 = interpolated_zoopValue * 1000
  )

#combine daylight and interpolated dateframes 
GenevaZoopsFull <- GenevaZoops_interp %>%
  left_join(resultsLight, by = "day_of_year")

#combine temp dataframe with others 
GenevaZoopsFull <- GenevaZoopsFull %>%
  left_join(Temp_interp_select, by = "date")


################# Daily fish growth simulation

## Daily daylight hours (changes every day of the year) 
daylight_hours <- GenevaZoopsFull$hoursdaylight
max_days <- 200
start_days <- which(
  month(GenevaZoopsFull$date) %in% 3:5
)

## Initial fish length (mm)
Pred_Length <- numeric(n_days + 1)
Pred_Length[1] <- 10 #assume fish starts off at 10mm

## Average Daphnia length (mm)
zoop_length_mm_df <- GenevaZoopLength %>%
  group_by(taxa_name_orig) %>%
  summarise(mean_length = mean(length_mm))

zoop_length_mm <- zoop_length_mm_df$mean_length[4] #select specifically the daphnia 

#ggplot(GenevaZoopLength,aes(x= as.factor(taxa_name_orig), y = length_mm)) + geom_point()


## Convert prey length to weight (µg)

#regression terms from Watkins paper
#Ln(W) = Ln(alpha) + Beta*Ln(L), W in ug and L in mm

#coefficients for copepods and cladocern:
#Copepods...Ln(alpha) = 1.953
#Copepods...Beta = 2.399
#Daphnia (Cladoceran)...Ln(alpha) = 1.468 
#Daphnia...Beta = 2.829

#W = alpha * L^2.829

Zoop_Weight_ug <- exp(1.468) * zoop_length_mm^2.829 

## Convert zooplankton density from density per L to density per m^3 (for encounter rate equation)
GenevaZoops_interp$density_m3 <-GenevaZoops_interp$interpolated_zoopValue * 1000 


## store results from simulation

Weight_g             <- numeric(n_days + 1) #weight of the fish after it has eaten a day's worth of zoops (calculated at end of sim and feeds back into beginning of next sim)
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


## Daily growth model
for(i in seq_len(n_days)) {
  
  ## Environmental conditions for this day
  daylight_today <- GenevaZoopsFull$hoursdaylight[i] #hours of daylight
  density_today  <- GenevaZoopsFull$density_m3[i] #density of daphnia 
  
  ## Swimming speed (m/s)
  Fish_Speed[i] <- Pred_Length[i] / 1000 #swimming speed is dependent on the current length of the fish
  
  ## Visual acuity - Letcher formula gives an angle in degrees
  acuityaverage[i] <- 
    0.0167 *
    exp(
      9.14 -
        2.4 * log(Pred_Length[i]) +
        0.229 * log(Pred_Length[i])^2
    )
  
  ## Convert visual angle from degrees to radians for tan() 
  acuity_rad <- acuityaverage[i] * pi / 180
  
  ## Reactive distance (mm)
  RD_mm[i] <- zoop_length_mm /
    (2 * tan(acuity_rad / 2))
  
  ## Reactive distance (m)
  RD_m[i] <- RD_mm[i] / 1000
  
  ## Search volume (m^3/s)
  SV[i] <- Fish_Speed[i] * pi * RD_m[i]^2 * 0.5
  
  ## Encounter rate (prey/s)
  ER_s[i] <- SV[i] * density_today
  
  ## Encounter rate during daylight
  ERdaylight[i] <- ER_s[i] * daylight_today * 3600
  
  ## Capture success - Anneville 2010 https://link.springer.com/article/10.1007/s10641-010-9755-1
  #Capture success of prey depends on the size of the fish, higher success with longer (older) fish
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
  
  ## Bioenergetics losses - Huuskonen et al. 1998 (based on coefficients from Table in paper)
  #Egestion and ecretion kept as constants and not temperature dependent 
  #"It is possible tomodel egestion and excretion as functions of temperature and ration (e.g. Elliott,1976; Cui & Wootton, 1988) 
  #but egestion and excretion have been observed to contribute only slightly to prediction errors in sensitivity analyses of bioenergetics models (Bartell et al., 1986).
  F_egest[i]  <- 0.19 * consumption_ugPerday[i]
  U_excret[i] <- 0.07 * consumption_ugPerday[i]
  SDA[i]      <- 0.17 * (consumption_ugPerday[i] - F_egest[i])

  
  ## Respiration
  #R = RA * W^RB * e^RQT* A
  #W is fish mass (g)
  #T is water temp
  #RA, RB, RQ are fitted coefficients 
  
  #Respiration Coefficients 
  RA = 0.00584
  RB = -0.05341
  RQ = 0.05060
  A <- 1 ## Activity coefficient
  
  
  #This was based on coefficients in the paper, but needs to be generalized 
  # R[i] <-
  #   0.00584 *
  #   Weight_g[i]^(-0.05341) *
  #   exp(0.5060)
  
  Temp_today <- GenevaZoopsFull$interpolated_temp[i]
  
  #Respiration modified to be temperature dependent 
  R[i] <-
    RA *
    Weight_g[i]^RB *
    exp(RQ * Temp_today) * A
  
  ## Daily growth (µg/day) 
  growth_ugPerday[i] <-
    consumption_ugPerday[i] -
    F_egest[i] -
    U_excret[i] -
    SDA[i] -
    (R[i] * A)
  
  ## Convert growth to grams ####Feels like this is where the issue is
  growth_g[i] <- growth_ugPerday[i] /1000000 #growth shouldn't be negative 
  
  ## Update fish weight
  Weight_g[i + 1] <- Weight_g[i] + growth_g[i]
  
  # Convert grams to µg
  Weight_ug <- Weight_g[i + 1] * 1e6
  
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
  
  # Convert weight back to length
  Pred_Length[i + 1] <-
    (Weight_ug / exp(1.468))^(1/2.829)
  
  if (Pred_Length[i + 1] >= 25) {
    message("Fish reached 25 mm (recruitment age) on day ", i)
    
    # Fill remaining values with NA
    Pred_Length[(i + 2):(n_days + 1)] <- NA
    Weight_g[(i + 2):(n_days + 1)] <- NA
    
    break #stop for loop
  }
}


## Results
results <- data.frame(
  Date = GenevaZoopsFull$date,
  Year = GenevaZoopsFull$year,
  Day_of_Year = GenevaZoopsFull$day_of_year,
  Daylight_hours = GenevaZoopsFull$hoursdaylight,
  Zooplankton_density_m3 = GenevaZoopsFull$density_m3,
  Temperature = GenevaZoopsFull$interpolated_temp,
  Fish_length_mm = Pred_Length[-(n_days + 1)],
  Fish_weight_g = Weight_g[-(n_days + 1)],
  Consumption_ug = consumption_ugPerday,
  Growth_ug = growth_ugPerday,
  Growth_g = growth_g
  
)


####Daily fish growth simulation
#Fish start at 10mm on every `hatch_day` from March 1 - May 31 from 1974 to 2024
#simulation runs until fish...
#reaches 25mm (recruitment size)
#hits `max_days` = 200, so doesn't run indefinitely

run_growth_sim <- function(hatch_day, env, zoop_length_mm, Zoop_Weight_ug,
                           start_length = 10, max_days = 200) {
  
  n_env    <- nrow(env)
  n_steps <- min(max_days, n_env - hatch_day + 1) 
  Pred_Length <- numeric(n_steps + 1)
  Weight_g    <- numeric(n_steps + 1)
  growth_g        <- numeric(n_steps)
  growth_ugPerday <- numeric(n_steps)
  
  Pred_Length[1] <- start_length
  Weight_g[1]    <- 0.00000188 * Pred_Length[1]^3.296
  
  recruited   <- FALSE #creates True/False for whether fish has reached recruitment (25 mm)
  days_elapsed    <- NA_integer_ #creates empty integer for # of days until recruitment 
  
  RA <- 0.00584; RB <- -0.05341; RQ <- 0.05060; A <- 1
  
  for (k in seq_len(n_steps)) {
    i <- hatch_day + k - 1  # index into environmental data 
    
    daylight_today <- env$hoursdaylight[i]
    density_today  <- env$density_m3[i]
    Temp_today     <- env$interpolated_temp[i]
    
    Fish_Speed <- Pred_Length[k] / 1000
    
    acuityaverage <- 0.0167 * exp(9.14 - 2.4 * log(Pred_Length[k]) +
                                    0.229 * log(Pred_Length[k])^2)
    acuity_rad <- acuityaverage * pi / 180
    RD_mm <- zoop_length_mm / (2 * tan(acuity_rad / 2))
    RD_m  <- RD_mm / 1000
    SV    <- Fish_Speed * pi * RD_m^2 * 0.5
    ER_s  <- SV * density_today
    ERdaylight <- ER_s * daylight_today * 3600
    
    CS <- if (Pred_Length[k] < 15) 0.15382444902005643 else 0.5467894455883474
    captured <- ERdaylight * CS
    consumption_ugPerday <- captured * Zoop_Weight_ug
    
    F_egest  <- 0.19 * consumption_ugPerday
    U_excret <- 0.07 * consumption_ugPerday
    SDA      <- 0.17 * (consumption_ugPerday - F_egest)
    
    R <- RA * Weight_g[k]^RB * exp(RQ * Temp_today) * A
    
    growth_ugPerday[k] <- consumption_ugPerday - F_egest - U_excret - SDA - (R * A)
    growth_g[k] <- growth_ugPerday[k] / 1e6
    
    new_weight_g <- Weight_g[k] + growth_g[k]
    
    Weight_g[k + 1] <- new_weight_g
    Weight_ug <- new_weight_g * 1e6
    Pred_Length[k + 1] <- (Weight_ug / exp(1.468))^(1 / 2.829)
    
    if (Pred_Length[k + 1] >= 25) {
      recruited <- TRUE
      days_elapsed  <- k
      break
    }
  }
  
  final_len <- max(Pred_Length)
  
  list(
    start_date          = env$date[hatch_day],
    hatch_day           = hatch_day,
    days_elapsed        = days_elapsed,
    recruited           = recruited,
    mean_growth_ug_day  = mean(growth_ugPerday[1:days_elapsed]), 
    mean_growth_mm_day  = (Pred_Length[days_elapsed + 1] - Pred_Length[1]) / days_elapsed,
    final_length        = final_len
  )
}

#Run simulartion starting only on days in the March 1 - May 31 hatch window (each year, 1974-2024)
n_days   <- nrow(GenevaZoopsFull)
max_days <- 200

start_days <- which(month(GenevaZoopsFull$date) %in% 3:5) #restrict to time of year when whitefish likely hatching

cohort_results <- vector("list", length(start_days))

for (j in seq_along(start_days)) {
  s <- start_days[j]
  cohort_results[[j]] <- run_growth_sim(
    hatch_day      = s,
    env            = GenevaZoopsFull,
    zoop_length_mm = zoop_length_mm,
    Zoop_Weight_ug = Zoop_Weight_ug,
    max_days       = max_days
  )
}

cohort_summary <- bind_rows(lapply(cohort_results, as.data.frame)) %>%
  mutate(start_year = year(start_date))

yearly_growth <- cohort_summary %>%
  filter(recruited) %>%
  group_by(start_year) %>%
  summarise(
    mean_days_to_recruitment    = mean(days_elapsed),
    mean_growth_mm_day         = mean(mean_growth_mm_day),
    mean_growth_ug_day         = mean(mean_growth_ug_day),
    n_cohorts                  = n()
  )

ggplot(yearly_growth, aes(x = start_year, y = mean_growth_mm_day)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Mean daily growth rate (mm/day)"
  ) +
  theme_minimal()

ggplot(yearly_growth, aes(x = start_year, y = mean_days_to_recruitment)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Mean days to recruitment (25mm)") +
  theme_minimal()

