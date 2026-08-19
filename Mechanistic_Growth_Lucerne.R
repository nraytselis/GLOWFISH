setwd("~/Desktop/Whitefish/LakeLucerne")
library(readr)
library(dplyr)
library(tidyr)
library(lubridate) #add rows for missing dates
library(zoo) #interpolate values 

Zoops = read_csv("GLOWFISH_Lucerne_Selz_zoop_abundance.csv")
WaterChem = read_csv("GLOWFISH_Lucerne_Selz_water_chem.csv")

###Filter Water Chemistry Data for just temperatures rows
#multiple parameters exist in the W_value column because of drop_downs in excel sheets

#There is a gap in temp data between 2004-2011
#data goes from 1976-2004 and 2011-2025
Temp = WaterChem %>% filter(W_parameter == "temperature")

#add date column 
Temp <- Temp %>% 
  rename(
    year = year_yyyy,
    month = month_mm,
    day = day_of_month_dd
  )

Temp <- Temp %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-")), day_of_year = yday(date)) 
#Take an average Temp across all depths before interpolating

#add rows for additional dates
Temp_interp <- Temp %>%
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


ggplot(Temp_interp_select, aes(x=date, y = interpolated_temp)) + geom_line() + theme_minimal() +
  labs(x="Year", y = "Temp (C)")


####Zoops####
unique(Zoops$taxa_name_orig)

copepods <- c(
  "Eudiaptomus gracilis",
  "Mixodiaptomus laciniatus",
  "Cyclops abyssorum",
  "Cyclops sp",
  "Megacyclops gigas",
  "Mesocyclops leuckarti",
  "Cyclops vicinus",
  "Cyclops bohater",
  "Diacyclops sp",
  "Mixodiaptomus sp",
  "Thermocyclops sp"
)


Zoops <- Zoops %>%
  filter(taxa_name_orig %in% copepods)

ggplot(Zoops, aes(x = taxa_name_orig, y = zoop_value)) + 
  geom_col() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + facet_wrap(~ year_yyyy) 

Zoopsfiltered = Zoops %>% filter(year_yyyy %in% c("1975","1976","1977","1978","1979","1980","1981","1982","1983"))

ggplot(Zoopsfiltered, aes(x = taxa_name_orig, y = zoop_value)) + 
  geom_col() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + facet_wrap(~ year_yyyy) 


#filter for dominant species, a lot are at very low levels
Zoops = Zoops %>% filter(taxa_name_orig == "Cyclops abyssorum" |
                           taxa_name_orig == "Cyclops sp" |
                           taxa_name_orig == "Cyclops vicinus" |
                           taxa_name_orig == "Eudiaptomus gracilis")

#add date column 
Zoops <- Zoops %>% 
  rename(
    year = year_yyyy,
    month = month_mm,
    day = day_of_month_dd
  )

Zoops <- Zoops %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-")), day_of_year = yday(date))


#convert to density/m3 from density/m2
#density/m2 / Depth in meters

Zoops$ave_depth_m = (Zoops$zoop_depth_max_m + Zoops$zoop_depth_min_m)/2

Zoops$density_m3 = Zoops$zoop_value/Zoops$ave_depth_m

#Zoop weight x density
#For now, let's assume any Cyclops labeled "cyclops sp." are "Cyclops abyssorum" as these are the most abundant of the cyclops species, aside from the unspecified cyclops sp.

#Take average density across stations for a given time point

Zoops <- Zoops %>% group_by(date,year,month,day_of_year,taxa_name_orig) %>%
  summarise(mean_density = mean(density_m3))


#Total density of all copepod species present at a given time point
ZoopsTotal <- Zoops %>% group_by(date,year,month,day_of_year) %>%
  summarise(TotalCopepodDensity = sum(as.numeric(mean_density))) 


ZoopsTotal_interp <- ZoopsTotal %>%
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

#Remove rows for 1975 because no data for spring or summer months

ZoopsTotal_interp = ZoopsTotal_interp %>% filter(year != "1975")

#interpolate densities - day to day changes are linear, not logistic
ZoopsTotal_interp <- ZoopsTotal_interp %>%
  group_by(year) %>%
  arrange(date) %>%
  mutate(
    interpolated_zoopValue = na.approx(
      TotalCopepodDensity,
      x = date,
      rule = 2,
      na.rm = FALSE
    ),
    day_of_year = yday(date)
  ) %>%
  ungroup()


ggplot(ZoopsTotal_interp, aes(x = date, y = interpolated_zoopValue)) + geom_line() +
  ylim(0,20000) + theme_minimal() + labs(x = "Year", y = "Interpolated Copepod Density (m3)") #Major decrease in zoops in late 1980s  


####Light####
#No lat long provided
#Lake Lucerne in central Switzerland is centered approximately at latitude 47.0093° N and longitude 8.4499° E

#DayLength = 2/15*arcos*(-tan(Lat)*tan(solardeclination))
#solardeclination based on day of the year
#arcos provides angle in degrees which is multiplied by 2/15 to get hours, since earth rotates 15 degrees per hour
#solardeclination = 23.45*sin(360/365*(284+d))
#full formula:    dayLength = 2/15*arcos*(-tan(Lat)*tan(23.45*sin(360/365*(284+day)))) 
# Set your date range for multiple years
start_date <- as.Date("1976-01-01")
end_date <- as.Date("2025-12-31")
date_seq <- seq.Date(start_date, end_date, by = "day")

lat <- 47.0093
deg2rad <- function(deg) deg * pi / 180
rad2deg <- function(rad) rad * 180 / pi

resultsLight <- data.frame(
  day_of_year = numeric(),
  hoursdaylight = numeric()
)

for (d in date_seq) {
  curr_date <- as.Date(d, origin = "1970-01-01")
  doy <- yday(curr_date)
  
  solardec <- 23.45 * sin(deg2rad(360 / 365 * (284 + doy)))
  
  # Protect against extreme latitudes/polar values where acos is undefined
  tan_val <- tan(deg2rad(lat)) * tan(deg2rad(solardec))
  tan_val <- max(min(tan_val, 1), -1) 
  
  dayLength_rad <- acos(tan_val)
  dayLength_deg <- rad2deg(dayLength_rad)
  hoursdaylight <- 2 * (1 / 15) * dayLength_deg
  
  resultsLight <- rbind(resultsLight, data.frame(
    day_of_year = doy,
    hoursdaylight = hoursdaylight
  ))
}

#combine daylight and interpolated dateframes 
ZoopsLight <- ZoopsTotal_interp %>%
  left_join(resultsLight, by = "day_of_year")

#combine temp dataframe with others 
ZoopsLightTemp <- ZoopsLight %>%
  left_join(Temp_interp_select, by = "date")



####### Mechanistic Whitefish Larval Bioenergetics Model

##### Parameters #####

#Energy densitiy of zooplankton - https://esapubs.org/archive/appl/A025/122/appendix-C.php
#Erik R. Schoen, David A. Beauchamp, Anna R. Buettner, and Nathanael C. Overman. 2015. Temperature and depth mediate resource competition and apparent competition between Mysis diluviana and kokanee. Ecological Applications 25:1962–1975. http://dx.doi.org/10.1890/14-1822.1
#Above paper cites Luecke and Brandt 1993
#Cladocerns - 1620 J / g
#Copepods - 2260 J / g

# Daphnia energy density
#PreyED <- 967.5 / 1e6        # J/ug wet weight #Hoefnagel et al. 2018, https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.3933 

PreyED <- 2260 / 1e6 # J/ug, copepods


# Assimilation efficiency - Huuskonen assumes a constant assimilation efficiency of 75.3% of ingested energy
Ae <- 0.753

# Oxygen conversion
#Oxy_cal is the constant converting grams of oxygen to Joules
#1 µmol Oxygen = 0.45 J
#molar mass of oyxgen is 32 g/mol
#0.45/32*1e-6 g O2 = 14,062.5 J/mg Oxygen

Oxy_cal <- 14.0625           # J per mg O2

# Energy density of fish tissue - Huuskonen
Fish_energy_density <- 4286  # J/g wet mass


####Zoop Weight####
#there are no data in the dataset for zoop length
#look up zoop lengths for the species mentioned
#https://www.st.nmfs.noaa.gov/copepedia/taxa/T4003215/html/biometricframe.html  

#"Eudiaptomus gracilis", 0.5–0.8 mm, ave is 0.65, https://pmc.ncbi.nlm.nih.gov/articles/PMC6084575/ 
exp(1.953) * (0.5+0.8)/2^2.399 #1.7 ug equation from Watkins 
#"Mixodiaptomus laciniatus",0.8 - 0.9 mm, https://www.researchgate.net/publication/286866943_Life_cycles_size_and_reproduction_of_the_two_coexisting_calanoid_copepods_Arctodiaptomus_alpinus_IMHOF_1885_and_Mixodiaptomus_laciniatus_LILLJEBORG_1889_in_a_small_high-altitude_lake
exp(1.953) * 0.85^2.399 #4.77 ug
#"Cyclops abyssorum", 1.4 mm, 12.60 ug
#"Cyclops sp", NA
#"Megacyclops gigas", 1.98 mm, 40.55 ug 
#"Mesocyclops leuckarti", 0.94 mm, 17.00 ug 
#"Cyclops vicinus", 0.154 - 0.176 mm, ave is 0.165, 18.30 ug
#"Cyclops bohater", 2.3 mm, 42.70 ug 
#"Diacyclops sp", 0.92 mm, 2.600 ug 
#"Mixodiaptomus sp", 14-26 ug (depends on sex, ave 27), https://www.researchgate.net/publication/236146469_Relationship_Between_NP_Ratio_and_Growth_Rate_During_the_Life_Cycle_of_Calanoid_Copepods_An_in_situ_Measurement 
#"Thermocyclops sp", 9.883 ug (oithonoides), 9.705 ug (crassus) https://pmc.ncbi.nlm.nih.gov/articles/PMC5546047/ (which thermos in the Alps)

#Average zoop weight across 4 dominant taxa
Zoop_Weight_ug = sum(1.7,12.60,12.60,18.30)/4 

zoop_length_mm = sum(0.65+1.4+1.4+0.165)/4

#Daylight Hours
daylight_hours <- ZoopsLightTemp$hoursdaylight

## Fish length-weight model - https://onlinelibrary.wiley.com/doi/full/10.1111/eff.12498 Pothoven 2019 
## W(g) = a * L(mm)^b

LW_a <- 0.00000188
LW_b <- 3.296


##### Initial conditions #####
ZoopsLightTemp = ZoopsLightTemp %>% filter(year <= 2004)

n_days <- nrow(ZoopsLightTemp) #number of days in df from 1976-2024; however there is a gap in temp data from 2005-2010

Weight_g <- numeric(n_days + 1)
Weight_mg <- numeric(n_days)

Fish_Length <- numeric(n_days + 1)

growth_g <- numeric(n_days)
growth_mg <- numeric(n_days)
growth_rate <- numeric(n_days)

consumption_J_day <- numeric(n_days)
assimilated_J_day <- numeric(n_days)
respiration_J_day <- numeric(n_days)
net_growth_J_day <- numeric(n_days)

##### Feeding function #####

feeding_model <- function(length_mm,
                          weight_g,
                          prey_density,
                          daylight_hours,
                          temp_C) {
  
  
  # Fish swimming speed
  speed <- length_mm / 1000   # m/s
  
  # Visual acuity
  acuity_deg <- 0.0167 * exp(9.14 -
                               2.4 * log(length_mm) +
                               0.229 * log(length_mm)^2)
  
  acuity_rad <- acuity_deg * pi / 180
  
  # Reactive distance
  RD_mm <- zoop_length_mm /(2 * tan(acuity_rad / 2))
  
  RD_m <- RD_mm / 1000
  
  
  # Search volume
  SV <- speed * pi * RD_m^2 * 0.5
  
  # Encounter rate
  encounter_rate <- SV * prey_density        # prey/s
  
  # Capture success - Anneville et al. 2010 (depends on fish predator length)
  capture_success <-
    ifelse(
      length_mm < 15,
      0.153824449,
      0.546789446
    )
  
  # Successful attack rate
  capture_rate <- encounter_rate *capture_success       # prey/s
  
  feeding_time_sec <- daylight_hours * 3600 #time in seconds
  
  captured_available <- capture_rate * feeding_time_sec #neglible handling_time, don't include Type II 
  
  
  ## Stomach capacity limitation - Brett 1971; Koski & Johnson 2002
  stomach_capacity_ug <-(weight_g * (14.1 - 4.95 * log(weight_g)) /100) * 1e6 
  
  evac_rate <- abs(-0.01958333 * temp_C + 0.03433333) #Karjalainen et al. 1990 (The gastric evacuation rate of vendace (Coregonus albula L. ) larvae predating on zooplankters in the laboratory) 
  
  stomach_biomass_limit <- stomach_capacity_ug * evac_rate * daylight_hours
  
  stomach_prey_limit <- stomach_biomass_limit / Zoop_Weight_ug
  
  ## Actual prey consumed - prey consumed cannot surpass stomach limits
  captured <-min(captured_available, stomach_prey_limit)
  
  # Convert prey to energy
  ug_captured <- captured * Zoop_Weight_ug
  
  J_consumed <- ug_captured * PreyED
  
  return(J_consumed)
}

##### Metabolism function #####

metabolism_model <- function(weight_mg,
                             weight_g,
                             temp_C){
  
  # Respiration coefficients - Huuskonen 1998
  RA <- 0.00584
  RB <- -0.05341
  RQ <- 0.05060
  
  
  # Oxygen consumption  - Huuskonen 1998 (how much oxygen a fish consumes per gram of body weight per day)
  R <- RA * weight_mg^RB * exp(RQ * temp_C)
  
  # Convert oxygen use to energy - energy expenditure (Joules consumed by metabolism per fish per day)
  respiration_J <- R * weight_g * Oxy_cal
  
  return(respiration_J)
}


##### Growth function #####

growth_model <- function(J_consumed,
                         respiration){
  
  assimilated <- J_consumed * Ae
  
  net_energy <- assimilated - respiration
  
  growth_g <- net_energy / Fish_energy_density
  
  return(list(
    assimilated = assimilated,
    net_energy = net_energy,
    growth_g = growth_g
  ))
}


######### Run daily simulation ########

# Initial fish length
Fish_Length[1] <- 10 #mm


for(i in seq_len(n_days)){
  
  #current fish size
  Weight_g[i] <- LW_a * Fish_Length[i]^LW_b
  
  Weight_mg[i] <- Weight_g[i] * 1000
  
  # Environmental conditions
  
  Temp_today <- ZoopsLightTemp$interpolated_temp[i]
  
  prey_density <- ZoopsLightTemp$interpolated_zoopValue[i] #already in correct units (density m3)
  
  daylight <-ZoopsLightTemp$hoursdaylight[i]
  
  
  # Feeding
  consumption_J_day[i] <-
    feeding_model(
      length_mm = Fish_Length[i],
      weight_g = Weight_g[i],
      prey_density = prey_density,
      daylight_hours = daylight,
      temp_C = Temp_today
    )
  
  
  # Respiration
  respiration_J_day[i] <-
    metabolism_model(
      weight_mg = Weight_mg[i],
      weight_g = Weight_g[i],
      temp_C = Temp_today
    )
  
  
  # Growth
  growth_output <-
    growth_model(
      J_consumed = consumption_J_day[i],
      respiration = respiration_J_day[i]
    )
  
  assimilated_J_day[i] <- growth_output$assimilated
  
  net_growth_J_day[i] <- growth_output$net_energy
  
  
  growth_g[i] <- growth_output$growth_g
  
  
  growth_mg[i] <- growth_g[i] * 1000
  
  
  growth_rate[i] <- (growth_g[i] / Weight_g[i]) * 100
  
  
  # Update size
  Weight_g[i+1] <- Weight_g[i] + growth_g[i]
  
  # Convert weight back to length
  Fish_Length[i+1] <-(Weight_g[i+1] / LW_a)^(1/LW_b)
  
  # Stop if fish exceeds 25mm
  if(Fish_Length[i+1] >= 25){
    Fish_Length[(i+2):(n_days+1)] <- NA
    Weight_g[(i+2):(n_days+1)] <- NA
    break
  }
  
}

## Results

results <- data.frame(
  
  Date = ZoopsLightTemp$date,
  
  Temperature = ZoopsLightTemp$interpolated_temp,
  
  Fish_length_mm = Fish_Length[-(n_days+1)],
  
  Fish_weight_g = Weight_g[-(n_days+1)],
  
  Food_energy_J = consumption_J_day,
  
  Assimilated_energy_J = assimilated_J_day,
  
  Respiration_J = respiration_J_day,
  
  Growth_energy_J = net_growth_J_day,
  
  Growth_g = growth_g,
  
  Growth_percent_mass_day = growth_rate
  
)


results


##################
#### Daily fish growth simulation
# Fish start at 10 mm on every hatch day from Feb 1 - Apr 30 (1974–2024)
# Uses the same mechanistic bioenergetics model as above
# Simulation ends when fish reaches 25 mm or after max_days

run_growth_sim <- function(hatch_day,
                           env,
                           start_length = 10,
                           max_days = 200) {
  
  n_env   <- nrow(env)
  n_steps <- min(max_days, n_env - hatch_day + 1)
  
  Pred_Length <- numeric(n_steps + 1)
  Weight_g    <- numeric(n_steps + 1)
  
  growth_g      <- rep(NA, n_steps)
  growth_rate   <- rep(NA, n_steps)
  consumption_J <- rep(NA, n_steps)
  respiration_J <- rep(NA, n_steps)
  
  Pred_Length[1] <- start_length
  Weight_g[1]    <- LW_a * start_length^LW_b
  
  recruited   <- FALSE
  days_elapsed <- NA_integer_
  
  for(k in seq_len(n_steps)) {
    
    i <- hatch_day + k - 1
    
    ## Current fish size
    Weight_g[k]  <- LW_a * Pred_Length[k]^LW_b
    Weight_mg    <- Weight_g[k] * 1000
    
    ## Environmental conditions
    
    #as a test, let's keep temp always constant and set it to the mean interpolated temp value
    #Temp_today <- mean(GenevaZoopsFull$interpolated_temp)
    Temp_today     <- env$interpolated_temp[i]
    #as a test, let's keep prey density always constant and set it to the mean prey density value
    #prey_density   <- mean(GenevaZoopsFull$density_m3)
    prey_density   <- env$interpolated_zoopValue[i]
    daylight_today <- env$hoursdaylight[i]
    
    ## Feeding
    consumption_J[k] <- feeding_model(
      length_mm      = Pred_Length[k],
      weight_g       = Weight_g[k],
      prey_density   = prey_density,
      daylight_hours = daylight_today,
      temp_C         = Temp_today
    )
    
    ## Respiration
    respiration_J[k] <- metabolism_model(
      weight_mg = Weight_mg,
      weight_g  = Weight_g[k],
      temp_C    = Temp_today
    )
    
    ## Growth
    growth_output <- growth_model(
      J_consumed  = consumption_J[k],
      respiration = respiration_J[k]
    )
    
    growth_g[k] <- growth_output$growth_g
    growth_rate[k] <- (growth_g[k] / Weight_g[k]) * 100
    
    ## Update weight
    Weight_g[k + 1] <- Weight_g[k] + growth_g[k]
    
    ## Prevent impossible negative weights
    if(Weight_g[k + 1] <= 0){
      break
    }
    
    ## Convert weight back to length
    Pred_Length[k + 1] <- (Weight_g[k + 1] / LW_a)^(1 / LW_b)
    
    ## Recruitment
    if(Pred_Length[k + 1] >= 25){
      recruited   <- TRUE
      days_elapsed <- k
      break
    }
  }
  
  if(is.na(days_elapsed)){
    days_elapsed <- sum(!is.na(growth_g))
  }
  
  
  list(
    start_date         = env$date[hatch_day],
    hatch_day          = hatch_day,
    days_elapsed       = days_elapsed,
    recruited          = recruited,
    mean_growth_g_day  = mean(growth_g[1:days_elapsed], na.rm = TRUE), 
    mean_growth_mg_day = 1000 * mean(growth_g[1:days_elapsed], na.rm = TRUE),#Average miligrams gained per day over the whole simulated period for a given hatch day
    mean_growth_mm_day = (Pred_Length[days_elapsed + 1] -
                            Pred_Length[1]) / days_elapsed,
    final_length       = max(Pred_Length, na.rm = TRUE),
    Temperature        = Temp_today
  )
}

###########################################################
## Run simulations for every hatch date from Feb 1–Apr 30

max_days <- 200

start_days <- which(
  month(ZoopsLightTemp$date) %in% c(3,4,5,6)
)

cohort_results <- vector("list", length(start_days))


for(j in seq_along(start_days)){
  
  cohort_results[[j]] <- run_growth_sim(
    hatch_day = start_days[j],
    env       = ZoopsLightTemp,
    max_days  = max_days
  )
  
}

###########################################################
cohort_summary <- bind_rows(lapply(cohort_results, as.data.frame)) %>%
  mutate(start_year = year(start_date), month = month(start_date), day = day(start_date)) 

cohort_summary$start_date <- as.character(cohort_summary$start_date) 

monthly_growth <- cohort_summary %>%
  group_by(month,start_year) %>%
  summarise(
    monthly_mean_days_to_recruitment = mean(days_elapsed),
    monthly_mean_growth_mm_day       = mean(mean_growth_mm_day),
    monthly_mean_growth_mg_day       = mean(mean_growth_mg_day),
    n_cohorts                = n()
  )


#Annual mean growth rate represents the average daily mass gain to recruitment across all recruited cohorts(hatch days) within a given year
yearly_growth <- cohort_summary %>%   
  filter(recruited) %>%   
  group_by(start_year) %>%   
  summarise(
    annual_mean_days_to_recruitment = mean(days_elapsed),
    annual_mean_growth_mm_day = mean(mean_growth_mm_day), 
    annual_mean_growth_mg_day = mean(mean_growth_mg_day),
    n_cohorts = n(),
    sd = sd(mean_growth_mm_day), 
    se_growth_mm_day = sd(mean_growth_mm_day) / sqrt(n_cohorts),
    sd_recruit = sd(days_elapsed),
    se_annual_days_to_recruitment = sd(days_elapsed) / sqrt(n_cohorts)
  )


ggplot(yearly_growth, aes(x = start_year, y = annual_mean_growth_mm_day)) + geom_line() + theme_minimal() +
  labs(x = "Year" , y = "Average Growth mm/day")

ggplot(yearly_growth, aes(x = start_year, y = annual_mean_days_to_recruitment)) + geom_line() + theme_minimal() +
  labs(x = "Year", y = "Average days to recruitment")

ggplot(monthly_growth,
       aes(x = month,
           y = monthly_mean_growth_mm_day, group = as.factor(start_year), color = as.factor(start_year))) +
  geom_line() +
  geom_point() +
  labs(
    x = "month",
    y = "Mean daily growth rate (mm/day)"
  ) +
  theme_minimal()

ggplot(monthly_growth,
       aes(x = month,
           y = monthly_mean_days_to_recruitment, group = as.factor(start_year), color = as.factor(start_year))) +
  geom_line() +
  geom_point() +
  labs(
    x = "month",
    y = "Mean days to recruitment (25 mm)"
  ) +
  theme_minimal()


