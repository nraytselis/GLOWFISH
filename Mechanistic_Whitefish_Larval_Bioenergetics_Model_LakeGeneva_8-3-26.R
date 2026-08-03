setwd("~/Desktop/Whitefish/LakeGeneva")
library(readr)
library(dplyr)
library(tidyr)
library(lubridate) #add rows for missing dates
library(zoo) #interpolate values 


#For each model day, we will use larval size and water temperature to model consumption, respiration (losses), and growth 
#For each lake and year, we will run the model until the larval fish is recruited (25mm metamorphasis)

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


## Average Daphnia length (mm)
zoop_length_mm_df <- GenevaZoopLength %>%
  group_by(taxa_name_orig) %>%
  summarise(mean_length = mean(length_mm))

zoop_length_mm <- zoop_length_mm_df$mean_length[4] #select specifically the daphnia 


# Gastric evacuation rate - Karjalainen et al. 1990
#The gastric evacuation rate of vendace (Coregonus albula L.) larvae predating on zooplankters in the laboratory
#Experimental data
df <- data.frame(
  Temp = c(6, 12, 18),
  R    = c(-0.045, -0.277, -0.280)
)

# Fit the linear models for relationship between a and Temp and R and Temp
model_R <- lm(R ~ Temp, data = df)

# 3. Extract the coefficients
m_R <- coef(model_R)["Temp"]      # Slope for R (-0.01958)
c_R <- coef(model_R)["(Intercept)"] # Intercept for R (0.03433)



####### Mechanistic Whitefish Larval Bioenergetics Model

##### Parameters #####

# Daphnia energy density
PreyED <- 967.5 / 1e6        # J/ug wet weight

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


#Zooplankton length-weight model (Daphnia) - Watkins et al. 2011
Zoop_Weight_ug <- exp(1.468) * zoop_length_mm^2.829 

#Daylight Hours
daylight_hours <- GenevaZoopsFull$hoursdaylight

## Fish length-weight model - https://onlinelibrary.wiley.com/doi/full/10.1111/eff.12498 Pothoven 2019 
## W(g) = a * L(mm)^b

LW_a <- 0.00000188
LW_b <- 3.296


##### Initial conditions #####

n_days <- nrow(GenevaZoopsFull) #number of days in df from 1974-2024

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
  
  
  # Oxygen consumption  - Huuskonen 1998
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

  Temp_today <- GenevaZoopsFull$interpolated_temp[i]
  
  prey_density <- GenevaZoopsFull$density_m3[i]
  
  daylight <-GenevaZoopsFull$hoursdaylight[i]
  
  
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
  
  Date = GenevaZoopsFull$date,
  
  Temperature = GenevaZoopsFull$interpolated_temp,
  
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
    Temp_today     <- env$interpolated_temp[i]
    prey_density   <- env$density_m3[i]
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
    mean_growth_mg_day = 1000 * mean(growth_g[1:days_elapsed], na.rm = TRUE),
    mean_growth_mm_day = (Pred_Length[days_elapsed + 1] -
                            Pred_Length[1]) / days_elapsed,
    final_length       = max(Pred_Length, na.rm = TRUE)
  )
}

###########################################################
## Run simulations for every hatch date from Feb 1–Apr 30

max_days <- 200

start_days <- which(
  month(GenevaZoopsFull$date) %in% c(2, 3, 4)
)

cohort_results <- vector("list", length(start_days))


for(j in seq_along(start_days)){
  
  cohort_results[[j]] <- run_growth_sim(
    hatch_day = start_days[j],
    env       = GenevaZoopsFull,
    max_days  = max_days
  )
  
}


###########################################################
cohort_summary <- bind_rows(lapply(cohort_results, as.data.frame)) %>%
  mutate(start_year = year(start_date), month = month(start_date), day = day(start_date)) 

cohort_summary$start_date <- as.character(cohort_summary$start_date) 


monthly_growth <- cohort_summary %>%
  filter(recruited) %>%
  group_by(month,start_year) %>%
  summarise(
    mean_days_to_recruitment = mean(days_elapsed),
    mean_growth_mm_day       = mean(mean_growth_mm_day),
    mean_growth_mg_day       = mean(mean_growth_mg_day),
    n_cohorts                = n()
  )

yearly_growth <- cohort_summary %>%
  filter(recruited) %>%
  group_by(start_year) %>%
  summarise(
    mean_days_to_recruitment = mean(days_elapsed),
    mean_growth_mm_day       = mean(mean_growth_mm_day),
    mean_growth_mg_day       = mean(mean_growth_mg_day),
    n_cohorts                = n()
  )

#select a few years to more easily compare - 1974, 1984, 1994, 2004, 2014, 2024

cohort_summary_select = cohort_summary %>% filter(start_year %in% c("1974", "1984", "1994", "2004", "2014", "2024"))

monthly_growth_select <- cohort_summary_select %>%
  filter(recruited) %>%
  group_by(month,start_year) %>%
  summarise(
    mean_days_to_recruitment = mean(days_elapsed),
    mean_growth_mm_day       = mean(mean_growth_mm_day),
    mean_growth_mg_day       = mean(mean_growth_mg_day),
    n_cohorts                = n()
  )

###########################################################
ggplot(monthly_growth,
       aes(x = month,
           y = mean_growth_mm_day, group = as.factor(start_year), color = as.factor(start_year))) +
  geom_line() +
  geom_point() +
  labs(
    x = "month",
    y = "Mean daily growth rate (mm/day)"
  ) +
  theme_minimal()

ggplot(monthly_growth,
       aes(x = month,
           y = mean_days_to_recruitment, group = as.factor(start_year), color = as.factor(start_year))) +
  geom_line() +
  geom_point() +
  labs(
    x = "month",
    y = "Mean days to recruitment (25 mm)"
  ) +
  theme_minimal()


ggplot(yearly_growth,
       aes(x = start_year,
           y = mean_growth_mm_day)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Mean daily growth rate (mm/day)"
  ) +
  theme_minimal()

ggplot(yearly_growth,
       aes(x = start_year,
           y = mean_days_to_recruitment)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Mean days to recruitment (25 mm)"
  ) +
  theme_minimal()


#select plots
ggplot(monthly_growth_select,
       aes(x = month,
           y = mean_growth_mm_day, group = as.factor(start_year), color = as.factor(start_year))) +
  geom_line() +
  geom_point() +
  labs(
    x = "month",
    y = "Mean daily growth rate (mm/day)"
  ) +
  theme_minimal()

ggplot(monthly_growth_select,
       aes(x = month,
           y = mean_days_to_recruitment, group = as.factor(start_year), color = as.factor(start_year))) +
  geom_line() +
  geom_point() +
  labs(
    x = "month",
    y = "Mean days to recruitment (25 mm)"
  ) +
  theme_minimal()


