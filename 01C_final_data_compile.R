#--------------------------------------------
# Final data compilation for tree ring gamms
#--------------------------------------------

library(tidyverse)
library(lubridate)
library(zoo) # for rollmean

# Simplifying dataset to only include variables of interest
# -- Choosing DayMet-only derived climate data
fulld <- read.csv("./data/ACAD_plot_core_climate_full_data.csv") |> filter(year >= 1980) |> 
  filter(!is.na(coreID))

names(fulld)

simpd <- fulld |> select(Plot_Name, 
                         Year = year, Unit = ParkSubUnit, coreID, species,
                         RRWmm, BAIcm2,
                         X = xCoordinate, Y = yCoordinate,
                         Physio = PhysiographySummary, Aspect, forest_type, elev_m, 
                         Pct_Rock, Pct_Bryophyte, PlotSlope, Stand_Structure, Pct_Crown_Closure,
                         Avg_Height_Codom, Avg_Height_Inter, live_stems_ha, liveBA_m2ha, 
                         horizon_depth, soilpH, pctTN, pctTC, Ca_Al,
                         C_N, BaseSat, AlSat, northiness, eastness, fire1947, 
                         Crown_Class, est_age, DBH = DBH_num, BA_pct_lg,
                         NO3, SO4, pH,
                         dm_ppt_sum_1:dm_ppt_sum_12, 
                         # pr_ppt_01:pr_tmax_12,
                         #dm_no_ppt_days_1:dm_no_ppt_days_12,
                         dm_tmax_mean_1:dm_tmax_mean_12, 
                         dm_tmin_mean_1:dm_tmin_mean_12,
                         dm_SPEI01_1:dm_SPEI03_12) 

# Clean up by dropping dm_, _sum and _mean  from daymet variables
names(simpd) <- gsub("dm_", "", names(simpd))
names(simpd) <- gsub("_sum", "", names(simpd))
names(simpd) <- gsub("_mean", "", names(simpd))

# Now to add rolling averages and lags
table(simpd$Crown_Class)

# Calculate minimum winter temperature- need to lead the months 10-12, b/c they're after the growing season
# that the core is measured.
names(simpd)
simpd2 <- simpd |> 
  group_by(Plot_Name, coreID) |> 
  mutate(tmin_10_lag = dplyr::lag(tmin_10, 1),
         tmin_11_lag = dplyr::lag(tmin_11, 1),
         tmin_12_lag = dplyr::lag(tmin_12, 1),
         # tmin_10_lag_pr = dplyr::lag(pr_tmin_10, 1),
         # tmin_11_lag_pr = dplyr::lag(pr_tmin_11, 1),
         # tmin_12_lag_pr = dplyr::lag(pr_tmin_12, 1),
         tmin_wint = pmin(tmin_1, tmin_2, tmin_3, 
                          tmin_10_lag, tmin_11_lag, tmin_12_lag),
         # tmin_wint_pr = pmin(pr_tmin_01, pr_tmin_02, pr_tmin_03, 
         #                     tmin_10_lag_pr, tmin_11_lag_pr, tmin_12_lag_pr),
         tmax_gs = pmax(tmax_4, tmax_5, tmax_6, tmax_7, tmax_8, tmax_9),
         #tmax_gs_pr = pmax(pr_tmax_04, pr_tmax_05, pr_tmax_06, pr_tmax_07, pr_tmax_08, pr_tmax_09),
         #ppt_gs_pr = pr_ppt_04 + pr_ppt_05 + pr_ppt_06 + pr_ppt_07 + pr_ppt_08 + pr_ppt_09),
         ppt_gs = ppt_4 + ppt_5 + ppt_6 + ppt_7 + ppt_8 + ppt_9) |> 
  ungroup() |> 
  select(Plot_Name:pH, ppt_4:ppt_9, 
         # ppt_4_pr = pr_ppt_04, ppt_5_pr = pr_ppt_05, ppt_6_pr = pr_ppt_06,
         # ppt_7_pr = pr_ppt_07, ppt_8_pr = pr_ppt_08, ppt_9_pr = pr_ppt_09,
         ppt_gs, tmax_gs, tmin_wint, SPEI01_1:SPEI03_12)
names(simpd2)
#----- Calculate growing season length. -----
comb_daymet <- function(plot){
  
  df <- read.table(paste0("./data/daymet/", plot, "_1980_2022.csv"), 
                   sep = ",", skip = 7, header = T) |> 
    mutate(Plot_Name = plot) |> 
    select(Plot_Name, everything())
  colnames(df) <- c("Plot_Name", "year", "yday", "dayls", 'ppt', 
                    'srad', 'swe', 'tmaxc', 'tminc', 'vp')
  return(df)
}

plots <- sort(unique(simpd2$Plot_Name)) 
daymet_comb <- map(plots, ~comb_daymet(.)) |> list_rbind()
daymet_comb <- daymet_comb |> mutate(date = as.Date(yday-1, origin = paste0(year, "-01-01")),
                                     month = month(date),
                                     tmean = (tmaxc + tminc)/2) 

# Start of growing season 5C is triggered by first 6-day span of daily mean temps > 5 C, 
# with the start the first day of that 6-day stretch. The end of the growing season is the 
# first day of the first stretch of 6-days with daily mean temp < 5 C after July 1 (doy = 183)
calc_gsl_5c <- function(dat, yr, plot){
  gs <- data.frame(rle(dat$tmean > 5)[1], rle(dat$tmean > 5)[2])
  gs$consdays <- cumsum(gs$lengths)
  gs_start <- gs$consdays[which(gs$values == TRUE & 
                                  gs$lengths >= 6)[1]-1] + 1 #takes the previous row then adds 1 day
  gs_end <- gs$consdays[which(gs$consdays > 183 & 
                              gs$lengths >= 6 &
                              gs$values == FALSE)[1]-1] + 1
  if(is.na(gs_end)){# backup if less than 6 consecutive cold days
      gs_end <- gs$consdays[which(gs$consdays > 183 &
                                  gs$lengths >= 2 &
                                  gs$values == FALSE)[1]-1] + 1}
  gsl <- (gs_end - 5) - (gs_start - 5) 
  gsldata <- data.frame(Plot_Name = plot, year = yr, gs_5c_length = gsl)
  return(gsldata)
}

years <- 1980:2022
plot_years <- expand.grid(plots, years) |> select(Plot_Name = Var1, year = Var2) |> 
  arrange(Plot_Name, year)

gsl_5c_data <- 
  map2_dfr(plot_years$year, plot_years$Plot_Name, 
           function(yr, plt){
                           data1 <- daymet_comb |> filter(year == yr) |> filter(Plot_Name == plt)
                           gsl <- calc_gsl_5c(dat = data1, yr = yr, plot = plt)
                     }
)

# Start of growing season FF is number of days between last killing freeze of -2.2C in spring and first
# killing freeze in fall, based on tmin.
calc_gsl_ff <- function(dat, yr, plot){
  gs <- data.frame(rle(dat$tminc > -2.22)[1], rle(dat$tminc > -2.22)[2]) # killing frost temp
  gs$consdays <- cumsum(gs$lengths)
  # Find at least 30 cons. days without killing frost, then take the previous rowsum then add 1 day
  gs_start <- gs$consdays[which(gs$values == TRUE & 
                                gs$lengths >= 30)[1]-1] + 1 
  # Find the first killing frost after July 1
  gs_end <- gs$consdays[which(gs$consdays > 183 & 
                              gs$values == FALSE)[1]-1] + 1
  gsl <- (gs_end - 5) - (gs_start - 5) 
  gsldata <- data.frame(Plot_Name = plot, year = yr, gs_ff_length = gsl)
  return(gsldata)
}

gsl_ff_data <- 
  map2_dfr(plot_years$year, plot_years$Plot_Name, 
           function(yr, plt){
             data1 <- daymet_comb |> filter(year == yr) |> filter(Plot_Name == plt)
             gsl <- calc_gsl_ff(dat = data1, yr = yr, plot = plt)
           }
  )

ggplot(gsl_ff_data |> filter(Plot_Name == "ACAD-003"), aes(x = year, y = gs_ff_length)) +
  geom_point() + geom_smooth() + forestNETN::theme_FHM()

# Calculating number of days where max temp is > 5c in spring and fall to capture shoulder growing seasons
calc_5c_days <- function(dat, yr, plot){
  spring <- dat |> filter(yday < 121) # May 1
  fall <- dat |> filter(yday > 273) # Sept 30
  spring_5c <- data.frame(rle(spring$tmax > 5)[1], rle(spring$tmax > 5)[2]) |> filter(values == "TRUE")
  spring_5c_sum <- sum(spring_5c$lengths)
  
  fall_5c <- data.frame(rle(fall$tmax > 5)[1], rle(fall$tmax > 5)[2]) |> filter(values == "TRUE")
  fall_5c_sum <- sum(fall_5c$lengths)
  
  gsldata <- data.frame(Plot_Name = plot, year = yr, 
                        spring_5c_days = spring_5c_sum, fall_5c_days = fall_5c_sum)
  return(gsldata)
}


seas_5c_data <- 
  map2_dfr(plot_years$year, plot_years$Plot_Name, 
           function(yr, plt){
             data1 <- daymet_comb |> filter(year == yr) |> filter(Plot_Name == plt)
             gsl <- calc_5c_days(dat = data1, yr = yr, plot = plt)
           }
  )

ggplot(seas_5c_data |> filter(Plot_Name == "ACAD-001"), aes(x = year, y = spring_5c_days)) +
  geom_point() + geom_smooth() + forestNETN::theme_FHM()
ggplot(seas_5c_data |> filter(Plot_Name == "ACAD-001"), aes(x = year, y = fall_5c_days)) +
  geom_point() + geom_smooth() + forestNETN::theme_FHM()

# Calculating number of days where min temp is -30C or below
  # -30c comes from Schaberg et al. 2011
tmin30c <- daymet_comb |> group_by(Plot_Name, year) |> 
           summarize(num_neg30c_days = sum(tminc <= -30), 
                                       .groups = 'drop')
table(tmin30c$num_neg30c_days) # only 0 or 1, not enough range to model. Won't include in dataset
  
# Calculate number of frost/thaw events  
calc_spring_ft_events <- function(dat, yr, plot){ 
  dat <- dat |> filter(yday < 183) 
  
  ft <- data.frame(rle(dat$tminc > 0)[1], rle(dat$tminc > 0)[2]) 
  ft$consdays <- cumsum(ft$lengths)

# Summarize number of spring thaw (defined as > 5 days with tmin > 0), freeze cycles
# By finding values = TRUE (>0C) for 5 days, that are followed by a freeze. 
# Subtract 1 to drop the final growing season return to > 0 temps
  ft_events <- sum(ft$values[which(ft$values == TRUE & 
                                 ft$lengths >= 5)])-1

  ftdata <- data.frame(Plot_Name = plot, year = yr, spring_ft_events = ft_events)  
  return(ftdata)

}

spring_ft_data <- 
  map2_dfr(plot_years$year, plot_years$Plot_Name, 
           function(yr, plt){
             data1 <- daymet_comb |> filter(year == yr) |> filter(Plot_Name == plt)
             ft <- calc_spring_ft_events(dat = data1, yr = yr, plot = plt)
           }
  )


ggplot(spring_ft_data |> filter(Plot_Name == "ACAD-001"), aes(x = year, y = spring_ft_events)) +
  geom_point() + geom_smooth() + forestNETN::theme_FHM()

daily_ft <- 
  daymet_comb |> filter(yday < 183) |> #July 1  
  group_by(Plot_Name, year) |> 
  mutate(daily_thaw = ifelse(tminc < 0 & tmaxc > 0, 1, 0)) |> 
  summarize(daily_thaws = sum(daily_thaw), .groups = 'drop')

ggplot(daily_ft |> filter(Plot_Name == "ACAD-001"), aes(x = year, y = daily_thaws)) +
  geom_point() + geom_smooth() + forestNETN::theme_FHM()

# Not enough range in daily ft events or seasonal ft events to model. Not going to include in 
# final dataset.

gs_dfs <- list(gsl_5c_data, gsl_ff_data, seas_5c_data)
gs_data <- reduce(gs_dfs, full_join, by = c("Plot_Name", "year"))
write.csv(gs_data, "./data/ACAD_growing_season_metrics.csv", row.names = F)

# Calculate a rolling average to capture cumulative effects of previous dry or wet years
vars <- c("SPEI01_4", "SPEI01_5", "SPEI01_6", "SPEI01_7", 
          "SPEI01_8", "SPEI01_9", "SPEI01_10", 
          "SPEI03_4", "SPEI03_5", "SPEI03_6", "SPEI03_7", 
          "SPEI03_8", "SPEI03_9", "SPEI03_10",
          "NO3", "SO4", "pH", "tmin_wint", "tmax_gs", "ppt_gs", 
          "ppt_4", "ppt_5", "ppt_6", "ppt_7", "ppt_8",
          "gs_5c_length", "gs_ff_length", "spring_5c_days", "fall_5c_days")

simpd3 <- full_join(simpd2, gs_data, by = c("Plot_Name", "Year" = "year")) |> 
  filter(!is.na(coreID))
names(simpd3)

core_rolls <- simpd3 %>% #select(Plot_Name, coreID, Year, ppt_4) |> 
  group_by(Plot_Name, coreID) %>%  
  mutate(across(all_of(vars), 
                ~dplyr::lag(., 1), 
                .names = "{.col}_lag1")) %>% 
  mutate(across(ends_with("lag1"), 
                ~zoo::rollmean(., k = 2, 
                               align = 'right', # calcs previous window 
                               fill = NA),
                .names = '{.col}_roll2')) %>%
  mutate(across(ends_with("lag1"), 
                ~zoo::rollmean(., k = 3, 
                               align = 'right', # calcs previous window 
                               fill = NA),
                .names = '{.col}_roll3')) %>%
  mutate(across(ends_with("lag1"), 
                ~zoo::rollmean(., k = 4, 
                               align = 'right', # calcs previous window 
                               fill = NA),
                .names = '{.col}_roll4')) %>%
  mutate(across(ends_with("lag1"), 
                ~zoo::rollmean(., k = 5, 
                               align = 'right', # calcs previous window 
                               fill = NA),
                .names = '{.col}_roll5')) %>%

 data.frame()

names(core_rolls)
old_names <- names(core_rolls[1:37])
new_order <- sort(names(core_rolls[,38:ncol(core_rolls)]))
core_rolls2 <- core_rolls[,c(old_names, new_order)]

names(core_rolls2) <- gsub("_lag1_roll", "_roll", names(core_rolls2))

table(core_rolls2$Crown_Class)
table(core_rolls2$Year) #413 from 1980 to 2022
write.csv(core_rolls2, "./data/ACAD_final_core_data_20240315.csv", row.names = F)


