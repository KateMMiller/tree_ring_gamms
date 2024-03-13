#--------------------------------------------
# Final data compilation for tree ring gamms
#--------------------------------------------

library(tidyverse)
library(zoo) # for rollmean

# Simplifying dataset to only include variables of interest
# -- Choosing DayMet-only derived climate data
fulld <- read.csv("./data/ACAD_plot_core_climate_full_data.csv") |> filter(Year >= 1980) |> 
  filter(!is.na(coreID))

fulld <- fulld |> mutate(forest_type = 
                           case_when(Group_1 == "Laurentian - Acadian Acidic Swamp" ~ "FSWP",
                                     Group_1 == "Northern Hardwood - Hemlock Hardwood Forest" ~ "NHWD",
                                     Group_1 == "North-Central Appalachian & Laurentian Rocky Outc*" ~ "OUTC",
                                     Group_1 == "Red Spruce - Fir Forest" ~ "SSF",
                                     Group_1 == "SSF- Aspen & Birch Phase" ~ "SSFA")) 
names(fulld)

simpd <- fulld |> select(Plot_Name, Year, Unit = ParkSubUnit, coreID, species,
                         RRWmm, BAIcm2,
                         X = xCoordinate, Y = yCoordinate,
                         Physio = PhysiographySummary, Aspect, forest_type, elev_m, 
                         Pct_Rock, Pct_Bryophyte, PlotSlope, Stand_Structure, Pct_Crown_Closure,
                         live_stems_ha, liveBA_m2ha, horizon_depth, soilpH, pctTN, pctTC, Ca_Al,
                         C_N, BaseSat, AlSat, northiness, eastness, fire1947, 
                         Crown_Class, est_age, DBH = DBH_num, BA_pct_lg,
                         NO3, SO4, pH,
                         dm_ppt_sum_1:dm_ppt_sum_12, 
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
         tmin_wint = pmin(tmin_1, tmin_2, tmin_3, 
                          tmin_10_lag, tmin_11_lag, tmin_12_lag),
         tmax_gs = pmax(tmax_4, tmax_5, tmax_6, tmax_7, tmax_8, tmax_9),
         ppt_gs = ppt_4 + ppt_5 + ppt_6 + ppt_7 + ppt_8 + ppt_9) |> 
  ungroup() |> 
  select(Plot_Name:pH, ppt_4:ppt_9, ppt_gs, tmax_gs, tmin_wint, 
         SPEI01_1:SPEI03_12)
head(simpd2)
names(simpd2)
?pmin

# Calculate a rolling average to capture cumulative effects of previous dry or wet years
var = "SPEI01_3"
roll = 5

roll_fun <- function(dat = simpd2, var, roll){
  #xcol <- paste0(var, "_roll", roll)
  dat1 <- dat |> select(Plot_Name, coreID, Year, all_of(var)) |> 
    arrange(Plot_Name, coreID, Year) |> 
    group_by(Plot_Name, coreID) |> 
    mutate(lag_x = dplyr::lag(!!sym(var), n = 1)) |> 
    mutate("{var}_roll{roll}" := 
             zoo::rollmean(lag_x, k = roll, 
                      align = 'right', # calcs previous window 
                      fill = NA)) |> 
    ungroup() |> select(-lag_x) |> 
    data.frame()
  col = c(paste0(var, "_roll", roll))
  dat2 <- data.frame(dat1[,ncol(dat1)])
  names(dat2) <- col
  return(dat2)
}

lag_fun <- function(dat, x, lag){
  dat1 <- dat |> select(Plot_Name, coreID, Year, all_of(x)) |> 
    arrange(Plot_Name, coreID, Year) |> 
    group_by(Plot_Name, coreID) |> 
    mutate("{x}_lag{lag}" := dplyr::lag(!!sym(x), n = lag)) |> 
    ungroup() |> select(-all_of(x)) |> data.frame()
  col = c(paste0(x, "_lag", lag))
  dat2 <- data.frame(dat1[,ncol(dat1)])
  names(dat2) <- col
  return(dat2)
}

# head(roll_fun(simpd2, "SPEI03_9", 2))
# head(lag_fun(simpd2, "SPEI03_9", 1))
# head(simpd2$SPEI03_9)

vars <- c("SPEI01_4", "SPEI01_5", "SPEI01_6", "SPEI01_7", 
          "SPEI01_8", "SPEI01_9", "SPEI01_10", 
          "SPEI03_4", "SPEI03_5", "SPEI03_6", "SPEI03_7", 
          "SPEI03_8", "SPEI03_9", "SPEI03_10",
          "NO3", "SO4", "pH", "tmin_wint", "tmax_gs", "ppt_gs", 
          "ppt_8", "ppt_4")
rolls <- c(2:5)

roll_df <- data.frame(vars = rep(vars, each = length(rolls)),
                      rolls = rep(rolls, length(vars)))


core_rolls <- cbind(simpd2, ##|> select(Plot_Name, coreID, Year),
                    map2(roll_df$vars, roll_df$rolls,
                       ~roll_fun(simpd2, var = .x, roll = .y)) |> list_cbind())

core_lagrolls <- cbind(core_rolls, 
                       map(vars, 
                           ~lag_fun(core_rolls, var = ., lag = 1)) |> 
                         list_cbind())

ggplot(core_lagrolls |> filter(Plot_Name == "ACAD-101"),
       aes(x = Year, y = SPEI03_8_roll5)) + geom_point() + geom_smooth()


table(core_lagrolls$Crown_Class)
table(core_lagrolls$Year) #1980 to 2022
hist(core_lagrolls$Year)
write.csv(core_lagrolls, "./data/ACAD_final_core_data_20240311.csv", row.names = F)


