#------------------------
# Code used to compile plot-level covariates for tree growth gams
#------------------------
library(tidyverse)
library(forestNETN)
library(readxl)
library(sf) # for extracting elevation from forest plot shapefile
library(mgcv) # for imputing missing deposition data via gams
library(zoo)

importData()

#++++++++++++ NEXT STEPS 20240124 +++++++++++#
# Figure out whether to do a rolling 5-year average for climate and atm dep metrics,
# or a 5-year lag
# Also lag all climate variables by 1 year
# Rolling average will take into account previous and current year

#---- Compile data that have 1 record for each plot ----
#------ Plot-level data, including XY coords ------
plot_info <- joinLocEvent(park = "ACAD", from = 2018, to = 2021) %>% 
  select(Plot_Name, ParkSubUnit, xCoordinate, yCoordinate, PhysiographySummary, Aspect)

#------ Stand data from most recent visit ------
stand_info <- joinStandData(park = "ACAD", from = 2018, to = 2021) %>% 
  select(Plot_Name, Pct_Rock, Txt_Rock, Pct_Bryophyte, Txt_Bryophyte,
         PlotSlope, Stand_Structure, Stand_Structure_Code,
         Pct_Crown_Closure)

#------ Soil data from cycle 2, the most complete sample of data. Have to pull from cycle 1 for 029. ------
soil_info <- joinSoilLabData(park = "ACAD", from = 2010, to = 2013) %>%  # cycle 2: best coverage of all plots
  filter(Horizon_QC == "O") %>% 
  select(Plot_Name, soilYear = SampleYear, horizon_depth, soilpH, pctTN, pctTC, 
         Ca_Al, C_N, BaseSat, AlSat)
  
soil_plots <- left_join(plot_info %>% select(Plot_Name), soil_info, by = "Plot_Name")
miss_soil_plots <- soil_plots$Plot_Name[is.na(soil_plots$soilYear)]
miss_soil_plots # All but 029 are because of A horizons instead of O. Will use A where O isn't available.
# Have to grab cycle 1 record for 029 because we're missing cycle 2 data.

miss_soil_info <- joinSoilLabData(park = "ACAD", from = 2007, to = 2019) %>%  # cycle 2: best coverage of all plots
  filter(Plot_Name %in% miss_soil_plots) %>% 
  select(Plot_Name, soilYear = SampleYear, Horizon_QC, horizon_depth, soilpH, pctTN, pctTC, 
         Ca_Al, C_N, BaseSat, AlSat) %>% group_by(Plot_Name) %>% 
  slice(1) %>% select(-Horizon_QC)

soil_info_fix <- rbind(soil_info, miss_soil_info)

#------ Pull in tree density and BA ------
trees <- joinTreeData(park = "ACAD", from = 2018, to = 2021, status = "live") %>% 
         group_by(Plot_Name) %>% summarize(live_stems_ha = sum(num_stems, na.rm = T)*(10000/225),
                                           liveBA_m2ha = sum(BA_cm2, na.rm = T)/225)

#------ Pull in elevation data, which came from https://coast.noaa.gov/dataviewer/ ------
forest_plots <- st_read("C:/NETN/Monitoring_Projects/Forest_Health/Forest_plots_2010/ACAD_forest_plots_2018_elev.shp")
forest_plots_df <- data.frame(forest_plots) %>% mutate(Plot_Name = paste0("ACAD-", Plot)) %>% 
                   select(Plot_Name, Group_1, elev_m = RASTERVALU, fire47 = X1947_fire)

#---- Combine plot-level data sets ----
full_plot_data <- purrr::reduce(list(plot_info, forest_plots_df, stand_info, trees, soil_info_fix), 
                                left_join, by = "Plot_Name") %>% 
             # Convert aspect to something usable
                  mutate(aspect_fac = case_when(Aspect == 0 ~ "NONE", 
                                                between(Aspect, 1, 90) ~ "NE",
                                                between(Aspect, 271, 360) ~ "NW",
                                                between(Aspect, 91, 180) ~ "SW",
                                                between(Aspect, 181, 270) ~ "SE",
                                                TRUE ~ NA_character_),
                         northiness = ifelse(aspect_fac == "NONE", 0, round(cos(Aspect*pi/180), 3)),
                         eastness = ifelse(aspect_fac == "NONE", 0, round(sin(Aspect*pi/180), 3)), 
                         fire1947 = ifelse(fire47 == 'burned', 1, 0)) %>% select(-fire47)

full_plot_data$fire1947[is.na(full_plot_data$fire1947)] <- 0

table(complete.cases(full_plot_data)) #176 TRUE
head(full_plot_data)

#---- Compile data with 1 row per year in each core ----
#------ Core meta data ------
acad_core_meta <- read_excel("C:/NETN/collaborators/Schaberg/Core_Data/ACAD_tree_core_data_2013-2021.xlsx") %>% 
  arrange(coreID) %>% select(Plot_Name, coreID, Species, Crown_Class, DBH, est_age = `estimated #_Rings_USFS`)

head(acad_core_meta)
length(unique(acad_core_meta$Plot_Name)) #149
length(unique(acad_core_meta$coreID)) # 352

table(acad_core_meta$Species) # all spelled correctly

acad_core_meta$DBH[acad_core_meta$DBH == "NA"] <- NA
acad_core_meta$DBH_num <- as.numeric(acad_core_meta$DBH)

#------ Read in and compile ring width data ------
path = "C:/NETN/collaborators/Schaberg/Acadia increment core data.xlsx"
BAI <- excel_sheets(path)[grep("BAI", excel_sheets(path))]
RRW <- excel_sheets(path)[grep("RRW", excel_sheets(path))]

BAI_data <- BAI %>% map_dfr(function(sheet){
                df = read_xlsx(path = path, sheet = sheet) %>% #mutate(species = substr(sheet, 1, 4)) %>% 
                     select(Year, everything())
                df_long <- df %>% pivot_longer(-Year, names_to = 'coreID', values_to = 'BAIcm2') %>% 
                                  mutate(species = substr(sheet, 1, 4))
}) %>% filter(!is.na(BAIcm2)) %>% arrange(coreID, Year)

length(unique(BAI_data$coreID)) #420

# link ring width data to tree core metadata
BAI_comb <- left_join(BAI_data, acad_core_meta |> select(-DBH), by = "coreID")
head(BAI_data)

BAI_check <- BAI_comb |> group_by(Plot_Name, coreID) |> 
  summarize(missing1 = sum(ifelse(is.na(DBH_num), 1, 0)), 
            missing = ifelse(missing1 > 0, 1, 0)) 

sum(BAI_check$missing) #7 cores missing from acad_core_meta

table(BAI_comb$Crown_Class, useNA = 'always')

# Include measure of competition following Salas-Eljatib & Weiskittel. 2020. 10.1016/j.foreco.2020.118369.
trees <- joinTreeData(from = 2015, to = 2018, status = 'live') |> 
  select(Plot_Name, BA_cm2, DBHcm)

BAI_tree1 <- left_join(BAI_comb, trees, by = "Plot_Name", relationship = 'many-to-many') 
BAI_tree2 <- BAI_tree1 |> 
  mutate(DBH_lg = ifelse(DBH_num <= DBHcm, DBHcm, 0),
         BAm2ha_lg = (((DBH_lg/2)^2)*pi)/225,
         BAm2ha_tot = (((DBHcm/2)^2)*pi)/225
)

BAI_tree <- BAI_tree2 |> 
  group_by(Plot_Name, coreID, Crown_Class, species, est_age) |> 
  summarize(DBH_orig = first(DBH_num),
            BA_pct_lg = (sum(BAm2ha_lg) / sum(BAm2ha_tot))
  )

BAI_comb2 <- left_join(BAI_comb, BAI_tree, by = c("coreID", "species", "Plot_Name", "Crown_Class", "est_age"))

write.csv(BAI_comb2, "./data/BAI_tree_data_long.csv", row.names = F)

RRW_data <- RRW %>% map_dfr(function(sheet){
  df = read_xlsx(path = path, sheet = sheet) %>% #mutate(species = substr(sheet, 1, 4)) %>% 
    select(Year, everything())
  df_long <- df %>% pivot_longer(-Year, names_to = 'coreID', values_to = 'RRWmm') %>% 
    mutate(species = substr(sheet, 1, 4))
}) %>% filter(!is.na(RRWmm)) %>% arrange(coreID, Year)

write.csv(RRW_data, "./data/RRW_tree_data_long.csv", row.names = F)

#------ Read in climate data ------
# PRISM
prism <- read.csv("C:/NETN/collaborators/Schaberg/PRISM_simp.csv") %>% 
  mutate(Year = as.numeric(substr(Date, 1, 4)),
         month = as.numeric(substr(Date, 6, 7))) %>% 
  select(Year, month, ppt_mm = ppt..mm., tminC = tmin..degrees.C., 
         tmeanC = tmean..degrees.C., tmaxC = tmax..degrees.C.) %>% 
  pivot_wider(names_from = month, values_from = c(ppt_mm, tminC, tmeanC, tmaxC))

# SPEI
spei <- read.csv("C:/NETN/collaborators/Schaberg/ACAD_SPEI_thru2018.csv") %>% 
  pivot_wider(names_from = Month, values_from = c(SPEI01, SPEI03))

climate_comb <- full_join(prism, spei, by = "Year")
table(climate_comb$Year)

write.csv(climate_comb, "./data/PRISM_SPEI_wide.csv", row.names = F)

#----- Read in deposition data -----
dep <- read.csv("C:/NETN/collaborators/Schaberg/NTN-me98-annual-mgl.csv") %>% select(Year = yr, NO3, SO4, pH)
dep$NO3[dep$NO3 == -9.0] <- NA
dep$SO4[dep$SO4 == -9.0] <- NA
dep$pH[dep$pH == -9.0] <- NA

# Use GAM to impute years with missing deposition data
modNO3 <- gam(NO3 ~ s(Year, k = 10), # thinplate spline is default
              data = dep, method = 'REML')
gam.check(modNO3) # edf < than k and p is high, looks okay

depNAs <- dep[is.na(dep$NO3),] %>% filter(Year > 1983)
depNAs <- depNAs %>% mutate(NO3 = predict(modNO3, depNAs))

ggplot(dep, aes(x = Year, y = NO3))+
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 10))+ 
  geom_point(data = depNAs, aes(x = Year, y = NO3), color = 'red')+
  geom_point(data = dep, aes(x = Year, y = NO3), color = 'blue') 

# There's a lot more variance in year to year dep, but there's at least
# a strong trend in the data that the imputed values are derived from. 

modSO4 <- gam(SO4 ~ s(Year, k = 10), # thinplate spline is default
              data = dep, method = 'REML')
gam.check(modSO4) # edf < than k and p is high, looks okay
plot(modSO4)

depNAs <- depNAs %>% mutate(SO4 = predict(modSO4, depNAs))

ggplot(dep, aes(x = Year, y = SO4))+
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 10))+ 
  geom_point(data = depNAs, aes(x = Year, y = SO4), color = 'red')+
  geom_point(data = dep, aes(x = Year, y = SO4), color = 'blue') 

# There's a lot more variance in year to year dep, but there's at least
# a strong trend in the data that the imputed values are derived from. 

modpH <- gam(pH ~ s(Year, k = 10), # thinplate spline is default
             data = dep, method = 'REML')
gam.check(modpH) # edf < than k and p is high, looks okay
plot(modpH)

depNAs <- depNAs %>% mutate(pH = predict(modpH, depNAs))

ggplot(dep, aes(x = Year, y = pH))+
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 10))+ 
  geom_point(data = depNAs, aes(x = Year, y = pH), color = 'red')+
  geom_point(data = dep, aes(x = Year, y = pH), color = 'blue') 

# There's less variance in year to year dep than the other metrics 

# Adding imputed values to the data, so more complete
dep_complete <- dep[complete.cases(dep),]
dep_final <- rbind(dep_complete, depNAs) %>% arrange(Year)

#---- Combine climate and core BAI data ----
core_comb <- BAI_comb2 %>% filter(!is.na(Year))
head(BAI_comb2)
length(unique(core_comb$coreID)) #420

core_clim <- full_join(core_comb, climate_comb, by = 'Year') %>% 
  filter(!is.na(ppt_mm_1)) %>%  # drops early years without climate data; time series starts at 1895
  filter(!is.na(Plot_Name)) # drops years with climate data but no cores

core_clim_dep <- full_join(core_clim, dep_final, by = "Year")

#---- Combine 9-month SPI 1895-present for Hancock County from drought.gov
drt <- read.csv("./data/Hancock_Cty_9mo_SPI-23009_1895-2023.csv")
drt <- drt |> mutate(year = as.numeric(substr(DATE, 3, 6)),
                     month = as.numeric(substr(DATE, 7, 8)))

head(drt)

# Data are April to October % of days in various drought or moisture states
gr_seas_drt <- drt |> filter(month %in% 4:10) |> 
  group_by(year) |> summarize(D0_sum = (sum(D0))/7,
                              D1_sum = (sum(D1))/7,
                              D2_sum = (sum(D2))/7,
                              D3_sum = (sum(D3))/7,
                              D4_sum = (sum(D4))/7,
                              W0_sum = (sum(W0))/7,
                              W1_sum = (sum(W1))/7,
                              W2_sum = (sum(W3))/7,
                              W3_sum = (sum(W4))/7,
                              W4_sum = (sum(W4))/7,
                              pct_D0_D4 = (D0_sum + D1_sum + D2_sum + D3_sum + D4_sum)/5,
                              pct_D1_D4 = (D1_sum + D2_sum + D3_sum + D4_sum)/4,
                              pct_D2_D4 = (D2_sum + D3_sum + D4_sum)/3,
                              pct_W0_W4 = (W0_sum + W1_sum + W2_sum + W3_sum + W4_sum)/5,
                              pct_W1_W4 = (W1_sum + W2_sum + W3_sum + W4_sum)/4,
                              pct_W2_W4 = (W2_sum + W3_sum + W4_sum)/3) |> 
  select(year, pct_D0_D4:pct_W2_W4)

#---- Combine all data ----
head(full_plot_data)
head(gr_seas_drt)
str(gr_seas_drt)
names(core_clim_dep)

core_clim_dep2 <- left_join(core_clim_dep, gr_seas_drt, by = c("Year" = "year"))

plot_core_clim <- full_join(full_plot_data, core_clim_dep2 |> select(-DBH_orig), by = "Plot_Name") %>% 
  filter(!is.na(coreID)) %>% data.frame()

table(complete.cases(plot_core_clim)) # 562 F, from Dep data that starts at 1984
head(plot_core_clim)

names(plot_core_clim)
write.csv(plot_core_clim, "./data/Full_plot_core_climate_deposition_dataset.csv", row.names = F)

plot_core_clim <- read.csv("./data/Full_plot_core_climate_deposition_dataset.csv")

# this is a working dataset for the gam as of 2022. 

# Now to add rolling averages and lags
names(plot_core_clim)

core_simp <- plot_core_clim |> 
  mutate(pcID = paste0(Plot_Name, "_", coreID)) |> 
  select(pcID, Plot_Name, ParkSubUnit, coreID, Year, 
         xCoordinate:fire1947, BAIcm2:BA_pct_lg,
         tmaxC_4:tmaxC_9,
         tminC_1:tminC_3, tminC_10:tminC_12,
         SPEI03_6:SPEI03_9, NO3:pct_W2_W4) |> 
  filter(Year > 1975) |> 
  arrange(Plot_Name, coreID, Year)

# Calculate minimum winter temperature- need to lead the months 10-12, b/c they're after the growing season
# that the core is measured.
core_simp <- core_simp |> 
  group_by(Plot_Name, coreID) |> 
  mutate(tminC_10_lag = dplyr::lag(tminC_10, 1),
         tminC_11_lag = dplyr::lag(tminC_11, 1),
         tminC_12_lag = dplyr::lag(tminC_12, 1),
         tmin_wint = pmin(tminC_1, tminC_2, tminC_3, tminC_10_lag, tminC_11_lag, tminC_12_lag),
         tmax_gs = pmax(tmaxC_4, tmaxC_5, tmaxC_6, tmaxC_7, tmaxC_8, tmaxC_9)) |> 
  ungroup() |> select(pcID:BA_pct_lg, SPEI03_6, SPEI03_9, NO3, SO4, pH, pct_D0_D4:pct_W2_W4,
                      tmin_wint, tmax_gs) |> 
  data.frame()

#dat$SPEI03_8_roll5 <- rollmean(stats::lag(dat$SPEI03_8, k = 5), k = 5, fill = NA)

# Calculate a rolling average to capture cumulative effects of previous dry or wet years
roll_fun <- function(dat, x, roll){
  dat1 <- dat |> select(Plot_Name, coreID, Year, all_of(x)) |> 
    arrange(Plot_Name, coreID, Year) |> 
    group_by(Plot_Name, coreID) |> 
    mutate("{x}_roll{roll}" := 
             rollmean(!!sym(x), k = roll, 
                      align = 'right', # calcs previous window 
                      fill = NA)) |> 
    ungroup() |> data.frame()
  col = c(paste0(x, "_roll", roll))
  dat2 <- data.frame(dat1[,ncol(dat1)])
  names(dat2) <- col
  return(dat2)
}

head(roll_fun(core_simp, "SPEI03_6", 5))

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

vars <- c("SPEI03_6", "SPEI03_9", "NO3", "SO4", "pH", "pct_D0_D4",
          "pct_D1_D4", "pct_D2_D4", "pct_W0_W4", "pct_W1_W4",
          "pct_W2_W4", "tmin_wint", "tmax_gs")
rolls <- c(2:5)
roll_df <- data.frame(vars = rep(vars, each = length(rolls)),
                      rolls = rep(rolls, length(vars)))
names(core_simp)
core_rolls <- cbind(core_simp,
                    purrr::map2_dfc(roll_df$vars, roll_df$rolls, ~roll_fun(dat = core_simp, x = .x, roll = .y))
              )

core_lagrolls <- cbind(core_rolls, purrr::map_dfc(vars, ~lag_fun(dat = core_rolls, ., 1)))

head(core_lagrolls)
write.csv(core_lagrolls, "./data/core_data_all_lagrolls_20240223.csv", row.names = F)

# Check that things roughly make rough sense (ie the joins were correct)
ggplot(plot_core_clim, aes(x = SPEI03_7, y = BAIcm2)) +
       geom_point()+
       geom_smooth()
       
ggplot(plot_core_clim, aes(x = BaseSat, y = BAIcm2)) +
  geom_point()+
  geom_smooth()

ggplot(plot_core_clim, aes(x = northiness, y = BAIcm2)) +
  geom_point()+
  geom_smooth()

ggplot(plot_core_clim, aes(x = soilpH, y = BAIcm2)) +
  geom_point()+
  geom_smooth()



