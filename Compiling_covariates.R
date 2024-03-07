#------------------------
# Code used to compile plot-level covariates for tree growth gams
#------------------------
library(tidyverse)
library(forestNETN)
library(readxl)
library(sf) # for extracting elevation from forest plot shapefile
library(mgcv) # for imputing missing deposition data via gams
library(zoo)
library(stars) # for extracting prism and spei data
library(starsExtra) # for nearest neighbor fixes
library(lubridate)
library(prism)
library(daymetr)
library(SPEI)
#devtools::install_github("WillemMaetens/standaRdized")
library(standaRdized) # for batch calculating SPEI

importData()

#---- Compile NETN forest data to have 1 record for each plot ----
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
forest_plots <- st_read("./data/ACAD_forest_plots_2018_elev.shp")
forest_plots_df <- data.frame(forest_plots) %>% mutate(Plot_Name = paste0("ACAD-", Plot)) %>% 
                   select(Plot_Name, Group_1, elev_m = RASTERVALU, fire47 = X1947_fire)

#---- Combine NETN plot-level data sets ----
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

#---- Compile data with 1 row per year in each core ----
#------ Core meta data ------
acad_core_meta1 <- read_excel(
  "./data/ACAD_tree_core_data_2013-2021.xlsx") %>% 
  arrange(coreID) %>% select(Plot_Name, coreID, Species, Crown_Class, DBH, est_age = `estimated #_Rings_USFS`)

head(acad_core_meta1)
length(unique(acad_core_meta1$Plot_Name)) #149 plots, missing ones appear to me M-A-R
length(unique(acad_core_meta1$coreID)) # 447
table(acad_core_meta1$Plot_Name)

table(acad_core_meta1$Species) # all spelled correctly

acad_core_meta1$DBH[acad_core_meta1$DBH == "NA"] <- NA
acad_core_meta1$DBH_num <- as.numeric(acad_core_meta1$DBH)

# Add Long/lat to dataset
acad_core_meta <- left_join(acad_core_meta1, 
                         VIEWS_NETN$Plots_NETN |> select(Plot_Name, Lat, Long), 
                         by = "Plot_Name")

#------ Read in and compile ring width data ------
path = "./data/Acadia increment core data 20240228.xlsx"
BAI <- excel_sheets(path)[grep("BAI", excel_sheets(path))]
RRW <- excel_sheets(path)[grep("RRW", excel_sheets(path))]

# Basal Area Increment Data
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
table(BAI_comb$Year)

# Raw ring width data
RRW_data <- RRW %>% map_dfr(function(sheet){
  df = read_xlsx(path = path, sheet = sheet) %>% #mutate(species = substr(sheet, 1, 4)) %>% 
    select(Year, everything())
  df_long <- df %>% pivot_longer(-Year, names_to = 'coreID', values_to = 'RRWmm') %>% 
    mutate(species = substr(sheet, 1, 4))
}) %>% filter(!is.na(RRWmm)) %>% arrange(coreID, Year)

length(unique(RRW_data$coreID)) #420

# link ring width data to tree core metadata
RRW_comb <- left_join(RRW_data, acad_core_meta |> select(-DBH), by = "coreID")
head(RRW_data)

RRW_check <- RRW_comb |> group_by(Plot_Name, coreID) |> 
  summarize(missing1 = sum(ifelse(is.na(DBH_num), 1, 0)), 
            missing = ifelse(missing1 > 0, 1, 0)) 

sum(RRW_check$missing) #9 cores missing from acad_core_meta

table(RRW_comb$Crown_Class, useNA = 'always')

write.csv(RRW_data, "./data/RRW_tree_data_long.csv", row.names = F)

# Include measure of competition following Salas-Eljatib & Weiskittel. 2020. 10.1016/j.foreco.2020.118369.
# Using 2015-2018 because that covers most cores
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
intersect(names(BAI_comb2), names(RRW_comb))
bai_rrw_comb <- left_join(BAI_comb2, RRW_comb, by = c("Year", "coreID", "species", "Plot_Name",
                                                      "Species", "Crown_Class", "est_age", "DBH_num"))

bai_rrw_comb$Crown_Class[bai_rrw_comb$Crown_Class == "Codom"] <- "3"
bai_rrw_comb$Crown_Class[bai_rrw_comb$Crown_Class == "Inter"] <- "4"
bai_rrw_comb$Crown_Class[bai_rrw_comb$Crown_Class == "Open Grown"] <- "1"

table(bai_rrw_comb$Crown_Class)

write.csv(bai_rrw_comb, "./data/BAI_RRW_data_long.csv", row.names = F)

#---- Read in and compile climate data ----
# list of ACAD plots
acadplots <- VIEWS_NETN$Plots_NETN |> filter(ParkUnit == "ACAD") |> 
  select(Plot_Name, long = Long, lat = Lat)

acad_sf <- st_as_sf(acadplots, coords = c("long", "lat"), crs = 4326)

bndbox <- st_as_sfc(st_bbox(acad_sf)) #xmin: -74.56284 ymin: 40.74542 xmax: -68.0704 ymax: 44.40644

#------ PRISM ------
# Annoyingly, PRISM makes you set the WD to download to the folder you want
# Commented out because only have to run it once to download
# orig_wd <- getwd()
# prismfile <- ("./data/prism")
# setwd(prismfile)
# options(prism.path = getwd())
# get_prism_monthlys("ppt", years = 1974:2022, mon = 1:12)
# get_prism_monthlys("tmin", years = 1974:2022, mon = 1:12)
# get_prism_monthlys("tmax", years = 1974:2022, mon = 1:12)
# setwd(orig_wd)

options(prism.path = "./data/prism/")
all_files <- list.files("./data/prism/") 
# Set up vectors to filter and iterate on
years <- 1974:2022
mos <- c("01", "02", "03", "04" ,"05", "06", 
         "07", "08", "09", "10", "11", "12")
yearmos_vec <- paste0(rep(years, each = 12), rep(mos, times = length(years)))
yearmos <- paste0(yearmos_vec, collapse = "|")
# yearmos then can filter out the year only .bil

# Function to compile/extract prism values for each forest plot
#type = 'ppt'; year = 1980; month = 4

comb_prism_year <- function(type, year, month){
  files <- all_files[!grepl(".zip", all_files) & grepl(type, all_files)]
  year_files <- files[grepl(year, files)]
  mos <- year_files[grepl(yearmos, year_files)]
  
  ras1 <- read_stars(paste0("./data/prism/", mos[month], "/", mos[month], ".bil"))
  bbox2 <- st_transform(bndbox, crs = st_crs(ras1))
  ras1b <- st_crop(ras1, bbox2)
  
  plots_ras <- cbind(acad_sf, st_extract(ras1b, st_coordinates(acad_sf)))
  col <-  paste0(type, "_", str_pad(month, 2, "0", side = 'left'))
  cols <- c("Plot_Name", col, "geometry")
  colnames(plots_ras) <- cols
  
  # Take nearest neighbor value for plots missing values
  miss <- st_transform(plots_ras[is.na(plots_ras[,col]),], crs = st_crs(acad_sf))
  acad_comp <- plots_ras |> filter(!Plot_Name %in% miss$Plot_Name)
  neigh <- st_nearest_feature(miss, acad_comp)
  neigh_plots <- st_join(miss |> select(Plot_Name_miss = Plot_Name, geometry), 
                         acad_comp, join = st_nearest_feature)
  
  
  plots_comb <- rbind(acad_comp |> st_drop_geometry(), 
                       neigh_plots |> select(Plot_Name = Plot_Name_miss, all_of(col)) |> 
                         st_drop_geometry()) |> 
                       arrange(Plot_Name)
  
  plots_final <- plots_comb[!is.na(plots_comb$Plot_Name),]
  
  return(data.frame(plots_final))
}

# For plots located on grid cells without prism data, take nearest neighbor value.

ppt <- 
purrr::map_dfr(1974:2022, function(yr){
  comb_prism_year(type = 'ppt', year = yr, month = 1:12) |> mutate(year = yr)
})

# test <- comb_prism_year(type = 'ppt', year = 1975, month = 1:12)
#length(1974:2022) * 176 

head(ppt)

write.csv(ppt, "./data/ppt_new_approach.csv", row.names = F)

tmin <- 
  purrr::map_dfr(1974:2022, function(yr){
    comb_prism_year(type = 'tmin', year = yr, month = 1:12) |> mutate(year =  yr) 
  })

write.csv(tmin, "./data/tmin_new_approach.csv", row.names = F)

tmax <- 
  purrr::map_dfr(1974:2022, function(yr){
    comb_prism_year(type = 'tmax', year = yr, month = 1:12) |> mutate(year =  yr) 
  })

write.csv(tmax, "./data/tmax_new_approach.csv", row.names = F)

clims <- list(ppt, tmin, tmax)
prism_comb <- reduce(clims, left_join, by = c("Plot_Name", "year"))

write.csv(prism_comb, "./data/prism_data_comb.csv", row.names = F)

# DayMet data
# create csv for batch download
acad_daymet_df <- acadplots |> select(site = Plot_Name, latitude = lat, longitude = long)
write.csv(acad_daymet_df, "./data/acadplots_for_daymet.csv", row.names = F)

# batch download data for each plot- commented out because only run once
# download_daymet_batch(file_location = "./data/acadplots_for_daymet.csv",
#                       start = 1980, 
#                       end = 2022, 
#                       path = "./data/daymet",
#                       simplify = TRUE,
#                       internal = FALSE)

# read in and combine the data into 1 large data frame
# need to iterate reading in the files, drop the deader text on the top, and make data wide by month
comb_daymet_by_mon <- function(plot){

  df <- read.table(paste0("./data/daymet/", plot, "_1980_2022.csv"), 
                   sep = ",", skip = 7, header = T) |> 
    mutate(Plot_Name = plot) |> select(Plot_Name, everything())

  colnames(df) <- c("Plot_Name", "year", "yday", "dayls", 'ppt', 'srad', 'swe', 'tmax', 'tmin', 'vp')
  df_sum <- df |> mutate(date = as.Date(yday-1, origin = "1980-01-01"),
                             month = month(date)) |> 
    group_by(Plot_Name, year, month) |> 
    summarise(dm_ppt_sum = sum(ppt, na.rm = T),
              dm_no_ppt_days = sum(ppt < 0.1),
              dm_srad_mjm2 = mean((srad*dayls)/1000000),
              dm_tmax = max(tmax, na.rm = T),
              dm_tmin = min(tmin, na.rm = T),
              dm_vpmax = max(vp, na.rm = T),
              dm_vpavg = mean(vp, na.rm = T), 
              .groups = 'drop')

 df_wide <- df_sum |> pivot_wider(names_from = month, 
                                  values_from = c(dm_ppt_sum, dm_no_ppt_days, dm_srad_mjm2,
                                                  dm_tmax, dm_tmin, dm_vpmax, dm_vpavg))
return(df_wide)
}

# test <- comb_daymet(plot = "ACAD-002")
# View(test)
plots <- sort(unique(acadplots$Plot_Name)) 
daymet_comb <- map_dfr(plots, ~comb_daymet_by_mon(.))
nrow(daymet_comb) == length(1980:2022) * 176 # 7586; TRUE
daymet_final <- inner_join(acadplots, data.frame(daymet_comb), by = "Plot_Name")
write.csv(daymet_final, "./data/daymet_climate_data_1980_2022.csv")

#------ SPEI ------
#-- SPEI by hand using dayment data --
comb_daymet_long <- function(plot){
  
  dat <- read.table(paste0("./data/daymet/", plot, "_1980_2022.csv"), 
                   sep = ",", skip = 7, header = T) |> 
    mutate(Plot_Name = plot,
           origin = as.Date(paste0(year, "-01-01"),tz = "UTC") -1 ,
           date = as.Date(yday, origin = origin, tz = "UTC"), 
           month = month(date)) |> 
    select(Plot_Name, everything())
  
  names(dat)
  colnames(dat) <- c("Plot_Name", "year", "yday", "dayls", 'ppt', 'srad', 'swe', 'tmax', 'tmin', 'vp',
                     'origin', 'date', 'month')
  dat$srad_mjm2 <- (dat$srad * dat$dayls/1000000)
  dat2 <- inner_join(acadplots, dat, by = "Plot_Name")
  return(dat2)
  }

dm_comb <- map_dfr(plots, ~comb_daymet_long(.))

dm_month <- dm_comb |> group_by(Plot_Name, year, month) |> 
  summarize(lat = first(lat),
            tmax = max(tmax),
            tmin = min(tmin),
            srad = mean(srad_mjm2),
            ppt = sum(ppt), 
            .groups = 'drop')

head(dm_month)
# data <- data.frame(tmax = tmax_mon[,1], tmin = tmin_mon[,1], ppt = ppt_mon[,1], srad = srad_mon[,1])
dm_month$PEThar <- hargreaves(Tmin = dm_month$tmin,
                              Tmax = dm_month$tmax,
                              lat = dm_month$lat,
                              Ra = dm_month$srad,
                              Pre = dm_month$ppt)

dm_month$BAL <- dm_month$ppt - dm_month$PEThar

# dm_ts <- ts(dm_month[, c("lat", "ppt", "tmax", "tmin", "srad_mjm2", "PEThar", "BAL")],
#             start = c(1980, 01, 01), end = c(2022, 12, 31), frequency = 365)
# summary(dm_ts)

SPEI_dat <- data.frame(dm_month, 
                       SPEI01 = spei(dm_month$BAL, scale = 1)$fitted,
                       SPEI03 = spei(dm_month$BAL, scale = 3)$fitted,
                       SPEI06 = spei(dm_month$BAL, scale = 6)$fitted,
                       SPEI12 = spei(dm_month$BAL, scale = 12)$fitted)
head(SPEI_dat)
SPEI_wide <- SPEI_dat |> select(Plot_Name, year, month, SPEI01:SPEI12) |> 
  pivot_wider(names_from = month, values_from = c(SPEI01, SPEI03, SPEI06, SPEI12))
head(SPEI_wide)

#-- Online SPEI --
# download nc datasets from web:
#   https://spei.csic.es/spei_database/#map_name=spei01#map_position=1463
# URL for SPEI01

# Function to crop and extract values for ACAD plots from the spei database
#spei = "spei01"
extract_spei <- function(spei){
  ras <- read_stars(paste0("./data/spei_nc/", spei, ".nc"))
  spei_crop <- st_crop(ras, bndbox)
  names(spei_crop) <- "spei"
  spei_ex <- st_extract(spei_crop, acad_sf)
  spei_df <- as.data.frame(spei_ex, xy = T)
  spei_df$long <- st_coordinates(spei_df$geometry)[,1]
  spei_df$lat <- st_coordinates(spei_df$geometry)[,2]
  spei_yr <- spei_df |> filter(time >= "1974-01-01") # make dataset smaller
  
  spei_plots <- left_join(spei_yr, acadplots, by = c("long", "lat")) |> 
    mutate(year = year(time), month = month(time),
           value = as.numeric(spei)) |> 
    select(Plot_Name, long, lat, value, year, month)  
  
  spei_wide <- spei_plots |> 
    pivot_wider(names_from = month, values_from = value,
                names_prefix = paste0(toupper(spei), "_"))
}

spei01 <- extract_spei("spei01")
spei02 <- extract_spei("spei02")
spei03 <- extract_spei("spei03")
spei06 <- extract_spei("spei06")
spei12 <- extract_spei("spei12")
spei36 <- extract_spei("spei36")

spei_list <- list(spei01, spei02, spei03, spei06, spei12, spei36)
spei_all <- purrr::reduce(spei_list, full_join, 
                          by = c("Plot_Name", "long", "lat", "year"))
head(SPEI_wide)
write.csv(spei_all, "./data/ACAD_SPEI_1974-2022.csv", row.names = F)



#---- Read in deposition data ----
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
core_comb <- bai_rrw_comb %>% filter(!is.na(Year))
head(BAI_comb2)
length(unique(core_comb$coreID)) #420

core_clim <- full_join(core_comb, climate_comb, by = 'Year') %>% 
  filter(!is.na(ppt_mm_1)) %>%  # drops early years without climate data; time series starts at 1895
  filter(!is.na(Plot_Name)) # drops years with climate data but no cores

core_clim_dep <- full_join(core_clim, dep_final, by = "Year")
names(core_clim_dep)
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
names(plot_core_clim)

write.csv(plot_core_clim, "./data/Full_plot_core_climate_deposition_dataset.csv", row.names = F)

#plot_core_clim <- read.csv("./data/Full_plot_core_climate_deposition_dataset.csv")

# this is a working dataset for the gam as of 2022. 

# Now to add rolling averages and lags
head(plot_core_clim)
table(plot_core_clim$Crown_Class)

core_simp <- plot_core_clim |> 
  mutate(pcID = paste0(Plot_Name, "_", coreID)) |> 
  select(pcID, Plot_Name, ParkSubUnit, coreID, Year, 
         BAIcm2, RRWmm,
         xCoordinate:fire1947, species:BA_pct_lg,
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

names(core_lagrolls)
table(core_lagrolls$Crown_Class)
table(core_lagrolls$Year) #1975 to 2020
write.csv(core_lagrolls, "./data/core_data_all_lagrolls_20240305.csv", row.names = F)

