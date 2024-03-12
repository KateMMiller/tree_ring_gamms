#------------------------
# Code used to compile climate covariates for tree growth gams
#------------------------

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         Run 01A_compile_plot_data.R first 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(sf) # for extracting elevation from forest plot shapefile
library(mgcv) # for imputing missing deposition data via gams
library(stars) # for extracting prism and spei data
library(starsExtra) # for nearest neighbor fixes
library(lubridate) # date wrangling
library(prism) # Download historic daily gridded climate data 
library(daymetr) # Download daily gridded climate data starting at 1980
library(SPEI) # calculate SPEI using DayMet data
library(readr) # for rcc-acis data compile
library(httr) # for rcc-acis data download
library(jsonlite) # for rcc-acis data compile

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
  col <-  paste0('pr_', type, "_", str_pad(month, 2, "0", side = 'left'))
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

write.csv(ppt, "./data/prism_ppt_1974-2022.csv", row.names = F)

tmin <- 
  purrr::map_dfr(1974:2022, function(yr){
    comb_prism_year(type = 'tmin', year = yr, month = 1:12) |> mutate(year =  yr) 
  })

write.csv(tmin, "./data/prism_tmin_1974-2022.csv", row.names = F)

tmax <- 
  purrr::map_dfr(1974:2022, function(yr){
    comb_prism_year(type = 'tmax', year = yr, month = 1:12) |> mutate(year =  yr) 
  })

write.csv(tmax, "./data/prism_tmax_1974-2022.csv", row.names = F)

clims <- list(ppt, tmin, tmax)
prism_comb <- reduce(clims, left_join, by = c("Plot_Name", "year"))

write.csv(prism_comb, "./data/prism_data_1974-2022.csv", row.names = F)

#------  DayMet data ------ 
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
  
  colnames(df) <- c("Plot_Name", "year", "yday", "dayls", 'ppt', 
                    'srad', 'swe', 'tmax', 'tmin', 'vp')
  df_sum <- df |> mutate(date = as.Date(yday-1, origin = "1980-01-01"),
                         month = month(date)) |> 
    group_by(Plot_Name, year, month) |> 
    summarise(dm_ppt_sum = sum(ppt, na.rm = T),
              dm_no_ppt_days = sum(ppt < 0.1),
              dm_srad_mjm2 = mean((srad*dayls)/1000000),
              dm_tmax = max(tmax, na.rm = T),
              dm_tmax_mean = mean(tmax, na.rm = T),
              dm_tmin = min(tmin, na.rm = T),
              dm_tmin_mean = mean(tmin, na.rm = T),
              dm_vpmax = max(vp, na.rm = T),
              dm_vpavg = mean(vp, na.rm = T), 
              .groups = 'drop')
  
  df_wide <- df_sum |> pivot_wider(names_from = month, 
                                   values_from = c(dm_ppt_sum, dm_no_ppt_days, dm_srad_mjm2,
                                                   dm_tmax, dm_tmax_mean, 
                                                   dm_tmin, dm_tmin_mean,
                                                   dm_vpmax, dm_vpavg))
  return(df_wide)
}

# test <- comb_daymet(plot = "ACAD-002")
# View(test)
plots <- sort(unique(acadplots$Plot_Name)) 
daymet_comb <- map_dfr(plots, ~comb_daymet_by_mon(.))
nrow(daymet_comb) == length(1980:2022) * 176 # 7586; TRUE
daymet_final <- inner_join(acadplots, data.frame(daymet_comb), by = "Plot_Name")
write.csv(daymet_final, "./data/DayMet_climate_data_1980-2022.csv")

#------ SPEI ------
#-- SPEI by hand using DayMet data --
comb_daymet_long <- function(plot){
  
  dat <- read.table(paste0("./data/daymet/", plot, "_1980_2022.csv"), 
                    sep = ",", skip = 7, header = T) |> 
    mutate(Plot_Name = plot,
           origin = as.Date(paste0(year, "-01-01"),tz = "UTC") -1 ,
           date = as.Date(yday, origin = origin, tz = "UTC"), 
           month = month(date)) |> 
    select(Plot_Name, everything())
  
  names(dat)
  colnames(dat) <- c("Plot_Name", "year", "yday", "dayls", 'ppt', 
                     'srad', 'swe', 'tmax', 'tmin', 'vp',
                     'origin', 'date', 'month')
  dat$srad_mjm2 <- (dat$srad * dat$dayls/1000000)
  dat2 <- inner_join(acadplots, dat, by = "Plot_Name")
  return(dat2)
}

dm_comb <- map_dfr(plots, ~comb_daymet_long(.))

dm_month <- dm_comb |> group_by(Plot_Name, year, month) |> 
  summarize(lat = first(lat),
            tmax = max(tmax),
            tmax_mean = mean(tmax),
            tmin = min(tmin),
            tmin_mean = mean(tmin),
            srad = mean(srad_mjm2),
            ppt = sum(ppt), 
            .groups = 'drop')

dm_month$PEThar <- hargreaves(Tmin = dm_month$tmin,
                              Tmax = dm_month$tmax,
                              lat = dm_month$lat,
                              Ra = dm_month$srad,
                              Pre = dm_month$ppt)

dm_month$BAL <- dm_month$ppt - dm_month$PEThar

dm_month$PEThar_mean <- hargreaves(Tmin = dm_month$tmin_mean,
                                   Tmax = dm_month$tmax_mean,
                                   lat = dm_month$lat,
                                   Ra = dm_month$srad,
                                   Pre = dm_month$ppt)

dm_month$BAL_mean <- dm_month$ppt - dm_month$PEThar_mean

SPEI_dat <- data.frame(dm_month, 
                       dm_SPEI01 = spei(dm_month$BAL, scale = 1)$fitted,
                       dm_SPEI03 = spei(dm_month$BAL, scale = 3)$fitted,
                       dm_SPEI06 = spei(dm_month$BAL, scale = 6)$fitted,
                       dm_SPEI12 = spei(dm_month$BAL, scale = 12)$fitted,
                       dm_SPEI01m = spei(dm_month$BAL_mean, scale = 1)$fitted,
                       dm_SPEI03m = spei(dm_month$BAL_mean, scale = 3)$fitted,
                       dm_SPEI06m = spei(dm_month$BAL_mean, scale = 6)$fitted,
                       dm_SPEI12m = spei(dm_month$BAL_mean, scale = 12)$fitted)

SPEI_wide <- SPEI_dat |> select(Plot_Name, year, month, dm_SPEI01:dm_SPEI12m) |> 
  pivot_wider(names_from = month, 
              values_from = c(dm_SPEI01, dm_SPEI03, dm_SPEI06, dm_SPEI12,
                              dm_SPEI01m, dm_SPEI03m, dm_SPEI06m, dm_SPEI12m))
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
  spei_df <- as.data.frame(spei_ex, xy = F)
  spei_df$long <- st_coordinates(spei_df$geometry)[,1]
  spei_df$lat <- st_coordinates(spei_df$geometry)[,2]
  spei_yr <- spei_df |> filter(time >= "1974-01-01") |> st_drop_geometry() # make dataset smaller
  
  spei_plots <- inner_join(spei_yr, acadplots, by = c("long", "lat")) |> 
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
names(spei_all)
colnames(spei_all) <- c("Plot_Name", "long", "lat", "year", paste0("nc_", names(spei_all[,5:ncol(spei_all)])))

spei_all_comb <- left_join(spei_all, SPEI_wide, by = c("Plot_Name", "year")) |> 
  filter(year >= 1980)

write.csv(spei_all_comb, "./data/SPEI_DayMet_CSIC_1980-2022.csv", row.names = F)

# Note SPEI 
#   1. Not Drought > -0.5;  
#   2. Mild = -0.5 to -1 ; 
#   3. Moderate = -1.5 to -1;
#   4. Severe = -2 to -1.5
#   5. Extreme = < -2

#------ Weather Station Data ------
# ACAD MARS WS GHCN ID: USR0000MMCF
# ACAD_station_IDs: 170100 (Coop); USC00170100 (GHCN); BHRM1 (NWS LI)
# Adapted from https://github.com/mhlinder/weather-data-now/blob/master/fetch_temp.R

get_wdat <- function(field, stn, reduce = c("sum", "mean")){
  urlbase <- "http://data.rcc-acis.org/StnData?params=%s"
  elem_template <- function(field) {
    list(name = field,
         interval = "mly",
         duration = "mly",
         reduce = reduce)}

  params <-
    list(sid   = stn,  
         sdate = "por", #"1980-01-01",
         edate = "por", #"2022-12-21",
         elems = lapply(field, elem_template))

  query <- sprintf(urlbase, URLencode(toJSON(params, auto_unbox = TRUE)))

  response <- GET(query) |> content(as = "text", encoding = "UTF-8") |> fromJSON()
 
  wdat <-
    response$data %>%
    as.data.frame(stringsAsFactors = FALSE)

  names(wdat) <- c("Date", field)
  wdat[,field][wdat[,field] == "M"] <- NA_real_ # replace M with NA for missing data
  wdat[,field] <- as.numeric(wdat[,field])
  return(wdat)
}

pcp_in <- get_wdat('pcpn', stn = "BHRM1", reduce = 'sum')
tmax_f <- get_wdat('maxt', stn = "USC00170100", reduce = 'mean')
tmin_f <- get_wdat('mint', stn = "USC00170100", reduce = 'mean')


ws_ls <- list(pcp_in, tmax_f, tmin_f)
ws_comb <- reduce(ws_ls, full_join, by = "Date")
str(ws_comb)

ws_comb$ws_pcpmm <- ws_comb$pcpn * 25.4
ws_comb$ws_tmaxc <- (ws_comb$maxt - 32) * (5/9)
ws_comb$ws_tminc <- (ws_comb$mint - 32) * (5/9)

head(ws_comb)
str(ws_comb)

ws_wide <- ws_comb |> mutate(month = as.numeric(substr(Date, 6, 7)), 
                             year = as.numeric(substr(Date, 1, 4))) |> 
  select(year, month, ws_pcpmm, ws_tmaxc, ws_tminc) |> 
  pivot_wider(names_from = month, values_from = c(ws_pcpmm, ws_tmaxc, ws_tminc))
table(ws_wide$year)

write.csv(ws_wide, "./data/NOAA_weather_station_data_1980-2014.csv", row.names = F)

#---- Read in deposition data ----
# annual_dep <- read.csv("https://nadp2.slh.wisc.edu/datalib/ntn/cydep/NTN-ME98-cydep.csv")

# using annual precipitation-weighted mean concentrations for dep values
dep <- read.csv("https://nadp2.slh.wisc.edu/datalib/ntn/cy/NTN-ME98-cy.csv") |> 
  select(year = yr, NO3, SO4, pH) # Only goes up through 2020
summary(dep) # No -9 values in min, so don't need to replace with NAs
table(complete.cases(dep)) # All T
# dep$NO3[dep$NO3 == -9.0] <- NA
# dep$SO4[dep$SO4 == -9.0] <- NA
# dep$pH[dep$pH == -9.0] <- NA

#----- Combine datasets ------
names(prism_comb) #PlotName year
names(daymet_final) #PlotName, year lat/long
names(spei_all_comb) #PlotName, year lat/long
names(ws_wide) #year
names(dep) #year

clim_comb1 <- full_join(daymet_final, prism_comb, by = c("Plot_Name", "year"))
clim_comb2 <- full_join(clim_comb1, spei_all_comb, by = c("Plot_Name", "year", "lat", "long"))
clim_comb3 <- full_join(clim_comb2, ws_wide, by = "year")
clim_comb4 <- full_join(clim_comb3, dep, by = "year")
write.csv(clim_comb4, "./data/ACAD_climate_data.csv", row.names = F)

core_clim <- full_join(plot_data_comb, clim_comb4, by = c("Plot_Name", "Year" = "year"))
write.csv(core_clim, "./data/ACAD_plot_core_climate_full_data.csv", row.names = F)
