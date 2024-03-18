#------------------------
# Code used to compile plot-level covariates for tree growth gams
#------------------------
library(tidyverse)
library(forestNETN)
library(readxl) # importing core data 
library(sf) # for extracting elevation from forest plot shapefile
#library(lubridate) 

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
trees_plot <- joinTreeData(park = "ACAD", from = 2018, to = 2021, status = "live") %>% 
  group_by(Plot_Name) %>% summarize(live_stems_ha = sum(num_stems, na.rm = T)*(10000/225),
                                    liveBA_m2ha = sum(BA_cm2, na.rm = T)/225)

#------ Pull in elevation data, which came from https://coast.noaa.gov/dataviewer/ ------
forest_plots <- st_read("./data/ACAD_forest_plots_2018_elev.shp")
forest_plots_df <- data.frame(forest_plots) %>% mutate(Plot_Name = paste0("ACAD-", Plot)) %>% 
  select(Plot_Name, Group_1, elev_m = RASTERVALU, fire47 = X1947_fire)

#---- Combine NETN plot-level data sets ----
full_plot_data <- purrr::reduce(list(plot_info, forest_plots_df, stand_info, trees_plot, 
                                     soil_info_fix), 
                                left_join, by = "Plot_Name") |> unique() |> 
  # Convert aspect to something usable
  mutate(aspect_fac = case_when(Aspect == 0 ~ "NONE", 
                                between(Aspect, 1, 90) ~ "NE",
                                between(Aspect, 271, 360) ~ "NW",
                                between(Aspect, 91, 180) ~ "SW",
                                between(Aspect, 181, 270) ~ "SE",
                                TRUE ~ NA_character_),
         northiness = ifelse(aspect_fac == "NONE", 0, round(cos(Aspect*pi/180), 3)),
         eastness = ifelse(aspect_fac == "NONE", 0, round(sin(Aspect*pi/180), 3)), 
         fire1947 = ifelse(fire47 == 'burned', 1, 0), 
         forest_type = case_when(Group_1 == "Laurentian - Acadian Acidic Swamp" ~ "FSWP",
                                     Group_1 == "Northern Hardwood - Hemlock Hardwood Forest" ~ "NHWD",
                                     Group_1 == "North-Central Appalachian & Laurentian Rocky Outc*" ~ "OUTC",
                                     Group_1 == "Red Spruce - Fir Forest" ~ "SSF",
                                     Group_1 == "SSF- Aspen & Birch Phase" ~ "SSF")) |> 
         select(-fire47) |> select(-Group_1)

full_plot_data$fire1947[is.na(full_plot_data$fire1947)] <- 0

# Expand so records for each year in analysis
plot_years <- expand.grid(Plot_Name = unique(full_plot_data$Plot_Name), Year = 1980:2020)
plot_data_comb <- full_join(plot_years, full_plot_data, by = c("Plot_Name"))

table(complete.cases(full_plot_data)) #176 TRUE
head(full_plot_data)
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
bai_rrw_comb <- left_join(BAI_comb2, RRW_comb, by = c("Year", "coreID", "species", "Plot_Name", "Lat", "Long",
                                                      "Species", "Crown_Class", "est_age", "DBH_num"))

bai_rrw_comb$Crown_Class[bai_rrw_comb$Crown_Class == "Codom"] <- "3"
bai_rrw_comb$Crown_Class[bai_rrw_comb$Crown_Class == "Inter"] <- "4"
bai_rrw_comb$Crown_Class[bai_rrw_comb$Crown_Class == "Open Grown"] <- "1"

head(bai_rrw_comb)
table(bai_rrw_comb$Crown_Class)
write.csv(bai_rrw_comb, "./data/BAI_RRW_data_long.csv", row.names = F)

plot_data_comb <- left_join(full_plot_data, bai_rrw_comb, by = c("Plot_Name")) |> filter(Year >= 1980)
# add missing years, for climate variables

write.csv(plot_data_comb, "ACAD_Plot_RRW_BAI_data_1980_2020.csv", row.names = F)
