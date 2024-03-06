
# Testing
names(acadplots)
acad1 <- acadplots[1,]
acad1
acad1 <- download_daymet(lat = 44.36945, 
                         lon = -68.0704,
                         start = 1980,
                         end = 2022, 
                         internal = T, 
                         simplify = T)
table(acad1$measurement)

plot = "ACAD-001"

acad1_df <- acad1 |> pivot_wider(names_from = measurement, values_from = value) |> 
  rename(dayl_s = dayl..s.,  prcp_mmday = prcp..mm.day., 
         srad_wm2 = srad..W.m.2.,  swe_kgm2 = swe..kg.m.2., 
         tmax_c = tmax..deg.c., tmin_c = tmin..deg.c.,  vp_pa = vp..Pa.) |> 
  mutate(phx_mean = (tmax_c + tmin_c)/2, 
         srad_mjm2 = dayl_s/1000000, 
         origin = as.Date(paste0(year, "-01-01"),tz = "UTC") -1 ,
         date = as.Date(yday, origin = origin, tz = "UTC"),
         yearmonth = as.Date(format(date, format = "%Y-%m-01")),
         site = plot)

acad_month <- acad1_df %>%
  group_by(yearmonth) %>%
  summarise(latitude = latitude[1],
            month_sum = sum(phx_mean),
            monthmt_max = max(tmax_c),
            monthmt_mean = mean(phx_mean),
            srad_mean = mean(srad_wm2),
            srad_mean_mjm2 = mean(srad_mjm2),
            monthmt_min = min(tmin_c),
            monthpre = sum(prcp_mmday)) %>%
  ungroup() 

# From: https://search.r-project.org/CRAN/refmans/SPEI/html/Potential-evapotranspiration.html and
# https://daymet.ornl.gov/overview need to use mjm2 for solar radiation metric in hargreaves
acad_1 <- acad_month[1,]

acad_month <- acad_month |> 
  mutate(
    PEThar = SPEI::hargreaves(
      Tmin = acad_month$monthmt_min, 
      Tmax = acad_month$monthmt_max, 
      lat = acad_month$latitude,
      Ra = acad_month$srad_mean_mjm2,
      Pre = acad_month$monthpre),
    BAL = monthpre - PEThar)

acadSPEI <- data.frame(acad_month, SPEI_1 = SPEI::spei(acad_month$BAL, 1)$fitted)
head(acadSPEI)

