#------------------------------------------
# Comparing climate datasets to determine most appropriate for ACAD
#------------------------------------------

library(tidyverse)
clim_dat <- read.csv('./data/temp/ACAD_climate_data.csv')

#---- Compare weather station, prism and daymet data for near and far plots ----
ws_comp_39 <- clim_dat |> filter(between(year, 1985, 2014)) |> filter(Plot_Name == "ACAD-039") # closest to WS
ws_comp_64 <- clim_dat |> filter(between(year, 1985, 2014)) |> filter(Plot_Name == "ACAD-064") # IAH- farthest
ws_comp_36 <- clim_dat |> filter(between(year, 1985, 2014)) |> filter(Plot_Name == "ACAD-036") # SCH- farthest

ggplot(ws_comp_39, aes(x = year, y = ws_tmaxc_7)) + forestNETN::theme_FHM() + 
  geom_point(color = 'grey') +
  geom_point(aes(y = pr_tmax_07), color = "steelblue") +
  geom_point(aes(y = dm_tmax_mean_7), color = 'springgreen') +
  labs(y = "Max. Temp. C : July", subtitle = "WS: grey; DayMet: green; prism: blue")

ggplot(ws_comp_64, aes(x = year, y = ws_tmaxc_7)) + forestNETN::theme_FHM() + 
  geom_point(color = 'grey') +
  geom_point(aes(y = pr_tmax_07), color = "steelblue") +
  geom_point(aes(y = dm_tmax_mean_7), color = 'green') +
  labs(y = "Max. Temp. C : July", subtitle = "WS: grey; DayMet: green; prism: blue")

ggplot(ws_comp_36, aes(x = year, y = ws_tmaxc_7)) + forestNETN::theme_FHM() + 
  geom_point(color = 'grey') +
  geom_point(aes(y = pr_tmax_07), color = "dodgerblue") +
  geom_point(aes(y = dm_tmax_mean_7), color = 'orange') +
  labs(y = "Max. Temp. C : July", subtitle = "WS: grey; DayMet: orange; prism: blue")

range(ws_comp_39$ws_tmaxc_7, na.rm = T)
range(ws_comp_39$pr_tmax_07, na.rm = T)

tmax_39 <- 
ggplot(ws_comp_39, aes(x = ws_tmaxc_7, y = pr_tmax_07)) + forestNETN::theme_FHM() +
  geom_point(color = 'dodgerblue')+
  geom_point(aes(y = dm_tmax_mean_7), color = 'orange') +
  geom_abline(slope = 1, intercept = 0) + xlim(18, 30) + ylim(18, 30) +
  labs(y = "Gridded TMax C July", x = "Weather Station TMax C July",
       subtitle = "Plot closest to WS: Prism = blue; Daymet = orange")

tmax_36 <- 
ggplot(ws_comp_36, aes(x = ws_tmaxc_7, y = pr_tmax_07)) + forestNETN::theme_FHM() +
  geom_point(color = 'dodgerblue')+
  geom_point(aes(y = dm_tmax_mean_7), color = 'orange') +
  geom_abline(slope = 1, intercept = 0) + xlim(18, 30) + ylim(18, 30) +
  labs(y = "Gridded TMax C July", x = "Weather Station TMax C July", 
       subtitle = "Plot1 farthest from WS: Prism = blue; Daymet = orange")

tmax_64 <- 
ggplot(ws_comp_64, aes(x = ws_tmaxc_7, y = pr_tmax_07)) + forestNETN::theme_FHM() +
  geom_point(color = 'dodgerblue')+
  geom_point(aes(y = dm_tmax_mean_7), color = 'orange') +
  geom_abline(slope = 1, intercept = 0) + xlim(18, 30) + ylim(18, 30) +
  labs(y = "Gridded TMax C July", x = "Weather Station TMax C July",
       subtitle = "Plot2 farthest from WS: Prism = blue; Daymet = orange")

cowplot::plot_grid(tmax_39, tmax_36, tmax_64, nrow = 3)
# Going with DayMet because less separation from weather station at farther distances.
# Not sure I trust the smoothing algorithm behind PRISM enough to believe the differences.

names(ws_comp)

pcp_39 <- 
ggplot(ws_comp_39, aes(x = year)) + forestNETN::theme_FHM() +
  geom_point(aes(y = ws_pcpmm_7), color = 'dimgrey', alpha = 0.5) +
  geom_smooth(aes(y = ws_pcpmm_7), se = F, color = 'dimgrey') + 
  
  geom_point(aes(y = dm_ppt_sum_7), color = 'orange', alpha = 0.5) +
  geom_smooth(aes(y = dm_ppt_sum_7, x = year), se = F, color = 'orange') + 
  
  geom_point(aes(y = pr_ppt_07), color = 'dodgerblue', alpha = 0.5) +
  geom_smooth(aes(y = pr_ppt_07, x = year), se = F, color = 'dodgerblue') +
  labs(y = "July Precip (mm)", 
       subtitle = "WS closest plot 39: ws = grey; prism = blue; daymet = orange")

pcp_36 <- 
ggplot(ws_comp_36, aes(x = year)) + forestNETN::theme_FHM() +
  geom_point(aes(y = ws_pcpmm_7), color = 'dimgrey', alpha = 0.5) +
  geom_smooth(aes(y = ws_pcpmm_7), se = F, color = 'dimgrey') + 
  
  geom_point(aes(y = dm_ppt_sum_7), color = 'orange', alpha = 0.5) +
  geom_smooth(aes(y = dm_ppt_sum_7, x = year), se = F, color = 'orange') + 
  
  geom_point(aes(y = pr_ppt_07), color = 'dodgerblue', alpha = 0.5) +
  geom_smooth(aes(y = pr_ppt_07, x = year), se = F, color = 'dodgerblue') +
  labs(y = "July Precip (mm)", 
       subtitle = "WS farthest plot 36: ws = grey; prism = blue; daymet = orange")

pcp_64 <-
ggplot(ws_comp_64, aes(x = year)) + forestNETN::theme_FHM() +
  geom_point(aes(y = ws_pcpmm_7), color = 'dimgrey', alpha = 0.5) +
  geom_smooth(aes(y = ws_pcpmm_7), se = F, color = 'dimgrey') + 
  
  geom_point(aes(y = dm_ppt_sum_7), color = 'orange', alpha = 0.5) +
  geom_smooth(aes(y = dm_ppt_sum_7, x = year), se = F, color = 'orange') + 
  
  geom_point(aes(y = pr_ppt_07), color = 'dodgerblue', alpha = 0.5) +
  geom_smooth(aes(y = pr_ppt_07, x = year), se = F, color = 'dodgerblue') +
  labs(y = "July Precip (mm)", 
       subtitle = "WS farthest plot 64: ws = grey; prism = blue; daymet = orange")

cowplot::plot_grid(pcp_39, pcp_36, pcp_64, ncol = 1)

ppt_39 <- 
  ggplot(ws_comp_39, aes(x = ws_pcpmm_7, y = pr_ppt_07)) + forestNETN::theme_FHM() +
  geom_point(color = 'dodgerblue')+
  geom_point(aes(y = dm_ppt_sum_7), color = 'orange') +
  geom_abline(slope = 1, intercept = 0) + 
  labs(y = "Gridded Precip July", x = "Weather Station Precip July",
       subtitle = "Plot closest to WS: Prism = blue; Daymet = orange")

ppt_36 <- 
  ggplot(ws_comp_36, aes(x = ws_pcpmm_7, y = pr_ppt_07)) + forestNETN::theme_FHM() +
  geom_point(color = 'dodgerblue')+
  geom_point(aes(y = dm_ppt_sum_7), color = 'orange') +
  geom_abline(slope = 1, intercept = 0) + 
  labs(y = "Gridded Precip July", x = "Weather Station Precip July", 
       subtitle = "Plot1 farthest from WS: Prism = blue; Daymet = orange")

ppt_64 <- 
  ggplot(ws_comp_64, aes(x = ws_pcpmm_7, y = pr_ppt_07)) + forestNETN::theme_FHM() +
  geom_point(color = 'dodgerblue')+
  geom_point(aes(y = dm_ppt_sum_7), color = 'orange') +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Gridded Precip July", x = "Weather Station Precip July", 
       subtitle = "Plot2 farthest from WS: Prism = blue; Daymet = orange")

cowplot::plot_grid(ppt_39, ppt_36, ppt_64, nrow = 3)

# Going with DayMet- less separation of values with distance from weather station than PRISM.