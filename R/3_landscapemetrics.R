# Get the stable version from CRAN

#install.packages("landscapemetrics")

library(landscapemetrics)
library(landscapetools)
library(tidyverse)
library(dplyr)
library(raster)
library(rgdal)
library(rgeos)


# Preparing species names for the loop

sp.names <- read.csv("./dispersal_ability.csv", head=T)
sp.names <- sp.names$spp
n <-length(sp.names)

# See the list of metrics from the package landscapemetrics

#aaa <- list_lsm()


# Loop

for(i in 1:n) {

  # Reading raster per species
  
  cu <-
    raster(paste0("./outputs_bin_mpc20/", sp.names[i], "_present_bin_wgs.tif"))
  fu <-
    raster(paste0("./outputs_bin_mpc20/", sp.names[i], "_future_bin_wgs.tif"))
  
  # Checking my rasters
  
  show_landscape(cu)
  check_landscape(cu)
  
  # Running the metrics
  
  landscape_cu <-
    calculate_lsm(cu, level = "landscape", full_name = TRUE)
  landscape_fu <-
    calculate_lsm(fu, level = "landscape", full_name = TRUE)
  
  # Writing results
  
  write.csv(
    landscape_cu,
    paste0(
      "./Results_2_landscapemetrics/",
      sp.names[i],
      "_landscape_cu_wgs84.csv"),
      row.names = F
    
  )
  write.csv(
    landscape_fu,
    paste0(
      "./Results_2_landscapemetrics/",
      sp.names[i],
      "_landscape_fu_wgs84.csv"),
      row.names = F
    
  )
  
}


# _present_bin_albers.tif
# cu[cu == 1] <- 1
# cu[cu <= 1] <- 0
# Error: The area of the circumscribing circle is currently only implemented for equal resolutions.
# > check_landscape(cu)
# layer       crs units   class n_classes OK
# 1     1 projected     m integer         1  ✓
# 
# wgs
# layer        crs   units   class n_classes OK
# 1     1 geographic degrees integer         1  x


for(i in 1:n) {
  
  # Reading raster per species
  
  cu2 <-
    raster(paste0("./outputs_bin_mpc20/", sp.names[i], "_present_bin_albers.tif"))
  fu2 <-
    raster(paste0("./outputs_bin_mpc20/", sp.names[i], "_future_bin_albers.tif"))
  # Checking my rasters
  
  cu2[cu2 > 0] <- 1
  cu2[cu2 <= 0] <- 0
  show_landscape(cu2)
  check_landscape(cu2)

  
  #Running the metrics
  
  
  landscape_cu2 <-
    calculate_lsm(cu2, level = "landscape", full_name = TRUE)
  landscape_fu2 <-
    calculate_lsm(fu2, level = "landscape", full_name = TRUE)
  
  # Writing results
  
  write.csv(
    landscape_cu2,
    paste0(
      "./Results_2_landscapemetrics/",
      sp.names[i],
      "_landscape_cu_albers.csv"),
    row.names = F
    
  )
  write.csv(
    landscape_fu2,
    paste0(
      "./Results_2_landscapemetrics/",
      sp.names[i],
      "_landscape_fu_albers.csv"),
    row.names = F
    
  )
  
}





#Show correlation

source("./Scripts/show_correlation.R")

show_correlation(
  landscape_cu,
  method = "pearson",
  diag = TRUE,
  labels = FALSE,
  vjust = 0,
  text_size = 15
)


correlations <- calculate_correlation(landscape_cu)
