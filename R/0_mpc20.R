library(raster)
library(rgdal)
library(rgeos)

sp.d <- read.csv("./dispersal_ability.csv")
sp.d$number <- 1:71
sp.names <- read.csv("./dispersal_ability.csv")
sp.names <- sp.names$spp

file.names <- list.files(paste0("./polygons/"), 
                         pattern="*.shp", full.names=T, recursive=FALSE)
allShapes = lapply(file.names, readOGR)


for(sp.n in sp.names[4:71]){
  
  cu <- raster(paste0("./outputs_bin/", sp.n, "_present_ensemble_0.5_consensus.tif"))
  fu <- raster(paste0("./outputs_bin/", sp.n, "_future_ensemble_0.5_consensus.tif"))
  cu_cont <- raster(paste0("./outputs_cont/", sp.n, "_present_TSSmax_ensemble_weighted_average.tif"))
  fu_cont1 <- raster(paste0("./outputs_cont/", sp.n, "_ac85bi50_TSSmax_ensemble_weighted_average.tif"))
  fu_cont2 <- raster(paste0("./outputs_cont/", sp.n, "_he85bi50_TSSmax_ensemble_weighted_average.tif"))
  fu_cont3 <- raster(paste0("./outputs_cont/", sp.n, "_mp85bi50_TSSmax_ensemble_weighted_average.tif"))
  fu_cont <- (fu_cont1 + fu_cont2 + fu_cont3)/3
  
  crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  proj4string(fu) <- crs.wgs84
  proj4string(fu_cont) <- crs.wgs84
  
  writeRaster(fu_cont, filename = paste0("./outputs_cont/", sp.n,
                                         "_future_TSSmax_ensemble_weighted_average.tif"),
              format="GTiff", overwrite=T)

  
  for (m in sp.d[grep(sp.n, sp.d$spp),4]){
    pol <- allShapes[[m]]
  }
  
  
  #crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
  crs.albers <- ("+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +datum=WGS84 +units=m +no_defs")
  proj4string(pol) <- crs.wgs84
  pol_20 <- spTransform(pol, crs.albers)
  b <- (gArea(pol_20)*2e-07)
  pol_20 <- gBuffer(pol_20, width = b)
  pol_20 <- spTransform(pol_20, crs.wgs84)
  
  
  # CURRENT
  
  cu_pol_20 <- crop(cu, pol_20)
  cu_pol_20_crop_mask <- mask(cu_pol_20, pol_20)
  
  cu_pol_20_crop_mask_albers <- cu_pol_20_crop_mask
  proj4string(cu_pol_20_crop_mask_albers) <- crs.wgs84
  cu_pol_20_crop_mask_albers <- projectRaster(cu_pol_20_crop_mask_albers, crs=crs.albers)
  
  writeRaster(cu_pol_20_crop_mask, filename = paste0("./outputs_bin_mpc20/", sp.n,
                                             "_present_bin_wgs.tif"), format="GTiff",
              overwrite=T)
  
  writeRaster(cu_pol_20_crop_mask_albers, filename = paste0("./outputs_bin_mpc20/", sp.n,
                                             "_present_bin_albers.tif"), format="GTiff",
              overwrite=T)
  
  
  cu_cont_pol_20 <- crop(cu_cont, pol_20)
  cu_cont_pol_20_crop_mask <- mask(cu_cont_pol_20, pol_20)
  
  cu_cont_pol_20_crop_mask_albers <- cu_cont_pol_20_crop_mask
  proj4string(cu_cont_pol_20_crop_mask_albers) <- crs.wgs84
  cu_cont_pol_20_crop_mask_albers <- projectRaster(cu_cont_pol_20_crop_mask_albers, crs=crs.albers)
  
  writeRaster(cu_cont_pol_20_crop_mask, filename = paste0("./outputs_cont_mpc20/", sp.n,
                                                     "_present_cont_wgs.tif"), format="GTiff",
              overwrite=T)
  
  writeRaster(cu_cont_pol_20_crop_mask_albers, filename = paste0("./outputs_cont_mpc20/", sp.n,
                                                            "_present_cont_albers.tif"), format="GTiff",
              overwrite=T)
  
  # FUTURE
  
  fu_pol_20 <- crop(fu, pol_20)
  fu_pol_20_crop_mask <- mask(fu_pol_20, pol_20)
  
  fu_pol_20_crop_mask_albers <- fu_pol_20_crop_mask
  proj4string(fu_pol_20_crop_mask_albers) <- crs.wgs84
  fu_pol_20_crop_mask_albers <- projectRaster(fu_pol_20_crop_mask_albers, crs=crs.albers)
  
  writeRaster(fu_pol_20_crop_mask, filename = paste0("./outputs_bin_mpc20/", sp.n,
                                                     "_future_bin_wgs.tif"), format="GTiff",
              overwrite=T)
  
  writeRaster(fu_pol_20_crop_mask_albers, filename = paste0("./outputs_bin_mpc20/", sp.n,
                                                            "_future_bin_albers.tif"), format="GTiff",
              overwrite=T)
  
  
  fu_cont_pol_20 <- crop(fu_cont, pol_20)
  fu_cont_pol_20_crop_mask <- mask(fu_cont_pol_20, pol_20)
  
  fu_cont_pol_20_crop_mask_albers <- fu_cont_pol_20_crop_mask
  proj4string(fu_cont_pol_20_crop_mask_albers) <- crs.wgs84
  fu_cont_pol_20_crop_mask_albers <- projectRaster(fu_cont_pol_20_crop_mask_albers, crs=crs.albers)
  
  writeRaster(fu_cont_pol_20_crop_mask, filename = paste0("./outputs_cont_mpc20/", sp.n,
                                                          "_future_cont_wgs.tif"), format="GTiff",
              overwrite=T)
  
  writeRaster(fu_cont_pol_20_crop_mask_albers, filename = paste0("./outputs_cont_mpc20/", sp.n,
                                                                 "_future_cont_albers.tif"), format="GTiff",
              overwrite=T)
  
}
