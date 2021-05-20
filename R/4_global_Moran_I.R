# https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/#:~:text=Moran's%20I%20is%20a%20measure,and%20then%20library(ape).
# https://stackoverflow.com/questions/15711013/spatial-autocorrelation-analysis-in-r


library(raster)
library(ape)
library(letsR)


# Preparing species names for the loop

sp.names <- read.csv("./dispersal_ability.csv", head=T)
sp.names <- sp.names$spp
n <-length(sp.names)


# Current -----------------------------------------------------------------


Moran_table <- matrix(nrow = n, ncol = 9)

colnames(Moran_table) <-
  c("Species",
    "moran_I_observed",
    "moran_I_expected",
    "moran_I_sd",
    "moran_I_p.value",
    "moran_I_observed_fu",
    "moran_I_expected_fu",
    "moran_I_sd_fu",
    "moran_I_p.value_fu")

for (i in 1:n) {
  cu_cont <-
    raster(paste0("./cont_mpc_wgs/", sp.names[i], "_present_cont_wgs.tif"))
  number_of_cells <-  1:ncell(cu_cont)
  values_of_cells <- values(cu_cont)
  coorden_of_cells <- xyFromCell(cu_cont, number_of_cells) #lonlat
  
  matrix_species_with_NA <-
    cbind(number_of_cells, values_of_cells, coorden_of_cells)
  matrix_species <-
    matrix_species_with_NA[!is.na(matrix_species_with_NA[, 2]),]
  
  matrix_pf_dist <- lets.distmat(matrix_species[, 3:4])
  
  dists.inv <- 1 / as.matrix(matrix_pf_dist)
  diag(dists.inv) <- 0
  
  dists.inv[1:5, 1:5]
  
  moran_I_cu <- Moran.I(matrix_species[, 2], dists.inv)
  moran_I_observed <- moran_I_cu$observed
  moran_I_expected <- moran_I_cu$expected
  moran_I_sd <- moran_I_cu$sd
  moran_I_p.value <- moran_I_cu$p.value
  
  tiff(
    paste0("correlogram_", sp.names[i], "_cu.tiff"),
    units = "in",
    width = 5,
    height = 4,
    res = 100,
    compression = 'lzw'
  )
  
  moran_correlogram <-
    lets.correl(matrix_species[, 2], matrix_pf_dist, 5, T, T)
  
  dev.off()
  
  write.csv(
    moran_correlogram,
    paste0("./Moran_table/", sp.names[i],
           "moran_correlogram_cu.csv")
  )
  

# Future ------------------------------------------------------------------


  fu_cont <-
    raster(paste0("./cont_mpc_wgs/", sp.names[i], "_future_cont_wgs.tif"))
  
  number_of_cells2 <-  1:ncell(fu_cont)
  values_of_cells2 <- values(fu_cont)
  coorden_of_cells2 <- xyFromCell(fu_cont, number_of_cells2) #lonlat
  
  matrix_species_with_NA2 <-
    cbind(number_of_cells2, values_of_cells2, coorden_of_cells2)
  matrix_species2 <-
    matrix_species_with_NA2[!is.na(matrix_species_with_NA2[, 2]),]
  
  matrix_pf_dist2 <- lets.distmat(matrix_species2[, 3:4])
  
  dists.inv2 <- 1 / as.matrix(matrix_pf_dist2)
  diag(dists.inv2) <- 0
  
  dists.inv2[1:5, 1:5]
  
  moran_I_fu <- Moran.I(matrix_species2[, 2], dists.inv2)
  moran_I_observed_fu <- moran_I_fu$observed
  moran_I_expected_fu <- moran_I_fu$expected
  moran_I_sd_fu <- moran_I_fu$sd
  moran_I_p.value_fu <- moran_I_fu$p.value
  
  tiff(
    paste0("correlogram_", sp.names[i], "_fu.tiff"),
    units = "in",
    width = 5,
    height = 4,
    res = 100,
    compression = 'lzw'
  )
  
  moran_correlogram <-
    lets.correl(matrix_species2[, 2], matrix_pf_dist2, 5, T, T)
  
  dev.off()
  
  write.csv(
    moran_correlogram,
    paste0("./Moran_table/", sp.names[i],
           "moran_correlogram_fu.csv")
  )
  
  Moran_table_[i,] <-
    c(sp.names[i],
      moran_I_observed,
      moran_I_expected,
      moran_I_sd,
      moran_I_p.value,
      moran_I_observed_fu,
      moran_I_expected_fu,
      moran_I_sd_fu,
      moran_I_p.value_fu)
  
  
}

write.csv(Moran_table, "./Moran_I.csv", row.names = F)