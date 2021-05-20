library(tidyverse)
library(dplyr)
library(raster)
library(rgdal)
library(rgeos)

setwd("~/Metricas_paisagem_no_clima")

sp.names <- read.csv("./dispersal_ability.csv", head=T)
sp.names <- sp.names$spp
n <-length(sp.names)


#sp.names <- as.character(unique(sp.names$spp))



# Building table to HM ----------------------------------------------------

Results_diff_area <- matrix(nrow = n, ncol = 8)
colnames(Results_diff_area) <- c("Species", "unchanged_porcent_bin", "loss_porcent_bin", 
                          "gain_porcent_bin", "mean_diff_fut_cur_porcent",
                          "ocorrencia_area_calc_pres", "ocorrencia_area_calc2_fut",
                          "ocorrencia_area_calc3_nogain")


for(i in 1:n){

cu <- raster(paste0("./outputs_bin_mpc20/", sp.names[i], "_present_bin_wgs.tif"))
fu <- raster(paste0("./outputs_bin_mpc20/", sp.names[i], "_future_bin_wgs.tif"))
cu_cont <- raster(paste0("./outputs_cont_mpc20/", sp.names[i], "_present_cont_wgs.tif"))
fu_cont <- raster(paste0("./outputs_cont_mpc20/", sp.names[i], "_future_cont_wgs.tif"))


#BINARIO

cu2 = cu
cu2[cu2 == 1] <- 2
# fut-pres
# 1-2 = -1 -> Aqui seria 1-1 = 0, i.e., nos dois n mudoram, estavel
# 1-0 = 1 -> ganho
# 0-2 = -2 ->  Aqui seria 0-1 = -1, i.e., perda
# 0-0 = 0 -> area inadequada em ambos cenarios, logo, proximo passo transformamos logo em NA

diff_bin <- fu - cu2
diff_bin[diff_bin == 0] <- NA
#unsuitable_both <- ncell(which(diff_bin[] == NA))
unchanged <- ncell(which(diff_bin[] == -1))
loss <- ncell(which(diff_bin[] == -2))
gain <- ncell(which(diff_bin[] == 1))
total_cells <- unchanged + loss + gain

unchanged_porcent_bin <- (unchanged/total_cells)*100
loss_porcent_bin <- (loss/total_cells)*100
gain_porcent_bin <- (gain/total_cells)*100
#extent_ecoph <- unchanged + gain


#reescalonando de zero a 1 para salvar o raster

diff_bin[diff_bin == -1] <- 0
diff_bin[diff_bin == -2] <- -1

writeRaster(diff_bin, filename = paste0("./Diff_raster/", sp.names[i],
                                       "_diff_bin.tif"),
            format="GTiff", overwrite=T)


#CONTINUO

testt <- t.test(fu_cont[], cu_cont[], paired=T)
#Mean of difference - Negative value means loss, positive value means gain
chars <- capture.output(print(testt))
mean_diff_estimate <- as.data.frame(testt$estimate)
mean_diff_fut_cur <- mean_diff_estimate[1,]
mean_diff_fut_cur_porcent <- mean_diff_fut_cur*100
mean_diff_conf.int <- as.data.frame(testt$conf.int)
diff_cont <- fu_cont - cu_cont
writeRaster(diff_cont, filename = paste0("./Diff_raster/", sp.names[i],
                                         "_diff_cont.tif"),
            format="GTiff", overwrite=T)



# AREA --------------------------------------------------------------------

#current
cu_area = cu
# B.torq_RMR_cur_bin #262, 151, 39562  (nrow, ncol, ncell)
cu_area[cu_area == 0] <- NA
# B.torq_RMR_fut_bin #262, 151, 39562  (nrow, ncol, ncell)
r = raster(nrow = cu_area@nrows, ncol = cu_area@ncols, xmn = cu_area@extent@xmin, 
           xmx = cu_area@extent@xmax, ymn = cu_area@extent@ymin, ymx = cu_area@extent@ymax) 
# calcula area de cada celula (a area muda ao longo das latitudes/longitudes)
x = raster::area(r) 
#plot(x)
ocorrencia_area = x * cu_area
ocorrencia_area_calc <- cellStats(ocorrencia_area, stat='sum', na.rm=TRUE, asSample=TRUE)
#ocorrencia_area_calc #EM KM2

#future
fu_area = fu 
fu_area[fu_area == 0] <- NA
ocorrencia_area2 = x * fu_area
ocorrencia_area_calc2 <- cellStats(ocorrencia_area2, stat='sum', na.rm=TRUE, asSample=TRUE)


# Non-dispersal scenario
# mudanca = perda

diff_bin_perda <- diff_bin
diff_bin_perda[diff_bin_perda == 1] <- NA
diff_bin_perda[diff_bin_perda == -1] <- NA
diff_bin_perda[diff_bin_perda == 0] <- 1
diff_bin_perda[diff_bin_perda == -2] <- NA
ocorrencia_area = x * diff_bin_perda
ocorrencia_area_calc3 <- cellStats(ocorrencia_area, stat='sum', na.rm=TRUE, asSample=TRUE)
#ocorrencia_area_calc3 #EM KM2


Results_diff_area[i, ] <- c(sp.names[i], unchanged_porcent_bin, loss_porcent_bin,
                     gain_porcent_bin, mean_diff_fut_cur_porcent,
                     ocorrencia_area_calc, ocorrencia_area_calc2,
                     ocorrencia_area_calc3)

}

write.csv(Results_diff_area, "./1_calculo_de_area.csv", row.names = F)






# testar ------------------------------------------------------------------

# . -----------------------------------------------------------------------
# Mapbiomas ---------------------------------------------------------------
# . -----------------------------------------------------------------------

# Essa parte abaixo ja fiz, ai ja chamo o raster pronto depois

# Organizando e resample os dados do Mapbiomas

# Mapbiomas_2000_raw <- raster("~/land_use/B.torq_Mapbiomas_2000.tif")
# Mapbiomas_2000 <- raster("~/land_use/forest2000.tif")
# 
# Mapbiomas_2000_crop <- crop(Mapbiomas_2000, B.torq_RMR_cur_bin)
# Mapbiomas_2000_resample <- resample(Mapbiomas_2000_crop, B.torq_RMR_cur_bin, method="bilinear")
#B.torq_RMR_cur_bin_resample <- resample(B.torq_RMR_cur_bin, Mapbiomas_2000_crop, method="bilinear")
# writeRaster(Mapbiomas_2000_resample, "~/land_use/B.torq_Mapbiomas_2000_resample.tif"
#             ,format="GTiff",overwrite=T)

# Mapbiomas_2000_resample2 = Mapbiomas_2000_resample
# Mapbiomas_2000_resample2[Mapbiomas_2000_resample2 >= 0.1] <- 1
# Mapbiomas_2000_resample2[Mapbiomas_2000_resample2 < 0.1] <- 0
# plot(Mapbiomas_2000_resample2)
# 
# Mapbiomas_2000_resample3 <-mask(Mapbiomas_2000_resample2,B.torq_RMR_cur_bin)
# 
# writeRaster(Mapbiomas_2000_resample2, "~/land_use/B.torq_Mapbiomas_2000_resample2.tif"
#             ,format="GTiff",overwrite=T)

Mapbiomas_2000_resample3 <- raster("~/land_use/B.torq_Mapbiomas_2000_resample3.tif")

B.torq_RMR_cur_bin_Mapbiomas <- Mapbiomas_2000_resample3 * B.torq_RMR_cur_bin
B.torq_RMR_fut_bin_Mapbiomas <- Mapbiomas_2000_resample3 * B.torq_RMR_fut_bin
B.torq_RMR_cur_cont_Mapbiomas <- Mapbiomas_2000_resample3 * B.torq_RMR_cur_cont
B.torq_RMR_fut_cont_Mapbiomas <- Mapbiomas_2000_resample3 * B.torq_RMR_fut_cont

writeRaster(B.torq_RMR_cur_bin_Mapbiomas, "~/Rasters_result_cap3/B.torq_RMR_cur_bin_Mapbiomas.tif"
            ,format="GTiff",overwrite=T)
writeRaster(B.torq_RMR_fut_bin_Mapbiomas, "~/Rasters_result_cap3/B.torq_RMR_fut_bin_Mapbiomas.tif"
            ,format="GTiff",overwrite=T)
writeRaster(B.torq_RMR_cur_cont_Mapbiomas, "~/Rasters_result_cap3/B.torq_RMR_cur_cont_Mapbiomas.tif"
            ,format="GTiff",overwrite=T)
writeRaster(B.torq_RMR_fut_cont_Mapbiomas, "~/Rasters_result_cap3/B.torq_RMR_fut_cont_Mapbiomas.tif"
            ,format="GTiff",overwrite=T)

#BINARIO

#Calcullando perda e ganho

B.torq_RMR_cur_bin_Mapbiomas_2 = B.torq_RMR_cur_bin_Mapbiomas
B.torq_RMR_cur_bin_Mapbiomas_2[B.torq_RMR_cur_bin_Mapbiomas_2 == 1] <- 2

diff_bin_Mapbiomas <- B.torq_RMR_fut_bin_Mapbiomas - B.torq_RMR_cur_bin_Mapbiomas_2
diff_bin_Mapbiomas[diff_bin_Mapbiomas == 0] <- NA
#unsuitable_both <- ncell(which(diff_bin_Mapbiomas[] == NA))
unchanged <- ncell(which(diff_bin_Mapbiomas[] == -1))
loss <- ncell(which(diff_bin_Mapbiomas[] == -2))
gain <- ncell(which(diff_bin_Mapbiomas[] == 1))
total_cells <- unchanged + loss + gain

unchanged_porcent_bin <- (unchanged/total_cells)*100
loss_porcent_bin <- (loss/total_cells)*100
gain_porcent_bin <- (gain/total_cells)*100
#extent_ecoph <- unchanged + gain


#reescalonando de zero a 1 para salvar o raster
diff_bin_Mapbiomas[diff_bin_Mapbiomas == -1] <- 0
diff_bin_Mapbiomas[diff_bin_Mapbiomas == -2] <- -1
writeRaster(diff_bin_Mapbiomas, "~/Rasters_result_cap3/B.torq_Mapbiomas_bin_diff.tif"
            ,format="GTiff",overwrite=T)

#CONTINUO

# B.torq_RMR_cur_cont_Mapbiomas <- Mapbiomas_2000_resample3 * B.torq_RMR_cur_cont
# B.torq_RMR_fut_cont_Mapbiomas <- Mapbiomas_2000_resample3 * B.torq_RMR_fut_cont

testt <- t.test(B.torq_RMR_fut_cont_Mapbiomas[], B.torq_RMR_cur_cont_Mapbiomas[], paired=T)
chars <- capture.output(print(testt))
mean_diff_estimate <- as.data.frame(testt$estimate)
mean_diff_fut_cur <- mean_diff_estimate[1,]
mean_diff_fut_cur_porcent <- mean_diff_fut_cur*100 
mean_diff_conf.int <- as.data.frame(testt$conf.int)
table_results$mean_diff_fut_cur_porcent_Mapbiomas <- mean_diff_fut_cur_porcent

B.torq_RMR_cont_diff_Mapbiomas <- B.torq_RMR_fut_cont_Mapbiomas - B.torq_RMR_cur_cont_Mapbiomas
writeRaster(B.torq_RMR_cont_diff_Mapbiomas, 
            "~/Rasters_result_cap3/B.torq_RMR_cont_diff_Mapbiomas.tif"
            ,format="GTiff",overwrite=T)


# Area --------------------------------------------------------------------
#current
# B.torq_RMR_cur_bin #262, 151, 39562  (nrow, ncol, ncell)
B.torq_RMR_cur_bin_Mapbiomasarea = B.torq_RMR_cur_bin_Mapbiomas
B.torq_RMR_cur_bin_Mapbiomasarea[B.torq_RMR_cur_bin_Mapbiomasarea == 0] <- NA
ocorrencia_area = x * B.torq_RMR_cur_bin_Mapbiomasarea
ocorrencia_area_calc <- cellStats(ocorrencia_area, stat='sum', na.rm=TRUE, asSample=TRUE)
ocorrencia_area_calc #EM KM2

#future
B.torq_RMR_fut_bin_Mapbiomasarea = B.torq_RMR_fut_bin_Mapbiomas 
B.torq_RMR_fut_bin_Mapbiomasarea[B.torq_RMR_fut_bin_Mapbiomasarea == 0] <- NA
ocorrencia_area2 = x * B.torq_RMR_fut_bin_Mapbiomasarea
ocorrencia_area_calc2 <- cellStats(ocorrencia_area2, stat='sum', na.rm=TRUE, asSample=TRUE)

table_results$ocorrencia_area_km2_cur_Mapbiomas <- ocorrencia_area_calc
table_results$ocorrencia_area_km2_fut_Mapbiomas <- ocorrencia_area_calc2

# Antes era a mudanca liquida = ganhos - perda, agora eh mudanca = perda
diff_bin_perda <- diff_bin_Mapbiomas
diff_bin_perda[diff_bin_perda == 1] <- NA
diff_bin_perda[diff_bin_perda == -1] <- NA
diff_bin_perda[diff_bin_perda == 0] <- 1
ocorrencia_area = x * diff_bin_perda
ocorrencia_area_calc <- cellStats(ocorrencia_area, stat='sum', na.rm=TRUE, asSample=TRUE)
ocorrencia_area_calc #EM KM2
table_results$area_less_loss_Mapbiomas <- ocorrencia_area_calc


#table_results$Species <- c("B. torquatus", "B. tridactylus", "B. variegatus")
#table_results <- table_results[,c(19,1:18)]
write.csv(table_results, "~/Rasters_result_cap3/table_results_cap3.csv", row.names = F)

save.image("~/Rasters_result_cap3/Btorq_diff_clima_LUH_Mapbiomas.rData")
load("~/Rasters_result_cap3/Btorq_diff_clima_LUH_Mapbiomas.rData")

