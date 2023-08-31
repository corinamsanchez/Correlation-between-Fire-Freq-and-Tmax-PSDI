
### Exploratory practice code, not to be used for other purposes ###


library(spatialEco)
library(raster)
library(sp)
library(terra)
library(sf)


#Import data --------------------------------------------------------
#max temp and fire rasters from cheyenne's analysis
load(file="~/Documents/spatial_temporal_model/Covariates28.rdata")

#remove all the large objects that we don't need
rm(annherb)
rm(bare)
rm(bigsage)
rm(burnlag1yr)
rm(burnlag2yr)
rm(burnlag4yr)
rm(burnlag8yr)
rm(distroad)
rm(distspr)
rm(disttowns)
rm(distwater)
rm(elev)
rm(herb)
rm(litter)
rm(month)
rm(month.sq)
rm(ndvi)
rm(pdsilag1yr)
rm(pdsilag2yr)
rm(pdsilag3yr)
rm(pdsilag4yr)
rm(pdsilag5yr)
rm(pdsilag8yr)
rm(ppt)
rm(ppt.previous)
rm(ppt.previous.full)
rm(ppt1992)
rm(pptlag1yr)
rm(pptlag1yr1992)
rm(pptlag2yr)
rm(pptlag2yr1992)
rm(pptlag3yr)
rm(pptlag3yr1992)
rm(pptlag4yr)
rm(pptlag4yr1992)
rm(pptlag5yr)
rm(pptlag5yr1992)
rm(pptlag8yr)
rm(pptlag8yr1992)
rm(sage)
rm(sageheight)
rm(shrub)
rm(slope)
rm(tmean)
rm(tmin)


#burn, tmax, and pdsi are rasterbricks with 228 layers (months)
#same resolution, extent, etc


# Test -------------------------------------------------------------
# test one month

#take month 100 from burn rasterbrick and from tmax rasterbrick, then combine the 2 layers

m100_bt <- c(burn[[100]], tmax[[100]])
m100_bt <- stack(m100_bt)

#got this code from stackexchange 
jnk=layerStats(m100_bt, 'pearson', na.rm=T)
corr_matrix=jnk$'pearson correlation coefficient'
c100 <- corr_matrix[2,1] #correlation coefficient for burn vs tmax

#Expand to all months and make output a vector of correlation coefficients by month

# Burn vs Tmax ---------------------------------------------------------

n_layers <- nlayers(burn)
burn_tmax_list <- list()
for(i in 1:n_layers){
  burn_tmax_list[[i]] <- c(burn[[i]], tmax[[i]]) 
}

burn_tmax_list2 <- list()
for(j in 1:n_layers){
  burn_tmax_list2[[j]] <- stack(burn_tmax_list[[j]]) 
}

rm(burn_tmax_list) #takes up too much memory to keep

burn_tmax_list3 <- list()
for(k in 1:n_layers){
  burn_tmax_list3[[k]] <- layerStats(burn_tmax_list2[[k]], "pearson", na.rm = T) 
}

rm(burn_tmax_list2) #takes up too much memory to keep

burn_tmax_list4 <- list() #going to use for list of matrices 
for(m in 1:n_layers){
  burn_tmax_list4[[m]] <- burn_tmax_list3[[m]]$`pearson correlation coefficient` 
}

rm(burn_tmax_list3) #takes up too much memory to keep

burn_tmax_correlation <- numeric() #going to use for vector of correlation coefficients
for(n in 1:n_layers){
  burn_tmax_correlation[[n]] <- burn_tmax_list4[[n]][2,1] 
}


summary(burn_tmax_correlation)
# 89 NA's represent 89 months with no fires...I think


#save output min, max, mean

# Burn vs PDSI --------------------------------------------------------

n_layers <- nlayers(burn)
burn_pdsi_list <- list()
for(i in 1:n_layers){
  burn_pdsi_list[[i]] <- c(burn[[i]], pdsi[[i]]) 
}

burn_pdsi_list2 <- list()
for(j in 1:n_layers){
  burn_pdsi_list2[[j]] <- stack(burn_pdsi_list[[j]]) 
}

rm(burn_pdsi_list) #takes up too much memory to keep

burn_pdsi_list3 <- list()
for(k in 1:n_layers){
  burn_pdsi_list3[[k]] <- layerStats(burn_pdsi_list2[[k]], "pearson", na.rm = T) 
}

rm(burn_pdsi_list2) #takes up too much memory to keep

burn_pdsi_list4 <- list() #going to use for list of matrices 
for(m in 1:n_layers){
  burn_pdsi_list4[[m]] <- burn_pdsi_list3[[m]]$`pearson correlation coefficient` 
}

rm(burn_pdsi_list3) #takes up too much memory to keep

burn_pdsi_correlation <- numeric() #going to use for vector of correlation coefficients
for(n in 1:n_layers){
  burn_pdsi_correlation[[n]] <- burn_pdsi_list4[[n]][2,1] 
}


summary(burn_pdsi_correlation)


# Scrap code ------------------------------------------------------
# This code used the "rasterCorrelation" function but failed to give
# appropriate outputs (i.e., coefficients -1 < or > 1)

#convert to terra SpatRasters before using rasterCorrelation()
burn_sr <- rast(burn)
tmax_sr <- rast(tmax)
pdsi_sr <- rast(pdsi)

#test
test1 <- rasterCorrelation(burn_sr[[100]], tmax_sr[[100]], type = "pearson")


### Fire freq ~ Tmax correlation

#forloop of correlation tests through each month (layer) of burn raster data and
#max temp raster data 
#the output raster values are correlation for each layer
burn_temp_cor <- rast() #create an empty raster for upcoming forloop
n_layers <- nlyr(burn_sr)
nlyr(burn_temp_cor) <- n_layers
nrow(burn_temp_cor) <- nrow(burn_sr[[1]])
ncol(burn_temp_cor) <- ncol(burn_sr[[1]])
ext(burn_temp_cor) <- ext(burn_sr[[1]])
crs(burn_temp_cor) <- crs(burn_sr)
for(k in 1:n_layers){
  burn_temp_cor[[k]] <- rasterCorrelation(burn_sr[[k]], tmax_sr[[k]], type = "pearson") 
}

temp_minmax <- minmax(burn_temp_cor)

### Fire freq ~ PDSI correlation

#forloop of correlation tests through each month (layer) of burn raster data and
#pdsi raster data 
#the output raster values are correlation for each layer
burn_pdsi_cor <- rast() #create an empty raster for upcoming forloop
n_layers <- nlyr(burn_sr)
nlyr(burn_pdsi_cor) <- n_layers
nrow(burn_pdsi_cor) <- nrow(burn_sr[[1]])
ncol(burn_pdsi_cor) <- ncol(burn_sr[[1]])
ext(burn_pdsi_cor) <- ext(burn_sr[[1]])
crs(burn_pdsi_cor) <- crs(burn_sr)
for(k in 1:n_layers){
  burn_pdsi_cor[[k]] <- rasterCorrelation(burn_sr[[k]], pdsi_sr[[k]], type = "pearson") 
}

min(burn_pdsi_cor[[100:105]], na.rm=TRUE)
max(burn_pdsi_cor, na.rm=TRUE)
mean(burn_pdsi_cor, na.rm=TRUE)



### *just to test* Tmax ~ PDSI correlation

#forloop of correlation tests through each month (layer) of burn raster data and
#pdsi raster data 
#the output raster values are correlation for each layer
tmax_pdsi_cor <- rast() #create an empty raster for upcoming forloop
n_layers <- nlyr(tmax_sr)
nlyr(tmax_pdsi_cor) <- n_layers
nrow(tmax_pdsi_cor) <- nrow(tmax_sr[[1]])
ncol(tmax_pdsi_cor) <- ncol(tmax_sr[[1]])
ext(tmax_pdsi_cor) <- ext(tmax_sr[[1]])
crs(tmax_pdsi_cor) <- crs(tmax_sr)
for(k in 1:n_layers){
  tmax_pdsi_cor[[k]] <- rasterCorrelation(tmax_sr[[k]], pdsi_sr[[k]], s=1, type = "pearson") 
}

### this code is giving me problems because some of the output correlation coefficients
# are > 1 which means there is something wrong with the test. I thought this was 
# maybe due to the the burn variable being binary data but that doesn't seem to be the problem
# I think the real problem is that this is a moving window analysis and idk 
# exactly how to use that appropriately 

