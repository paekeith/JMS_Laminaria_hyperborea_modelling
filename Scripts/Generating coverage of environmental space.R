################################
# Code for the estimation of coverage of the environmental space as part of Eskuche-Keith et al. "Abundance-based distribution modelling reveals climate-driven loss of functionally important kelp forest habitat".

# This follows the approach of https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13035
# steps involved:
#1. create a 1km buffer around each of the PA records. This will be used to remove cells from the predictors before sampling background values for the environmental coverage. The reason for this is because we suspect our coverage predictions are being influenced by the granularity of the rasters. This was done in QGIS
#2. convert buffers to raster
#3. crop the predictor raster by the buffer raster. 
#4. then extract the cell IDs for the remaining cells and cut extra ones from the predictor df 
#5. extract the environmental data from all the presences
#6. randomly sample the predictors

# Last edited 22/05/2026
################################


# Load packages and data --------------------------------
rm(list=ls())
gc()
library(terra)
library(raster)
library(dismo)
library(gbm)
library(cowplot)

#Loading and prepping the species data
data <- readRDS("Data/Refined_data/full_species_data_thinned.rds")

#extract relevant variables
data <- data[,c("cell","Origin","Lat","Lon","SACFORN","SACFORN_int","temp_mean","nitrate_mean","cur_vel_dep_mean","bathymetry","wave_fetch","PAR_depth")]

## Prepping the sample data ##
data$sample_site <- 1 #identifying where samples are present

#splitting the records by origin, as the background cells need to be selected for UK and Ireland seperately
data_UK <- data[data$Origin=="UK",]
data_IRE <- data[data$Origin=="IRE",]

#preparing the predictor data
records_buffer <- vect("Data/Refined_data/full_species_data_thinned_1km_buffer.shp")

#loading predictor layers
# NOTE: The predictor .rds files used for fitting the original SDM can't be used because the IDs no longer match those of the species data (due to the tiling process). Instead, the cell IDs need to be re-linked using the associated xy coordinates. This will need to be done separately for the Irish and UK data, as they were seperate rasters and might therefore have duplicated cell values

predictors_UK_pres_rast <- rast("Data/Refined_data/Predictors_UK_present.tif")
predictors_UK_pres_df <- readRDS("Data/Refined_data/Predictors_UK_present.rds")
predictors_UK_245_rast <- rast("Data/Refined_data/Predictors_UK_ssp245.tif")
predictors_UK_245_df <- readRDS("Data/Refined_data/Predictors_UK_ssp245.rds")
predictors_UK_585_rast <- rast("Data/Refined_data/Predictors_UK_ssp585.tif")
predictors_UK_585_df <- readRDS("Data/Refined_data/Predictors_UK_ssp585.rds")
predictors_IRE_pres_rast <- rast("Data/Refined_data/Predictors_IRE_present.tif")
predictors_IRE_pres_df <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")
predictors_IRE_245_rast <- rast("Data/Refined_data/Predictors_IRE_ssp245.tif")
predictors_IRE_245_df <- readRDS("Data/Refined_data/Predictors_IRE_ssp245.rds")
predictors_IRE_585_rast <- rast("Data/Refined_data/Predictors_IRE_ssp585.tif")
predictors_IRE_585_df <- readRDS("Data/Refined_data/Predictors_IRE_ssp585.rds")

#extracting cell IDs
predictors_UK_pres_df$cell <- cellFromXY(predictors_UK_pres_rast, predictors_UK_pres_df[, c("x", "y")])
predictors_UK_245_df$cell <- cellFromXY(predictors_UK_245_rast, predictors_UK_245_df[, c("x", "y")])
predictors_UK_585_df$cell <- cellFromXY(predictors_UK_585_rast, predictors_UK_585_df[, c("x", "y")])

predictors_IRE_pres_df$cell <- cellFromXY(predictors_IRE_pres_rast, predictors_IRE_pres_df[, c("x", "y")])
predictors_IRE_245_df$cell <- cellFromXY(predictors_IRE_245_rast, predictors_IRE_245_df[, c("x", "y")])
predictors_IRE_585_df$cell <- cellFromXY(predictors_IRE_585_rast, predictors_IRE_585_df[, c("x", "y")])

#setting bathymetry values above 0 to 0
predictors_UK_pres_df$bathymetry[predictors_UK_pres_df$bathymetry<0] <- 0
predictors_UK_245_df$bathymetry[predictors_UK_245_df$bathymetry<0] <- 0
predictors_UK_585_df$bathymetry[predictors_UK_585_df$bathymetry<0] <- 0
predictors_IRE_pres_df$bathymetry[predictors_IRE_pres_df$bathymetry<0] <- 0
predictors_IRE_245_df$bathymetry[predictors_IRE_245_df$bathymetry<0] <- 0
predictors_IRE_585_df$bathymetry[predictors_IRE_585_df$bathymetry<0] <- 0

# Prepping the predictor layers --------------------
# Rasterising the buffer layer
values(records_buffer) <- 1

buffer_rast_UK <- rasterize(records_buffer, predictors_UK_pres_rast$bathymetry, field="value", fun="mean", background=NA)
buffer_rast_IRE <- rasterize(records_buffer, predictors_IRE_pres_rast$bathymetry, field="value", fun="mean", background=NA)

#3. cut the buffer layer from the predictor rasters for training the model with
preds_pres_UK_crop <- terra::mask(predictors_UK_pres_rast$wave_fetch, buffer_rast_UK,inverse=TRUE)
preds_pres_IRE_crop <- terra::mask(predictors_IRE_pres_rast$wave_fetch, buffer_rast_IRE,inverse=TRUE)

#4. extract cell IDs for the remaining cells after cutting by buffer layer
preds_pres_UK_crop_df <- as.data.frame(preds_pres_UK_crop,cells=TRUE) #produces a dataframe with the non-NA cells 
preds_pres_IRE_crop_df <- as.data.frame(preds_pres_IRE_crop,cells=TRUE) #produces a dataframe with the non-NA cells 

#5. removing rows from predictor_dfs which do not match the cell IDs in the predictor layer cropped by the buffer layer
predictors_UK_pres_df_reduced <- predictors_UK_pres_df[which(predictors_UK_pres_df$cell %in% preds_pres_UK_crop_df$cell),]

predictors_IRE_pres_df_reduced <- predictors_IRE_pres_df[which(predictors_IRE_pres_df$cell %in% preds_pres_IRE_crop_df$cell),]

imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))

#sampling the reduced predictor rasters for background points ---- FIRST CUT RASTER TO REMOVE LAT
train_ind_pres_UK <- sample(seq_len(nrow(predictors_UK_pres_df_reduced)), size = nrow(data_UK)*2)
preddat_pres_UK <- predictors_UK_pres_df_reduced[train_ind_pres_UK, ]
preddat_pres_UK$sample_site <- 0 #setting up the sample site absences
preddat_pres_UK <- preddat_pres_UK[,c("cell",imp.var,"sample_site")]

train_ind_pres_IRE <- sample(seq_len(nrow(predictors_IRE_pres_df_reduced)), size = nrow(data_IRE)*2)
preddat_pres_IRE <- predictors_IRE_pres_df_reduced[train_ind_pres_IRE, ]
preddat_pres_IRE$sample_site <- 0 #setting up the sample site absences
preddat_pres_IRE <- preddat_pres_IRE[,c("cell",imp.var,"sample_site")]

#combining the species records and the background samples into one dataset for modelling
data <- data[,c("cell",imp.var,"sample_site")]
preddat_pres_comb <- rbind(preddat_pres_UK,preddat_pres_IRE)
data_preddat_pres_comb <- rbind(preddat_pres_comb,data)
summary(data_preddat_pres_comb)

## fitting the model for environmental coverage  -----------------
model_pres <- dismo::gbm.step(data=data_preddat_pres_comb, gbm.x = imp.var, # environmental variable columns
                              gbm.y = "sample_site", # presence absence of samples
                              family = "bernoulli", tree.complexity = 2,
                              learning.rate = 0.1, bag.fraction = 0.6, n.folds=10, 
                              max.trees = 5000, plot.main = F, tolerance.method = "fixed",
                              tolerance = 0.01, verbose = T)

#saving the model to disk
saveRDS(model_pres,"Outputs/Env_cover_model.rds")

pred.map.UK <- predict.gbm(model_pres, predictors_UK_pres_df, 
                           n.trees = model_pres$gbm.call$best.trees, type = "response")
pred.map.IRE <- predict.gbm(model_pres, predictors_IRE_pres_df, 
                           n.trees = model_pres$gbm.call$best.trees, type = "response")

pred.map.245.UK <- predict.gbm(model_pres, predictors_UK_245_df, 
                           n.trees = model_pres$gbm.call$best.trees, type = "response")
pred.map.245.IRE <- predict.gbm(model_pres, predictors_IRE_245_df, 
                            n.trees = model_pres$gbm.call$best.trees, type = "response")

pred.map.585.UK <- predict.gbm(model_pres, predictors_UK_585_df, 
                           n.trees = model_pres$gbm.call$best.trees, type = "response")
pred.map.585.IRE <- predict.gbm(model_pres, predictors_IRE_585_df, 
                            n.trees = model_pres$gbm.call$best.trees, type = "response")

pred.map.UK.rast <- rasterFromXYZ(data.frame(x = predictors_UK_pres_df[,"x"],y = predictors_UK_pres_df[,"y"],z = pred.map.UK),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
writeRaster(pred.map.UK.rast,"Outputs/Present_UK_env_cover.tif",overwrite=TRUE)

pred.map.245.UK.rast <- rasterFromXYZ(data.frame(x = predictors_UK_245_df[,"x"],y = predictors_UK_245_df[,"y"],z = pred.map.245.UK),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
writeRaster(pred.map.245.UK.rast,"Outputs/ssp245_UK_env_cover.tif",overwrite=TRUE)

pred.map.585.UK.rast <- rasterFromXYZ(data.frame(x = predictors_UK_585_df[,"x"],y = predictors_UK_585_df[,"y"],z = pred.map.585.UK),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
writeRaster(pred.map.585.UK.rast,"Outputs/ssp585_UK_env_cover.tif",overwrite=TRUE)

pred.map.IRE.rast <- rasterFromXYZ(data.frame(x = predictors_IRE_pres_df[,"x"],y = predictors_IRE_pres_df[,"y"],z = pred.map.IRE),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
writeRaster(pred.map.IRE.rast,"Outputs/Present_IRE_env_cover.tif",overwrite=TRUE)

pred.map.245.IRE.rast <- rasterFromXYZ(data.frame(x = predictors_IRE_245_df[,"x"],y = predictors_IRE_245_df[,"y"],z = pred.map.245.IRE),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
writeRaster(pred.map.245.IRE.rast,"Outputs/ssp245_IRE_env_cover.tif",overwrite=TRUE)

pred.map.585.IRE.rast <- rasterFromXYZ(data.frame(x = predictors_IRE_585_df[,"x"],y = predictors_IRE_585_df[,"y"],z = pred.map.585.IRE),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
writeRaster(pred.map.585.IRE.rast,"Outputs/ssp585_IRE_env_cover.tif",overwrite=TRUE)

## Further steps (not provided here) include clipping the coverage layers to the lowest astronomical tide, and masking cells with soft sediment. The necessary LAT and substrate files are provided in the raw data folder

