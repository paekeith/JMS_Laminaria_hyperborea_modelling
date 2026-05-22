################################
# Code for the processing of input data (species records and environmental datasets) for use in the SDM
# as part of Eskuche-Keith et al. "Abundance-based distribution modelling reveals climate-driven loss of functionally important kelp forest habitat"
# Last edited 22/05/2026
################################

####========== 1. PREPARING SPECIES RECORDS  ================================####

## Loading packages and data ##
rm(list=ls())
library(terra)
library(tidyverse)
library(biooracler)

#loading the records obtained from the UK Archive for Marine Species and Habitats data (DASSH)
hyperborea_records_DASSH <- read.csv("Data/Raw_data/Kelp_records/L_hyperborea_DASSH.csv")
#loading records from the JNCC MNCR (supplied by Mike Burrows)
hyperborea_records_MNCR <- read.csv("Data/Raw_data/Kelp_records/L_hyperborea_MNCR.csv")

## processing the DASSH data ##
# removing all DASSH records that were part of the JNCC MNCR, to avoid possible duplicates with the MNCR data
hyperborea_records_DASSH <- hyperborea_records_DASSH[-grep("JNCCMNCR",hyperborea_records_DASSH$survey_code),]

#subsetting the DASSH data to exclude irrelevant columns
names(hyperborea_records_DASSH)
hyperborea_records_DASSH <- hyperborea_records_DASSH%>%
  dplyr::select(originator,latitude,longitude,date,abundance,abundance_type,depth)

#generating a presence/absence column based on the abundance column
hyperborea_records_DASSH$Presence <- NA
hyperborea_records_DASSH$Presence[which(hyperborea_records_DASSH$abundance=="N")] <- 0
hyperborea_records_DASSH$Presence[which(is.na(hyperborea_records_DASSH$Presence)==TRUE)] <- 1

#converting values with percent cover to SACFORN, following https://mhc.jncc.gov.uk/media/1009/sacfor.pdf:
# S = >80%
# A = 40-80%
# C = 20-40%
# F = 10 - 20%
# O = 5 - 10%
# R = 1 - 5%
# N = <1%

percent_cover <- hyperborea_records_DASSH%>%
  filter(abundance_type %in% c("percentage","Percentage","percent"))

percent_cover <- percent_cover %>%
  mutate(abundance = str_remove(abundance, "%"),
         abundance = as.numeric(abundance))

percent_cover$abundance[which(percent_cover$abundance==0)] <- "N" #converting all the definite numerical absences (from count and percent cover data) to binary absence
percent_cover$abundance[which(percent_cover$abundance >=80)] <- "S"
percent_cover$abundance[which(percent_cover$abundance >=40 & percent_cover$abundance < 80)] <- "A"
percent_cover$abundance[which(percent_cover$abundance >=20 & percent_cover$abundance < 40)] <- "C"
percent_cover$abundance[which(percent_cover$abundance >=10 & percent_cover$abundance < 20)] <- "F"
percent_cover$abundance[which(percent_cover$abundance >=5 & percent_cover$abundance < 10)] <- "O"
percent_cover$abundance[which(percent_cover$abundance >=1 & percent_cover$abundance < 5)] <- "R"
percent_cover$abundance[which(percent_cover$abundance <1)] <- "N"

table(percent_cover$abundance)

#adding back to main dataset, first by removing these rows and then by adding the new transformed rows back in
hyperborea_records_DASSH <- hyperborea_records_DASSH%>%
  filter(!abundance_type %in% c("percentage","Percentage","percent"))
hyperborea_records_DASSH <- rbind(hyperborea_records_DASSH,percent_cover)

#extracting the SACFORN values to a dedicated column but excluding absences (N) because it isn't clear whether these were based on a consistent methodology - additional absences will be added later
hyperborea_records_DASSH$SACFORN <- NA
sacfor_vals <- c("O","C","F","A","R","S")

hyperborea_records_DASSH$SACFORN[which((hyperborea_records_DASSH$abundance %in% sacfor_vals)==TRUE)] <- hyperborea_records_DASSH$abundance[which((hyperborea_records_DASSH$abundance %in% sacfor_vals)==TRUE)]

#generating ordered integer SACFORN scale
hyperborea_records_DASSH$SACFORN_int <- as.integer(factor(hyperborea_records_DASSH$SACFORN,levels = c("R","O","F","C","A","S"),ordered = TRUE))

#filtering the data to exclude now-irrelevant columns
hyperborea_records_DASSH <- hyperborea_records_DASSH%>%
  dplyr::select(originator,latitude,longitude,date,Presence,SACFORN,SACFORN_int,depth)
#renaming
names(hyperborea_records_DASSH) <- c("Origin","Lat","Lon","Date","Presence","SACFORN","SACFORN_int","depth")

#subsetting to remove all DASSH records that do not have a SACFORN value associated:
hyperborea_records_DASSH <- hyperborea_records_DASSH%>%
  drop_na(SACFORN_int)

hyperborea_records_DASSH$Origin <- "DASSH"

## Processing the MNCR data ##
#selecting relevant columns
hyperborea_records_MNCR <- hyperborea_records_MNCR%>%
  dplyr::select(STARTDATE,LAT,LONG,nLaminaria_hyperb,avgd)

#renaming
names(hyperborea_records_MNCR) <- c("Date","Lat","Lon","SACFORN_int","depth")

hyperborea_records_MNCR$SACFORN <- NA
hyperborea_records_MNCR$SACFORN[which(hyperborea_records_MNCR$SACFORN_int==6)] <- "S"
hyperborea_records_MNCR$SACFORN[which(hyperborea_records_MNCR$SACFORN_int==5)] <- "A"
hyperborea_records_MNCR$SACFORN[which(hyperborea_records_MNCR$SACFORN_int==4)] <- "C"
hyperborea_records_MNCR$SACFORN[which(hyperborea_records_MNCR$SACFORN_int==3)] <- "F"
hyperborea_records_MNCR$SACFORN[which(hyperborea_records_MNCR$SACFORN_int==2)] <- "O"
hyperborea_records_MNCR$SACFORN[which(hyperborea_records_MNCR$SACFORN_int==1)] <- "R"
hyperborea_records_MNCR$SACFORN[which(hyperborea_records_MNCR$SACFORN_int==0)] <- "N"

## Combining the datasets 
hyperborea_records_MNCR$Origin <- "MNCR"

#flipping the depths to align with DASSH
hyperborea_records_MNCR$depth <- hyperborea_records_MNCR$depth*(-1)

dat.combined <- rbind(hyperborea_records_MNCR[c("Date","Lat","Lon","Origin","SACFORN","SACFORN_int","depth")],hyperborea_records_DASSH[c("Origin","Date","Lat","Lon","SACFORN","SACFORN_int","depth")])

#dropping cells with missing depth values
dat.combined <- dat.combined %>%
  drop_na(depth)

## there are some rows which have duplicate Date, coordinates and depth. For these locations, we combine the duplicates by keeping the maximum SACFORN value
dat.combined.refin <- dat.combined%>%
  dplyr::group_by(Date,Lat,Lon,depth)%>% #including depth keeps many values from same coordinates with unique depths
  dplyr::slice_max(
    SACFORN_int, 
    n = 1, # Ensures only one row is selected per group
    with_ties = FALSE # Optional: If two rows have the same max SACFORN_int, this keeps only one. 
  ) %>%  dplyr::ungroup()

## Now adding in the absences derived from DASSH for areas of infralittoral and sublittoral hard substrate ##
absence_data <- read.csv("Data/Raw_data/Kelp_records/L_hyperborea_absences_DASSH.csv")

#remove all locations which contained kelp taxa to get full kelp absences
absence_data <- absence_data %>%
  group_by(latitude,longitude) %>%
  mutate(group_contains_string = any(grepl(paste(c("Laminaria","Saccharina","Saccorhiza"),collapse = "|"), taxon))) %>% filter(group_contains_string == FALSE)%>%
  ungroup()

# remove all rows which have 0 or N for SACFORN
absence_data <- absence_data %>%
  filter(!taxon %in% c(""))%>%
  filter(!abundance %in% c(0,"N",""))

#dropping all MNCR records
absence_data <- absence_data %>%
  filter(!str_detect(survey_code, "MNCR"))

#now keeping only the unique locations
absence_data <- absence_data%>%
  distinct(across(all_of(c("latitude","longitude"))), .keep_all = TRUE)

#adding a presence/absence column and setting all to absent
absence_data$PA <- 0

#refining the data
absence_data <- absence_data%>%
  select(latitude,longitude,depth,date,PA,biotope_name)

absence_data$Date <- absence_data$date
absence_data$Lat <- absence_data$latitude
absence_data$Lon <- absence_data$longitude
absence_data$SACFORN_int <- 0
absence_data$SACFORN <- "N"
absence_data$Origin <- "DASSH_absences"

#dropping records with NA values for depth
summary(absence_data$depth)
absence_data <- absence_data%>%
  drop_na(depth)

#combining the absences with the other kelp data
full_dataset <- rbind(dat.combined.refin,absence_data[,c("Date","Lat","Lon","Origin","SACFORN","SACFORN_int","depth")])

#restricting the data to depth values between 0 and 30m
summary(full_dataset$depth)
full_dataset <- full_dataset%>%
  filter(depth<30 & depth >= 0)

saveRDS(full_dataset,"Data/Refined_data/full_species_data.rds")
write.csv(full_dataset,"Data/Refined_data/full_species_data.csv",row.names = FALSE)

#### ==========    2. PROCESSING PRESENT-DAY PREDICTOR DATASETS  ================================####
# In this step all of the predictor datasets are loaded and processed to match extents and resolutions
# This must be done seperately for the UK and Ireland due to the fact that a seperate bathymetry dataset is in use for Ireland
## Loading the predictor layers ##

# Reading in the bathymetry and wave fetch datasets:
wave_fetch <- terra::rast("Data/Raw_data/Predictors/wave_fetch.tif")
bathy_data_UK <- terra::rast("Data/Raw_data/Predictors/UK_marine_DEM_1_arcsecond_6nm.tif")
bathy_data_IRE <- terra::rast("Data/Raw_data/Predictors/Ireland_25m_bathymetry.tif")

#reading in a shapefile with the UK and Ireland nautical limits, simply for quickly cutting down the size of the rasters to reduce runtimes
UK_6nm <- terra::vect("Data/Raw_data/Geographical_layers/UK_6nm.shp")
IRE_6nm <- terra::vect("Data/Raw_data/Geographical_layers/Ireland_6nm.shp")

# Loading the other environmental variables from Bio-oracle:
# First for the UK:
# setting the general desired extent of the rasters
latitude = c(49.10065,61.10212)
longitude = c(-10.38137,3.88545)

# Temperature (at depth)
# biooracler::list_layers("temperature") #can be used to identify suitable datasets
# info_layer("thetao_baseline_2000_2019_depthmean") #investigate the individual layers within the dataset

dataset_id <- c("thetao_baseline_2000_2019_depthmean")
variables = c("thetao_mean") #getting the average value for the coldest/warmest month, and the overall decadal annual averages
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
temperature_layers <- download_layers(dataset_id, variables, constraints)

#getting the average value across both time periods
temperature_layers$temp_mean <- mean(c(temperature_layers$thetao_mean_1,temperature_layers$thetao_mean_2))

# Nitrate (at depth)
dataset_id <- c("no3_baseline_2000_2018_depthmean")
variables = c("no3_mean") #getting the mean and min
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
no3_layers <- download_layers(dataset_id, variables, constraints)

no3_layers$no3_mean <- mean(c(no3_layers$no3_mean_1,no3_layers$no3_mean_2))

# Kd_PAR (surface)
dataset_id <- c("kdpar_mean_baseline_2000_2020_depthsurf")
variables = c("kdpar_mean_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
kd_PAR_layers <- download_layers(dataset_id, variables, constraints)
kd_PAR_layers$kdpar_mean <- mean(c(kd_PAR_layers$kdpar_mean_mean_1,kd_PAR_layers$kdpar_mean_mean_2))

#PAR (surface)
dataset_id <- c("par_mean_baseline_2000_2020_depthsurf")
variables = c("par_mean_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
PAR_layers <- download_layers(dataset_id, variables, constraints)
PAR_layers$PAR_mean <- mean(c(PAR_layers$par_mean_mean_1,PAR_layers$par_mean_mean_2))

#seawater velocity (depth)
dataset_id <- c("sws_baseline_2000_2019_depthmean")
variables = c("sws_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
sea_vel_depth_layers <- download_layers(dataset_id, variables, constraints)

sea_vel_depth_layers$veloc_mean <- mean(c(sea_vel_depth_layers$sws_mean_1,sea_vel_depth_layers$sws_mean_2))

# Selecting predictors 
predictors_UK <- c(temperature_layers$temp_mean,no3_layers$no3_mean,kd_PAR_layers$kdpar_mean,PAR_layers$PAR_mean,sea_vel_depth_layers$veloc_mean)

# Now for Ireland:
latitude = c(50,56)
longitude = c(-11,-5.5)

# Temperature (at depth)
dataset_id <- c("thetao_baseline_2000_2019_depthmean")
variables = c("thetao_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
temperature_layers <- download_layers(dataset_id, variables, constraints)

#getting the average value across both time periods
temperature_layers$temp_mean <- mean(c(temperature_layers$thetao_mean_1,temperature_layers$thetao_mean_2))

# Nitrate (at depth)
dataset_id <- c("no3_baseline_2000_2018_depthmean")
# dataset_id <- c("no3_baseline_2000_2018_depthsurf")
variables = c("no3_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
no3_layers <- download_layers(dataset_id, variables, constraints)

no3_layers$no3_mean <- mean(c(no3_layers$no3_mean_1,no3_layers$no3_mean_2))

# Kd_PAR (surface)
dataset_id <- c("kdpar_mean_baseline_2000_2020_depthsurf")
variables = c("kdpar_mean_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
kd_PAR_layers <- download_layers(dataset_id, variables, constraints)

kd_PAR_layers$kdpar_mean <- mean(c(kd_PAR_layers$kdpar_mean_mean_1,kd_PAR_layers$kdpar_mean_mean_2))

#PAR (surface)
dataset_id <- c("par_mean_baseline_2000_2020_depthsurf")
variables = c("par_mean_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
PAR_layers <- download_layers(dataset_id, variables, constraints)

PAR_layers$PAR_mean <- mean(c(PAR_layers$par_mean_mean_1,PAR_layers$par_mean_mean_2))

#seawater velocity (depth)
dataset_id <- c("sws_baseline_2000_2019_depthmean")
variables = c("sws_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
sea_vel_depth_layers <- download_layers(dataset_id, variables, constraints)

sea_vel_depth_layers$veloc_mean <- mean(c(sea_vel_depth_layers$sws_mean_1,sea_vel_depth_layers$sws_mean_2))

# Selecting predictors 
predictors_IRE <- c(temperature_layers$temp_mean,no3_layers$no3_mean,kd_PAR_layers$kdpar_mean,PAR_layers$PAR_mean,sea_vel_depth_layers$veloc_mean)

## Processing predictor layers ##

#masking  the predictors to the 6nm limit to speed up computation (not needed for wave fetch as this is already within the limit)
predictors_UK <- (terra::crop(predictors_UK, UK_6nm, snap = "out", mask = TRUE))
bathy_data_UK <- (terra::crop(bathy_data_UK, UK_6nm, snap = "out", mask = TRUE))

predictors_IRE <- (terra::crop(predictors_IRE, Ire_6nm, snap = "out", mask = TRUE))
bathy_data_IRE <- (terra::crop(bathy_data_IRE, project(Ire_6nm,crs(bathy_data_IRE)), snap = "out", mask = TRUE))

#reprojecting layers to British National Grid for consistency:
predictors_UK_reproj <- project(predictors_UK,crs(wave_fetch))
bathy_UK_reproj <- project(bathy_data_UK,crs(wave_fetch))

predictors_IRE_reproj <- project(predictors_IRE,crs(wave_fetch))
bathy_IRE_reproj <- project(bathy_data_IRE,crs(wave_fetch))

# Filling in missing raster data for the Bio-Oracle predictors by imputing data for empty cells
predictors_UK_focal <- terra::focal(predictors_UK_reproj, 3, "mean",na.policy="only", na.rm=TRUE)
predictors_IRE_focal <- terra::focal(predictors_IRE_reproj, 3, "mean",na.policy="only", na.rm=TRUE)

## Resampling environmental layers to bathymetry resolution ##
predictors_UK_focal_resampled <- terra::resample(predictors_UK_focal,bathy_UK_reproj, method="bilinear")
wave_fetch_UK_resampled <- terra::resample(wave_fetch, bathy_UK_reproj, method="bilinear")

predictors_IRE_focal_resampled <- terra::resample(predictors_IRE_focal,bathy_IRE_reproj, method="bilinear")
wave_fetch_IRE_resampled <- terra::resample(wave_fetch, bathy_IRE_reproj, method="bilinear")

## combining the bathymetry, wave fetch and other data
predictor_layers_UK_resampled <- c(predictors_UK_focal_resampled,bathy_UK_reproj,wave_fetch_UK_resampled)
predictor_layers_IRE_resampled <- c(predictors_IRE_focal_resampled,bathy_IRE_reproj,wave_fetch_IRE_resampled)

#double check order before renaming!
names(predictor_layers_UK_resampled) <- c("temp_mean","nitrate_mean","kdpar_mean","PAR_mean","cur_vel_dep_mean","bathymetry","wave_fetch")

names(predictor_layers_IRE_resampled) <-c("temp_mean","nitrate_mean","kdpar_mean","PAR_mean","cur_vel_dep_mean","bathymetry","wave_fetch")

## masking the predictors ##
## cropping the data to the wave fetch coverage to align raster extents
predictor_layers_UK_resampled_crop <- (terra::crop(predictor_layers_UK_resampled, predictor_layers_UK_resampled$wave_fetch, snap = "out", mask = TRUE))

predictor_layers_IRE_resampled_crop <- (terra::crop(predictor_layers_IRE_resampled, predictor_layers_IRE_resampled$wave_fetch, snap = "out", mask = TRUE))

## Cutting the bathymetry data to only include values above -50m and below 5m. The -50m is so that later regridding does not exclude records which are within cells with a mean depth above 30m, although the actual depth values from the samples will be used for fitting models and rasters down to 30m will be used for predictions. Using +5m is because the raster depths aren't always completely accurate and so ensures we don't lose many records which are just beyond the shallow cutoff.

predictor_layers_UK_resampled_crop$bathymetry[predictor_layers_UK_resampled_crop$bathymetry>(5) | predictor_layers_UK_resampled_crop$bathymetry<=(-50)] <- NA

predictor_layers_IRE_resampled_crop$bathymetry[predictor_layers_IRE_resampled_crop$bathymetry>(5) | predictor_layers_IRE_resampled_crop$bathymetry<=(-50)] <- NA

# now masking the other predictors to the bathymetry layer to ensure all are aligned
predictor_layers_UK_resampled_crop <- (terra::crop(predictor_layers_UK_resampled_crop, predictor_layers_UK_resampled_crop$bathymetry, snap = "out", mask = TRUE))

predictor_layers_IRE_resampled_crop <- (terra::crop(predictor_layers_IRE_resampled_crop, predictor_layers_IRE_resampled_crop$bathymetry, snap = "out", mask = TRUE))

## creating PAR at depth variable by combining the mean depth in each cell with the PAR and kPAR. following the formula:
# PAR(z) = PAR(0) * e^(-KPAR * z)
predictor_layers_UK_resampled_crop$bathymetry <- predictor_layers_UK_resampled_crop$bathymetry*(-1) #converting depths to positive values
par_data <- predictor_layers_UK_resampled_crop$PAR_mean
kpar_data <- predictor_layers_UK_resampled_crop$kdpar_mean
depth <- predictor_layers_UK_resampled_crop$bathymetry
PAR_depth <- par_data * exp(-kpar_data * depth)
names(PAR_depth) <- "PAR_depth"

#adding this to the predictors and removing the now-defunct kdPAR and PAR_mean variables
predictor_layers_UK_resampled_crop <- c(predictor_layers_UK_resampled_crop,PAR_depth)
predictor_layers_UK_resampled_crop <- predictor_layers_UK_resampled_crop[[-which(names(predictor_layers_UK_resampled_crop)=="kdpar_mean"|names(predictor_layers_UK_resampled_crop)=="PAR_mean")]] #removing unnecessary variables

# writeRaster(predictor_layers_UK_resampled_crop, "../Data/Refined_data/Kelp_SDM/Predictors_UK_present_depth50m.tif",overwrite=TRUE)

#now cropping to 30m for modelling
predictor_layers_UK_resampled_crop$bathymetry[predictor_layers_UK_resampled_crop$bathymetry>(30)] <- NA

predictor_layers_UK_resampled_crop <- (terra::crop(predictor_layers_UK_resampled_crop, predictor_layers_UK_resampled_crop$bathymetry, snap = "out", mask = TRUE))

writeRaster(predictor_layers_UK_resampled_crop, "Data/Refined_data/Predictors_UK_present.tif",overwrite=TRUE)

## BE AWARE THAT THE PREDICTOR RASTER NEEDS TO BE CONVERTED TO .RDS FOR THE MODELLING. THE UK RASTER IS TOO LARGE, SO THIS WAS ACHIEVED BY TILING IT AND THEN CONVERTING TILES INDIVIDUALLY BEFORE RECOMBINING. THE CODE IS NOT PROVIDED HERE.

predictor_layers_IRE_resampled_crop$bathymetry <- predictor_layers_IRE_resampled_crop$bathymetry*(-1)
par_data <- predictor_layers_IRE_resampled_crop$PAR_mean
kpar_data <- predictor_layers_IRE_resampled_crop$kdpar_mean
depth <- predictor_layers_IRE_resampled_crop$bathymetry
PAR_depth <- par_data * exp(-kpar_data * depth)
names(PAR_depth) <- "PAR_depth"

#adding this to the predictors and removing the now-defunct kdPAR and PAR_mean variables
predictor_layers_IRE_resampled_crop <- c(predictor_layers_IRE_resampled_crop,PAR_depth)
predictor_layers_IRE_resampled_crop <- predictor_layers_IRE_resampled_crop[[-which(names(predictor_layers_IRE_resampled_crop)=="kdpar_mean"|names(predictor_layers_IRE_resampled_crop)=="PAR_mean")]] #removing unnecessary variables

# writeRaster(predictor_layers_IRE_resampled_crop, "../Data/Refined_data/Kelp_SDM/Predictors_ireland_present_depth50m.tif",overwrite=TRUE)

#now cropping to 30m for modelling
predictor_layers_IRE_resampled_crop$bathymetry[predictor_layers_IRE_resampled_crop$bathymetry>(30)] <- NA

predictor_layers_IRE_resampled_crop <- (terra::crop(predictor_layers_IRE_resampled_crop, predictor_layers_IRE_resampled_crop$bathymetry, snap = "out", mask = TRUE))

writeRaster(predictor_layers_IRE_resampled_crop, "Data/Refined_data/Predictors_IRE_present.tif",overwrite=TRUE)

# REMEMBER, THE IRE PREDICTOR RASTER NEEDS TO BE CONVERTED TO .RDS FOR MODELLING. 

#### ========== 3. PROCESSING FUTURE PREDICTOR DATASETS ============ #########
#Here we extract raster layers for the same predictors but under an SSP2-4.5 and SSP5-8.5 scenaRio
#loading the bathymetry, PAR and wave fetch rasters for present day as these are unchanged for future modelling
predictors_UK_present <- terra::rast("Data/Refined_data/Predictors_UK_present.tif")
predictors_IRE_present <- terra::rast("Data/Refined_data/Predictors_IRE_present.tif")

#extracting the bathymetry, PAR and wave fetch layers
PAR_wave_bathy_UK <- predictors_UK_present[[c("bathymetry","wave_fetch","PAR_depth")]] 
PAR_wave_bathy_IRE <- predictors_IRE_present[[c("bathymetry","wave_fetch","PAR_depth")]] 

## Processing predictors for UK:
# setting the general desired extent of the rasters
latitude = c(49.10065,61.10212)
longitude = c(-10.38137,3.88545)

# Temperature (at depth)
dataset_id <- c("thetao_ssp585_2020_2100_depthmean")
variables = c("thetao_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
temperature_585_layers <- download_layers(dataset_id, variables, constraints)

dataset_id <- c("thetao_ssp245_2020_2100_depthmean")
variables = c("thetao_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
temperature_245_layers <- download_layers(dataset_id, variables, constraints)

# Nitrate (at depth)
dataset_id <- c("no3_ssp585_2020_2100_depthmean")
# dataset_id <- c("no3_baseline_2000_2018_depthsurf")
variables = c("no3_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
no3_585_layers <- download_layers(dataset_id, variables, constraints)

dataset_id <- c("no3_ssp245_2020_2100_depthmean")
# dataset_id <- c("no3_baseline_2000_2018_depthsurf")
variables = c("no3_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
no3_245_layers <- download_layers(dataset_id, variables, constraints)

#seawater velocity (depth)
dataset_id <- c("sws_ssp585_2020_2100_depthmean")
variables = c("sws_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
sea_vel_depth_585_layers <- download_layers(dataset_id, variables, constraints)

dataset_id <- c("sws_ssp245_2020_2100_depthmean")
variables = c("sws_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
sea_vel_depth_245_layers <- download_layers(dataset_id, variables, constraints)

# selecting the specific time periods of interest (2090-2100)
predictors585 <- c(temperature_585_layers$thetao_mean_8,no3_585_layers$no3_mean_8,sea_vel_depth_585_layers$sws_mean_8)

predictors245 <- c(temperature_245_layers$thetao_mean_8,no3_245_layers$no3_mean_8,sea_vel_depth_245_layers$sws_mean_8)

names(predictors585) <- c("temp_mean_585","nitrate_mean_585","velocity_mean_585")
names(predictors245) <- c("temp_mean_245","nitrate_mean_245","velocity_mean_245")

predictors_UK_future <- c(predictors585,predictors245)

#masking  the predictors to the spatial limits to speed up computation (not needed for wave fetch as this is already within the limit)
predictors_UK_future <- (terra::crop(predictors_UK_future, UK_6nm, snap = "out", mask = TRUE))
#reprojecting layers to BNG for consistency:
predictors_UK_future_reproj <- project(predictors_UK_future,crs(PAR_wave_bathy_UK))

predictors_UK_future_focal <- terra::focal(predictors_UK_future_reproj, 3, "mean",na.policy="only", na.rm=TRUE)

## Resampling environmental layers to bathymetry resolution
predictors_UK_future_focal_resampled <- terra::resample(predictors_UK_future_focal,PAR_wave_bathy_UK$bathymetry, method="bilinear")

## combining the bathymetry, wave fetch and other data
predictors_UK_future_focal_resampled <- c(predictors_UK_future_focal_resampled,PAR_wave_bathy_UK)

## masking the predictors ##
#masking the other predictors to the bathymetry layer
predictors_UK_future_focal_resampled2 <- (terra::crop(predictors_UK_future_focal_resampled, predictors_UK_future_focal_resampled$bathymetry, snap = "out", mask = TRUE))

predictors_UK_245 <- predictors_UK_future_focal_resampled2[[c("temp_mean_245","nitrate_mean_245","velocity_mean_245","bathymetry","wave_fetch","PAR_depth")]]
predictors_UK_585 <- predictors_UK_future_focal_resampled2[[c("temp_mean_585","nitrate_mean_585","velocity_mean_585","bathymetry","wave_fetch","PAR_depth")]]

names(predictors_UK_245) <- c("temp_mean","nitrate_mean","cur_vel_dep_mean","bathymetry","wave_fetch","PAR_depth")
names(predictors_UK_585) <- c("temp_mean","nitrate_mean","cur_vel_dep_mean","bathymetry","wave_fetch","PAR_depth")

writeRaster(predictors_UK_245, "Data/Refined_data/Predictors_UK_ssp245.tif", overwrite=TRUE)
writeRaster(predictors_UK_585, "Data/Refined_data/Predictors_UK_ssp585.tif", overwrite=TRUE)

#REMEMBER, THE PREDICTOR RASTERS MUST BE CONVERTED TO .RDS FOR MODELLING

## Processing predictors for IRE:
## Loading the other environmental variables from Bio-oracle:
# setting the general desired extent of the rasters
latitude = c(50,56)
longitude = c(-11,-5.5)

dataset_id <- c("thetao_ssp585_2020_2100_depthmean")
variables = c("thetao_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
temperature_585_layers <- download_layers(dataset_id, variables, constraints)

dataset_id <- c("thetao_ssp245_2020_2100_depthmean")
variables = c("thetao_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
temperature_245_layers <- download_layers(dataset_id, variables, constraints)

# Nitrate (at depth)
biooracler::list_layers("nitrate")
dataset_id <- c("no3_ssp585_2020_2100_depthmean")
# dataset_id <- c("no3_baseline_2000_2018_depthsurf")
variables = c("no3_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
no3_585_layers <- download_layers(dataset_id, variables, constraints)

dataset_id <- c("no3_ssp245_2020_2100_depthmean")
# dataset_id <- c("no3_baseline_2000_2018_depthsurf")
variables = c("no3_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
no3_245_layers <- download_layers(dataset_id, variables, constraints)

#seawater velocity (depth)
dataset_id <- c("sws_ssp585_2020_2100_depthmean")
variables = c("sws_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
sea_vel_depth_585_layers <- download_layers(dataset_id, variables, constraints)

dataset_id <- c("sws_ssp245_2020_2100_depthmean")
variables = c("sws_mean")
constraints = list(latitude, longitude)
names(constraints) = c("latitude", "longitude")
sea_vel_depth_245_layers <- download_layers(dataset_id, variables, constraints)

# selecting the specific time periods of interest (2090-2100)
predictors585 <- c(temperature_585_layers$thetao_mean_8,no3_585_layers$no3_mean_8,sea_vel_depth_585_layers$sws_mean_8)

predictors245 <- c(temperature_245_layers$thetao_mean_8,no3_245_layers$no3_mean_8,sea_vel_depth_245_layers$sws_mean_8)

names(predictors585) <- c("temp_mean_585","nitrate_mean_585","velocity_mean_585")
names(predictors245) <- c("temp_mean_245","nitrate_mean_245","velocity_mean_245")

predictors_IRE_future <- c(predictors585,predictors245)

#masking  the predictors to the spatial limits to speed up computation (not needed for wave fetch as this is already within the limit)
predictors_IRE_future <- (terra::crop(predictors_IRE_future, IRE_6nm, snap = "out", mask = TRUE))
#reprojecting layers to BNG for consistency:
predictors_IRE_future_reproj <- project(predictors_IRE_future,crs(PAR_wave_bathy_IRE))

predictors_IRE_future_focal <- terra::focal(predictors_IRE_future_reproj, 3, "mean",na.policy="only", na.rm=TRUE)

## Resampling environmental layers to bathymetry resolution
predictors_IRE_future_focal_resampled <- terra::resample(predictors_IRE_future_focal,PAR_wave_bathy_IRE$bathymetry, method="bilinear")

## combining the bathymetry, wave fetch and other data
predictors_IRE_future_focal_resampled <- c(predictors_IRE_future_focal_resampled,PAR_wave_bathy_IRE)

## masking the predictors ##
#masking the other predictors to the bathymetry layer
predictors_IRE_future_focal_resampled2 <- (terra::crop(predictors_IRE_future_focal_resampled, predictors_IRE_future_focal_resampled$bathymetry, snap = "out", mask = TRUE))

predictors_IRE_245 <- predictors_IRE_future_focal_resampled2[[c("temp_mean_245","nitrate_mean_245","velocity_mean_245","bathymetry","wave_fetch","PAR_depth")]]
predictors_IRE_585 <- predictors_IRE_future_focal_resampled2[[c("temp_mean_585","nitrate_mean_585","velocity_mean_585","bathymetry","wave_fetch","PAR_depth")]]

names(predictors_IRE_245) <- c("temp_mean","nitrate_mean","cur_vel_dep_mean","bathymetry","wave_fetch","PAR_depth")
names(predictors_IRE_585) <- c("temp_mean","nitrate_mean","cur_vel_dep_mean","bathymetry","wave_fetch","PAR_depth")

writeRaster(predictors_IRE_245, "Data/Refined_data/Predictors_IRE_ssp245.tif", overwrite=TRUE)
writeRaster(predictors_IRE_585, "Data/Refined_data/Predictors_IRE_ssp585.tif", overwrite=TRUE)

#REMEMBER, THE PREDICTOR RASTERS MUST BE CONVERTED TO .RDS FOR MODELLING

#### ==========  4. REGRIDDING AND SPATIALLY THINNING SPECIES RECORDS  ============= ####
rm(list=ls())
library(terra)  # mapping
library(tidyverse) # for cleaning

#load the full species dataset
species_data <- readRDS("Data/Refined_data/full_species_data.rds")

species_data$Year <- year(parse_date_time(species_data$Date, orders = c("dmy", "ymd")))

## load the full extent predictor layers - for this using the layers down to 50m depth to make sure that presence records are not cut unneccesarily (for fitting the models the actual depth values are being used anyway)
predictors_UK <- rast("Data/Refined_data/Predictors_UK_present.tif")
predictors_IRE <- rast("Data/Refined_data/Predictors_IRE_present.tif")

## Regridding ##

#converting the species data to spatial object
species_data_spat <- terra::vect(species_data, geom = c("Lon", "Lat"), keepgeom = TRUE, crs = "EPSG:4326")
species_data_spat <- project(species_data_spat,crs(predictors_UK))

## Regridding the records to the scale of the bathymetry data ##
## Regidding the species records by manually extracting the cell number for the bathymetry raster for each species record and extracting the max SACFORN value in each cell. 
# This needs to be done seperately for UK and IRE

#extracting the raster cell values including cell ID and adding to species data
raster_extract_UK <- terra::extract(predictors_UK,species_data_spat,cells=TRUE,xy=TRUE)
species_data_comb_UK <- cbind(species_data,raster_extract_UK)

raster_extract_IRE <- terra::extract(predictors_IRE,species_data_spat,cells=TRUE,xy=TRUE)
species_data_comb_IRE <- cbind(species_data,raster_extract_IRE)

#now for each grid cell, take the maximum SACFOR value observed.
species_data_agg_UK <- species_data_comb_UK%>%
  group_by(cell)%>%
  slice_max(order_by = SACFORN_int,n=1,with_ties = FALSE) %>%
  ungroup()

species_data_agg_IRE <- species_data_comb_IRE%>%
  group_by(cell)%>%
  slice_max(order_by = SACFORN_int,n=1,with_ties = FALSE) %>%
  ungroup()

#replacing the location coordinates for the species with the central values for each cell
species_data_agg_UK$Lon <- xyFromCell(predictors_UK, species_data_agg_UK$cell)[,1]
species_data_agg_UK$Lat <- xyFromCell(predictors_UK, species_data_agg_UK$cell)[,2]

species_data_agg_IRE$Lon <- xyFromCell(predictors_IRE, species_data_agg_IRE$cell)[,1]
species_data_agg_IRE$Lat <- xyFromCell(predictors_IRE, species_data_agg_IRE$cell)[,2]

# A lot of the records occur in areas that have missing values for one or more predictors. Much of this will be due to points from beyond the focal study region and raster (e.g. isle of Man, Guernsey etc.). 
summary(species_data_agg_UK)
summary(species_data_agg_IRE)

#These points with missing values are removed
species_data_agg_filt_UK <- species_data_agg_UK%>%
  drop_na(temp_mean,nitrate_mean,cur_vel_dep_mean,PAR_depth,bathymetry,wave_fetch)

summary(species_data_agg_filt_UK)

species_data_agg_filt_IRE <- species_data_agg_IRE%>%
  drop_na(temp_mean,nitrate_mean,cur_vel_dep_mean,PAR_depth,bathymetry,wave_fetch)
summary(species_data_agg_filt_IRE)

table(species_data_agg_filt_UK$SACFORN_int)
table(species_data_agg_filt_IRE$SACFORN_int)

## Spatial thinning ##
# Many of the records are highly clustered, which could introduce biases. Will manually thin the data, keeping the maximum SACFOR values. This will be done for the presences and absences separately to ensure the dataset remains balanced

# extracting presences
presences_UK <- species_data_agg_filt_UK %>%
  filter(SACFORN_int>0)
presences_UK_spat <- terra::vect(presences_UK, geom = c("x", "y"), keepgeom = TRUE, crs = "EPSG:27700")

presences_IRE <- species_data_agg_filt_IRE %>%
  filter(SACFORN_int>0)
presences_IRE_spat <- terra::vect(presences_IRE, geom = c("x", "y"), keepgeom = TRUE, crs = "EPSG:27700")

#extracting absences
absences_UK <- species_data_agg_filt_UK %>%
  filter(SACFORN_int==0)
absences_UK_spat <- terra::vect(absences_UK, geom = c("Lon", "Lat"), keepgeom = TRUE, crs = "EPSG:27700")

absences_IRE <- species_data_agg_filt_IRE %>%
  filter(SACFORN_int==0)
absences_IRE_spat <- terra::vect(absences_IRE, geom = c("Lon", "Lat"), keepgeom = TRUE, crs = "EPSG:27700")

#Will use a 100m grid to thin the data
#creating the grid by resampling the predictor layers
raster_UK_100m <- rast(ext(predictors_UK$bathymetry), res = 100)
crs(raster_UK_100m) <- "EPSG:27700"
values(raster_UK_100m) <- 1 #generating a mask by setting all non-NA values to 1

raster_IRE_100m <- rast(ext(predictors_IRE$bathymetry), res = 100)
crs(raster_IRE_100m) <- "EPSG:27700"
values(raster_IRE_100m) <- 1 #generating a mask by setting all non-NA values to 1

#extracting the cell IDs for each presence location and assigning to the presences
presences_UK_extract_100 <- terra::extract(raster_UK_100m,presences_UK_spat,cells=TRUE)
presences_UK$cell100 <- presences_UK_extract_100$cell

absences_UK_extract_100 <- terra::extract(raster_UK_100m,absences_UK_spat,cells=TRUE)
absences_UK$cell100 <- absences_UK_extract_100$cell

presences_IRE_extract_100 <- terra::extract(raster_IRE_100m,presences_IRE_spat,cells=TRUE)
presences_IRE$cell100 <- presences_IRE_extract_100$cell

absences_IRE_extract_100 <- terra::extract(raster_IRE_100m,absences_IRE_spat,cells=TRUE)
absences_IRE$cell100 <- absences_IRE_extract_100$cell

#keeping only one record with the maximum SACFORN within each of the 100m cells
presences_UK <- presences_UK%>%
  group_by(cell100)%>%
  slice_max(SACFORN_int, n = 1,with_ties = FALSE)%>%
  ungroup()

absences_UK <- absences_UK%>%
  group_by(cell100)%>%
  slice_max(SACFORN_int, n = 1,with_ties = FALSE)%>%
  ungroup()

presences_IRE <- presences_IRE%>%
  group_by(cell100)%>%
  slice_max(SACFORN_int, n = 1,with_ties = FALSE)%>%
  ungroup()

absences_IRE <- absences_IRE%>%
  group_by(cell100)%>%
  slice_max(SACFORN_int, n = 1,with_ties = FALSE)%>%
  ungroup()

#combining the data
thinned_data_UK <- rbind(presences_UK,absences_UK)
thinned_data_IRE <-  rbind(presences_IRE,absences_IRE)

thinned_data_UK$Origin <- "UK"
thinned_data_IRE$Origin <- "IRE"

combined_data <- rbind(thinned_data_UK,thinned_data_IRE)
names(combined_data)

combined_data <- combined_data%>%
  dplyr::select(c(Date,Lat,Lon,SACFORN,SACFORN_int,depth,temp_mean,nitrate_mean,cur_vel_dep_mean,bathymetry,wave_fetch,PAR_depth,cell))

#changing column names for sampled depth and bathymetry from raster to ensure matching names for modelling and prediction
names(combined_data)[which(names(combined_data)=="bathymetry")] <- "bathymetry_rast"
names(combined_data)[which(names(combined_data)=="depth")] <- "bathymetry"

nrow(combined_data[combined_data$SACFORN_int>0,])

saveRDS(combined_data,"Data/Refined_data/full_species_data_thinned.rds")
write.csv(combined_data,"Data/Refined_data/full_species_data_thinned.csv")
