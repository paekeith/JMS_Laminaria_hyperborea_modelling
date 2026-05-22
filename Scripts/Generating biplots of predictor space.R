################################
# Code for generating biplots of the distribution of samples within the predictor space
# As part of Eskuche-Keith et al. "Abundance-based distribution modelling reveals climate-driven loss of functionally important kelp forest habitat"
# Last edited 22/05/2026
################################

library(terra)
library(ggplot2)
library(tidyverse)
library(ggExtra) #for adding histograms to the scatterplot margins if desired
library(cowplot)

#loading predictors
predictors_UK_present <- readRDS("Data/Refined_data/Predictors_UK_present.rds")
predictors_IRE_present <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")
predictors_present <- rbind(predictors_UK_present,predictors_IRE_present)
predictors_present <- predictors_present %>%
  na.omit()

predictors_UK_sp245 <- readRDS("Data/Refined_data/Predictors_UK_ssp245.rds")
predictors_IRE_sp245 <- readRDS("Data/Refined_data/Predictors_IRE_ssp245.rds")
predictors_sp245 <- rbind(predictors_UK_sp245,predictors_IRE_sp245)
predictors_sp245 <- predictors_sp245 %>%
  na.omit()

predictors_UK_sp585 <- readRDS("Data/Refined_data/Predictors_UK_ssp585.rds")
predictors_IRE_sp585 <- readRDS("Data/Refined_data/Predictors_IRE_ssp585.rds")
predictors_sp585 <- rbind(predictors_UK_sp585,predictors_IRE_sp585)
predictors_sp585 <- predictors_sp585 %>%
  na.omit()

#identifying the ranges of each predictor for plotting
summary(predictors_present$temp_mean)
summary(predictors_sp245$temp_mean)
summary(predictors_sp585$temp_mean)

summary(predictors_present$nitrate_mean)
summary(predictors_sp245$nitrate_mean)
summary(predictors_sp585$nitrate_mean)

summary(predictors_present$cur_vel_dep_mean)
summary(predictors_sp245$cur_vel_dep_mean)
summary(predictors_sp585$cur_vel_dep_mean)


#generating the 2d surface for the predictors
row_num <- 100
biplot_data <- as.data.frame(matrix(nrow=row_num,ncol=0))
biplot_data$temp_mean <- seq(min(predictors_present$temp_mean),max(predictors_sp585$temp_mean),length.out=row_num)
biplot_data$nitrate_mean <- seq(min(predictors_sp585$nitrate_mean),max(predictors_present$nitrate_mean),length.out=row_num)
biplot_data$current_vel <- seq(min(predictors_sp245$cur_vel_dep_mean),max(predictors_present$cur_vel_dep_mean),length.out=row_num)

#loading the kelp data and extracting the future conditions at the points
data <- as.data.frame(readRDS("Data/Refined_data/full_species_data_thinned.rds"))

# splitting into just the presences
# i dont think i need to do this for forest, as I'm interested in how the conditions are changing across all of the datapoints we used.
# data_park <- data[which(data$SACFORN_int>=3),] 
# data_forest <- data[which(data$SACFORN_int>=5),] 

#loading the predictor rasters for future conditions
predictors_UK_245_rast <- rast("Data/Refined_data/Predictors_UK_ssp245.tif")
predictors_IRE_245_rast <- rast("Data/Refined_data/Predictors_IRE_ssp245.tif")
predictors_UK_585_rast <- rast("Data/Refined_data/Predictors_UK_ssp585.tif")
predictors_IRE_585_rast <- rast("Data/Refined_data/Predictors_IRE_ssp585.tif")

data_spat <- terra::vect(data, geom = c("Lon", "Lat"), keepgeom = TRUE, crs = "EPSG:27700")
# data_forest_spat <- terra::vect(data_forest, geom = c("Lon", "Lat"), keepgeom = TRUE, crs = "EPSG:27700")

raster_extract_245_UK <- terra::extract(predictors_UK_245_rast,data_spat,cells=TRUE,xy=TRUE)
raster_extract_245_IRE <- terra::extract(predictors_IRE_245_rast,data_spat,cells=TRUE,xy=TRUE)
raster_extract_585_UK <- terra::extract(predictors_UK_585_rast,data_spat,cells=TRUE,xy=TRUE)
raster_extract_585_IRE <- terra::extract(predictors_IRE_585_rast,data_spat,cells=TRUE,xy=TRUE)

raster_extract_245_UK <- raster_extract_245_UK %>%
  na.omit()
raster_extract_245_IRE <- raster_extract_245_IRE %>%
  na.omit()
raster_extract_585_UK <- raster_extract_585_UK %>%
  na.omit()
raster_extract_585_IRE <- raster_extract_585_IRE %>%
  na.omit()

data_sp245 <- rbind(raster_extract_245_UK,raster_extract_245_IRE)
data_sp585 <- rbind(raster_extract_585_UK,raster_extract_585_IRE)

data$group <- "Present"
data_sp245$group <- "SSP2-4.5"
data_sp585$group <- "SSP5-8.5"

data <- data[,c("temp_mean","nitrate_mean","cur_vel_dep_mean","group")]
data_sp245 <- data_sp245[,c("temp_mean","nitrate_mean","cur_vel_dep_mean","group")]
data_sp585 <- data_sp585[,c("temp_mean","nitrate_mean","cur_vel_dep_mean","group")]

data_full <- rbind(data,data_sp245,data_sp585)

summary(data_sp585$temp_mean)

#plotting temperature against nitrates
biplot_data_exp1 <- expand.grid(
  temp_mean = biplot_data$temp_mean,
  nitrate_mean = biplot_data$nitrate_mean
)

temp_nitr_plot <- ggplot(data=biplot_data_exp1, aes(x = temp_mean, y = nitrate_mean)) +
  geom_raster(fill = "white")+labs(x="Temperature (°C)",y="Nitrate conc. (mmol.m-3)")+
  theme_bw()+
  theme(legend.position = c(0.8,0.7,1,1),legend.direction = "vertical",axis.text = element_text(size=15,colour="black",family="serif"),axis.title=element_text(size=17,family="serif"),legend.title = element_text(size=12,family="serif"),legend.text = element_text(size=12,family="serif"),plot.margin = unit(c(0.5, 0.5, 0, 0), "cm"))+scale_x_continuous(expand=c(0,0),breaks=seq(8.5,15.5,1))+scale_y_continuous(breaks=seq(0,9.5,1),expand=c(0,0))+
  geom_point(data=data_full,aes(x=temp_mean,y=nitrate_mean,colour=group),size=1,alpha=0.3,pch=16)+
  scale_color_manual(name = "",values = c("#1B9E77", "#D95F02", "#E7298A"))+
  annotate("rect", xmin = min(predictors_present$temp_mean), xmax = max(predictors_present$temp_mean), ymin = min(predictors_present$nitrate_mean), ymax = max(predictors_present$nitrate_mean),
           alpha = 0,color = "black", linetype = "dashed")+
  guides(color = guide_legend(override.aes = list(size = 4)))

temp_nitr_plot <- ggMarginal(temp_nitr_plot, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.3)

#plotting temperature against current
biplot_data_exp2 <- expand.grid(
  temp_mean = biplot_data$temp_mean,
  current_vel = biplot_data$current_vel
)

temp_curr_plot <- ggplot(data=biplot_data_exp2, aes(x = temp_mean, y = current_vel)) +
  geom_raster(fill = "white")+labs(x="Temperature (°C)",y="Current velocity (m.s-1)")+
  # scale_fill_gradientn(colors = hcl.colors(20, "viridis"))+
  theme_bw()+
  theme(legend.position = "none",axis.text = element_text(size=15,colour="black",family="serif"),axis.title=element_text(size=17,family="serif"),plot.margin = unit(c(0.5, 0.5, 0, 0), "cm"))+scale_x_continuous(expand=c(0,0),breaks=seq(8.5,15.5,1))+scale_y_continuous(breaks=seq(0,0.18,0.05),expand=c(0,0))+
  geom_point(data=data_full,aes(x=temp_mean,y=cur_vel_dep_mean,colour=group),size=1,alpha=0.3,pch=16)+
  scale_color_manual(name = "Scenario",values = c("#1B9E77", "#D95F02", "#E7298A"))+
  annotate("rect", xmin = min(predictors_present$temp_mean), xmax = max(predictors_present$temp_mean), ymin = min(predictors_present$cur_vel_dep_mean), ymax = max(predictors_present$cur_vel_dep_mean),
           alpha = 0,color = "black", linetype = "dashed")

temp_curr_plot <- ggMarginal(temp_curr_plot, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.3)

#plotting nitrate against current
biplot_data_exp3 <- expand.grid(
  nitrate_mean = biplot_data$nitrate_mean,
  current_vel = biplot_data$current_vel
)

nitr_curr_plot <- ggplot(data=biplot_data_exp3, aes(x = nitrate_mean, y = current_vel)) +
  geom_raster(fill = "white")+
  labs(x="Nitrate conc. (mmol.m-3)",y="Current velocity (m.s-1)")+
  # scale_fill_gradientn(colors = hcl.colors(20, "viridis"))+
  theme_bw()+
  theme(legend.position = "none",axis.text = element_text(size=15,colour="black",family="serif"),axis.title=element_text(size=17,family="serif"),plot.margin = unit(c(0.5, 0.5, 0, 0), "cm"))+scale_x_continuous(expand=c(0,0),breaks=seq(0,9.5,1))+scale_y_continuous(limits=c(0,max(biplot_data_exp3$current_vel)+0.001),breaks=seq(0,0.18,0.05),expand=c(0,0))+
  geom_point(data=data_full,aes(x=nitrate_mean,y=cur_vel_dep_mean,colour=group),size=1,alpha=0.3,pch=16)+
  scale_color_manual(name = "Scenario",values = c("#1B9E77", "#D95F02", "#E7298A"))+
  annotate("rect", xmin = min(predictors_present$nitrate_mean), xmax = max(predictors_present$nitrate_mean), ymin = min(predictors_present$cur_vel_dep_mean), ymax = max(predictors_present$cur_vel_dep_mean),
           alpha = 0,color = "black", linetype = "dashed")

nitr_curr_plot <- ggMarginal(nitr_curr_plot, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.3)

jpeg(".jpg",width=10,height=8,units="in",res=300)
plot_grid(temp_nitr_plot,temp_curr_plot,nitr_curr_plot)
dev.off()
