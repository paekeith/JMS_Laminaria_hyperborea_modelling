## This script is used to generate 3d partial dependence plots for the kelp park and forest models, using the full range of temperature and nitrate across the climate projections

#Steps required:
#generated the PD predictions for nitrate and temperature across the full ranges
# generate a grid of combinations for the two conditions, then extract the change grid cells with the closest match - will need to do this seperately for the IRE and UK data and find a way to combine without duplicates
#overlay the selected change cells as points colour coded by change category onto the 3d plot
rm(list=ls())
library(tidyverse)
library(randomForest) 
library(terra)
library(scales)

## Generating 3d kelp Park PD plots --------------------
#loading the kelp park model ensemble
model_list <- readRDS("Outputs/RF_models_park_GB.rds") #reading in the model list

## creating pd plot using current-day ranges ---------
#loading predictors
predictors_UK_pres <- readRDS("Data/Refined_data/Predictors_UK_present.rds")
predictors_IRE_pres <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")

predictors_pres <- rbind(predictors_UK_pres,predictors_IRE_pres)
predictors_pres <- predictors_pres %>%
  na.omit()

predictors_pres <- predictors_pres%>%
  select(temp_mean,nitrate_mean,cur_vel_dep_mean,bathymetry,wave_fetch,PAR_depth)

#need to run this for every model bootstrap and then get the average to plot
n.grid <- 75
predictor.df <- as.data.frame(matrix(nrow=n.grid,ncol=2))
names(predictor.df) <- c("temp_mean","nitrate_mean")

predictor.df$temp_mean <- seq(min(predictors_pres$temp_mean), max(predictors_pres$temp_mean), length.out = n.grid)

predictor.df$nitrate_mean <- seq(min(predictors_pres$nitrate_mean), max(predictors_pres$nitrate_mean), length.out = n.grid)

#creating an expanded version of the dataframe
df_combinations <- expand.grid(
  temp_mean = predictor.df$temp_mean,
  nitrate_mean = predictor.df$nitrate_mean
)
names(predictors_pres)
#adding in the other predictors at their averages
df_combinations$cur_vel_dep_mean <- mean(predictors_pres$cur_vel_dep_mean)
df_combinations$bathymetry <- mean(predictors_pres$bathymetry)
df_combinations$wave_fetch <- mean(predictors_pres$wave_fetch)
df_combinations$PAR_depth <- mean(predictors_pres$PAR_depth)

#predicting for all the models and then averaging
predictions <- as.data.frame(matrix(nrow=nrow(df_combinations),ncol=100))
for(i in 1:length(model_list)){
  predictions[,i] <- predict(model_list[[i]], newdata = df_combinations, type = "prob")[,2]
  print(i)
}

df_combinations$prob <- rowMeans(predictions)

#converting the predictions to wide format
wide_df <- df_combinations %>%
  select(temp_mean,nitrate_mean,prob)%>%
  pivot_wider(names_from = nitrate_mean, values_from = prob)

#remove the temperature column and the column names and convert to matrix
Z_matrix <- as.matrix(wide_df[, -1]) 
colnames(Z_matrix) <- NULL

X_temp_mean <- sort(unique(df_combinations$temp_mean))
Y_nitrate_mean <- sort(unique(df_combinations$nitrate_mean))

# generate colour map #
facet_dim <- n.grid - 1

# Define the color palette (adjust colors as desired)
ncolors <- n.grid
color_palette <- colorRampPalette(c("blue", "cyan", "yellow", "red"))(ncolors)

# Map Z-values (Probability) to the color palette indices
Z_flat <- as.vector(Z_matrix) 
color_index <- round(scales::rescale(Z_flat, 
                                     to = c(1, ncolors), 
                                     from = c(0, 1)))

z_min <- 0 
z_max <- 1

color_matrix_full <- matrix(color_palette[color_index], nrow = n.grid, ncol = n.grid)

# Select the N-1 x N-1 portion for the facets (required by persp)
facet_col <- color_matrix_full[1:facet_dim, 1:facet_dim]

# tiff(".tif",width=18,height=18,units = "in",res=300)
par(fig = c(0.02, 0.87, 0, 1), new = FALSE, mar = c(0, 4, 0, 0) + 0.1)
p_transform <- persp(x = X_temp_mean, y = Y_nitrate_mean, z = Z_matrix,
                     xlab = "Temperature (°C)",
                     ylab = "Nitrate conc (mmol.m-3)",
                     zlab = "Probability",
                     theta = 40, # Angle for viewing x-axis
                     phi = 30,   # Angle for viewing y-axis
                     d=5,
                     col = facet_col, # Color of the surface
                     shade = 0.3, # Shading for depth perception
                     ticktype = "detailed", # Detailed ticks on axes
                     nticks = 5, # Number of ticks on the axes
                     cex.lab = 3.5,
                     font.lab= 6,
                     cex.axis = 2.5,
                     font.axis= 6,
                     border="#00000033"
)

par(fig = c(0.85, 0.95, 0, 1), new = TRUE, mar = c(5, 0, 3, 4) + 0.1,family = "serif")

plot(1, 1, type = "n", 
     xlim = c(0, 1), ylim = c(0, n.grid), 
     xaxt = 'n', yaxt = 'n', ylab = "", xlab = "", bty = "n")

points(
  rep(0.5, n.grid), # X-coordinates centered at 0.5
  1:n.grid,         # Y-coordinates from 1 to N
  col = color_palette, 
  pch = 15,    # Solid squares
  cex = 10     # Make them large enough to fill the space
)

probability_labels <- round(seq(z_min, z_max, length.out = 5), 2)
tick_positions <- seq(0-1.2, n.grid+2.2, length.out = 5) # 5 evenly spaced positions from 1 to N

axis(side = 4, 
     at = tick_positions, 
     labels = probability_labels,
     cex.axis = 2.5,
     las = 1 
)

mtext("Probability", 
      side = 4, 
      line = 6, 
      cex = 3.5, 
      font = 1,
      las = 0) 

par(mfrow = c(1, 1), fig = c(0, 1, 0, 1), new = FALSE, mar = c(5, 4, 4, 2) + 0.1)

# dev.off()

## Generating 3d kelp forest PD plot --------------------
#loading the kelp forest model ensemble

model_list <- readRDS("Outputs/RF_models_forest_GB.rds") #reading in the model list

## creating pd plot using current-day ranges ---------
#loading predictors
predictors_UK_pres <- readRDS("Data/Refined_data/Predictors_UK_present.rds")
predictors_IRE_pres <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")

predictors_pres <- rbind(predictors_UK_pres,predictors_IRE_pres)
predictors_pres <- predictors_pres %>%
  na.omit()

predictors_pres <- predictors_pres%>%
  select(temp_mean,nitrate_mean,cur_vel_dep_mean,bathymetry,wave_fetch,PAR_depth)

#need to run this for every model bootstrap and then get the average to plot
#need to run this for every model bootstrap and then get the average to plot
n.grid <- 75
predictor.df <- as.data.frame(matrix(nrow=n.grid,ncol=2))
names(predictor.df) <- c("temp_mean","nitrate_mean")

predictor.df$temp_mean <- seq(min(predictors_pres$temp_mean), max(predictors_pres$temp_mean), length.out = n.grid)

predictor.df$nitrate_mean <- seq(min(predictors_pres$nitrate_mean), max(predictors_pres$nitrate_mean), length.out = n.grid)

#creating an expanded version of the dataframe
df_combinations <- expand.grid(
  temp_mean = predictor.df$temp_mean,
  nitrate_mean = predictor.df$nitrate_mean
)
names(predictors_pres)
#adding in the other predictors at their averages
df_combinations$cur_vel_dep_mean <- mean(predictors_pres$cur_vel_dep_mean)
df_combinations$bathymetry <- mean(predictors_pres$bathymetry)
df_combinations$wave_fetch <- mean(predictors_pres$wave_fetch)
df_combinations$PAR_depth <- mean(predictors_pres$PAR_depth)

#predicting for all the models and then averaging
predictions <- as.data.frame(matrix(nrow=nrow(df_combinations),ncol=100))
for(i in 1:length(model_list)){
  predictions[,i] <- predict(model_list[[i]], newdata = df_combinations, type = "prob")[,2]
  print(i)
}

df_combinations$prob <- rowMeans(predictions)

#converting the predictions to wide format
wide_df <- df_combinations %>%
  select(temp_mean,nitrate_mean,prob)%>%
  pivot_wider(names_from = nitrate_mean, values_from = prob)

#remove the temperature column and the column names and convert to matrix
Z_matrix <- as.matrix(wide_df[, -1]) 
colnames(Z_matrix) <- NULL

X_temp_mean <- sort(unique(df_combinations$temp_mean))
Y_nitrate_mean <- sort(unique(df_combinations$nitrate_mean))

# generate colour map #
facet_dim <- n.grid - 1

# Define the color palette (adjust colors as desired)
ncolors <- n.grid
color_palette <- colorRampPalette(c("blue", "cyan", "yellow", "red"))(ncolors)

# Map Z-values (Probability) to the color palette indices
Z_flat <- as.vector(Z_matrix) 
color_index <- round(scales::rescale(Z_flat, 
                                     to = c(1, ncolors), 
                                     from = c(0, 1)))

z_min <- 0 
z_max <- 1

color_matrix_full <- matrix(color_palette[color_index], nrow = n.grid, ncol = n.grid)

# Select the N-1 x N-1 portion for the facets (required by persp)
facet_col <- color_matrix_full[1:facet_dim, 1:facet_dim]

windows(record=T)

# tiff(".tif",width=18,height=18,units = "in",res=300)
par(fig = c(0.02, 0.87, 0, 1), new = FALSE, mar = c(0, 4, 0, 0) + 0.1)
p_transform <- persp(x = X_temp_mean, y = Y_nitrate_mean, z = Z_matrix,
                     xlab = "Temperature (°C)",
                     ylab = "Nitrate conc (mmol.m-3)",
                     zlab = "Probability",
                     theta = 40, # Angle for viewing x-axis
                     phi = 30,   # Angle for viewing y-axis
                     d=5,
                     col = facet_col, # Color of the surface
                     shade = 0.3, # Shading for depth perception
                     ticktype = "detailed", # Detailed ticks on axes
                     nticks = 5, # Number of ticks on the axes
                     cex.lab = 3.5,
                     font.lab= 6,
                     cex.axis = 2.5,
                     font.axis= 6,
                     border="#00000033"
)

par(fig = c(0.85, 0.95, 0, 1), new = TRUE, mar = c(5, 0, 3, 4) + 0.1,family = "serif")

plot(1, 1, type = "n", 
     xlim = c(0, 1), ylim = c(0, n.grid), 
     xaxt = 'n', yaxt = 'n', ylab = "", xlab = "", bty = "n")

points(
  rep(0.5, n.grid), # X-coordinates centered at 0.5
  1:n.grid,         # Y-coordinates from 1 to N
  col = color_palette, 
  pch = 15,    # Solid squares
  cex = 10     # Make them large enough to fill the space
)

probability_labels <- round(seq(z_min, z_max, length.out = 5), 2)
tick_positions <- seq(0-1.2, n.grid+2.2, length.out = 5) # 5 evenly spaced positions from 1 to N

axis(side = 4, 
     at = tick_positions, 
     labels = probability_labels,
     cex.axis = 2.5,
     las = 1 
)

# 5. Add a title/label to the legend
mtext("Probability", 
      side = 4, 
      line = 6, 
      cex = 3.5, 
      font = 1,
      las = 0) 

par(mfrow = c(1, 1), fig = c(0, 1, 0, 1), new = FALSE, mar = c(5, 4, 4, 2) + 0.1)

# dev.off()