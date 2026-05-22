################################
# Code for development of Random Forest models as part of Eskuche-Keith et al. "Abundance-based distribution modelling reveals climate-driven loss of functionally important kelp forest habitat"

# Due to memory allocation issues, the structure of this code is as follows:
# run bootstraps fitting models and calculating model fits, variable importance and plotting partial dependence. Save models to disk. Then independently load models and separately fit to each scenario (present day, SSP245, SSP585), generate average and SD rasters, save to disk, wiping memory between each section.
# Last edited 22/05/2026
################################

####==========    1. LOAD  PACKAGES  ================================####
rm(list = ls())
gc()
library(randomForest) 
library(terra)
library(ggplot2)
library(raster)

# library(tidyverse)

# Kelp park - Fitting models -----------------------------------------------------

#reading the full dataset of presences and background points
data <- as.data.frame(readRDS("Data/Refined_data/full_species_data_thinned.rds"))

data_park <- data
data_park$presence <- NA #setting up empty column for binary P/A

#now converting SACFORN to presence absence
data_park$presence[which(data_park$SACFORN_int>=3)] <- 1
data_park$presence[which(data_park$SACFORN_int<3)] <- 0

#setting the presence data to factor for the random forest
data_park$presence <- as.factor(data_park$presence)
table(data_park$presence)

#seperating the presences and absences for splitting train and test data
presences <- data_park[data_park$presence=="1",]
absences <- data_park[data_park$presence=="0",]

### -----------------------------------BOOTSTRAP------------------------------------###
n.boot <- 100
#setting up a data frame for storing the model evaluation parameters
mod_eval <- data.frame(matrix(nrow=4,ncol=n.boot))
row.names(mod_eval) <- c("AUC","TSS","Threshold_train","Threshold_test")

#setting up an object with the names of the predictors to be used in the model
imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))
plot_names <- sort(c("Temperature (°C)","Nitrate conc. (mmol.m-3)","Current velocity (m.s-1)","PAR at depth (E.m-2.day-1)","Bathymetry (m)","Wave fetch (log10 km*5)"))

influence_mat <- array(0,c(length(imp.var),n.boot)) # for saving information on predictors 
rownames(influence_mat) <- imp.var

# create environmental gradient files for each taxa
PredMins <- apply(data_park[,imp.var],2,min)
PredMaxs <- apply(data_park[,imp.var],2,max)

EnvRanges <- as.data.frame(array(0,c(100,length(imp.var))))
names(EnvRanges) <- imp.var

for (i in c(1:length(imp.var))) {
  EnvRanges[,i] <- seq(PredMins[imp.var[i]],PredMaxs[imp.var[i]],length = 100)
}

# create 3D array for saving env preds from boots
PD <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot))
dimnames(PD)[[2]] <- imp.var
dimnames(PD)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")

model_list <- list() #list for saving model object to
start_time <- Sys.time()

for (i in 1:n.boot){
  
  print(paste(paste("Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  rnd.train_P <- presences[sample(nrow(presences),nrow(presences),replace = T),] #   bootstrap sampling
  rnd.train_A <- absences[sample(nrow(absences),nrow(absences),replace = T),] 
  train_PA <- rbind(rnd.train_P,rnd.train_A)
  
  rnd.eval_P <- presences[!presences$cell %in% rnd.train_P$cell, ]
  rnd.eval_A <- absences[!absences$cell %in% rnd.train_A$cell, ]
  eval_PA <- rbind(rnd.eval_P,rnd.eval_A)
  
  model_PA <- tuneRF(x = train_PA[ ,imp.var],
                     y = train_PA[ , "presence"],
                     mtryStart = 2,
                     ntreeTry = 1000,
                     stepFactor = 2,
                     improve = 0.005,
                     trace = FALSE, plot = F, doBest = T)
  
  model_list[[i]] <- model_PA
  
  # Calculating model evaluation metrics
  train_pred <- predict(model_PA, newdata = train_PA, type = "prob")[,2]
  train_roc <- pROC::roc(response = train_PA$presence, predictor = train_pred, quiet = TRUE)
  train_coords <- pROC::coords(train_roc, x = "best",input = "threshold",best.method = "youden",transpose = FALSE)
  mod_eval["Threshold_train", i] <- train_coords["threshold"]
  
  eval_pred <- predict(model_PA, newdata = eval_PA, type = "prob")[,2]
  eval_roc <- pROC::roc(response = eval_PA$presence, predictor = eval_pred, quiet = TRUE)
  mod_eval["AUC",i] <- eval_roc$auc
  eval_coords <- pROC::coords(eval_roc, x = "best",input = "threshold",best.method = "youden",transpose = FALSE)
  mod_eval["TSS", i] <- eval_coords["sensitivity"] + eval_coords["specificity"] - 1 
  # Calculate AUC thresh
  mod_eval["Threshold_test", i] <- eval_coords["threshold"]
  
  # Environmental influence
  M1_contrib <- as.data.frame(importance(model_PA, type=2))
  M1_contrib$norm <- 0
  for (j in 1:nrow(M1_contrib)){M1_contrib$norm[j] <- (M1_contrib[j,1] * 100)/sum(M1_contrib[,1])}
  influence_mat[,i] <- M1_contrib[,2]
  
  # Generating partial dependence plots
  for (k in 1:length(imp.var)){
    # predict pd plot for each var
    temp.var <- as.data.frame(array(0,c(100,length(imp.var))))
    colnames(temp.var) <- imp.var
    temp.var[,k] <- seq(PredMins[imp.var[k]],PredMaxs[imp.var[k]],length = 100)
    for (m in 1:nrow(temp.var)){temp.var[m,-c(k)] <- apply(data_park[,imp.var[-k]],2,mean)}
    y <- predict(model_PA, newdata = temp.var, type = "prob")[,2]
    PD[,imp.var[[k]],i] <- loess(y ~ EnvRanges[,imp.var[[k]]])$y
  }
  
}
Sys.time()-start_time 

#save models
saveRDS(model_list,"Outputs/RF_models_park_GB.rds") 

# model fits
M.Fits <- round(apply(mod_eval, 1, function(x) c(Mean = mean(x), SD = sd(x))), 3)
write.csv(M.Fits, file =  "Outputs/GB_SDM/RF_fits_GB_park.csv")

# variable importance +- SD
pred_inf <-t(round(apply(influence_mat, 1, function(x) c(Mean = mean(x), SD = sd(x))), 1))
write.csv(pred_inf, file =  "Outputs/GB_SDM/RF_var_imp_GB_park.csv")

# plot PDs + 95 PI and save to file
svg(file = "Outputs/GB_SDM/RF_GB_park_pdplots.svg",width=15,height=10)
par(mar=c(4.2, 4, 1, 1))
par(mfrow=c(2,3))

for (i in c(1:length(imp.var))){
  plot(EnvRanges[,i],apply(PD[,i,],1,mean), col = "black",type='l',
       xlab = paste(plot_names[i], " (",pred_inf[imp.var[i],1], " ± ",pred_inf[imp.var[i],2],")", sep = ""), 
       ylab = '',
       ylim = c(min(PD[,,]), max((PD[,,]))),cex.axis=2,cex.lab=2) #max(boot_array_EnvTran[,,])))
  
  # 95% PI
  UC <- na.omit(cbind(EnvRanges[,imp.var[i]],
                      apply(PD[,i,],1, quantile, probs= c(0.05)),
                      apply(PD[,i,],1, quantile, probs= c(0.95))))
  polygon(c(UC[,1], rev(UC[,1])),c(UC[,2], rev(UC[,3])), col = rgb(0,0,0, 0.25), border = NA)
  rug(quantile(data_park[,imp.var[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
  # if (i == 1) title('Sample size - 1000')
}
dev.off()



# Kelp park - Predicting to current conditions ----------------------------
rm(list=ls())
gc()
model_list <- readRDS("Outputs/RF_models_park_GB.rds") #reading in the model list

#loading predictors
predictors_UK_present <- readRDS("Data/Refined_data/Predictors_UK_present.rds")
predictors_IRE_present <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")

#setting all depth values above 0m to 0, to deal with issues around chart datum etc.
predictors_UK_present$bathymetry[predictors_UK_present$bathymetry<0] <- 0
predictors_IRE_present$bathymetry[predictors_IRE_present$bathymetry<0] <- 0

n.boot <- length(model_list)

length.map.UK.pres <- nrow(predictors_UK_present)
length.map.IRE.pres <- nrow(predictors_IRE_present)

boot_matrix_UK_pres <- array(0, c(length.map.UK.pres,n.boot))
boot_matrix_IRE_pres <- array(0, c(length.map.IRE.pres,n.boot))

imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))

start_time <- Sys.time()
for (i in 1:n.boot){
  print(paste(paste("Predicting present, Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  model_PA <- model_list[[i]]  
  
  boot_matrix_UK_pres[,i] <- predict(model_PA, newdata = predictors_UK_present[,imp.var], type = "prob")[,2]
  boot_matrix_IRE_pres[,i] <- predict(model_PA, newdata = predictors_IRE_present[,imp.var], type = "prob")[,2]
  
}
Sys.time()-start_time 

#Saving mean and SD
boot.mean.pres.UK <-apply(boot_matrix_UK_pres,1,mean)
boot.sd.pres.UK <-apply(boot_matrix_UK_pres,1,sd)
saveRDS(boot.mean.pres.UK,"Outputs/RF_UK_preds_mean_park.rds")
saveRDS(boot.sd.pres.UK,"Outputs/RF_UK_preds_sd_park.rds")

boot.mean.pres.IRE <-apply(boot_matrix_IRE_pres,1,mean)
boot.sd.pres.IRE <-apply(boot_matrix_IRE_pres,1,sd)
saveRDS(boot.mean.pres.IRE,"Outputs/GB_SDM/RF_IRE_preds_mean_park.rds")
saveRDS(boot.sd.pres.IRE,"Outputs/GB_SDM/RF_IRE_preds_sd_park.rds")

# Kelp park - Predicting to SSP245 conditions ----------------------------
rm(list=ls())
gc()
model_list <- readRDS("Outputs/RF_models_park_GB.rds") #reading in the model list

#loading predictors
predictors_UK_sp245 <- readRDS("Data/Refined_data/Predictors_UK_ssp245.rds")
predictors_IRE_sp245 <- readRDS("Data/Refined_data/Predictors_IRE_ssp245.rds")

#setting all depth values above 0m to 0, to deal with issues around chart datum etc.
predictors_UK_sp245$bathymetry[predictors_UK_sp245$bathymetry<0] <- 0
predictors_IRE_sp245$bathymetry[predictors_IRE_sp245$bathymetry<0] <- 0

n.boot <- length(model_list)
length.map.UK.sp245 <- nrow(predictors_UK_sp245)
length.map.IRE.sp245 <- nrow(predictors_IRE_sp245)

boot_matrix_UK_sp245 <- array(0, c(length.map.UK.sp245,n.boot))
boot_matrix_IRE_sp245 <- array(0, c(length.map.IRE.sp245,n.boot))

imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))

start_time <- Sys.time()
for (i in 1:n.boot){
  print(paste(paste("Predicting ssp245, Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  model_PA <- model_list[[i]]  
  
  boot_matrix_UK_sp245[,i] <- predict(model_PA, newdata = predictors_UK_sp245[,imp.var], type = "prob")[,2]
  boot_matrix_IRE_sp245[,i] <- predict(model_PA, newdata = predictors_IRE_sp245[,imp.var], type = "prob")[,2]
  
}
Sys.time()-start_time 

#Saving mean and SD
boot.mean.245.UK <-apply(boot_matrix_UK_sp245,1,mean)
boot.sd.245.UK <-apply(boot_matrix_UK_sp245,1,sd)
saveRDS(boot.mean.245.UK,"Outputs/RF_UK_preds_mean_park_ssp245.rds")
saveRDS(boot.sd.245.UK,"Outputs/RF_UK_preds_sd_park_ssp245.rds")

#Saving mean and SD
boot.mean.245.IRE <-apply(boot_matrix_IRE_sp245,1,mean)
boot.sd.245.IRE <-apply(boot_matrix_IRE_sp245,1,sd)
saveRDS(boot.mean.245.IRE,"Outputs/RF_IRE_preds_mean_park_ssp245.rds")
saveRDS(boot.sd.245.IRE,"Outputs/RF_IRE_preds_sd_park_ssp245.rds")

# Kelp park - Predicting to SSP585 conditions ----------------------------
rm(list=ls())
gc()
model_list <- readRDS("Outputs/RF_models_park_GB.rds") #reading in the model list

#loading predictors
predictors_UK_sp585 <- readRDS("Data/Refined_data/Predictors_UK_ssp585.rds")
predictors_IRE_sp585 <- readRDS("Data/Refined_data/Predictors_IRE_ssp585.rds")

#setting all depth values above 0m to 0, to deal with issues around chart datum etc.
predictors_UK_sp585$bathymetry[predictors_UK_sp585$bathymetry<0] <- 0
predictors_IRE_sp585$bathymetry[predictors_IRE_sp585$bathymetry<0] <- 0

n.boot <- length(model_list)
length.map.UK.sp585 <- nrow(predictors_UK_sp585)
length.map.IRE.sp585 <- nrow(predictors_IRE_sp585)

boot_matrix_UK_sp585 <- array(0, c(length.map.UK.sp585,n.boot))
boot_matrix_IRE_sp585 <- array(0, c(length.map.IRE.sp585,n.boot))

imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))

start_time <- Sys.time()
for (i in 1:n.boot){
  print(paste(paste("Predicting ssp585, Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  model_PA <- model_list[[i]]  
  
  boot_matrix_UK_sp585[,i] <- predict(model_PA, newdata = predictors_UK_sp585[,imp.var], type = "prob")[,2]
  boot_matrix_IRE_sp585[,i] <- predict(model_PA, newdata = predictors_IRE_sp585[,imp.var], type = "prob")[,2]
  
}
Sys.time()-start_time 

#Saving mean and SD
boot.mean.585.UK <-apply(boot_matrix_UK_sp585,1,mean)
boot.sd.585.UK <-apply(boot_matrix_UK_sp585,1,sd)
saveRDS(boot.mean.585.UK,"Outputs/RF_UK_preds_mean_park_ssp585.rds")
saveRDS(boot.sd.585.UK,"Outputs/RF_UK_preds_sd_park_ssp585.rds")

boot.mean.585.IRE <-apply(boot_matrix_IRE_sp585,1,mean)
boot.sd.585.IRE <-apply(boot_matrix_IRE_sp585,1,sd)
saveRDS(boot.mean.585.IRE,"Outputs/RF_IRE_preds_mean_park_ssp585.rds")
saveRDS(boot.sd.585.IRE,"Outputs/RF_IRE_preds_sd_park_ssp585.rds")

# Kelp forest - Fitting models -----------------------------------------------------
rm(list=ls())
gc()
#reading the full dataset of presences and background points
data <- as.data.frame(readRDS("Data/Refined_data/full_species_data_thinned.rds"))

data_forest <- data
data_forest$presence <- NA #setting up empty column for binary P/A

#now converting SACFORN to presence absence
data_forest$presence[which(data_forest$SACFORN_int>=5)] <- 1
data_forest$presence[which(data_forest$SACFORN_int<5)] <- 0

#setting the presence data to factor for the random forest
data_forest$presence <- as.factor(data_forest$presence)
table(data_forest$presence)

#seperating the presences and absences for splitting train and test data
presences <- data_forest[data_forest$presence=="1",]
absences <- data_forest[data_forest$presence=="0",]

### -----------------------------------BOOTSTRAP------------------------------------###
n.boot <- 100
#setting up a data frame for storing the model evaluation parameters
mod_eval <- data.frame(matrix(nrow=4,ncol=n.boot))
row.names(mod_eval) <- c("AUC","TSS","Threshold_train","Threshold_test")

#setting up an object with the names of the predictors to be used in the model
imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))
plot_names <- sort(c("Temperature (°C)","Nitrate conc. (mmol.m-3)","Current velocity (m.s-1)","PAR at depth (E.m-2.day-1)","Bathymetry (m)","Wave fetch (log10 km*5)"))

influence_mat <- array(0,c(length(imp.var),n.boot)) # for saving information on predictors 
rownames(influence_mat) <- imp.var

# create environmental gradient files for each taxa
PredMins <- apply(data_forest[,imp.var],2,min)
PredMaxs <- apply(data_forest[,imp.var],2,max)

EnvRanges <- as.data.frame(array(0,c(100,length(imp.var))))
names(EnvRanges) <- imp.var

for (i in c(1:length(imp.var))) {
  EnvRanges[,i] <- seq(PredMins[imp.var[i]],PredMaxs[imp.var[i]],length = 100)
}

# create 3D array for saving env preds from boots
PD <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot))
dimnames(PD)[[2]] <- imp.var
dimnames(PD)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")

model_list <- list() #list for saving model object to
start_time <- Sys.time()

for (i in 1:n.boot){
  
  print(paste(paste("Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  rnd.train_P <- presences[sample(nrow(presences),nrow(presences),replace = T),] #   bootstrap sampling
  rnd.train_A <- absences[sample(nrow(absences),nrow(absences),replace = T),] 
  train_PA <- rbind(rnd.train_P,rnd.train_A)
  
  rnd.eval_P <- presences[!presences$cell %in% rnd.train_P$cell, ]
  rnd.eval_A <- absences[!absences$cell %in% rnd.train_A$cell, ]
  eval_PA <- rbind(rnd.eval_P,rnd.eval_A)
  
  model_PA <- tuneRF(x = train_PA[ ,imp.var],
                     y = train_PA[ , "presence"],
                     mtryStart = 2,
                     ntreeTry = 1000,
                     stepFactor = 2,
                     improve = 0.005,
                     trace = FALSE, plot = F, doBest = T)
  
  model_list[[i]] <- model_PA
  
  # Calculating model evaluation metrics
  train_pred <- predict(model_PA, newdata = train_PA, type = "prob")[,2]
  train_roc <- pROC::roc(response = train_PA$presence, predictor = train_pred, quiet = TRUE)
  train_coords <- pROC::coords(train_roc, x = "best",input = "threshold",best.method = "youden",transpose = FALSE)
  mod_eval["Threshold_train", i] <- train_coords["threshold"]
  
  eval_pred <- predict(model_PA, newdata = eval_PA, type = "prob")[,2]
  eval_roc <- pROC::roc(response = eval_PA$presence, predictor = eval_pred, quiet = TRUE)
  mod_eval["AUC",i] <- eval_roc$auc
  eval_coords <- pROC::coords(eval_roc, x = "best",input = "threshold",best.method = "youden",transpose = FALSE)
  mod_eval["TSS", i] <- eval_coords["sensitivity"] + eval_coords["specificity"] - 1 
  # Calculate AUC thresh
  mod_eval["Threshold_test", i] <- eval_coords["threshold"]
  
  # Environmental influence
  M1_contrib <- as.data.frame(importance(model_PA, type=2))
  M1_contrib$norm <- 0
  for (j in 1:nrow(M1_contrib)){M1_contrib$norm[j] <- (M1_contrib[j,1] * 100)/sum(M1_contrib[,1])}
  influence_mat[,i] <- M1_contrib[,2]
  
  # Generating partial dependence plots
  for (k in 1:length(imp.var)){
    # predict pd plot for each var
    temp.var <- as.data.frame(array(0,c(100,length(imp.var))))
    colnames(temp.var) <- imp.var
    temp.var[,k] <- seq(PredMins[imp.var[k]],PredMaxs[imp.var[k]],length = 100)
    for (m in 1:nrow(temp.var)){temp.var[m,-c(k)] <- apply(data_forest[,imp.var[-k]],2,mean)}
    y <- predict(model_PA, newdata = temp.var, type = "prob")[,2]
    PD[,imp.var[[k]],i] <- loess(y ~ EnvRanges[,imp.var[[k]]])$y
  }
  
}
Sys.time()-start_time 

#save models
saveRDS(model_list,"Outputs/RF_models_forest_GB.rds") 

# model fits
M.Fits <- round(apply(mod_eval, 1, function(x) c(Mean = mean(x), SD = sd(x))), 3)
write.csv(M.Fits, file =  "Outputs/RF_fits_GB_forest.csv")

# variable importance +- SD
pred_inf <-t(round(apply(influence_mat, 1, function(x) c(Mean = mean(x), SD = sd(x))), 1))
write.csv(pred_inf, file =  "Outputs/RF_var_imp_GB_forest.csv")

# plot PDs + 95 PI and save to file
svg(file = "Outputs/RF_GB_forest_pdplots.svg",width=15,height=10)
par(mar=c(4.2, 4, 1, 1))
par(mfrow=c(2,3))

for (i in c(1:length(imp.var))){
  plot(EnvRanges[,i],apply(PD[,i,],1,mean), col = "black",type='l',
       xlab = paste(plot_names[i], " (",pred_inf[imp.var[i],1], " ± ",pred_inf[imp.var[i],2],")", sep = ""), 
       ylab = '',
       ylim = c(min(PD[,,]), max((PD[,,]))),cex.axis=2,cex.lab=2) #max(boot_array_EnvTran[,,])))
  
  # 95% PI
  UC <- na.omit(cbind(EnvRanges[,imp.var[i]],
                      apply(PD[,i,],1, quantile, probs= c(0.05)),
                      apply(PD[,i,],1, quantile, probs= c(0.95))))
  polygon(c(UC[,1], rev(UC[,1])),c(UC[,2], rev(UC[,3])), col = rgb(0,0,0, 0.25), border = NA)
  rug(quantile(data_forest[,imp.var[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
  # if (i == 1) title('Sample size - 1000')
}
dev.off()



# Kelp forest - Predicting to current conditions ----------------------------
rm(list=ls())
gc()
model_list <- readRDS("Outputs/RF_models_forest_GB.rds") #reading in the model list

#loading predictors
predictors_UK_present <- readRDS("Data/Refined_data/Predictors_UK_present.rds")
predictors_IRE_present <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")

#setting all depth values above 0m to 0, to deal with issues around chart datum etc.
predictors_UK_present$bathymetry[predictors_UK_present$bathymetry<0] <- 0
predictors_IRE_present$bathymetry[predictors_IRE_present$bathymetry<0] <- 0

n.boot <- length(model_list)
length.map.UK.pres <- nrow(predictors_UK_present)
length.map.IRE.pres <- nrow(predictors_IRE_present)

boot_matrix_UK_pres <- array(0, c(length.map.UK.pres,n.boot))
boot_matrix_IRE_pres <- array(0, c(length.map.IRE.pres,n.boot))

imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))

start_time <- Sys.time()
for (i in 1:n.boot){
  print(paste(paste("Predicting present, Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  model_PA <- model_list[[i]]  
  
  boot_matrix_UK_pres[,i] <- predict(model_PA, newdata = predictors_UK_present[,imp.var], type = "prob")[,2]
  boot_matrix_IRE_pres[,i] <- predict(model_PA, newdata = predictors_IRE_present[,imp.var], type = "prob")[,2]
  
}
Sys.time()-start_time 

#Saving mean and SD
boot.mean.pres.UK <-apply(boot_matrix_UK_pres,1,mean)
boot.sd.pres.UK <-apply(boot_matrix_UK_pres,1,sd)
saveRDS(boot.mean.pres.UK,"Outputs/RF_UK_preds_mean_forest.rds")
saveRDS(boot.sd.pres.UK,"Outputs/RF_UK_preds_sd_forest.rds")

boot.mean.pres.IRE <-apply(boot_matrix_IRE_pres,1,mean)
boot.sd.pres.IRE <-apply(boot_matrix_IRE_pres,1,sd)
saveRDS(boot.mean.pres.IRE,"Outputs/RF_IRE_preds_mean_forest.rds")
saveRDS(boot.sd.pres.IRE,"Outputs/RF_IRE_preds_sd_forest.rds")

# Kelp forest - Predicting to SSP245 conditions ----------------------------
rm(list=ls())
gc()
model_list <- readRDS("Outputs/RF_models_forest_GB.rds") #reading in the model list

#loading predictors
predictors_UK_sp245 <- readRDS("Data/Refined_data/Predictors_UK_ssp245.rds")
predictors_IRE_sp245 <- readRDS("Data/Refined_data/Predictors_IRE_ssp245.rds")

#setting all depth values above 0m to 0, to deal with issues around chart datum etc.
predictors_UK_sp245$bathymetry[predictors_UK_sp245$bathymetry<0] <- 0
predictors_IRE_sp245$bathymetry[predictors_IRE_sp245$bathymetry<0] <- 0

n.boot <- length(model_list)
length.map.UK.sp245 <- nrow(predictors_UK_sp245)
length.map.IRE.sp245 <- nrow(predictors_IRE_sp245)

boot_matrix_UK_sp245 <- array(0, c(length.map.UK.sp245,n.boot))
boot_matrix_IRE_sp245 <- array(0, c(length.map.IRE.sp245,n.boot))

imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))

start_time <- Sys.time()
for (i in 1:n.boot){
  print(paste(paste("Predicting ssp245, Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  model_PA <- model_list[[i]]  
  
  boot_matrix_UK_sp245[,i] <- predict(model_PA, newdata = predictors_UK_sp245[,imp.var], type = "prob")[,2]
  boot_matrix_IRE_sp245[,i] <- predict(model_PA, newdata = predictors_IRE_sp245[,imp.var], type = "prob")[,2]
  
}
Sys.time()-start_time 

#Saving mean and SD
boot.mean.245.UK <-apply(boot_matrix_UK_sp245,1,mean)
boot.sd.245.UK <-apply(boot_matrix_UK_sp245,1,sd)
saveRDS(boot.mean.245.UK,"Outputs/RF_UK_preds_mean_forest_ssp245.rds")
saveRDS(boot.sd.245.UK,"Outputs/RF_UK_preds_sd_forest_ssp245.rds")

gc()
boot.mean.245.IRE <-apply(boot_matrix_IRE_sp245,1,mean)
boot.sd.245.IRE <-apply(boot_matrix_IRE_sp245,1,sd)
saveRDS(boot.mean.245.IRE,"Outputs/RF_IRE_preds_mean_forest_ssp245.rds")
saveRDS(boot.sd.245.IRE,"Outputs/RF_IRE_preds_sd_forest_ssp245.rds")

# Kelp forest - Predicting to SSP585 conditions ----------------------------
rm(list=ls())
gc()
model_list <- readRDS("Outputs/RF_models_forest_GB.rds") #reading in the model list

#loading predictors
predictors_UK_sp585 <- readRDS("Data/Refined_data/Predictors_UK_ssp585.rds")
predictors_IRE_sp585 <- readRDS("Data/Refined_data/Predictors_IRE_ssp585.rds")

#setting all depth values above 0m to 0, to deal with issues around chart datum etc.
predictors_UK_sp585$bathymetry[predictors_UK_sp585$bathymetry<0] <- 0
predictors_IRE_sp585$bathymetry[predictors_IRE_sp585$bathymetry<0] <- 0

n.boot <- length(model_list)
length.map.UK.sp585 <- nrow(predictors_UK_sp585)
length.map.IRE.sp585 <- nrow(predictors_IRE_sp585)

boot_matrix_UK_sp585 <- array(0, c(length.map.UK.sp585,n.boot))
boot_matrix_IRE_sp585 <- array(0, c(length.map.IRE.sp585,n.boot))

imp.var <- sort(c("temp_mean","nitrate_mean","cur_vel_dep_mean","PAR_depth","bathymetry","wave_fetch"))

start_time <- Sys.time()
for (i in 1:n.boot){
  print(paste(paste("Predicting ssp585, Run",i),paste("start:",format(Sys.time(), "%H:%M:%S"))))
  
  model_PA <- model_list[[i]]  
  
  boot_matrix_UK_sp585[,i] <- predict(model_PA, newdata = predictors_UK_sp585[,imp.var], type = "prob")[,2]
  boot_matrix_IRE_sp585[,i] <- predict(model_PA, newdata = predictors_IRE_sp585[,imp.var], type = "prob")[,2]
  
}
Sys.time()-start_time 

#Saving mean and SD
boot.mean.585.UK <-apply(boot_matrix_UK_sp585,1,mean)
boot.sd.585.UK <-apply(boot_matrix_UK_sp585,1,sd)
saveRDS(boot.mean.585.UK,"Outputs/RF_UK_preds_mean_forest_ssp585.rds")
saveRDS(boot.sd.585.UK,"Outputs/RF_UK_preds_sd_forest_ssp585.rds")

boot.mean.585.IRE <-apply(boot_matrix_IRE_sp585,1,mean)
boot.sd.585.IRE <-apply(boot_matrix_IRE_sp585,1,sd)
saveRDS(boot.mean.585.IRE,"Outputs/RF_IRE_preds_mean_forest_ssp585.rds")
saveRDS(boot.sd.585.IRE,"Outputs/RF_IRE_preds_sd_forest_ssp585.rds")

# Converting prediction dataframes to rasters ------------------------------------
# Present Day Kelp Park --------------------------------------
## UK Mean ##
rm(list=ls())
gc()
#loading predictors
predictors_UK_present <- readRDS("Data/Refined_data/Predictors_UK_present.rds")

#loading the predictions
boot.mean.pres.UK <- readRDS("Outputs/RF_UK_preds_mean_park.rds")

R.mean.pres.UK <- rasterFromXYZ(data.frame(x = predictors_UK_present[,"x"],y = predictors_UK_present[,"y"],z = boot.mean.pres.UK),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files

terra::writeRaster(R.mean.pres.UK, "Outputs/RF_UK_preds_mean_park.tif", overwrite=TRUE)

## UK SD ##
gc()

boot.sd.pres.UK <- readRDS("Outputs/RF_UK_preds_sd_park.rds")

R.sd.pres.UK <- rasterFromXYZ(data.frame(x = predictors_UK_present[,"x"],y = predictors_UK_present[,"y"],z = boot.sd.pres.UK),crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files

terra::writeRaster(R.sd.pres.UK, "Outputs/RF_UK_preds_sd_park.tif", overwrite=TRUE)

## IRE Mean ##
rm(list=ls())
gc()

predictors_IRE_present <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")

boot.mean.pres.IRE <- readRDS("Outputs/RF_IRE_preds_mean_park.rds")

R.mean.pres.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_present[,"x"],y = predictors_IRE_present[,"y"],z = boot.mean.pres.IRE),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files

terra::writeRaster(R.mean.pres.IRE, "Outputs/RF_IRE_preds_mean_park.tif", overwrite=TRUE)

## IRE SD ##
boot.sd.pres.IRE <- readRDS("Outputs/RF_IRE_preds_sd_park.rds")

R.sd.pres.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_present[,"x"],y = predictors_IRE_present[,"y"],z = boot.sd.pres.IRE),crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files

terra::writeRaster(R.sd.pres.IRE, "Outputs/RF_IRE_preds_sd_park.tif", overwrite=TRUE)

# SSP245 Kelp Park --------------------------------------
rm(list=ls())
gc()

## UK Mean ##
#loading predictors
predictors_UK_sp245 <- readRDS("Data/Refined_data/Predictors_UK_ssp245.rds")
#loading the predictions
boot.mean.245.UK <- readRDS("Outputs/RF_UK_preds_mean_park_ssp245.rds")

#creating raster
R.mean.245.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp245[,"x"], 
                                          y = predictors_UK_sp245[,"y"],
                                          z = boot.mean.245.UK),
                               crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
#writing to disk
terra::writeRaster(R.mean.245.UK, "Outputs/RF_UK_preds_mean_park_ssp245.tif", overwrite=TRUE)

## UK SD ##
gc()

boot.sd.245.UK <- readRDS("Outputs/RF_UK_preds_sd_park_ssp245.rds")

R.sd.245.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp245[,"x"], 
                                        y = predictors_UK_sp245[,"y"],
                                        z = boot.sd.245.UK),
                             crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files
terra::writeRaster(R.sd.245.UK, "Outputs/RF_UK_preds_sd_park_ssp245.tif", overwrite=TRUE)

## IRE Mean ##
rm(list=ls())
gc()
predictors_IRE_sp245 <- readRDS("Data/Refined_data/Predictors_IRE_ssp245.rds")

boot.mean.245.IRE <- readRDS("Outputs/RF_IRE_preds_mean_park_ssp245.rds")
R.mean.245.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp245[,"x"],
                                           y = predictors_IRE_sp245[,"y"],
                                           z = boot.mean.245.IRE),
                                crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
terra::writeRaster(R.mean.245.IRE, "Outputs/RF_IRE_preds_mean_park_ssp245.tif", overwrite=TRUE)

## IRE SD ##
gc()
boot.sd.245.IRE <- readRDS("Outputs/RF_IRE_preds_sd_park_ssp245.rds")
R.sd.245.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp245[,"x"], 
                                         y = predictors_IRE_sp245[,"y"],
                                         z = boot.sd.245.IRE),
                              crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files
terra::writeRaster(R.sd.245.IRE, "Outputs/RF_IRE_preds_sd_park_ssp245.tif", overwrite=TRUE)


# SSP585 Kelp Park --------------------------------------
## UK Mean ##
rm(list=ls())
gc()
predictors_UK_sp585 <- readRDS("Data/Refined_data/Predictors_UK_ssp585.rds")

boot.mean.585.UK <- readRDS("Outputs/RF_UK_preds_mean_park_ssp585.rds")
R.mean.585.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp585[,"x"], 
                                          y = predictors_UK_sp585[,"y"],
                                          z = boot.mean.585.UK),
                               crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
terra::writeRaster(R.mean.585.UK, "Outputs/RF_UK_preds_mean_park_ssp585.tif", overwrite=TRUE)

## UK SD ##
gc()
boot.sd.585.UK <- readRDS("Outputs/RF_UK_preds_sd_park_ssp585.rds")

R.sd.585.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp585[,"x"], 
                                        y = predictors_UK_sp585[,"y"],
                                        z = boot.sd.585.UK),
                             crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files

terra::writeRaster(R.sd.585.UK, "Outputs/RF_UK_preds_sd_park_ssp585.tif", overwrite=TRUE)

## IRE Mean ##
rm(list=ls())
gc()

predictors_IRE_sp585 <- readRDS("Data/Refined_data/Predictors_IRE_ssp585.rds")

boot.mean.585.IRE <- readRDS("Outputs/RF_IRE_preds_mean_park_ssp585.rds")
#converting to rasters and saving
R.mean.585.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp585[,"x"], 
                                           y = predictors_IRE_sp585[,"y"],
                                           z = boot.mean.585.IRE),
                                crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
terra::writeRaster(R.mean.585.IRE, "Outputs/RF_IRE_preds_mean_park_ssp585.tif", overwrite=TRUE)

## IRE SD ##
gc()
boot.sd.585.IRE <- readRDS("Outputs/RF_IRE_preds_sd_park_ssp585.rds")
R.sd.585.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp585[,"x"], 
                                         y = predictors_IRE_sp585[,"y"],
                                         z = boot.sd.585.IRE),
                              crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files
terra::writeRaster(R.sd.585.IRE, "Outputs/RF_IRE_preds_sd_park_ssp585.tif", overwrite=TRUE)

# Present Day Kelp forest --------------------------------------
## UK Mean ##
rm(list=ls())
gc()
#loading predictors
predictors_UK_present <- readRDS("Data/Refined_data/Predictors_UK_present.rds")

#loading the predictions
boot.mean.pres.UK <- readRDS("Outputs/RF_UK_preds_mean_forest.rds")

R.mean.pres.UK <- rasterFromXYZ(data.frame(x = predictors_UK_present[,"x"],y = predictors_UK_present[,"y"],z = boot.mean.pres.UK),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files

terra::writeRaster(R.mean.pres.UK, "Outputs/RF_UK_preds_mean_forest.tif", overwrite=TRUE)

## UK SD ##
gc()

boot.sd.pres.UK <- readRDS("Outputs/RF_UK_preds_sd_forest.rds")

R.sd.pres.UK <- rasterFromXYZ(data.frame(x = predictors_UK_present[,"x"],y = predictors_UK_present[,"y"],z = boot.sd.pres.UK),crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files

terra::writeRaster(R.sd.pres.UK, "Outputs/RF_UK_preds_sd_forest.tif", overwrite=TRUE)

## IRE Mean ##
rm(list=ls())
gc()

predictors_IRE_present <- readRDS("Data/Refined_data/Predictors_IRE_present.rds")

boot.mean.pres.IRE <- readRDS("Outputs/RF_IRE_preds_mean_forest.rds")

R.mean.pres.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_present[,"x"],y = predictors_IRE_present[,"y"],z = boot.mean.pres.IRE),crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files

terra::writeRaster(R.mean.pres.IRE, "Outputs/RF_IRE_preds_mean_forest.tif", overwrite=TRUE)

## IRE SD ##
boot.sd.pres.IRE <- readRDS("Outputs/RF_IRE_preds_sd_forest.rds")

R.sd.pres.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_present[,"x"],y = predictors_IRE_present[,"y"],z = boot.sd.pres.IRE),crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files

terra::writeRaster(R.sd.pres.IRE, "Outputs/RF_IRE_preds_sd_forest.tif", overwrite=TRUE)

# SSP245 Kelp forest --------------------------------------
rm(list=ls())
gc()

## UK Mean ##
#loading predictors
predictors_UK_sp245 <- readRDS("Data/Refined_data/Predictors_UK_ssp245.rds")
#loading the predictions
boot.mean.245.UK <- readRDS("Outputs/RF_UK_preds_mean_forest_ssp245.rds")

#creating raster
R.mean.245.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp245[,"x"], 
                                          y = predictors_UK_sp245[,"y"],
                                          z = boot.mean.245.UK),
                               crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
#writing to disk
terra::writeRaster(R.mean.245.UK, "Outputs/RF_UK_preds_mean_forest_ssp245.tif", overwrite=TRUE)

## UK SD ##
gc()

boot.sd.245.UK <- readRDS("Outputs/RF_UK_preds_sd_forest_ssp245.rds")

R.sd.245.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp245[,"x"], 
                                        y = predictors_UK_sp245[,"y"],
                                        z = boot.sd.245.UK),
                             crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files
terra::writeRaster(R.sd.245.UK, "Outputs/RF_UK_preds_sd_forest_ssp245.tif", overwrite=TRUE)

## IRE Mean ##
rm(list=ls())
gc()
predictors_IRE_sp245 <- readRDS("Data/Refined_data/Predictors_IRE_ssp245.rds")

boot.mean.245.IRE <- readRDS("Outputs/RF_IRE_preds_mean_forest_ssp245.rds")
R.mean.245.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp245[,"x"],
                                           y = predictors_IRE_sp245[,"y"],
                                           z = boot.mean.245.IRE),
                                crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
terra::writeRaster(R.mean.245.IRE, "Outputs/RF_IRE_preds_mean_forest_ssp245.tif", overwrite=TRUE)

## IRE SD ##
gc()
boot.sd.245.IRE <- readRDS("Outputs/RF_IRE_preds_sd_forest_ssp245.rds")
R.sd.245.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp245[,"x"], 
                                         y = predictors_IRE_sp245[,"y"],
                                         z = boot.sd.245.IRE),
                              crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files
terra::writeRaster(R.sd.245.IRE, "Outputs/RF_IRE_preds_sd_forest_ssp245.tif", overwrite=TRUE)


# SSP585 Kelp forest --------------------------------------
## UK Mean ##
rm(list=ls())
gc()
predictors_UK_sp585 <- readRDS("Data/Refined_data/Predictors_UK_ssp585.rds")

boot.mean.585.UK <- readRDS("Outputs/RF_UK_preds_mean_forest_ssp585.rds")
R.mean.585.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp585[,"x"], 
                                          y = predictors_UK_sp585[,"y"],
                                          z = boot.mean.585.UK),
                               crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
terra::writeRaster(R.mean.585.UK, "Outputs/RF_UK_preds_mean_forest_ssp585.tif", overwrite=TRUE)

## UK SD ##
gc()
boot.sd.585.UK <- readRDS("Outputs/RF_UK_preds_sd_forest_ssp585.rds")

R.sd.585.UK <- rasterFromXYZ(data.frame(x = predictors_UK_sp585[,"x"], 
                                        y = predictors_UK_sp585[,"y"],
                                        z = boot.sd.585.UK),
                             crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files

terra::writeRaster(R.sd.585.UK, "Outputs/RF_UK_preds_sd_forest_ssp585.tif", overwrite=TRUE)

## IRE Mean ##
rm(list=ls())
gc()

predictors_IRE_sp585 <- readRDS("Data/Refined_data/Predictors_IRE_ssp585.rds")

boot.mean.585.IRE <- readRDS("Outputs/RF_IRE_preds_mean_forest_ssp585.rds")
#converting to rasters and saving
R.mean.585.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp585[,"x"], 
                                           y = predictors_IRE_sp585[,"y"],
                                           z = boot.mean.585.IRE),
                                crs = crs("+init=epsg:27700")) # same proj as orginal env variables tiff files
terra::writeRaster(R.mean.585.IRE, "Outputs/RF_IRE_preds_mean_forest_ssp585.tif", overwrite=TRUE)

## IRE SD ##
gc()
boot.sd.585.IRE <- readRDS("Outputs/RF_IRE_preds_sd_forest_ssp585.rds")
R.sd.585.IRE <- rasterFromXYZ(data.frame(x = predictors_IRE_sp585[,"x"], 
                                         y = predictors_IRE_sp585[,"y"],
                                         z = boot.sd.585.IRE),
                              crs = crs("+init=epsg:27700"))  # same proj as orginal env variables tiff files
terra::writeRaster(R.sd.585.IRE, "Outputs/RF_IRE_preds_sd_forest_ssp585.tif", overwrite=TRUE)

#### ADDITIONAL STEPS ####
#Additional steps (not provided here) include clipping the output rasters to lowest astronomical tide and to areas of hard substrate, and thresholding them to identify areas where the habitats are likely to be present. The relevant layers for substrate and LAT are provided in the raw data folder, and the threshold values for both models are available in the manuscript and from the fits .csv files in the outputs folder