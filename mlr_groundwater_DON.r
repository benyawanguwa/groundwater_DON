
## load the library 
library(raster)
library(sp)
library(randomForest)
library(magrittr)
library(rgdal)
library(gstat)
library(ggplot2)
library(mlr)
library(SemiPar)
library(Hmisc)
library(foreign)
library(maptools)
library(prettymapr)
library(mlrMBO)
library(parallelMap)
library(caret)
library(automap)
library(reshape2)
library(rminer)

## start the parallel 
parallelStartSocket(16)

## set the geo-reference 
WGS84 <- CRS("+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

## load the shapefiles
study_area <- shapefile("study_area.shp")
water <- shapefile("water.shp")

## load the veg, soil, land use
Soil <- raster("soil1.ovr")
Veg <- raster("vegetation.ovr")
Land_use <- raster("landuse1.ovr")

## define the function for preprocess 
study_area <- spTransform(study_area, WGS84)
extent <- c(study_area@bbox[1, 1:2], study_area@bbox[2, 1:2])

water <- spTransform(water, WGS84)

pre <- function(x) {
  projection(x) <- WGS84
  extent(x) <- extent
  x <- raster::mask(x, study_area)
  return(x)
}

read_points <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SP <- spTransform(SP, WGS84)
  SP@bbox <- study_area@bbox
  if (length(zerodist(SP)) >= 1) {
    SP <- SP[-(zerodist(SP)[, 1]),]
  }
  return(SP)
}

read_pointDataframes <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SPD <- SpatialPointsDataFrame(SP, read_data)
  SPD <- spTransform(SPD, WGS84)
  SPD@bbox <- study_area@bbox
  if (length(zerodist(SPD)) >= 1) {
    SPD <- SPD[-(zerodist(SPD)[, 1]),]
  }
  return(SPD)
}

reclass <- function(df, i, j) {
  df[, "DON"][df[, "DON"] <= i] <- "Low"
  df[, "DON"][df[, "DON"] < j] <- "Medium"
  df[, "DON"][(df[, "DON"] != "Low") & (df[, "DON"] != "Medium")] <- "High"
  df[, "DON"] <- factor(df[, "DON"], levels = c("Low", "Medium", "High"))
  return(df)
}

reclass4<-function(df,i,j){
  for (t in c(1,2)){
    df[, t][df[, t] <=i] <- "Low"
    df[, t][df[, t] < j] <- "Medium"
    df[, t][(df[, t] != "Low") & (df[, t] != "Medium")] <- "High"
    df[, t] <- factor(df[, t], levels = c("Low", "Medium", "High"))
  }
  return(df)
}

# Add X and Y to training 
add_S1S2 <- function(dataset) {
  dataset$s1 <- coordinates(dataset)[, 1]
  dataset$s2 <- coordinates(dataset)[, 2]
  return(dataset)
}

get_landscape<-function(df){
  landscape_all<-data.frame()
  for (ii in seq(1,length(df))){
    aa<-as.data.frame(df[[ii]])
    aa<-subset(aa,aa$Soil!="NA")
    Soil=tail(names(sort(table(aa[,1]))),1)
    Veg=tail(names(sort(table(aa[,2]))),1)
    Landuse=tail(names(sort(table(aa[,3]))),1)
    Catchment=tail(names(sort(table(aa[,4]))),1)
    GW_depth=mean(aa[,5])
    Distance=mean(aa[,6])
    Distance_GWC=mean(aa[,7])
    sing_land<-data.frame(Soil,Veg,Landuse,Catchment,GW_depth,Distance,Distance_GWC)
    landscape_all<-rbind(landscape_all,sing_land)
  }
  return(landscape_all)
}

## preprocess the landscape raster
Soil <- pre(Soil)
Veg <- pre(Veg)
Land_use <- pre(Land_use)

## combine similar vegetation and landuse types 
v_Veg<-values(Veg)
v_Veg[v_Veg %in% c(2,3,4)]=1
v_Veg[v_Veg %in% c(8,9)]=8
v_Veg[v_Veg %in% c(12,13)]=12
v_Veg[v_Veg %in% c(18,19,20)]=18
values(Veg)<-v_Veg

v_land<-values(Land_use)
v_land[v_land %in% c(1,2,5,6,7,11,12,13)]=1
v_land[v_land %in% c(3,4)]=3
v_land[v_land %in% c(8,10)]=8
values(Land_use)<-v_land

# Create an empty grid where n is the total number of cells
r <- raster(study_area)
res(r) <- res(Soil) 
base_grid <- as(r, 'SpatialGrid')
plot(base_grid)

## groundwater depth 
depth <- read.csv("depth_to_groundwater.csv",header=T) %>% read_pointDataframes(.)

## kriging interpotation for groundwater depth 
f_depth <- as.formula(sampling_d ~ 1)
depth<-add_S1S2(depth)
var.depth <- variogram(f_depth, depth)
dat.fit_depth <- fit.variogram(var.depth,vgm(c("Sph","Exp")))
depth_k <- krige(f_depth, depth, base_grid, dat.fit_depth) %>% raster(.) %>% raster::mask(., study_area)
depth_k@data@names<-"GW_depth"

## distance
water <- raster::rasterize(water, depth_k)
water_distance <- raster::mask(distance(water),study_area)
water_distance@data@names<-"Distance_to_water"

GW_center<-data.frame(Latitude=c(6495000,6475000,6460000,6448000,6403000),Longitude=rep(402000,5),values=1)
GW_center <- SpatialPoints(GW_center[, c(2:1)], proj4string = WGS84)
GW_center@bbox <- study_area@bbox
base_GWC<-water 
values(base_GWC)<-1
Distance_GWC<-distanceFromPoints(base_GWC,GW_center)
Distance_GWC@data@names<-"Distance_GWC"

## combine the landscapes  
landscapes<-stack(Soil,Veg,Land_use,depth_k,water_distance,Distance_GWC)
names(landscapes) <- c("Soil", "Veg", "Landuse","GW_depth", "Distance","Distance_GWC")

## load the nutrient data 
all_results<-data.frame()
all_data<-read.csv("groundwater_nutrient.csv",header = T)

## set the threshold for groundwater DON
a1=0.5
a2=2.5

## split the dataset 
set.seed(91)
trainIndex <- createDataPartition(all_data$DON, p = .75, list = FALSE, times = 1)

training<-all_data[ trainIndex,]
testing<-all_data[-trainIndex,]

## 10-fold cross-validation
rdesc = makeResampleDesc("CV", iters = 10,stratify = TRUE)
## Classification tree, set it up for predicting probabilities
classif.lrn = makeLearner("classif.randomForest", predict.type = "prob", fix.factors.prediction = TRUE)

## set the model parameters for random forest
para_rf = makeParamSet(
  makeIntegerParam("ntree", lower = 400, upper = 1000),
  makeIntegerParam("mtry", lower = 2, upper = 7),
  makeIntegerParam("nodesize", lower = 1, upper = 6)
)

## define the search stratgy
ctrl = makeTuneControlIrace(maxExperiments = 200L)

## function to build the model 
model_build <- function(dataset, n_target) {
  WP3_target = makeClassifTask(id = "WP3_target", data = dataset, target = n_target)
  rin = makeResampleInstance(rdesc, task = WP3_target)
  res_rf = mlr::tuneParams(classif.lrn, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                           show.info = FALSE)
  lrn_rf = setHyperPars(classif.lrn, par.vals = res_rf$x)
  ## train the final model 
  #set.seed(719)
  rf <- mlr::train(lrn_rf, WP3_target)
  return(rf)
}

  ## load the point data 
  training_df <- read_pointDataframes(training)
  testing_df <-  read_pointDataframes(testing) 
  
  training_points<- read_points(training)
  testing_points <- read_points(testing)
  
  ## map1, using kringing for DON interpolation
  # Add X and Y to training 
  training_df<-add_S1S2(training_df)
  testing_df<-add_S1S2(testing_df)
  
  ## define the areas to extract the landscape for points 
  a=100
  b=200
  capture_zone_land<-function(df){
    num<-nrow(df)
    landscape_data<-data.frame()
    for (r in seq(1,num)){
      p1_long<-df@coords[r,1]
      p1_lat<-df@coords[r,2]
      pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long+a,p1_lat+b),c(p1_long+2*a,p1_lat+b),
                           c(p1_long+2*a,p1_lat-b),c(p1_long+a,p1_lat-b),c(p1_long,p1_lat)))  
      projection(pg)<- WGS84
      p1_landscape<-raster::extract(landscapes,pg)
      p1_landscape<-get_landscape(p1_landscape)
      landscape_data<-rbind(landscape_data,p1_landscape)
    }
    return(landscape_data)
  }
  
  landscape_train <- capture_zone_land(training_df)
  landscape_test <- capture_zone_land(testing_df)
  
  ## combine the nutrient data and landscapes
  M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","Collect_Month","date_")])
  M2_test <- cbind(as.data.frame(landscape_test), testing_df@data[c("DON","Collect_Month","date_")])
  
  names(M2_train) <- colnames(M2_test)
  
  common_landscape<-function(land){
    land_dataset<-data.frame(table(M2_train[,land]))
    land_common<-subset(land_dataset,land_dataset[,2]==max(land_dataset[,2]))[1]
    return(as.matrix(land_common))
  }
  
  soil_max = common_landscape("Soil")[1]
  veg_max=common_landscape("Veg")[1]
  landuse_max = common_landscape("Landuse")[1]
  
  max_list<-list(soil_max,veg_max,landuse_max)
  
  for (ii in 1:3){
    M2_train[,ii]<-factor(M2_train[,ii],levels = unique(values(landscapes[[ii]]))[-1])
    M2_test[,ii]<-factor(M2_test[,ii],levels=unique(values(landscapes[[ii]]))[-1])
    M2_test [(which(!(M2_test[,ii] %in% M2_train[,ii]))),ii]<-as.numeric(max_list[[ii]])
     
    M2_train[,ii]<-droplevels(M2_train[,ii])
    M2_test[,ii]<-factor(M2_test[,ii],levels = levels(M2_train[,ii]))
    }
   ## reclassify the dataset 
   M2_train<-reclass(M2_train,a1,a2)
   M2_test<-reclass(M2_test,a1,a2)

## build the model   
  rf_DON<- model_build(WP2Train,"DON")

## test in testing set
test_rf = predict(rf_DON, newdata = WP2Test)
## ConfusionMatrix
print(calculateConfusionMatrix(test_rf))
## get the prediction performance
train_rf = predict(rf_DON, newdata = WP2Train)

acc_test<-performance(test_rf,measures=acc)[1]
acc_train<-performance(train_rf,measures=acc)[1]

sing_acc<-data.frame(acc_test,acc_train)
all_results<-rbind(all_results,sing_acc)


## apply the sensitivity analysis for rf_DON model  
predrf=function(M,data) # the PRED function
{ return (as.matrix(predict(M,newdata=data)$data[,2:4]))}

## 1D-SA
I=Importance(rf_DON_m2,WP2Train,method="1D-SA",PRED=predrf,outindex=7) 
print(round(I$imp,digits=2))
L=list(runs=1,sen=t(I$imp),sresponses=I$sresponses)
mgraph(L,graph="IMP",leg=names(WP2Train),col="gray",Grid=10)

## SOIL
vecplot(I,graph="VEC",xval=1,Grid=10,TC=1)
vecplot(I,graph="VEC",xval=1,Grid=10,TC=2)
vecplot(I,graph="VEC",xval=1,Grid=10,TC=3)

## VEG

vecplot(I,graph="VEC",xval=2,Grid=10,TC=1)
vecplot(I,graph="VEC",xval=2,Grid=10,TC=2)
vecplot(I,graph="VEC",xval=2,Grid=10,TC=3)

## LAND
vecplot(I,graph="VEC",xval=3,Grid=10,TC=1)
vecplot(I,graph="VEC",xval=3,Grid=10,TC=2)
vecplot(I,graph="VEC",xval=3,Grid=10,TC=3)

## DEPTH
vecplot(I,graph="VEC",xval=4,Grid=10,TC=1,pch=1)
vecplot(I,graph="VEC",xval=4,Grid=10,TC=2,pch=2)
vecplot(I,graph="VEC",xval=4,Grid=10,TC=3,pch=3)

##GWC 
vecplot(I,graph="VEC",xval=6,Grid=10,TC=1,pch=1)
vecplot(I,graph="VEC",xval=6,Grid=10,TC=2,pch=2)
vecplot(I,graph="VEC",xval=6,Grid=10,TC=3,pch=3)

## DATE
vecplot(I,graph="VEC",xval=9,Grid=10,TC=1,pch=1)
vecplot(I,graph="VEC",xval=9,Grid=10,TC=2,pch=2)
vecplot(I,graph="VEC",xval=9,Grid=10,TC=3,pch=3)

## variable interaction in DSA 
I2=Importance(rf_DON_m2,WP2Train,method="DSA",PRED=predrf,outindex=7)

cm2=agg_matrix_imp(I2)
fcm2=cmatrixplot(cm2,threshold=c(0.01,0.01))
vecplot(I2,graph="VECC",xval=c(4,6))
