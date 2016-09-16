#---------------------------------------------------------------------------
# Programs to produce regression random forest models for vegetation metrics
#---------------------------------------------------------------------------

library(randomForest)
library(vegan)
library(forestFloor)

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Final_Metrics_Hattah_FP_with rich.csv") # load data with corrections to diversity, richness and abundance
image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Plots/";

data.matrix=data.matrix[data.matrix$total_abund>0,]

envdata=data.matrix[,c("d30", "d90", "d180", "d365", "Inundated","Flood_frequency","TSLF", "Easting", "Northing","Site.ID")]
envdata$TSLF[is.na(envdata$TSLF)] <- 10000 # woops! still had NA values in this where the TSLF exceeded the available flow record so I have used an arbitrary 10000 for this 
envdata_backup=envdata
rownames(envdata)=(data.matrix$Row.names)
envdata$Inundated = as.factor(envdata$Inundated)
envdata$Row.names=rownames(envdata)

Metrics=data.matrix[,c("Richness","exotic","T", "Tdr", "Tda", "A", "Atw", "Atl", "Ate", "Arp", "F", "S")]
Metrics$Terrestrial=Metrics$T+Metrics$Tda+Metrics$Tdr
Metrics$Aquatic=Metrics$A+Metrics$Atw+Metrics$Atl+Metrics$Ate+Metrics$Arp+Metrics$F+Metrics$S
Metrics$Aquatic_prop =(Metrics$Aquatic/(Metrics$Terrestrial+Metrics$Aquatic))*100
Metrics$Terrestrial_prop =(Metrics$Terrestrial/(Metrics$Terrestrial+Metrics$Aquatic))*100
Metrics$Total_abund=Metrics$Terrestrial+Metrics$Aquatic

envs.std<-envdata[,c("d30", "d90", "d180", "d365", "Flood_frequency","TSLF")]
metrics.rng<-Metrics

daten <-cbind(envs.std, metrics.rng)
daten<-cbind(daten, envdata[,c("Easting", "Northing")])

n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

envs.std<-datenSmall[,1:6]
Metrics<-datenSmall[,c("Richness", "exotic","Terrestrial", "Aquatic", "Total_abund")]

### ===============================================
## Random Forests

# Exotic count model

Exotic.mtry<-tuneRF(envs.std, Metrics$exotic, ntreeTry=1000, plot=TRUE)

Exotic.rf<-randomForest(Metrics$exotic~ ., data=envs.std, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1000, mtry=4, oob.prox=TRUE, 
      keep.forest=TRUE)

varImpPlot(Exotic.rf, sort=TRUE, type=1)
Exotic.rf

testdatav2=testdata
testdatav2=testdatav2[,c("exotic","d30", "d90", "d180","d365", "Flood_frequency", "TSLF")]
prediction <- predict(Exotic.rf, testdatav2)
qqnorm((prediction - testdatav2$exotic)/sd(prediction-testdatav2$exotic))
qqline((prediction - testdatav2$exotic)/sd(prediction-testdatav2$exotic))
RMSE.forest <- sqrt(mean((prediction-testdatav2$exotic)^2))
print(RMSE.forest/mean(testdatav2$exotic)) 

### ===============================================
# Checking for spatial autocorrelation in residuals

library(ape) 

resid_rf<- Metrics$exotic-predict(Exotic.rf,newdata=cbind(Metrics$exotic,envs.std))
resid_rf<- Metrics$exotic-predict(Exotic.rf)
resid_rf=as.data.frame(resid_rf)
daten.dists <- as.matrix(dist(cbind(datenSmall$Easting, datenSmall$Northing)))
daten.dists.inv <- 1/daten.dists
diag(daten.dists.inv) <- 0
daten.dists.inv[is.infinite(daten.dists.inv)] <- 0
Moran.I(resid_rf$resid_rf, daten.dists.inv)

### ===============================================
# Various plots

png(paste(image.dir,'Plots for aquatic count random forest fitted versus observed no autocovar.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
Fitted <- predict(AQUATIC.rf,newdata=cbind(Metrics$Aquatic,envs.std.auto))
RMSE.forest <- sqrt(mean((Fitted-Metrics$Aquatic)^2))
print(RMSE.forest/mean(Metrics$Aquatic)) 
plot(Fitted,Metrics$Aquatic, xlab="Fitted values", ylab="Observed values")
abline(0,1)
testdatav2=testdata
testdatav2=testdatav2[,c("Aquatic","d30", "d90", "d180","d365", "Flood_frequency", "TSLF")]
prediction <- predict(AQUATIC.rf, testdatav2)
#qqnorm((prediction - testdatav2$Aquatic)/sd(prediction-testdatav2$Aquatic))
RMSE.forest <- sqrt(mean((prediction-testdatav2$Aquatic)^2))
print(RMSE.forest/mean(testdatav2$Aquatic)) 
plot(prediction,testdatav2$Aquatic,xlab="Predicted values", ylab="Observed values")
abline(0,1)
dev.off()


png(paste(image.dir,'Plots for aquatic counts random forest variable importance.png',sep=''), width=2000, height=2000, units="px", res=500)
par(mfrow=c(1,1))
par(mar=c(5,4,1,1))
varImpPlot(AQUATIC.rf, sort=TRUE, type=1, main="Aquatic counts")
dev.off()


png(paste(image.dir,'Partial plots for aquatic count random forest.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
partialPlot(AQUATIC.rf,envs.std, TSLF, main=" ", ylab="Richness", xlab="TSLF (days)")
partialPlot(AQUATIC.rf,envs.std, Flood_frequency, main=" ", ylab="Richness",xlab="Flood frequency")
partialPlot(AQUATIC.rf,envs.std, d30, main=" ", ylab="Richness", xlab="d30 (mm)")
partialPlot(AQUATIC.rf,envs.std, d365, main=" ", ylab="Richness", xlab="d365 (mm)")
dev.off()


## ---------------------------------------------
## Develop autocorrelation model


library(raster)

yoi=c("_08", "_09", "_10", "_11", "_12", "_13","_14","_16")

for(y in yoi){

envdata_sub<-datenSmall[grepl(y, rownames(datenSmall)), ]
xy<-cbind(envdata_sub$Easting, envdata_sub$Northing)
colnames(xy)=c("Easting", "Northing")
xy=as.data.frame(xy)
grid_east=(max(xy$Easting))-(min(xy$Easting))
grid_north=(max(xy$Northing))-(min(xy$Northing))
rast_rf <-raster(ncol=grid_east+1000, nrow = grid_north+1000, ymn = (min(xy$Northing))-500, xmn = (min(xy$Easting))-500,xmx= (max(xy$Easting))+500, ymx=(max(xy$Northing))+500)
res(rast_rf) <-100
resid_rf<- datenSmall$exotic-predict(Exotic.rf,newdata=cbind(Metrics$exotic,envs.std)) # which residuals should i be using
resid_rf<-as.data.frame(resid_rf)
resid_rf_sub<-resid_rf[grepl(y, rownames(resid_rf)), ]
xy_res_rf <-cbind(xy,(resid_rf_sub))
rast_rf[cellFromXY(rast_rf, xy_res_rf)] <-xy_res_rf[,3]
w=matrix(1,nrow=3,ncol=3)
focal_rac_rast_rf_100 <-focal(rast_rf, w=w, fun = mean, na.rm = TRUE) 
focal_rac_vect_100 <-extract(focal_rac_rast_rf_100, xy)
tada=cbind(envdata_sub,focal_rac_vect_100)
if(y=="_08") {output=tada
} else {output=rbind(output,tada)
}
}

envs.std.auto = output[,c("d30","d90","d180","d365","Flood_frequency","TSLF", "focal_rac_vect_100")]

envs.std.auto =envs.std.auto [order(rownames(envs.std.auto)), ]
Metrics =Metrics [order(rownames(Metrics)), ]
## ---------------------------------------------
## rerun model with autocovariate included

# Exotic model


Exotic.mtry.auto<-tuneRF(envs.std.auto, Metrics$exotic, ntreeTry=1000, plot=TRUE)

Exotic.rf.auto<-randomForest(Metrics$exotic ~ ., data=envs.std.auto, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1001, mtry=2, oob.prox=TRUE, 
      keep.inbag=TRUE)

varImpPlot(Exotic.rf.auto, sort=TRUE, type=1)
Exotic.rf.auto


# various plots

png(paste(image.dir,'Plots for exotic random forest variable importance.png',sep=''), width=2000, height=2000, units="px", res=500)
par(mfrow=c(1,1))
par(mar=c(5,4,1,1))
varImpPlot(Exotic.rf.auto, sort=TRUE, type=1, main="Species richness")
dev.off()


png(paste(image.dir,'Partial plots for Exotic random forest.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
partialPlot(Exotic.rf.auto,envs.std.auto, TSLF, main=" ", ylab="Exotic count", xlab="TSLF (days)")
partialPlot(Exotic.rf.auto,envs.std.auto, Flood_frequency, main=" ", ylab="Exotic count",xlab="Flood frequency")
partialPlot(Exotic.rf.auto,envs.std.auto, d30, main=" ", ylab="Exotic count", xlab="d30 (mm)")
partialPlot(Exotic.rf.auto,envs.std.auto, d365, main=" ", ylab="Exotic count", xlab="d365 (mm)")
dev.off()


plot_dat <- cbind(as.data.frame(Metrics$exotic), envs.std.auto)
colnames(plot_dat)[1]="exotic"
ff = forestFloor(Exotic.rf.auto,envs.std.auto)

png(paste(image.dir,'Plots for exotic random forest fit.png',sep=''), width=2000, height=2000, units="px", res=300)
plot(ff,1:7,col=fcol(ff))
dev.off()



testdatav2=testdata
testdatav2$focal_rac_vect_100 <-0
testdatav2=testdatav2[,c("exotic","d30", "d90", "d180","d365", "Flood_frequency", "TSLF", "focal_rac_vect_100")]

prediction <- predict(Exotic.rf.auto, testdatav2)
qqnorm((prediction - testdatav2$exotic)/sd(prediction-testdatav2$exotic))
qqline((prediction - testdatav2$exotic)/sd(prediction-testdatav2$exotic))
RMSE.forest <- sqrt(mean((prediction-testdatav2$exotic)^2))
print(RMSE.forest/mean(testdatav2$exotic))


png(paste(image.dir,'Plots for exotic random forest fitted versus observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
tada=cbind(Metrics$exotic,envs.std.auto)
colnames(tada)[1]="exotic"
Fitted <- predict(Exotic.rf.auto, tada)
plot(Fitted,tada$exotic, xlab="Fitted values", ylab="Observed values")
abline(0,1)
prediction <- predict(Exotic.rf.auto, testdatav2)
plot(prediction,testdatav2$exotic,xlab="Predicted values", ylab="Observed values")
abline(0,1)
dev.off()

# Plots of different focal resolutions

png(paste(image.dir,'Plots for richness showing different focal distances.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,3))
plot(focal_rac_rast_rf_5)
plot(focal_rac_rast_rf_10)
plot(focal_rac_rast_rf_100)
plot(focal_rac_rast_rf_200)
plot(focal_rac_rast_rf_500)
plot(focal_rac_rast_rf_1000)
dev.off()


### ===============================================
# Checking for spatial autocorrelation in residuals

resid_rf<- Metrics$exotic-predict(Exotic.rf.auto,newdata=cbind(Metrics$exotic,envs.std.auto) )
resid_rf=as.data.frame(resid_rf)
daten.dists <- as.matrix(dist(cbind(datenSmall$Easting, datenSmall$Northing)))
daten.dists.inv <- 1/daten.dists
diag(daten.dists.inv) <- 0
daten.dists.inv[is.infinite(daten.dists.inv)] <- 0
Moran.I(resid_rf$resid_rf, daten.dists.inv)


