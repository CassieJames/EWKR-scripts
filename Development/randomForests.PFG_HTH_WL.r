#---------------------------------------------------------------------------
# Programs to produce regression random forest models for vegetation metrics
#---------------------------------------------------------------------------

library(randomForest)
library(vegan)

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


## --------------------------------------------
## Random Forests

# Richness model

RICH.mtry<-tuneRF(envs.std, Metrics$Richness, ntreeTry=1000, plot=TRUE)

RICH.rf<-randomForest(Metrics$Richness ~ ., data=envs.std, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1000, mtry=4, oob.prox=TRUE, 
      keep.forest=TRUE)

varImpPlot(RICH.rf, sort=TRUE, type=1)
RICH.rf

par(mfrow=c(3,3))
par(mar=c(4,4,1,1))
partialPlot(RICH.rf,envs.std, TSLF, main=" ", ylab="Richness")
partialPlot(RICH.rf,envs.std, Flood_frequency, main=" ", ylab="Richness")
partialPlot(RICH.rf,envs.std, d30, main=" ", ylab="Richness")
partialPlot(RICH.rf,envs.std, d365, main=" ", ylab="Richness")

## ---------------------------------------------
## Develop autocorrelation model-  done for each time period separately

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
resid_rf<- datenSmall$Richness-predict(RICH.rf) # this stage needs to be checked as I am not sure we want the OOB generalised error
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

envs.std.auto = output[,c("d30","d90","d180","d365","Flood_frequency","TSLF", "focal_rac_vect_100")] # new predictor set with autocovariate

envs.std.auto =envs.std.auto [order(rownames(envs.std.auto)), ] # Make sure everything is in the same order
Metrics =Metrics [order(rownames(Metrics)), ]

## ---------------------------------------------
## rerun model with autocovariate included

# Richness model


RICH.mtry.auto<-tuneRF(envs.std.auto, Metrics$Richness, ntreeTry=1000, plot=TRUE)

RICH.rf.auto<-randomForest(Metrics$Richness ~ ., data=envs.std.auto, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1001, mtry=2, oob.prox=TRUE, 
      keep.inbag=TRUE)

varImpPlot(RICH.rf.auto, sort=TRUE, type=1)
RICH.rf.auto

testdatav2=testdata
testdatav2$focal_rac_vect_100 <-0  
testdatav2=testdatav2[,c("Richness","d30", "d90", "d180","d365", "Flood_frequency", "TSLF", "focal_rac_vect_100")]

prediction <- predict(RICH.rf.auto, testdatav2)
qqnorm((prediction - testdatav2$Richness)/sd(prediction-testdatav2$Richness))
qqline((prediction - testdatav2$Richness)/sd(prediction-testdatav2$Richness))
RMSE.forest <- sqrt(mean((prediction-testdatav2$Richness)^2))
print(RMSE.forest/mean(testdatav2$Richness))

# various plots

png(paste(image.dir,'Plots for richness random forest variable importance.png',sep=''), width=2000, height=2000, units="px", res=500)
par(mfrow=c(1,1))
par(mar=c(5,4,1,1))
varImpPlot(RICH.rf.auto, sort=TRUE, type=1, main="Species richness")
dev.off()


png(paste(image.dir,'Partial plots for richness random forest.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
partialPlot(RICH.rf.auto,envs.std.auto, TSLF, main=" ", ylab="Richness", xlab="TSLF (days)")
partialPlot(RICH.rf.auto,envs.std.auto, Flood_frequency, main=" ", ylab="Richness",xlab="Flood frequency")
partialPlot(RICH.rf.auto,envs.std.auto, d30, main=" ", ylab="Richness", xlab="d30 (mm)")
partialPlot(RICH.rf.auto,envs.std.auto, d365, main=" ", ylab="Richness", xlab="d365 (mm)")
dev.off()


plot_dat <- cbind(as.data.frame(Metrics$Richness), envs.std.auto)
colnames(plot_dat)[1]="Richness"
ff = forestFloor(RICH.rf.auto,envs.std.auto)
png(paste(image.dir,'Plots for richness random forest fit.png',sep=''), width=2000, height=2000, units="px", res=300)
plot(ff,1:7,col=fcol(ff))
dev.off()



png(paste(image.dir,'Plots for richness random forest fitted versus observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
tada=cbind( Metrics$Richness,envs.std.auto)
colnames(tada)[1]="Richness"
Fitted <- predict(RICH.rf.auto, tada)
plot(Fitted,tada$Richness, xlab="Fitted values", ylab="Observed values")
abline(0,1)
prediction <- predict(RICH.rf.auto, testdatav2)
plot(prediction,testdatav2$Richness,xlab="Predicted values", ylab="Observed values")
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

resid_rf.auto<- datenSmall$Richness-predict(RICH.rf.auto)
resid_rf.auto<-as.data.frame(resid_rf.auto)

daten.dists <- as.matrix(dist(cbind(datenSmall$Easting, datenSmall$Northing)))
daten.dists.inv <- 1/daten.dists
diag(daten.dists.inv) <- 0
daten.dists.inv[is.infinite(daten.dists.inv)] <- 0
Moran.I(resid_rf.auto$resid_rf.auto, daten.dists.inv)

# result suggests there is spatial autocorrelation in the residuals of the environment only model(p = 0.001389218, Morans I = 0.0396, sd=0.013364)
# for 100 resolution morans I is 0.0103644, sd =0.013392 and p = 0.3152148