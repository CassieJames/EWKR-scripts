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


## --------------------------------------------
## Random Forests

# Aquatic count model

AQUATIC.mtry<-tuneRF(envs.std, Metrics$Aquatic, ntreeTry=1000, plot=TRUE)

AQUATIC.rf<-randomForest(Metrics$Aquatic~ ., data=envs.std, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1000, mtry=2, oob.prox=TRUE, 
      keep.forest=TRUE)

varImpPlot(AQUATIC.rf, sort=TRUE, type=1)
AQUATIC.rf


testdatav2=testdata
testdatav2=testdatav2[,c("Aquatic","d30", "d90", "d180","d365", "Flood_frequency", "TSLF")]
prediction <- predict(AQUATIC.rf, testdatav2)
qqnorm((prediction - testdatav2$Aquatic)/sd(prediction-testdatav2$Aquatic))
qqline((prediction - testdatav2$Aquatic)/sd(prediction-testdatav2$Aquatic))
RMSE.forest <- sqrt(mean((prediction-testdatav2$Aquatic)^2))
print(RMSE.forest/mean(testdatav2$Aquatic)) 



### ===============================================
# Checking for spatial autocorrelation in residuals

library(ape) 

resid_rf<- Metrics$Aquatic-predict(AQUATIC.rf,newdata=cbind(Metrics$Aquatic,envs.std) )
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
## Develop autocorrelation model not required because residuals were found not to be covary in space
