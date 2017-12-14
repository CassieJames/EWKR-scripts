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

# determine aquatic species richness

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
Spp_Env_matrix_HTH_FP=read.csv("Spp_Env_matrix_HTH_FP.csv") # read data
spp.matrix=Spp_Env_matrix_HTH_FP[,c(2:267)]
rownames(spp.matrix)=spp.matrix$Row.names
spp.matrix=spp.matrix[,-c(1)]
spp.matrix.pa=spp.matrix
spp.matrix.pa[ spp.matrix.pa> 0]<- 1 # make a copy of the spp matrix and turn into p/a data

mydata=data.frame(read.csv("HTH_FP.csv"))
mycodes=data.frame(read.csv("HL_FP_sp_codes.csv"))
mydata=merge(mydata,mycodes,by.x="Scientific_name", by.y="Scientific.name")

spp.matrix.pa=as.data.frame(spp.matrix.pa)
Outputt=t(spp.matrix.pa) # transpose matrix as my brain works better in this direction!
Outputt=as.data.frame(Outputt)
Outputt$Tdr <- NA # create empty columns
Outputt$Tda <- NA
Outputt$T <- NA
Outputt$A <- NA
Outputt$Atw <- NA
Outputt$Atl <- NA
Outputt$Ate <- NA
Outputt$Arp <- NA
Outputt$F <- NA
Outputt$S <- NA
specieslist=colnames(spp.matrix.pa)

FGs=c("T", "Tda", "Tdr", "A", "Arp","Atw", "Atl", "Ate", "Arp", "F", "S") # identify functional groups of interest - remove leaf litter and bare ground
for(spp in specieslist) { # loops through species list
Fung.group=unique(mydata[mydata$sp_code==spp,c("Functional.group")]) # for each species identify its functional group
if(Fung.group %in% FGs){Outputt[agrep(spp, rownames(Outputt)),match(Fung.group, colnames(Outputt))] <- 1} # is the functional group of the species is in FGs list then identify the species and match the correct functional group
} # agrep gives a lazy match - it seems to work where the exact match was falling over for a couple of species

T=colSums(Outputt[which(Outputt$T==1),]) # Sum all rows for which the functional group equals T and so on...
Tdr=colSums(Outputt[which(Outputt$Tdr==1),])
Tda=colSums(Outputt[which(Outputt$Tda==1),])
A=colSums(Outputt[which(Outputt$A==1),])
Atw=colSums(Outputt[which(Outputt$Atw==1),])
Atl=colSums(Outputt[which(Outputt$Atl==1),])
Ate=colSums(Outputt[which(Outputt$Ate==1),])
Arp=colSums(Outputt[which(Outputt$Arp==1),])
F=colSums(Outputt[which(Outputt$F==1),])
S=colSums(Outputt[which(Outputt$S==1),])


ttdata=rbind(T,Tdr) # lazy binding!
ttdata=rbind(ttdata,Tda)
ttdata=rbind(ttdata,A)
ttdata=rbind(ttdata,Atw)
ttdata=rbind(ttdata,Atl)
ttdata=rbind(ttdata,Ate)
ttdata=rbind(ttdata,Arp)
ttdata=rbind(ttdata,F)
ttdata=rbind(ttdata,S)

ttdata<-t(ttdata)
ttdata=ttdata[1:564,]
ttdata=as.data.frame(ttdata)
ttdata$Aquatic_Rich <-rowSums(ttdata[,4:10])
ttdata$Terrestrial_Rich <-rowSums(ttdata[,1:3])


envs.std<-envdata[,c("d30", "d90", "d180", "d365", "Flood_frequency","TSLF","Easting", "Northing")]
daten <-merge(envs.std,ttdata, by="row.names")


n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

rownames(datenSmall)=datenSmall$Row.names
envs.std<-datenSmall[,2:7]
Metrics<-datenSmall[,c("Aquatic_Rich", "Terrestrial_Rich")]


## --------------------------------------------
## Random Forests

# Richness model

AQUA_RICH.mtry<-tuneRF(envs.std, Metrics$Aquatic_Rich, ntreeTry=1000, plot=TRUE)

AQUA_RICH.rf<-randomForest(Metrics$Aquatic_Rich ~ ., data=envs.std, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1000, mtry=4, oob.prox=TRUE, 
      keep.forest=TRUE)

varImpPlot(AQUA_RICH.rf, sort=TRUE, type=1)
AQUA_RICH.rf

par(mfrow=c(3,3))
par(mar=c(4,4,1,1))
partialPlot(RICH.rf,envs.std, TSLF, main=" ", ylab="Richness")
partialPlot(RICH.rf,envs.std, Flood_frequency, main=" ", ylab="Richness")
partialPlot(RICH.rf,envs.std, d30, main=" ", ylab="Richness")
partialPlot(RICH.rf,envs.std, d365, main=" ", ylab="Richness")

### ===============================================
# Checking for spatial autocorrelation in residuals

library(ape) 

resid_rf<- Metrics$Aquatic_Rich-predict(AQUA_RICH.rf,newdata=cbind(Metrics$Aquatic_Rich,envs.std))
resid_rf=as.data.frame(resid_rf)
daten.dists <- as.matrix(dist(cbind(datenSmall$Easting, datenSmall$Northing)))
daten.dists.inv <- 1/daten.dists
diag(daten.dists.inv) <- 0
daten.dists.inv[is.infinite(daten.dists.inv)] <- 0
Moran.I(resid_rf$resid_rf, daten.dists.inv)
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
resid_rf<- datenSmall$Aquatic_Rich-predict(AQUA_RICH.rf, newdata=cbind(Metrics$Aquatic_Rich,envs.std)) # this stage needs to be checked as I am not sure we want the OOB generalised error
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


AQUA_RICH.mtry.auto<-tuneRF(envs.std.auto, Metrics$Aquatic_Rich, ntreeTry=1000, plot=TRUE)

AQUA_RICH.rf.auto<-randomForest(Metrics$Aquatic_Rich ~ ., data=envs.std.auto, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1001, mtry=2, oob.prox=TRUE, 
      keep.inbag=TRUE)

varImpPlot(RICH.rf.auto, sort=TRUE, type=1)
RICH.rf.auto

testdatav2=testdata
testdatav2$focal_rac_vect_100 <-0  
testdatav2=testdatav2[,c("Aquatic_Rich","d30", "d90", "d180","d365", "Flood_frequency", "TSLF", "focal_rac_vect_100")]

prediction <- predict(AQUA_RICH.rf.auto, testdatav2)
qqnorm((prediction - testdatav2$Aquatic_Rich)/sd(prediction-testdatav2$Aquatic_Rich))
qqline((prediction - testdatav2$Aquatic_Rich)/sd(prediction-testdatav2$Aquatic_Rich))
RMSE.forest <- sqrt(mean((prediction-testdatav2$Aquatic_Rich)^2))
print(RMSE.forest/mean(testdatav2$Aquatic_Rich))

# various plots

png(paste(image.dir,'Plots for aquatic richness random forest variable importance.png',sep=''), width=2000, height=2000, units="px", res=500)
par(mfrow=c(1,1))
par(mar=c(5,4,1,1))
varImpPlot(AQUA_RICH.rf.auto, sort=TRUE, type=1, main="Species richness")
dev.off()


png(paste(image.dir,'Partial plots for aquatic richness random forest.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
partialPlot(AQUA_RICH.rf.auto,envs.std.auto, TSLF, main=" ", ylab="Aquatic Richness", xlab="TSLF (days)")
partialPlot(AQUA_RICH.rf.auto,envs.std.auto, Flood_frequency, main=" ", ylab="Aquatic Richness",xlab="Flood frequency")
partialPlot(AQUA_RICH.rf.auto,envs.std.auto, d30, main=" ", ylab="Aquatic Richness", xlab="d30 (mm)")
partialPlot(AQUA_RICH.rf.auto,envs.std.auto, d365, main=" ", ylab="Aquatic Richness", xlab="d365 (mm)")
dev.off()


plot_dat <- cbind(as.data.frame(Metrics$Aquatic_Rich), envs.std.auto)
colnames(plot_dat)[1]="Aquatic_Rich"
ff = forestFloor(AQUA_RICH.rf.auto,envs.std.auto)
png(paste(image.dir,'Plots for aquatic richness random forest fit.png',sep=''), width=2000, height=2000, units="px", res=300)
plot(ff,1:7,col=fcol(ff))
dev.off()



png(paste(image.dir,'Plots for aquatic richness random forest fitted versus observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
tada=cbind( Metrics$Aquatic_Rich,envs.std.auto)
colnames(tada)[1]="Aquatic_Rich"
Fitted <- predict(AQUA_RICH.rf.auto, tada)
plot(Fitted,tada$Aquatic_Rich, xlab="Fitted values", ylab="Observed values")
abline(0,1)
prediction <- predict(AQUA_RICH.rf.auto, testdatav2)
plot(prediction,testdatav2$Aquatic_Rich,xlab="Predicted values", ylab="Observed values")
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
newdata=cbind(Metrics$Aquatic_Rich,envs.std.auto)
colnames(newdata)[1]="Aquatic_Rich"
resid_rf.auto<- datenSmall$Aquatic_Rich -predict(AQUA_RICH.rf.auto, newdata=newdata)
resid_rf.auto<-as.data.frame(resid_rf.auto)
daten.dists <- as.matrix(dist(cbind(datenSmall$Easting, datenSmall$Northing)))
daten.dists.inv <- 1/daten.dists
diag(daten.dists.inv) <- 0
daten.dists.inv[is.infinite(daten.dists.inv)] <- 0
Moran.I(resid_rf.auto$resid_rf.auto, daten.dists.inv)

# result suggests there is spatial autocorrelation in the residuals of the environment only model(p = 0.001389218, Morans I = 0.0396, sd=0.013364)
# for 100 resolution morans I is 0.0103644, sd =0.013392 and p = 0.3152148