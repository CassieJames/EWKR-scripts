##############################################################################
#### Model richness and abundance of wetland and dryland species 
library(corrplot)
library(mboost)
library(countreg)


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"

envdata=subset(envdata, envdata$Wet_Natives>0) # only model positive occurrences

# sort out environmental data - centre and log

envdata=envdata[!is.na(envdata$TSLW),] # Remove sites where extent of TSLW is beyond hydrology record

envdata$Inundated = as.factor(envdata$Inundated)
# highly skewed variables were log10 transformed before analysis
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
envdata$TSLW=as.numeric(scale(log10(envdata$TSLW+1),center=TRUE, scale=FALSE)) #This is a very strongly skewed variable 
envdata$MRI_length=as.numeric(scale(log10(envdata$MRI_length+1),center=TRUE, scale=FALSE)) #
envdata$MRI_meandepth=as.numeric(scale(log10(envdata$MRI_meandepth+1),center=TRUE, scale=FALSE))
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable
envdata$MeanTemp90=as.numeric(scale((envdata$MeanTemp90),center=TRUE, scale=FALSE))
envdata$Freq_d1=as.numeric(scale((envdata$Freq_d1),center=TRUE, scale=FALSE))
envdata$Freq_d3=as.numeric(scale((envdata$Freq_d3),center=TRUE, scale=FALSE))



#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Wet_Natives ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,center=TRUE, df=1, differences=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=1)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+
bols(MRI_meandepth, intercept=FALSE)+bbs(MRI_meandepth, center=TRUE, df=1)+
bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=1)+
bols(MRI_length, intercept=FALSE)+bbs(MRI_length, center=TRUE, df=1)+
bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE)

daten=envdata
# Create test and training data sets

n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

# seems to be a problem with random intercept so using weights to apply to model on full data to get predictions for OFB
weight0=rownames(testdata)
daten$weights<-1
daten$weights[which(rownames(daten) %in% weight0)]<-0


# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01), weights=daten$weights) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)

predicted<-list() # creates a whole bunch of empty lists
predicted.insample <-list()
nuvec<-list()
null.model.coef<-list()
null.model.nu<-list()


# Predictions for out-of-bootstrap data
predictions<-predict(Full.model,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
predictions_ofb <- predictions[which(rownames(predictions) %in% weight0),]
plot(predictions_ofb,daten$Wet_Natives[which(rownames(daten) %in% weight0)])
abline(0,1)


######################################################################
# Plots of fitted and predicted against observed

png(paste(image.dir,'Wet_Natives_fitted vs observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
plot(fitted(Full.model, type = "response"),daten$Wet_Natives, xlab='Fitted values', ylab='Observed values')
abline(0,1)
plot(predictions_ofb,daten$Wet_Natives[which(rownames(daten) %in% weight0)])
abline(0,1)
dev.off ()

################################################################################################################################################
# Model of wetland species counts for positive count results only 
#

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"

# sort out environmental data - centre and log
envdata=data.matrix.env.data
envdata$MRI_wet=envdata$MRI_shallow+envdata$MRI_deep
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond record

envdata=subset(envdata, envdata$Wet_Natives>0)
envdata$Inundated = as.factor(envdata$Inundated)
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
envdata$d90=as.numeric(scale(log10(envdata$d90+1),center=TRUE, scale=FALSE))
envdata$TSLW=as.numeric(scale((envdata$TSLW),center=TRUE, scale=FALSE))# this is a very strongly skewed variable 
envdata$d3Mon_wet=as.numeric(scale(log10(envdata$d3Mon_wet+1),center=TRUE, scale=FALSE)) #
envdata$d3Mon_shallow=as.numeric(scale(log10(envdata$d3Mon_shallow+1),center=TRUE, scale=FALSE)) # Highly skewed - log helps a bit but its heavily weighted to the left
envdata$d3Mon_deep=as.numeric(scale(log10(envdata$MRI_deep+1),center=TRUE, scale=FALSE))
envdata$MRI_wet=as.numeric(scale(log10(envdata$MRI_wet+1),center=TRUE, scale=FALSE)) #
envdata$MRI_shallow=as.numeric(scale(log10(envdata$MRI_shallow+1),center=TRUE, scale=FALSE)) # Highly skewed - log helps a bit but its heavily weighted to the left
envdata$MRI_deep=as.numeric(scale(log10(envdata$MRI_deep+1),center=TRUE, scale=FALSE))
envdata$d1yrs_wet=as.numeric(scale((envdata$d1yrs_wet),center=TRUE, scale=FALSE)) # quite strong peaks at low and high end of range
envdata$d1yrs_shallow=as.numeric(scale(log10(envdata$d1yrs_shallow+1),center=TRUE, scale=FALSE)) #this is a very strongly skewed variable as there are lots of low values
envdata$d1yrs_deep=as.numeric(scale(log10(envdata$d1yrs_deep+1),center=TRUE, scale=FALSE))
envdata$d3yrs_wet=as.numeric(scale((envdata$d3yrs_wet),center=TRUE, scale=FALSE))
envdata$d3yrs_shallow=as.numeric(scale(envdata$d3yrs_shallow,center=TRUE, scale=FALSE)) #this is a very strongly skewed variable as there are lots of low values
envdata$d3yrs_deep=as.numeric(scale(envdata$d3yrs_deep,center=TRUE, scale=FALSE))
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable
envdata$MaxTemp90=as.numeric(scale((envdata$MaxTemp90),center=TRUE, scale=FALSE))
envdata$MinTemp90=as.numeric(scale((envdata$MinTemp90),center=TRUE, scale=FALSE))
envdata$MeanTemp90=as.numeric(scale((envdata$MeanTemp90),center=TRUE, scale=FALSE))
envdata$Freq_d1=as.numeric(scale((envdata$Freq_d1),center=TRUE, scale=FALSE))
envdata$Freq_d3=as.numeric(scale((envdata$Freq_d3),center=TRUE, scale=FALSE))
envdata$Freq_d5=as.numeric(scale((envdata$Freq_d5),center=TRUE, scale=FALSE))

daten=envdata
# Create test and training data sets

n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

# seems to be a problem with random intercept so using weights to apply to model on full data to get predictions for OFB
weight0=rownames(testdata)
daten$weights<-1
daten$weights[which(rownames(daten) %in% weight0)]<-0


# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(), data = daten, control=boost_control(nu=0.01,mstop=200,trace=TRUE), weights=daten$weights) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)


# Predictions for out-of-bootstrap data
predictions<-predict(Full.model,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
predictions_ofb <- predictions[which(rownames(predictions) %in% weight0),]
plot(predictions_ofb,daten$Wet_Natives[which(rownames(daten) %in% weight0)])
abline(0,1)

# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)

######################################################################
# Plots of fitted and predicted against observed

png(paste(image.dir,'Wet_Natives_fitted vs observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
plot(fitted(Full.model, type = "response"),daten$Wet_Natives, xlab='Fitted values', ylab='Observed values')
abline(0,1)
plot(predictions_ofb,daten$Wet_Natives[which(rownames(daten) %in% weight0)])
abline(0,1)
dev.off ()

######################################################################
# Marginal plots

png(paste(image.dir,'Wet_Natives_positive_marginal_plots.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))

# plot using d90(accumulated rainfall in 90 days prior to sampling)
md90<-mean(log10(data.matrix.env.data$d90+1))
xmatLin <- extract(Full.model, which=5)
xmatSmooth <- extract(Full.model, which=6)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d90, df = 1, center = TRUE)`
plot(sort(daten$d90+md90),yvalues[order(daten$d90+md90)], type="l",xlab='log(d90+1)', ylab='f(log(d90+1))')
rug(sort(daten$d90+md90))


#  plot using Time Since Last Wet (TSLW)
mTSLF<-mean(log10(data.matrix.env.data$TSLW+1))
xmatSmooth <- extract(Full.model, which=12)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(TSLF, df = 1, center = TRUE)`
plot(sort(datenSmall$TSLF+mTSLF),yvalues[order(datenSmall$TSLF+mTSLF)], type="l",xlab='log(TSLF+1)', ylab='f(log(TSLF+1))')
rug(sort(datenSmall$TSLF+mTSLF))



################################################################################################################################################
# model of wetland species diversity for positive results only 

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"

# sort out environmental data - centre and log
envdata=data.matrix.env.data
envdata$MRI_wet=envdata$MRI_shallow+envdata$MRI_deep
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond record

envdata=subset(envdata, envdata$Wet_diversity>0)


formulaB <- Wet_diversity ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,center=TRUE, df=1, differences=1)+
bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=1)+
bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=1)+
bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=1)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+
bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=1)+
bols(d3Mon_shallow, intercept=FALSE)+bbs(d3Mon_shallow, center=TRUE, df=1)+
bols(d3Mon_deep, intercept=FALSE)+bbs(d3Mon_deep, center=TRUE, df=1)+
bols(MRI_wet, intercept=FALSE)+bbs(MRI_wet, center=TRUE, df=1)+
bols(MRI_shallow, intercept=FALSE)+bbs(MRI_shallow, center=TRUE, df=1)+
bols(MRI_deep, intercept=FALSE)+bbs(MRI_deep, center=TRUE, df=1)+
bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=1)+
bols(d1yrs_shallow, intercept=FALSE)+bbs(d1yrs_shallow, center=TRUE, df=1)+
bols(d1yrs_deep, intercept=FALSE)+bbs(d1yrs_deep, center=TRUE, df=1)+
bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=1)+
bols(d3yrs_shallow, intercept=FALSE)+bbs(d3yrs_shallow, center=TRUE, df=1)+
bols(d3yrs_deep, intercept=FALSE)+bbs(d3yrs_deep, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE)

daten=envdata
# Create test and training data sets

n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

# seems to be a problem with random intercept so using weights to apply to model on full data to get predictions for OFB
weight0=rownames(testdata)
daten$weights<-1
daten$weights[which(rownames(daten) %in% weight0)]<-0


# Run model on training set
Full.model <-gamboost(formulaB,family = Gaussian(), data = daten, control=boost_control(mstop=1000,trace=TRUE), weights=daten$weights) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)

predicted<-list() # creates a whole bunch of empty lists
predicted.insample <-list()
nuvec<-list()
null.model.coef<-list()
null.model.nu<-list()


# Predictions for out-of-bootstrap data
predictions<-predict(Full.model,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
predictions_ofb <- predictions[which(rownames(predictions) %in% weight0),]
plot(predictions_ofb,daten$Wet_diversity[which(rownames(daten) %in% weight0)])
abline(0,1)


######################################################################
# Plots of fitted and predicted against observed

png(paste(image.dir,'Wet_Natives_fitted vs observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
plot(fitted(Full.model, type = "response"),daten$Wet_Natives, xlab='Fitted values', ylab='Observed values')
abline(0,1)
plot(predictions_ofb,daten$Wet_Natives[which(rownames(daten) %in% weight0)])
abline(0,1)
dev.off ()

########################################################################
