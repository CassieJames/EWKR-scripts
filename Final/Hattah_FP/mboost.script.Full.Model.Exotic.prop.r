# Script to undertake analysis of Hattah Floodplains vegetation data using boosted generalized additive models
# Adapted by C James from sample code provided by Maloney et al. 2012 (Applying additive modelling, Methods in Ecology and Evolution vol 3, 116-128, Appendix E)
# 19th August 2016 

library(mboost)
library(gamboostLSS)
library(gamlss.dist)
library(MASS)

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Final_Metrics_Hattah_FP.csv") # load data with environmental variables attached
image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Plots/";

# sort out environmental data - centre and log
envdata=data.matrix.env.data[,c("d30", "d90", "d180", "d365", "Inundated","Flood_frequency","TSLF", "Easting", "Northing","Site.ID", "exotic")]
envdata$TSLF[is.na(envdata$TSLF)] <- 10000 # woops! still had NA values in this where the TSLF exceeded the available flow record so I have used an arbitrary 10000 for this 
envdata_backup=envdata
rownames(envdata)=(data.matrix$Row.names)
envdata$Inundated = as.factor(envdata$Inundated)
envdata$Row.names=rownames(envdata)

# Note: rainfall variables are (predictably) highly correlated so I am dropping the most correlated (d90 and d180) and am left with d30 and d365 are the least correlated
# highly skewed distributions were log10 transformed before analysis
envdata$d30=as.numeric(scale(log10(envdata$d30+1),center=TRUE, scale=FALSE)) # added one to avoid return of infinity values due to low mean rainfall over shorter time periods
envdata$d90=as.numeric(scale(log10(envdata$d90+1),center=TRUE, scale=FALSE))
envdata$d180=as.numeric(scale(log10(envdata$d180+1),center=TRUE, scale=FALSE))
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
envdata$Flood_frequency=as.numeric(scale(envdata$Flood_frequency,center=TRUE, scale=FALSE))
envdata$TSLF=as.numeric(scale(log10(envdata$TSLF+1),center=TRUE, scale=FALSE))
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable

####################################################    Total Abundance    #########################################################################################################
daten=envdata

n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

#specify model formula turning exotics into p/a- kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- exotic ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE) 

predicted<-list() # creates a whole bunch of empty lists
predicted.insample <-list()
nuvec<-list()
null.model.coef<-list()
null.model.nu<-list()

# specify training set

set.seed(1)
n<-nrow(datenSmall)
indvec<-sample(1:n,n,replace=TRUE)
traindata<-datenSmall[indvec,]

# Run Full model using zero inflatted model
Full.model <-gamboostLSS(formulaB,data=datenSmall, families=ZIPoLSS(), control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, Northing and d30

# Carry out  cross validation to determine optimal stopping iteration 

grid <- make.grid(max = c(mu = 1000, sigma = 1000), min = 20, length.out = 10, dense_mu_grid = FALSE)
plot(grid, pch = 20, cex = 0.2)
abline(0,1)
points(grid, pch = 20, col = "red")

densegrid <- make.grid(max = c(mu = 1000, sigma = 1000), min = 20, length.out = 10, dense_mu_grid = TRUE)
plot(densegrid, pch = 20, cex = 0.2)
abline(0,1)
points(grid, pch = 20, col = "red")

folds <- cv(model.weights(Full.model), type = "subsampling",B=5)
cvr <- cvrisk(Full.model, grid = densegrid, folds = folds)

plot(cvr)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)

# Predictions for out-of-bootstrap data
predictions_ofb<-predict(Full.model,newdata=datenSmall)

mu_predicts=predictions_ofb[[1]]
plot(exp(mu_predicts),datenSmall$exotic)
abline(0,1)

fitted_mu<-predict(Full.model,parameter="mu", newdata=datenSmall) # trying to extract fitted values for parameter mu
fitted_sigma<-predict(Full.model,parameter="sigma", newdata=datenSmall) # trying to extract fitted values for parameter sigma

# ZERO INFLATTED MODEL FOR POISSON WORKED BEST BUT MODEL STILL PERFORMED VERY poorly

##############################################################################################################################
# Second approach is to use just positive counts - checks revealed that data looked pretty good for a standard poisson so am using the mboost package as well


data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Final_Metrics_Hattah_FP.csv") # load data with environmental variables attached
image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Plots/";

# sort out environmental data - centre and log
envdata=data.matrix.env.data[,c("d30", "d90", "d180", "d365", "Inundated","Flood_frequency","TSLF", "Easting", "Northing","Site.ID", "exotic")]
envdata$TSLF[is.na(envdata$TSLF)] <- 10000 # woops! still had NA values in this where the TSLF exceeded the available flow record so I have used an arbitrary 10000 for this 
envdata_backup=envdata
rownames(envdata)=(data.matrix$Row.names)
envdata$Inundated = as.factor(envdata$Inundated)
envdata$Row.names=rownames(envdata)

# Note: rainfall variables are (predictably) highly correlated so I am dropping the most correlated (d90 and d180) and am left with d30 and d365 are the least correlated
# highly skewed distributions were log10 transformed before analysis
envdata$d30=as.numeric(scale(log10(envdata$d30+1),center=TRUE, scale=FALSE)) # added one to avoid return of infinity values due to low mean rainfall over shorter time periods
envdata$d90=as.numeric(scale(log10(envdata$d90+1),center=TRUE, scale=FALSE))
envdata$d180=as.numeric(scale(log10(envdata$d180+1),center=TRUE, scale=FALSE))
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
envdata$Flood_frequency=as.numeric(scale(envdata$Flood_frequency,center=TRUE, scale=FALSE))
envdata$TSLF=as.numeric(scale(log10(envdata$TSLF+1),center=TRUE, scale=FALSE))
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable

daten=envdata

datenExotics=daten[daten$exotic>0,]
n<- nrow(datenExotics) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenExoticsSmall <- datenExotics[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdataExotics <- datenExotics[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

formulaB <- exotic ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)

# Run model on training set
Full.model <-gamboost(formulaB,data=datenExoticsSmall, family=Poisson(), control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, Northing and d30

# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=25)
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


######################################################################
# Plots

png(paste(image.dir,'ExoticNo_fitted vs observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
plot(exp(fitted(Full.model)),datenExoticsSmall$exotic, xlab='Fitted values', ylab='Observed values')
abline(0,1)
plot(exp(predict(Full.model,newdata=testdataExotics)),testdataExotics$exotic, xlab='Predicted values', ylab='Observed values')
abline(0,1)
dev.off ()



png(paste(image.dir,'Exotic_poisson_marginal_plots.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

# plot using d30(accumulated rainfall in 30 days prior to sampling)
md30<-mean(log10(data.matrix.env.data$d30+1))
xmatLin <- extract(Full.model, which=5)
xmatSmooth <- extract(Full.model, which=6)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d30, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(d30, intercept = FALSE)`
plot(sort(datenExoticsSmall$d30+md30),yvalues[order(datenExoticsSmall$d30+md30)], type="l",xlab='log(d30+1)', ylab='f(log(d30+1))')
rug(sort(datenExoticsSmall$d30+md30))


# plot using d365(accumulated rainfall in 90 days prior to sampling)
md365<-mean(log10(data.matrix.env.data$d365+1))
xmatLin <- extract(Full.model, which=7)
xmatSmooth <- extract(Full.model, which=8)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d365, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(d365, intercept = FALSE)`
plot(sort(datenExoticsSmall$d365+md365),yvalues[order(datenExoticsSmall$d365+md365)], type="l",xlab='log(d365+1)', ylab='f(log(d365+1))')
rug(sort(datenExoticsSmall$d365+md365))


mFf<-mean(envdata_backup$Flood_frequency)
xmatSmooth <- extract(Full.model, which=10)
# the below line had to adapted from the paper to correct for the type of speech marks, the order or the components and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(Flood_frequency, df = 1, center = TRUE)`
plot(sort(datenExoticsSmall$Flood_frequency+mFf),yvalues[order(datenExoticsSmall$Flood_frequency+mFf)], type="l",xlab='Flood frequency', ylab='f(Flood frequency)')
rug(sort(datenExoticsSmall$Flood_frequency+mFf))

#  plot using Time Since Last Flood (TSLF)
mTSLF<-mean(log10(envdata_backup$TSLF+1))
xmatLin <- extract(Full.model, which=11)
xmatSmooth <- extract(Full.model, which=12)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(TSLF, df = 1, center = TRUE)`
plot(sort(datenExoticsSmall$TSLF+mTSLF),yvalues[order(datenExoticsSmall$TSLF+mTSLF)], type="l",xlab='log(TSLF+1)', ylab='f(log(TSLF+1))')
rug(sort(datenExoticsSmall$TSLF+mTSLF))


dev.off()


# Pretty terrible model!

# Predictions for out-of-bootstrap data
predictions<-predict(Full.model,newdata=testdataExotics)
plot(exp(predictions),datenExoticsSmall$exotic)
abline(0,1)

# Pseudo r2
m1 <- glm(exotic~1, data=datenExoticsSmall, family=poisson)
null.model.coef<-coef(m1)
y<-testdataExotics$exotic
lambda<-exp(predictions)
L1<-y*log(lambda)-lgamma(y+1)-lambda
c0<-exp(null.model.coef)
L0<-y*log(c0)-lgamma(y+1)-c0
n<-length(y)
r2<-(1-exp(-2/n*(sum(L1)-sum(L0))))/(1-exp(sum(L0))^{2/n})


 


