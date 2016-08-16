# Script to undertake analysis of Hattah Floodplains vegetation data using boosted generalized additive models
# 
# Adapted by C James from sample code provided by Maloney et al. 2012 (Applying additive modelling, Methods in Ecology and Evolution vol 3, 116-128, Appendix E)
# 
# 1st August 2016 
#
# Load data and libraries
library(mboost)
library(MASS)
library(PerformanceAnalytics)
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Spp_Env_matrix_HTH_FP.csv") # load data

# create diversity metric (just to work with a single response variable for the time being)
spp.matrix=data.matrix[,c(3:270)]
rownames(spp.matrix)=data.matrix$Row.names
spp.matrix.pa=spp.matrix
spp.matrix.pa[ spp.matrix.pa> 0]<- 1 # make a copy of the spp matrix and turn into p/a data
Richness=as.data.frame(rowSums(spp.matrix.pa))
colnames(Richness)=c("Richness")

H.index=(as.data.frame(diversity(spp.matrix, index="shannon")))
colnames(H.index)=c("H_index")

# sort out environmental data - centre and logarithmize
envdata=data.matrix[,c("d30", "d90", "d180", "d365", "Inundated.y","Flood_frequency","TSLF", "Easting", "Northing" )]
envdata_backup=envdata
rownames(envdata)=(data.matrix$Row.names)
envdata$Inundated.y=as.numeric(envdata$Inundated.y) # for some reason categoric variable turned into numeric value of 15?!
envdata$Inundated.y[envdata$Inundated.y==15] <- "YES"#envdata$Inundated.y is a categoric variable which describes whether the site was inundated/damp at the time of the survey

# Note: rainfall variables are (predictably) highly correlated so I am dropping the most correlated (d90 and d180) and am left with d30 and d365 which still have correlation coeff of 0.88! 
# 

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

daten=merge(envdata,H.index,by="row.names") # merge diversity estimate with environmental predictors
daten.rich=merge(envdata,Richness,by="row.names") # merge diversity estimate with environmental predictors

####################################################DIVERSITY analysis

n<- nrow(daten)
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:(n/2),] # create a subset of half the original data using the random numbers

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- H_index ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)+
bols(Inundated.y, intercept=FALSE) 

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
traindatav2=traindata[,-1] # remove site names from data frame

# Run model

Full.model <-gamboost(formulaB,data=traindatav2, family=Gaussian(), control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, Northing and d30
mopt <- mstop(aic <- AIC(Full.model)) # also suggests that mstop is 10000 during initial run

# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='kfold', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]

# Example of marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data
mFf<-mean(envdata_backup$Flood_frequency)

xmatLin <- extract(Full.model, which=8)
xmatSmooth <- extract(Full.model, which=9)

par(mfrow = c(3,2))

# the below line had to adapted from the paper to correct for the type of speech marks, the order or the components and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(Flood_frequency, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(Flood_frequency, intercept = FALSE)`
plot(sort(traindata$Flood_frequency+mFf),yvalues[order(traindata$Flood_frequency+mFf)], type="l",xlab='Flood frequency', ylab='f(Flood frequency)')
rug(sort(traindata$Flood_frequency+mFf))

#  plot using Time Since Last Flood (TSLF)
mTSLF<-mean(log10(envdata_backup$TSLF+1))
xmatSmooth <- extract(Full.model, which=11)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(TSLF, df = 1, center = TRUE)`
plot(sort(traindata$TSLF+mTSLF),yvalues[order(traindata$TSLF+mTSLF)], type="l",xlab='TSLF', ylab='f(log(TSLF+1))')
rug(sort(traindata$TSLF+mTSLF))

# plot using d365(accumulated rainfall in 90 days prior to sampling)
md365<-mean(log10(envdata_backup$d365+1))
xmatSmooth <- extract(Full.model, which=7)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d365, df = 1, center = TRUE)` 
plot(sort(traindata$d365+md365),yvalues[order(traindata$d365+md365)], type="l",xlab='d365', ylab='f(log(d365+1))')
rug(sort(traindata$d365+md365))

# Predictions for out-of-bootstrap data

testdata <- datenSmall[-indvec,]
predictions<-predict(Full.model,newdata=testdata)
plot(exp(predictions),testdata$H_index)
abline(1,0)

# Compute pseudo r^2

m1 <-glm(H_index~1,data=traindata, family=gaussian)
null.model.coef <- coef(m1)

y<-testdata$H_indexlamba<-exp(predictions)
lambda<-exp(predictions)

L1<-y*log(lambda)-lgamma(y+1)-lambda
c0<-exp(null.model.coef)
L0<-y*log(c0)-lgamma(y+1)-c0

n<-length(y)
r2<-(1-exp(-2/n*(sum(L1)-sum(L0))))/(1-exp(sum(L0))^{2/n})


# Compute out-of-bootstrap r^2 - NOT SURE I have done this correctly - I want to measure the predictive accuracy but I am not sure what I am doing here!

m1 <-glm(H_index~1,data=testdata, family=gaussian)
null.model.coef <- coef(m1)

y<-testdata$H_indexlamba<-exp(predictions)
lambda<-exp(predictions)

L1<-y*log(lambda)-lgamma(y+1)-lambda
c0<-exp(null.model.coef)
L0<-y*log(c0)-lgamma(y+1)-c0

n<-length(y)
r2<-(1-exp(-2/n*(sum(L1)-sum(L0))))/(1-exp(sum(L0))^{2/n})

####################################################RICHNESS analysis
n<- nrow(daten.rich)
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall.rich <- daten.rich[indvecL,][1:(n/2),] # create a subset of half the original data using the random numbers
datenSmall.rich$Inundated.y=as.factor(datenSmall.rich$Inundated.y)

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Richness ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)+
bols(Inundated.y, intercept=FALSE) 

predicted<-list() # creates a whole bunch of empty lists
predicted.insample <-list()
nuvec<-list()
null.model.coef<-list()
null.model.nu<-list()

# specify training set

set.seed(1)
n<-nrow(datenSmall.rich)
indvec<-sample(1:n,n,replace=TRUE)
traindata.rich<-datenSmall.rich[indvec,]
traindatav2.rich=traindata.rich[,-1] # remove site names from training data set

# Run model - FAILS

Full.model.rich <-gamboost(formulaB,data=traindatav2.rich, family=Poisson(), control=boost_control(mstop=2000,trace=TRUE)) # 
# This model FAILS with an error around singularity which I think means that there are responses that are NOT unique combinations of the predictors - this did not happen for the diversity index
# because this index takes into count abundance as well as p/a
mopt <- mstop(aic <- AIC(Full.model.rich)) 



#########################Further plotting of model and partial effects
# Plot partial effects all together
# plot just smoothed effects
par(mfrow = c(3,2))
plot(Full.model, which = "bbs")
# plot just linear effects
par(mfrow = c(3,3))
plot(Full.model, which = "bols")
# Plot spatial effects
plot(Full.model, which = "bspatial")