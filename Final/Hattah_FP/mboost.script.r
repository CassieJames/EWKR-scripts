# Script to undertake analysis of Hattah Floodplains vegetation data using mboost
# Adapted by C James from sample code provided in Maloney et al. 2012 (Applying additive modelling, Methods in Ecology and Evolution vol 3, 116-128)
# 
# 1st August 2016 
#
# Load data and libraries
library(mboost)
library(MASS)
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Spp_Env_matrix_HTH_FP.csv") # load data

# create diversity metric (just to work with a single response variable for the time being)

spp.matrix=data.matrix[,c(3:270)]
rownames(spp.matrix)=data.matrix$Row.names
H.index=(as.data.frame(diversity(spp.matrix, index="shannon")))
colnames(H.index)=c("H_index")
# sort out environmental data - centre and logarithmize
envdata=data.matrix[,c("d30", "d90", "d180", "d365", "Inundated.y","Flood_frequency","TSLF", "Easting", "Northing" )]
rownames(envdata)=(data.matrix$Row.names)
envdata$Inundated.y=as.numeric(envdata$Inundated.y)
envdata$Inundated.y[envdata$Inundated.y==15] <- "YES"

# highly skewed distributions were log10 transformed before analysis
envdata$d30=as.numeric(scale(log10(envdata$d30+1),center=TRUE, scale=FALSE)) # added one to avoid return of infinity values due to low mean rainfall over shorter time periods
envdata$d90=as.numeric(scale(log10(envdata$d90+1),center=TRUE, scale=FALSE))
envdata$d180=as.numeric(scale(log10(envdata$d180+1),center=TRUE, scale=FALSE))
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
envdata$Flood_frequency=as.numeric(scale(envdata$Flood_frequency,center=TRUE, scale=FALSE))
envdata$TSLF=as.numeric(scale(log10(envdata$TSLF+1),center=TRUE, scale=FALSE))
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
#envdata$Inundated.y=as.numeric(scale((envdata$Inundated.y),center=TRUE, scale=FALSE))

envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable

daten=merge(envdata,H.index,by="row.names") # merge diversity estimate with environmental predictors

n<- nrow(daten)
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:(n/2),] # create a subset of half the original data using the random numbers

#specify model formula - kept simple of the time being
formulaB <- H_index ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+
bols(d180, intercept=FALSE)+bbs(d180, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)
#+bols(Inundated.y, intercept=FALSE)+bbs(Inundated.y, centre=TRUE, df=1) # dropped for now as i am trying to work out how to incorporate factors

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

Full.model <-gamboost(formulaB,data=traindatav2, family=Gaussian(), control=boost_control(mstop=10000,trace=TRUE)) # might want to reduce down the mstop value for trial runs

# Carry out 5 fold cross validation to determine optimal stopping iteration
cv5f <- cv(model.weights(Full.model), type='kfold', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
st<-(mstop(cvm))
Full.model[st]

# example plot
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data
mFf<-mean(envdata$Flood_frequency)

xmatLin <- extract(Full.model, which=12)
xmatSmooth <- extract(Full.model, which=13)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(Flood_frequency, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(Flood_frequency, intercept = FALSE)`
plot(sort(traindata$Flood_frequency+mFf),yvalues[order(traindata$Flood_frequency+mFf)], type="l",xlab='Flood frequency', ylab='f(Flood frequency)')
rug(sort(traindata$Flood_frequency+mFf))


