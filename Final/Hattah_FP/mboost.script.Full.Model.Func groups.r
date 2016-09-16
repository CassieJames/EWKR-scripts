# Script to undertake analysis of Hattah Floodplains vegetation data using boosted generalized additive models
# Adapted by C James from sample code provided by Maloney et al. 2012 (Applying additive modelling, Methods in Ecology and Evolution vol 3, 116-128, Appendix E)
# 19th August 2016 


library(mboost)
library(MASS)
library(PerformanceAnalytics)

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Final_Metrics_Hattah_FP_with rich.csv") # load data with corrections to diversity, richness and abundance
image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Plots/";

# sort out environmental data - centre and log
envdata=data.matrix[,c("d30", "d90", "d180", "d365", "Inundated","Flood_frequency","TSLF", "Easting", "Northing","Site.ID", "T", "Tdr", "Tda", "A", "Atw", "Atl", "Ate", "Arp", "F", "S")]
envdata$TSLF[is.na(envdata$TSLF)] <- 10000 # woops! still had NA values in this where the TSLF exceeded the available flow record so I have used an arbitrary 10000 for this 
envdata_backup=envdata
rownames(envdata)=(data.matrix$Row.names)
envdata$Inundated = as.factor(envdata$Inundated)
envdata$Row.names=rownames(envdata)

envdata$Terrestrial=envdata$T+envdata$Tda+envdata$Tdr
envdata$Aquatic=envdata$A+envdata$Atw+envdata$Atl+envdata$Ate+envdata$Arp+envdata$F+envdata$S
envdata$Aquatic_prop =(envdata$Aquatic/(envdata$Terrestrial+envdata$Aquatic))*100
envdata$Terrestrial_prop =(envdata$Terrestrial/(envdata$Terrestrial+envdata$Aquatic))*100
envdata$Total_abund=envdata$Terrestrial+envdata$Aquatic
metrics = envdata[,c()]

# Note: rainfall variables are (predictably) highly correlated so I am dropping the most correlated (d90 and d180) and am left with d30 and d365 which still have correlation coeff of 0.88! 
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

####################################################    set up models   #########################################################################################################

 trial_transform=log10((datenSmall$Terrestrial_prop/100)+(1-0.99137931)/((1-(datenSmall$Terrestrial_prop/100))+(1-0.99137931)))

daten=envdata
n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set



#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- log10(Total_abund+1) ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
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

# Run Full model
Full.model <-gamboost(formulaB,data=traindata, family=Gaussian(), control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, Northing and d30

# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='kfold', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)

