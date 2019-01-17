# Script to undertake analysis of Hattah Floodplains vegetation data using boosted generalized additive models
# Adapted by C James from sample code provided by Maloney et al. 2012 (Applying additive modelling, Methods in Ecology and Evolution vol 3, 116-128, Appendix E)
# 19th August 2016 
#
# Load data and libraries
library(mboost)
library(MASS)
library(PerformanceAnalytics)
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Final_Metrics_Hattah_FP_with rich.csv") # load data
image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Plots/";

# sort out environmental data - centre and log
envdata=data.matrix.env.data[,c("d30", "d90", "d180", "d365", "Inundated","Flood_frequency","TSLF", "Easting", "Northing","Site.ID", "WRC", "MeanTemp30", "MaxTemp30","MinTemp30","MeanTemp365","MaxTemp365","MinTemp365", "H_index")]
envdata$TSLF[is.na(envdata$TSLF)] <- 10000 # woops! still had NA values in this where the TSLF exceeded the available flow record so I have used an arbitrary 10000 for this 
envdata_backup=envdata
rownames(envdata)=(data.matrix$Row.names)
envdata$Inundated = as.factor(envdata$Inundated )
envdata$WRC = as.factor(envdata$WRC)
# Note: rainfall variables are (predictably) highly correlated so I am dropping the most correlated (d90 and d180) and am left with d30 and d365 which still have correlation coeff of 0.88! 
# highly skewed variables were log10 transformed before analysis
envdata$d30=as.numeric(scale(log10(envdata$d30+1),center=TRUE, scale=FALSE)) # added one to avoid return of infinity values due to low mean rainfall over shorter time periods
envdata$d90=as.numeric(scale(log10(envdata$d90+1),center=TRUE, scale=FALSE))
envdata$d180=as.numeric(scale(log10(envdata$d180+1),center=TRUE, scale=FALSE))
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
envdata$Flood_frequency=as.numeric(scale(envdata$Flood_frequency,center=TRUE, scale=FALSE))
envdata$TSLF=as.numeric(scale(log10(envdata$TSLF+1),center=TRUE, scale=FALSE))
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable
envdata$MeanTemp30=as.numeric(scale((envdata$MeanTemp30),center=TRUE, scale=FALSE))
envdata$MaxTemp30=as.numeric(scale((envdata$MaxTemp30),center=TRUE, scale=FALSE))
envdata$MinTemp30=as.numeric(scale((envdata$MinTemp30),center=TRUE, scale=FALSE))
envdata$MeanTemp365=as.numeric(scale((envdata$MeanTemp365),center=TRUE, scale=FALSE))
envdata$MaxTemp365=as.numeric(scale((envdata$MaxTemp365),center=TRUE, scale=FALSE))
envdata$MinTemp365=as.numeric(scale((envdata$MinTemp365),center=TRUE, scale=FALSE))

####################################################     DIVERSITY analysis      #########################################################################################################

daten=envdata

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- H_index ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE)

# Create test and training data sets

n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

# Run model on training set
Full.model <-gamboost(formulaB,data=datenSmall, family=Gaussian(), control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
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

# Predictions for out-of-bootstrap data
predictions_ofb<-predict(Full.model,newdata=testdata)
plot((predictions_ofb),testdata$H_index)
abline(0,1)


# SCRIPT FOR EXTRACTING RESIDUALS
# create new data frame with residuals for plotting and also for use in t+1 - This is done on the entire dataset not the data split into training and testing datasets
newdat <- cbind(daten$Site.ID, as.data.frame(predict(Full.model,type="response")-daten$Wet_Natives)) # extract residuals



newdat=cbind(newdat, daten$Unique_site_year) # apend residuals to data
colnames(newdat)=c("Site.ID", "resid", "Site.year")

substrRight <- function(x, n){ # script to grab year off end of row names
  substr(x, nchar(x)-n+1, nchar(x))
}

Site.year=as.character(newdat$Site.year)
year=sapply(Site.year, function (x) substrRight(x, 2))
year=as.data.frame(year)
newdat=cbind(newdat, year)
newdat$year=as.numeric(newdat$year)
sitelist=unique(newdat$Site.ID)
newdat$newresid=NA # creates new column into which the residuals from the last time period will be added


for(s in sitelist) { # loop through each site - because some sites don't have year 4 I have had to create two sets of rules
yoi=unique(newdat[which(newdat$Site.ID==s),c("year")])
if("4" %in% yoi){
newdat[which(newdat$Site.ID==s & newdat$year==1),c("newresid")] = 0
roi_t2=newdat[which(newdat$Site.ID==s & newdat$year==1),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==2),c("newresid")] = roi_t2
roi_t3=newdat[which(newdat$Site.ID==s & newdat$year==2),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==3),c("newresid")] = roi_t3
roi_t4=newdat[which(newdat$Site.ID==s & newdat$year==3),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==4),c("newresid")] = roi_t4
roi_t5=newdat[which(newdat$Site.ID==s & newdat$year==4),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==5),c("newresid")] = roi_t5
roi_t6=newdat[which(newdat$Site.ID==s & newdat$year==5),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==6),c("newresid")] = roi_t6
roi_t7=newdat[which(newdat$Site.ID==s & newdat$year==6),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==7),c("newresid")] = roi_t7
roi_t8=newdat[which(newdat$Site.ID==s & newdat$year==7),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==8),c("newresid")] = roi_t8}

if(!"4" %in% yoi){
newdat[which(newdat$Site.ID==s & newdat$year==1),c("newresid")] = 0
roi_t2=newdat[which(newdat$Site.ID==s & newdat$year==1),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==2),c("newresid")] = roi_t2
roi_t3=newdat[which(newdat$Site.ID==s & newdat$year==2),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==3),c("newresid")] = roi_t3
roi_t5=newdat[which(newdat$Site.ID==s & newdat$year==3),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==5),c("newresid")] = roi_t5
roi_t6=newdat[which(newdat$Site.ID==s & newdat$year==5),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==6),c("newresid")] = roi_t6
roi_t7=newdat[which(newdat$Site.ID==s & newdat$year==6),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==7),c("newresid")] = roi_t7
roi_t8=newdat[which(newdat$Site.ID==s & newdat$year==7),c("resid")] 
newdat[which(newdat$Site.ID==s & newdat$year==8),c("newresid")] = roi_t8}

}

newdat$Row.names=rownames(newdat)
daten_resid=cbind(daten,newdat[,c("newresid", "Row.names")])
#daten_resid$newresid=as.numeric(scale(daten_resid$newresid, scale=TRUE))


# Plot of residuals for each site

png(paste(data.dir,'Full model_Abund_residuals_by year and site.png',sep=''), width=4500, height=3000, units="px", res=300)
xyplot(newdat$resid ~ newdat$year | newdat$Site.ID,
  panel=function(x, y){
    panel.xyplot(x, y)
    panel.loess(x, y, span = 0.75)
    panel.lmline(x, y, lty = 2)  # Least squares broken line
  }
)
dev.off()

png(paste(data.dir,'Full model_residuals_by year.png',sep=''), width=1500, height=1000, units="px", res=300)
xyplot(newdat$resid ~ newdat$year ,
  panel=function(x, y){
    panel.xyplot(x, y)
    panel.loess(x, y, span = 0.75)
    panel.lmline(x, y, lty = 2)  # Least squares broken line
  }
)
dev.off()

# create acf plots for each site

sitelist=unique(newdat$Site.ID)
sitelist_1=sitelist[1:12]

par(mfrow = c(3,4), mar=c(1,3,0,0))
for(s in sitelist_1) { # loop through each site 
yoi=unique(newdat[which(newdat$Site.ID==s),c("resid")])
acfdat=acf(yoi)
text(1,0.8,s)
}

sitelist_1=sitelist[13:24]

par(mfrow = c(3,4), mar=c(1,3,0,0))
for(s in sitelist_1) { # loop through each site 
yoi=unique(newdat[which(newdat$Site.ID==s),c("resid")])
acfdat=acf(yoi)
text(1,0.8,s)
}

sitelist_1=sitelist[25:36]

par(mfrow = c(3,4), mar=c(1,3,0,0))
for(s in sitelist_1) { # loop through each site 
yoi=unique(newdat[which(newdat$Site.ID==s),c("resid")])
acfdat=acf(yoi)
text(1,0.8,s)
}

sitelist_1=sitelist[37:48]

par(mfrow = c(3,4), mar=c(1,3,0,0))
for(s in sitelist_1) { # loop through each site 
yoi=unique(newdat[which(newdat$Site.ID==s),c("resid")])
acfdat=acf(yoi)
text(1,0.8,s)
}

sitelist_1=sitelist[49:60]

par(mfrow = c(3,4), mar=c(1,3,0,0))
for(s in sitelist_1) { # loop through each site 
yoi=unique(newdat[which(newdat$Site.ID==s),c("resid")])
acfdat=acf(yoi)
text(1,0.8,s)
}

sitelist_1=sitelist[61:72]

par(mfrow = c(3,4), mar=c(1,3,0,0))
for(s in sitelist_1) { # loop through each site 
yoi=unique(newdat[which(newdat$Site.ID==s),c("resid")])
acfdat=acf(yoi)
text(1,0.8,s)
}


# Run Final Model for diversity 

n<- nrow(daten) # this picks up the new data set with the residuals as predictors
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten_resid[indvecL,][1:(n/2),] # create a subset of half the original data using the random numbers

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- H_index ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)+



######################################################################
# Plots

png(paste(image.dir,'H_index_fitted vs observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
plot(fitted(Full.model,type="response"),daten$Wet_Natives, xlab='Fitted values', ylab='Observed values')
abline(0,1)
plot(predict(Full.model,newdata=testdata),testdata$H_index, xlab='Predicted values', ylab='Observed values')
abline(0,1)
dev.off ()



png(paste(image.dir,'H_index_marginal_plots.png',sep=''), width=2000, height=2000, units="px", res=300)
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
plot(sort(datenSmall$d30+md30),yvalues[order(datenSmall$d30+md30)], type="l",xlab='log(d30+1)', ylab='f(log(d30+1))')
rug(sort(datenSmall$d30+md30))


# plot using d365(accumulated rainfall in 90 days prior to sampling)
md365<-mean(log10(data.matrix.env.data$d365+1))
xmatSmooth <- extract(Full.model, which=8)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d365, df = 1, center = TRUE)` 
plot(sort(datenSmall$d365+md365),yvalues[order(datenSmall$d365+md365)], type="l",xlab='log(d365+1)', ylab='f(log(d365+1))')
rug(sort(datenSmall$d365+md365))


mFf<-mean(envdata_backup$Flood_frequency)
xmatLin <- extract(Full.model, which=9)
xmatSmooth <- extract(Full.model, which=10)
# the below line had to adapted from the paper to correct for the type of speech marks, the order or the components and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(Flood_frequency, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(Flood_frequency, intercept = FALSE)`
plot(sort(datenSmall$Flood_frequency+mFf),yvalues[order(datenSmall$Flood_frequency+mFf)], type="l",xlab='Flood frequency', ylab='f(Flood frequency)')
rug(sort(datenSmall$Flood_frequency+mFf))

#  plot using Time Since Last Flood (TSLF)
mTSLF<-mean(log10(envdata_backup$TSLF+1))
xmatSmooth <- extract(Full.model, which=12)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(TSLF, df = 1, center = TRUE)`
plot(sort(datenSmall$TSLF+mTSLF),yvalues[order(datenSmall$TSLF+mTSLF)], type="l",xlab='log(TSLF+1)', ylab='f(log(TSLF+1))')
rug(sort(datenSmall$TSLF+mTSLF))



dev.off()


par(mfrow = c(3,4))
plot(Final.model, which = "bbs")
# plot just linear effects
par(mfrow = c(3,3))
plot(Full.model, which = "bols")
# Plot spatial effects



Calculate pseudo R2 using Nagelkerke

library(lme4)

m1 <-lmer(H_index~1+(1|Site.ID),data=datenSmall) # null model
#m1 <-lm(H_index~1,data=datenSmall) # null model

logLik.lmer.intercept=logLik(m1)
logLik.lmer.full=logLik(Full.model)
N.lmer.full=nrow(datenSmall)
coxsnell<- 1 - exp((logLik.lmer.intercept - logLik.lmer.full) * (2/N.lmer.full))
 nagelkerke<- coxsnell / (1 - exp(logLik.lmer.intercept * (2/N.lmer.full)))
 


R2gauss<- function(y,model){
    moy<-mean(y)
    N<- length(y)
    p<-length(coef(model))-1
    SSres<- sum((y-predict(model))^2)
    SStot<-sum((y-moy)^2)
    R2<-1-(SSres/SStot)
    Rajust<-1-(((1-R2)*(N-1))/(N-p-1))
    return(data.frame(R2,Rajust,SSres,SStot))
}