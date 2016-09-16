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
library(vegan)
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Final_Metrics_Hattah_FP.csv") # load data

# sort out environmental data - centre and logarithmize
envdata=data.matrix[,c("d30", "d90", "d180", "d365", "Inundated.y","Flood_frequency","TSLF", "Easting", "Northing","Site.ID" )]
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

####################################################DIVERSITY analysis - first step was to create a full model with all the data and look for autocorrelation between successive samples at same site

n<- nrow(daten)
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:(n/2),] # create a subset of half the original data using the random numbers

#specify 1st model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- H_index ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(Flood_frequency, by=d30, intercept=FALSE)+bbs(Flood_frequency, by=d30, center=TRUE, df=1)+
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

Full.model <-gamboost(formulaB,data=daten, family=Gaussian(), control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, Northing and d30
mopt <- mstop(aic <- AIC(Full.model)) # also suggests that mstop is 10000 during initial run

# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
# cv5f <- cv(model.weights(Full.model), type='kfold', B=5)
# cvm <- cvrisk(Full.model, folds=cv5f)
# #plot(cvm)
# st<-(mstop(cvm))
# Full.model[st]


mTSLF<-mean(log10(envdata_backup$TSLF+1))
xmatSmooth <- extract(Full.model, which=14) # NOTE - make sure that the correct coefficients are being extracted by which =???
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(TSLF, df = 1, center = TRUE)`
plot(sort(daten$TSLF+mTSLF),yvalues[order(daten$TSLF+mTSLF)], type="l",xlab='TSLF', ylab='f(log(TSLF+1))')
rug(sort(daten$TSLF+mTSLF))


mFf<-mean(envdata_backup$Flood_frequency)
xmatLin <- extract(Full.model, which=9)
xmatSmooth <- extract(Full.model, which=10)
# the below line had to adapted from the paper to correct for the type of speech marks, the order or the components and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(Flood_frequency, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(Flood_frequency, intercept = FALSE)`
plot(sort(daten$Flood_frequency+mFf),yvalues[order(daten$Flood_frequency+mFf)], type="l",xlab='Flood frequency', ylab='f(Flood frequency)')
rug(sort(daten$Flood_frequency+mFf))


md30<-mean(log10(envdata_backup$d30+1))
xmatLin <- extract(Full.model, which=5)
xmatSmooth <- extract(Full.model, which=6)
# the below line had to adapted from the paper to correct for the type of speech marks, the order or the components and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d30, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(d30, intercept = FALSE)`
plot(sort(daten$d30+md30),yvalues[order(daten$d30+md30)], type="l",xlab='d30', ylab='f(d30)')
rug(sort(daten$d30+md30))


md365<-mean(log10(envdata_backup$d365+1))
xmatSmooth <- extract(Full.model, which=8) # NOTE - make sure that the correct coefficients are being extracted by which =???
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d365, df = 1, center = TRUE)`
plot(sort(daten$d365+md365),yvalues[order(daten$d365+md365)], type="l",xlab='d365', ylab='f(log(d365+1))')
rug(sort(daten$d365+md365))

# residual data extraction plots
newdat <- cbind(daten$Site.ID, as.data.frame(residuals(Full.model)))
newdat=cbind(newdat, daten$Row.names)
colnames(newdat)=c("Site.ID", "resid", "Site.year")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

Site.year=newdat$Site.year
year=sapply(Site.year, function (x) substrRight(x, 2))

year=as.data.frame(year)

newdat=cbind(newdat, year)
newdat$year=as.numeric(newdat$year)

# Plot of residuals for each site

png(paste(data.dir,'Full model_residuals_by year and site.png',sep=''), width=4500, height=3000, units="px", res=300)
xyplot(residuals(Full.model) ~ newdat$year | newdat$Site.ID,
  panel=function(x, y){
    panel.xyplot(x, y)
    panel.loess(x, y, span = 0.75)
    panel.lmline(x, y, lty = 2)  # Least squares broken line
  }
)
dev.off()

png(paste(data.dir,'Full model_residuals_by year.png',sep=''), width=1500, height=1000, units="px", res=300)
xyplot(residuals(Full.model) ~ newdat$year ,
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

##### Full model indicates some correlation amongst residuals but it appears to be fairly simple and linear
# so the idea is to remove the correlation by using the residuals in the previous year as a predictor in the next
# first thing to do is to loop through each Site.ID

sitelist=unique(newdat$Site.ID)
newdat$newresid=NA

for(s in sitelist) { # loop through each site - because some sites don't have year 4 I have had to create two sets of rules!
yoi=unique(newdat[which(newdat$Site.ID==s),c("year")])
if("4" %in% yoi){
newdat[which(newdat$Site.ID==s & newdat$year==1),c("newresid")] = 0.1
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
newdat[which(newdat$Site.ID==s & newdat$year==1),c("newresid")] = 0.1
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

daten_resid=merge(daten,newdat[,c("newresid", "Row.names")],by="Row.names") # merge data with new residuals

daten_resid$newresid=as.numeric(scale(daten_resid$newresid, scale=TRUE))



#specify 2nd model formula - with residuals of previous time frame included
formulaB <- H_index ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bbs(newresid, centre=TRUE,df=1)+
bspatial(Easting, Northing, knots=20, center=TRUE, df=1, differences=1)+
bols(d30, intercept=FALSE)+bbs(d30, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(Flood_frequency, intercept=FALSE)+bbs(Flood_frequency, center=TRUE, df=1)+
bols(Flood_frequency, by=d30, intercept=FALSE)+bbs(Flood_frequency, by=d30, center=TRUE, df=1)+
bols(TSLF, intercept=FALSE)+bbs(TSLF, center=TRUE, df=1)+
bols(Inundated.y, intercept=FALSE) 


Full.model <-gamboost(formulaB,data=daten_resid, family=Gaussian(), control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, Northing and d30
mopt <- mstop(aic <- AIC(Full.model)) # also suggests that mstop is 10000 during initial run

# examine residuals again over time and site

png(paste(data.dir,'Final model_residuals_by year and site v2.png',sep=''), width=4500, height=3000, units="px", res=300)
xyplot(residuals(Full.model) ~ newdat$year | newdat$Site.ID,
  panel=function(x, y){
    panel.xyplot(x, y)
    panel.loess(x, y, span = 0.75)
    panel.lmline(x, y, lty = 2)  # Least squares broken line
  }
)
dev.off()

png(paste(data.dir,'Final model_residuals_by year.png',sep=''), width=1500, height=1000, units="px", res=300)
xyplot(residuals(Full.model) ~ newdat$year ,
  panel=function(x, y){
    panel.xyplot(x, y)
    panel.loess(x, y, span = 0.75)
    panel.lmline(x, y, lty = 2)  # Least squares broken line
  }
)
dev.off()

acfdat<- acf(residuals(Full.model) )

png(paste(data.dir,'Final model_acf plots by site.png',sep=''), width=4500, height=3000, units="px", res=300)
xyplot(acfdat ~ newdat$year | newdat$Site.ID,
  panel=function(x, y){
    panel.xyplot(x, y)
    panel.loess(x, y, span = 0.75)
    panel.lmline(x, y, lty = 2)  # Least squares broken line
  }
)
dev.off()
