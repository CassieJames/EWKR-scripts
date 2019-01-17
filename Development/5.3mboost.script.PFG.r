##############################################################################
#### Model richness and abundance of wetland species at Hattah Lakes wetlands
#### C James October 2018


library(corrplot)
library(mboost)
library(countreg)
library(party)


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"
envdata=data.matrix.env.data # copy data

envdata=subset(envdata, envdata$Wet_Natives>0) # just model positive counts


# sort out environmental data - centre and log
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond record
envdata$TSLW[envdata$Innundated==TRUE]<-1 # if site is recorded as inundated change TSLW to 1 day
envdata$d3Mon_wet[is.na(envdata$d3Mon_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d1yrs_wet[is.na(envdata$d1yrs_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d3yrs_wet[is.na(envdata$d3yrs_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d3Mon_meandepth[is.na(envdata$d3Mon_meandepth)]<-0 # if site us never wet during time period replace with 0
envdata$d1yrs_meandepth[is.na(envdata$d1yrs_meandepth)]<-0 # if site us never wet during time period replace with 0
envdata$d3yrs_meandepth[is.na(envdata$d3yrs_meandepth)]<-0 # if site us never wet during time period replace with 0

envdata$Inundated = as.factor(envdata$Inundated)
# highly skewed variables were log10 transformed before analysis
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
envdata$d90=as.numeric(scale(log10(envdata$d90+1),center=TRUE, scale=FALSE))
envdata$TSLW=as.numeric(scale(log10(envdata$TSLW+1),center=TRUE, scale=FALSE)) #This is a very strongly skewed variable 
envdata$d3Mon_wet=as.numeric(scale(log10(envdata$d3Mon_wet+1),center=TRUE, scale=FALSE)) #
envdata$d3Mon_meandepth=as.numeric(scale(log10(envdata$d3Mon_meandepth+1),center=TRUE, scale=FALSE)) # Highly skewed - log helps a bit but its heavily weighted to the left
envdata$d1yrs_wet=as.numeric(scale((envdata$d1yrs_wet),center=TRUE, scale=FALSE)) # quite strong peaks at low and high end of range
envdata$d1yrs_meandepth=as.numeric(scale(log10(envdata$d1yrs_meandepth+1),center=TRUE, scale=FALSE)) #this is a very strongly skewed variable as there are lots of low values
envdata$d3yrs_wet=as.numeric(scale((envdata$d3yrs_wet),center=TRUE, scale=FALSE))
envdata$d3yrs_meandepth=as.numeric(scale(envdata$d3yrs_meandepth,center=TRUE, scale=FALSE)) #this is a very strongly skewed variable as there are lots of low values
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable
envdata$MeanTemp90=as.numeric(scale((envdata$MeanTemp90),center=TRUE, scale=FALSE))
envdata$Freq_d1=as.numeric(scale((envdata$Freq_d1),center=TRUE, scale=FALSE))
envdata$Freq_d3=as.numeric(scale((envdata$Freq_d3),center=TRUE, scale=FALSE))




#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Wet_Natives ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,center=TRUE, df=1, differences=1)+
bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+
bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=1)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+bbs(TSLW, by=d90,center=TRUE, df=1)+
bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=1)+bbs(d3Mon_meandepth, by=d90,center=TRUE, df=1)+
bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=1)+bbs(d3Mon_meandepth,by=d1yrs_meandepth, center=TRUE, df=1)+
bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=1)+bbs(d3Mon_meandepth,by=d3yrs_meandepth, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE)+
bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=1)+
bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=1)



# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))brandom(Site.ID, df = 1) 
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
plot(predictions[,1],daten$Wet_Natives)
abline(0,1)


# Plots of fitted and predicted against observed

png(paste(image.dir,'Wet_Natives_fitted vs observed October 2018.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
plot(fitted(Full.model, type = "response"),daten$Wet_Natives, xlab='Fitted values', ylab='Observed values')
abline(0,1)
plot(predictions_ofb,daten$Wet_Natives[which(rownames(daten) %in% weight0)])
abline(0,1)
dev.off ()


########################################################################################

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Wet_Natives ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+brandom(Site.ID,df=1)+
bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,center=TRUE, df=1, differences=1)+
bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+bbs(d90, by=TSLW,center=TRUE, df=1)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+
bbs(TSLW,by=d3Mon_meandepth, center=TRUE, df=1)+bbs(TSLW,by=d1yrs_meandepth, center=TRUE, df=1)+bbs(TSLW,by=d3yrs_meandepth, center=TRUE, df=1)+
bols(TSLW,by=d3Mon_meandepth, intercept=FALSE)+bols(TSLW,by=d1yrs_meandepth, intercept=FALSE)+bols(TSLW,by=d3yrs_meandepth, intercept=FALSE)


# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


