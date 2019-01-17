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
 
envdata$AT_ALL =envdata$ATl+envdata$ATe+envdata$ATw
envdata$AR_ALL =envdata$ARp+envdata$ARf
envdata$Amph=envdata$AT_ALL+envdata$AR_ALL


#envdata=subset(envdata, envdata$AT_ALL>0) # just model positive counts
#envdata=subset(envdata, envdata$AR_ALL>0) # just model positive counts
#envdata=subset(envdata, envdata$Tda>0) # just model positive counts
#envdata=subset(envdata, envdata$Tdr>0) # just model positive counts
envdata=subset(envdata, envdata$Amph>0) # just model positive counts


# sort out environmental data - centre and log
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond record
envdata$TSLW[envdata$Innundated==TRUE]<-1 # if site is recorded as inundated change TSLW to 1 day
envdata$d3Mon_wet[is.na(envdata$d3Mon_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d1yrs_wet[is.na(envdata$d1yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3yrs_wet[is.na(envdata$d3yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3Mon_meandepth[is.na(envdata$d3Mon_meandepth)]<-0 # if site us never wet during time period replace with 0
envdata$d1yrs_meandepth[is.na(envdata$d1yrs_meandepth)]<-0 # if site us never wet during time period replace with 0
envdata$d3yrs_meandepth[is.na(envdata$d3yrs_meandepth)]<-0 # if site us never wet during time period replace with 0

sorteddata = envdata
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
envdata$MinTemp90=as.numeric(scale((envdata$MinTemp90),center=TRUE, scale=FALSE))
envdata$MaxTemp90=as.numeric(scale((envdata$MaxTemp90),center=TRUE, scale=FALSE))
envdata$MeanTemp365=as.numeric(scale((envdata$MeanTemp365),center=TRUE, scale=FALSE))
envdata$MinTemp365=as.numeric(scale((envdata$MinTemp365),center=TRUE, scale=FALSE))
envdata$MaxTemp365=as.numeric(scale((envdata$MaxTemp365),center=TRUE, scale=FALSE))
envdata$Freq_d1=as.numeric(scale((envdata$Freq_d1),center=TRUE, scale=FALSE))
envdata$Freq_d3=as.numeric(scale((envdata$Freq_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d1=as.numeric(scale((envdata$CTF_prop_d1),center=TRUE, scale=FALSE))
envdata$CTF_prop_d3=as.numeric(scale((envdata$CTF_prop_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d5=as.numeric(scale((envdata$CTF_prop_d5),center=TRUE, scale=FALSE))


# Amphibious species abundance look at basic relationship between flood characteristics and response metrics

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Amph ~ brandom(Site.ID,df=1)+bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
bols(Inundated, intercept=FALSE)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+
bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=1)+
bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=1)+
bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=1)+
bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=1)+
bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=1)+
bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=1)+
bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=1)+
bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=1)+
bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=1)+
bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=1)+
bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=1)+
bols(MeanTemp365, intercept=FALSE)+bbs(MeanTemp365, center=TRUE, df=1)+
bols(MinTemp365, intercept=FALSE)+bbs(MinTemp365, center=TRUE, df=1)+
bols(MaxTemp365, intercept=FALSE)+bbs(MaxTemp365, center=TRUE, df=1)



daten=envdata # take a copy of the data


# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=5)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)


#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Wet_Natives ~ brandom(Site.ID,df=1)+bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
bols(Inundated, intercept=FALSE)+
bbs(TSLW, center=TRUE, df=1)+
bbs(d3Mon_wet, center=TRUE, df=1)+
bbs(d3Mon_meandepth, center=TRUE, df=1)+
bbs(d1yrs_wet, center=TRUE, df=1)+
bbs(d1yrs_meandepth, center=TRUE, df=1)+
bbs(d3yrs_wet, center=TRUE, df=1)+
bbs(d3yrs_meandepth, center=TRUE, df=1)+
bbs(Freq_d1, center=TRUE, df=1)+
bbs(Freq_d3, center=TRUE, df=1)+
bbs(d90, center=TRUE, df=1)+
bbs(d365, center=TRUE, df=1)+
bbs(MeanTemp90, center=TRUE, df=1)+
bbs(MinTemp90, center=TRUE, df=1)+
bbs(MaxTemp90, center=TRUE, df=1)+
bbs(MeanTemp365, center=TRUE, df=1)+
bbs(MinTemp365, center=TRUE, df=1)+
bbs(MaxTemp365, center=TRUE, df=1)

# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Wet_Natives ~ brandom(Site.ID,df=3)+bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
bols(Inundated, intercept=FALSE)+
bbs(TSLW, center=TRUE, df=3)+
bbs(d3Mon_wet, center=TRUE, df=3)+
bbs(d3Mon_meandepth, center=TRUE, df=3)+
bbs(d1yrs_wet, center=TRUE, df=3)+
bbs(d1yrs_meandepth, center=TRUE, df=3)+
bbs(d3yrs_wet, center=TRUE, df=3)+
bbs(d3yrs_meandepth, center=TRUE, df=3)+
bbs(Freq_d1, center=TRUE, df=3)+
bbs(Freq_d3, center=TRUE, df=3)+
bbs(d90, center=TRUE, df=3)+
bbs(d365, center=TRUE, df=3)+
bbs(MeanTemp90, center=TRUE, df=3)+
bbs(MinTemp90, center=TRUE, df=3)+
bbs(MaxTemp90, center=TRUE, df=3)+
bbs(MeanTemp365, center=TRUE, df=3)+
bbs(MinTemp365, center=TRUE, df=3)+
bbs(MaxTemp365, center=TRUE, df=3)

# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North

ctrl = party::ctree_control(stump = FALSE) 

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- Wet_Natives ~ brandom(Site.ID,df=1)+bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
btree(Inundated, tree_controls=ctrl)+
btree(TSLW, tree_controls=ctrl)+
btree(d3Mon_wet, tree_controls=ctrl)+
btree(d3Mon_meandepth, tree_controls=ctrl)+
btree(d1yrs_wet, tree_controls=ctrl)+
btree(d1yrs_meandepth, tree_controls=ctrl)+
btree(d3yrs_wet, tree_controls=ctrl)+
btree(d3yrs_meandepth, tree_controls=ctrl)+
btree(Freq_d1, tree_controls=ctrl)+
btree(Freq_d3, tree_controls=ctrl)+
btree(d90, tree_controls=ctrl)+
btree(d365, tree_controls=ctrl)+
btree(MeanTemp90, tree_controls=ctrl)+
btree(MinTemp90, tree_controls=ctrl)+
btree(MaxTemp90, tree_controls=ctrl)+
btree(MeanTemp365, tree_controls=ctrl)+
btree(MinTemp365, tree_controls=ctrl)+
btree(MaxTemp365, tree_controls=ctrl)


# Run model on training set
Full.model <-gamboost(formulaB,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North





png(paste(image.dir,'Amphibious_marginal_plots.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

# plot using d1yrs_wet
md1yrs_wet<-mean((sorteddata$d1yrs_wet))
xmatSmooth <- extract(Full.model,which=11)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d1yrs_wet, df = 1, center = TRUE)`
plot(sort(envdata$d1yrs_wet+md1yrs_wet),yvalues[order(envdata$d1yrs_wet+md1yrs_wet)], type="l",xlab='d1yrs_wet', ylab='f(d1yrs_wet)')
rug(sort(envdata$d1yrs_wet+md1yrs_wet))


md1yrs_meandepth<-mean(log10(sorteddata$d1yrs_meandepth+1))
xmatSmooth <- extract(Full.model,which=13)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d1yrs_meandepth, df = 1, center = TRUE)`
plot(sort(envdata$d1yrs_meandepth+md1yrs_meandepth),yvalues[order(envdata$d1yrs_meandepth+md1yrs_meandepth)], type="l",xlab='d1yrs_meandepth (log10+1)', ylab='f(d1yrs_meandepth)')
rug(sort(envdata$d1yrs_meandepth+md1yrs_meandepth))

md3Mon_wet<-mean((sorteddata$d3Mon_wet))
xmatSmooth <- extract(Full.model,which=7)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d3Mon_wet, df = 1, center = TRUE)`
plot(sort(envdata$d3Mon_wet+md3Mon_wet),yvalues[order(envdata$d3Mon_wet+md3Mon_wet)], type="l",xlab='dMon_wet', ylab='f(dMon_wet)')
rug(sort(envdata$dMon_wet+md3Mon_wet))