###########################################################################################################################
#### Model abundance of wetland species at Hattah Lakes wetlands
#### C James December 2018 TropWATER, JCU, Townsville, Australia
###########################################################################################################################

# Load libraries and data and set up output locations
library(corrplot)
library(mboost)
library(countreg)
library(partykit)
library(fitdistrplus)
library(rlist)
library(scales)

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"
envdata=data.matrix.env.data # copy data

###########################################################################################################################
# Have a look at response metric. Note that Wet_Natives is occurrence of species in 1:15 sub-quadrats but this value is then summed over all the 
# wetland species in that 1X15m transect so it tends to have a strange peak at 15 because there were lots of occasions with a single species that 
# occurred across all the 1m quads. This is an artifact of the method and might be a good reason to argue against the field approach?!
# Therefore the fit is generally good but underestimates in that interval if you look at the histogram.
# Envdata=subset(envdata, envdata$Wet_Natives>0) # I explored just modelling positive counts but the nbiom model looked okay with all the data
# negative binomial looks like a sensible distribution

envdata$WetNatRich=envdata$WetNatRich+envdata$TerrNatRich

fitdist(envdata$WetNatRich, "nbinom")
fitD <- dnbinom(0:100, size=0.4847895, mu=13.8976039)
hist(envdata$Wet_Natives,prob=TRUE)
lines(fitD, lwd="3", col="blue")

###########################################################################################################################
# I explored including previous vegetation abundance (T-1) as potential predictor in current year (T)
# This step could be removed if it was felt that the spatial and spatio-temporal terms delt with any potential autocorrelation 

sitelist=unique(envdata$Site.ID.x)
envdata$WetNatRich_T1 <-NA

for(s in sitelist) { # loop through each site - because some sites don't have year 4 I have had to create two sets of rules

mydata=envdata[which(envdata$Site.ID.x==s),]

yoi=unique(mydata$WaterYr)
yoi=yoi[yoi != "2016"]

for(y in yoi){

if(y==yoi[1]) next # if y is the first year of the list then skip 
envdata[which(envdata$WaterYr==y & envdata$Site.ID.x==s),c("WetNatRich_T1")]<-envdata[which(envdata$WaterYr==y-1& envdata$Site.ID.x==s),c("Wet_Natives")]
}

# Do 2016 separately as will have to use 2014 result as predictor as no 2015 data
envdata[which(envdata$WaterYr==2016 & envdata$Site.ID.x==s),c("Wet_Natives_T1")]<-envdata[which(envdata$WaterYr==2014 & envdata$Site.ID.x==s),c("Wet_Natives")]
}

envdata=envdata[!is.na(envdata$WetNatRich_T1),] # remove sites where there is no T-1 data

###########################################################################################################################
# Prepare environmental data - centering is strongly recommended for this method

# Create flood frequency by adding the pre 2005 to the post 2005 estimates from the different models
envdata$FF=envdata$Flood_frequency+envdata$Freq_ALL
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond hydrological record
envdata$TSLW[envdata$Innundated==TRUE]<-1 # if site is recorded as inundated TRUE change TSLW to 1 day
envdata$d3Mon_wet[is.na(envdata$d3Mon_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d1yrs_wet[is.na(envdata$d1yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3yrs_wet[is.na(envdata$d3yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3Mon_meandepth[is.na(envdata$d3Mon_meandepth)]<-0 # if site is never wet during time period replace with 0
envdata$d1yrs_meandepth[is.na(envdata$d1yrs_meandepth)]<-0 # if site is never wet during time period replace with 0
envdata$d3yrs_meandepth[is.na(envdata$d3yrs_meandepth)]<-0 # if site is never wet during time period replace with 0

sorteddata = envdata # take copy prior to centering data (this is for partial plots later on as I need to decenter data for the plots to be easily interpretable)
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
envdata$MeanTemp90=as.numeric(scale((envdata$MeanTemp90),center=TRUE, scale=FALSE))
envdata$MinTemp90=as.numeric(scale((envdata$MinTemp90),center=TRUE, scale=FALSE))
envdata$MaxTemp90=as.numeric(scale((envdata$MaxTemp90),center=TRUE, scale=FALSE))
envdata$Freq_d1=as.numeric(scale((envdata$Freq_d1),center=TRUE, scale=FALSE))
envdata$Freq_d3=as.numeric(scale((envdata$Freq_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d1=as.numeric(scale((envdata$CTF_prop_d1),center=TRUE, scale=FALSE))
envdata$CTF_prop_d3=as.numeric(scale((envdata$CTF_prop_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d5=as.numeric(scale((envdata$CTF_prop_d5),center=TRUE, scale=FALSE))
envdata$FF=as.numeric(scale((envdata$FF),center=TRUE, scale=FALSE))
envdata$WaterYr=as.numeric(envdata$WaterYr)
envdata <- as.data.frame(cbind(interc=1, envdata)) # and a column vector of 1's for the intercept
envdata$WetNatRichT1=as.numeric(scale((envdata$WetNatRich_T1),center=TRUE, scale=FALSE))


###########################################################################################################################				
# Candidate model formulas

dfd=1 # set degrees of freedom to 1

form_legacy <- WetNatRich~bols(WetNatRich_T1, intercept=FALSE)+bbs(WetNatRich_T1, center=TRUE, df=dfd)+bols(interc,intercept=FALSE)


form_spatial <- WetNatRich~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+bols(interc,intercept=FALSE)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) 

			
form_covary <- 	WetNatRich~
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=dfd)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=dfd)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=dfd)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=dfd)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=dfd)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=dfd)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=dfd)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=dfd)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=dfd)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=dfd)+
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=dfd)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=dfd)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=dfd)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=dfd)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)
				
form_covary_flow <- 	WetNatRich~
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=dfd)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=dfd)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=dfd)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=dfd)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=dfd)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=dfd)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=dfd)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=dfd)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=dfd)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=dfd)+bols(interc,intercept=FALSE)

				
form_covary_climate <- 	WetNatRich~
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=dfd)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=dfd)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=dfd)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=dfd)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)
				
				
form_covary_spatial <- WetNatRich~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) +
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=dfd)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=dfd)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=dfd)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=dfd)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=dfd)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=dfd)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=dfd)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=dfd)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=dfd)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=dfd)+
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=dfd)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=dfd)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=dfd)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=dfd)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)			
							

tctrl = partykit::ctree_control(stump = FALSE) 

				
form_tree_spatial <- WetNatRich~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) +
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree       <- WetNatRich~
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				
### interaction model				
				
form_covary_interact <- 	WetNatRich~
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=dfd)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=dfd)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=dfd)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=dfd)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=dfd)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=dfd)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=dfd)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=dfd)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=dfd)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=dfd)+
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=dfd)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=dfd)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=dfd)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=dfd)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)+
				bbs(d3Mon_meandepth, by = d90,center=TRUE, df=dfd)+
				bbs(FF, by = MeanTemp90,center=TRUE, df=dfd)+
				bbs(d1yrs_wet, by = FF,center=TRUE, df=dfd)+
				bbs(d1yrs_meandepth, by = FF,center=TRUE, df=dfd)+
				bbs(Freq_d3, by = d90,center=TRUE, df=dfd)
				

###########################################################################################################################
# Run models

daten=envdata

Model_legacy_WetNatRich <-mboost(form_legacy,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_spatial_WetNatRich <-mboost(form_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_WetNatRich <-gamboost(form_covary,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_WetNatRich<-gamboost(form_covary_climate,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow_WetNatRich<-gamboost(form_covary_flow,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial_WetNatRich <-gamboost(form_covary_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_WetNatRich <-mboost(form_tree_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
Model_tree_WetNatRich <-mboost(form_tree,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
Model_covary_interact_WetNatRich <-gamboost(form_covary_interact,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) 

cv5f <- cv(model.weights(Model_legacy_WetNatRich), type='subsampling', B=25)
cv_legacy_WetNatRich <- cvrisk(Model_legacy_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_spatial_WetNatRich), type='subsampling', B=25)
cv_spatial_WetNatRich <- cvrisk(Model_spatial_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_WetNatRich), type='subsampling', B=25)
cv_covar_WetNatRich <- cvrisk(Model_covary_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_climate_WetNatRich), type='subsampling', B=25)
cv_covar_climate_WetNatRich <- cvrisk(Model_covary_climate_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_flow_WetNatRich), type='subsampling', B=25)
cv_covar_flow_WetNatRich<- cvrisk(Model_covary_flow_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_spatial_WetNatRich), type='subsampling', B=25)
cv_covarspatial_WetNatRich<- cvrisk(Model_covary_spatial_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial_WetNatRich), type='subsampling', B=25)
cv_tree_spatial_WetNatRich<- cvrisk(Model_tree_spatial_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_WetNatRich), type='subsampling', B=25)
cv_tree_WetNatRich<- cvrisk(Model_tree_WetNatRich, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_interact_WetNatRich), type='subsampling', B=25)
cv_Model_covary_interact_WetNatRich<- cvrisk(Model_covary_interact_WetNatRich, folds=cv5f)

st<-(mstop(cv_legacy_WetNatRich)) # 343
Model_legacy_WetNatRich[st]

st<-(mstop(cv_spatial_WetNatRich)) # 352
Model_spatial_WetNatRichs[st]

st<-(mstop(cv_covar_WetNatRich)) # 10000
Model_covary_WetNatRich[st]

st<-(mstop(cv_covar_climate_WetNatRich)) # 10000
Model_covary_climate_WetNatRich[st]

st<-(mstop(cv_covar_flow_WetNatRich)) #10000
Model_covary_flow_WetNatRich[st]

st<-(mstop(cv_covarspatial_WetNatRich)) #10000
Model_covary_spatial_WetNatRich[st]

st<-(mstop(cv_tree_spatial_WetNatRich )) # 1195
Model_tree_spatial_WetNatRich[st]

st<-(mstop(cv_tree_WetNatRich )) # 687
Model_tree_WetNatRich[st]

st<-(mstop(cv_Model_covary_interact_WetNatRich)) #10000
Model_covary_interact_WetNatRich[st]


###########################################################################################################################
### Best model

dfd=1

form_covary_flow <- 	Wet_Natives~
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=dfd)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=dfd)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=dfd)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=dfd)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=dfd)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=dfd)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=dfd)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=dfd)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=dfd)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=dfd)+bols(interc,intercept=FALSE)

Model_covary_flow_Wet_Natives_final<-gamboost(form_covary_flow,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig

cv5f <- cv(model.weights(Model_covary_flow_Wet_Natives_final), type='subsampling', B=25)
cv_covar_flow_Wet_Natives_final<- cvrisk(Model_covary_flow_Wet_Natives_final, folds=cv5f)

st<-(mstop(cv_covar_flow_Wet_Natives_final)) #10000
Model_covary_flow_Wet_Natives_final[st]

###########################################################################################################################
# Evaluation of different models using multiplicity adjusted all-pairwise comparisons

library(multcomp)
library(lme4)
library(rlist)

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(Legacy)","(Spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro+spatial)","(Tree)","(Tree+spatial)")
tmp.abund.wet <- data.frame(cv = c(extrb(cv_legacy_WetNatRich),extrb(cv_spatial_WetNatRich), extrb(cv_covar_WetNatRich), extrb(cv_covar_climate_WetNatRich),extrb(cv_covar_flow_WetNatRich),extrb(cv_covarspatial_WetNatRich),
				 extrb(cv_tree_WetNatRich),extrb(cv_tree_spatial_WetNatRich)),
				model = gl(8, length(extrb(cv_spatial_WetNatRich))),b = factor(rep(1:length(extrb(cv_spatial_WetNatRich)), 8)))
				  

levels(tmp.abund.wet$model) <- nm

tmp.abund.rich=tmp.abund.wet

save(tmp.abund.wet, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.rich.RData") 
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.rich.RData")

png("WetNatRich NegLL comparing different models df1.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.abund.rich, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:7
		, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.rich$model)), par("usr")[3] - 0.008,srt = 30, adj = 1, label = levels(tmp.abund.rich$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.rich), tmp.abund.rich$b, function(x) lines(1:8, tmp.abund.rich[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

# Formal model comparison
sgmod_rich_Natives <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.abund.rich), mcp(model = "Tukey")))

save(sgmod_rich_Natives, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_rich_Natives.RData")
load( "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_rich_Natives.RData")

# stability selection of important parameters in additive model
sel_Rich_Natives <- stabsel(Model_covary_flow_WetNatRich,cutoff=0.75, q = 5)

save(sel_Rich_Natives, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Rich_Natives.RData") 
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Rich_Natives.RData")

# Pseudo r2 but see below for resampled estimates

		m1 <- glm.nb(WetNatRich ~ 1, data=daten) # null model	
		y <- daten$WetNatRich
		n <- length(y)
		r2 <- ( 1 - exp( - 2/n * (logLik(Model_covary_flow_WetNatRich) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n})


############################################################################################################################################################
# Evaluate fit of the 'best' model through boot strapping

datenSmall=daten # create copy of data
n=nrow(datenSmall)
set.seed(806)
predicted.in.rich <- list()
model.coefs.rich <-list()
datenSmall.rich <-list()
extracts.rich.TSLW <-list()
extracts.rich.d3Mon_meandepth <-list()
extracts.rich.d3Mon_wet <-list()

N = 100

for(i in 1:N) {

cat("iteration:", i, " \n") 

# Randomly permute the observations

indvecL <- sample(1:n, n, replace=FALSE) # random reordering
datenSmall <- daten[indvecL,][1:(n/2),] # subset to half data for test run

#### Run model
		Model_covary_flow_WetNatRichboot <-mboost(form_covary_flow,family = NBinomial(),data = datenSmall, control=boost_control(mstop=10000,trace=FALSE,nu=0.01)) # 
		cv5f <- cv(model.weights(Model_covary_flow_WetNatRichboot), type='subsampling', B=25)
		cv_covar_flow_WetNatRichboot<- cvrisk(Model_covary_flow_WetNatRichboot, folds=cv5f)
		st<-(mstop(cv_covar_flow_WetNatRichboot))
		Model_covary_flow_WetNatRichboot[st]
		mycoefs= coef(Model_covary_flow_WetNatRichboot) # store each model run coefficients incase we want to plot partial dependency plots for each test run to estimate ci
		extracts.rich.TSLW[[i]]=extract(Model_covary_flow_WetNatRichboot,which="TSLW")
		extracts.rich.d3Mon_meandepth[[i]]=extract(Model_covary_flow_WetNatRichboot,which="d3Mon_meandepth")
		extracts.rich.d3Mon_wet[[i]]=extract(Model_covary_flow_WetNatRichboot,which="d3Mon_wet")
		model.coefs.rich[[i]] = mycoefs
		datenSmall.rich[[i]] = datenSmall
#### goodness of fit in-bootstrap
		m1 <- glm.nb(WetNatRich ~ 1, data=datenSmall) # null model	
		y <- datenSmall$WetNatRich
		n.y <- length(y)
		r2 <- ( 1 - exp( - 2/n.y * (logLik(Model_covary_flow_WetNatRichboot) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n.y})
		predicted.in.rich[[i]] = r2
}

# pseudo R2 values with bootsrapped confidence interval		
r2.data=list.rbind(predicted.in.rich)
mean(r2.data)
quantile(r2.data, c(.1, .5, .9)) 

# Save out lists in order to plot validation runs on partial plots
save(extracts.rich.TSLW, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.TSLW.RData")
save(extracts.rich.d3Mon_meandepth, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.d3Mon_meandepth.RData")
save(extracts.rich.d3Mon_wet, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.d3Mon_wet")

save(predicted.in.rich, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.rich.RData")
save(model.coefs.rich, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.rich.RData") 
save(datenSmall.rich, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.rich.RData") 

# load R data files later if required
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.wet.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.TSLW.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.d3Mon_meandepth.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.d3Mon_wet.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.rich.RData")

############################################################################################################################################################
# Plot observed versus fitted values

png("WetNatRich predict versus observed df=1.png",width=12, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))

# Predictions for out-of-bootstrap data
predictions<-predict(Model_covary_flow_WetNatRich,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
plot(predictions[,1],daten$WetNatRich)
abline(0,1)

dev.off()

############################################################################################################################################################
# Plot partial dependency plots from best additive model
N= 100
png(paste(image.dir,'Wetland_Rich_marginal_plots_Feb2019.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

#################
# plot using TSLW
mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatLin <- extract(Model_covary_flow_WetNatRich,which=1)
xmatSmooth <- extract(Model_covary_flow_WetNatRich,which=2)
yvalues=xmatSmooth[[1]]%*%coef(Model_covary_flow_WetNatRich)$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin[[1]] * coef(Model_covary_flow_WetNatRich)$'bols(TSLW, intercept = FALSE)'
plot(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l",xlab='log10(TSLW+1)', ylab='f(TSLW)', ylim=c(-1.5,1.5))
rug(sort(envdata$TSLW+mTSLW))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.rich.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.TSLW.RData")

for(i in 1:N) {
mycoefs=model.coefs.rich[[i]]
datenSmall=datenSmall.rich[[i]]
xmats <- extracts.rich.TSLW[[i]]
xmatLin=xmats$'bols(TSLW, intercept = FALSE)'
xmatSmooth=xmats$`bbs(TSLW, df = dfd, center = TRUE)`
yvalues=xmatSmooth%*%mycoefs$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin * mycoefs$'bols(TSLW, intercept = FALSE)'
lines(sort(datenSmall$TSLW+mTSLW),yvalues[order(datenSmall$TSLW+mTSLW)], type="l", col=alpha("grey", 0.2))
}

mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatLin <- extract(Model_covary_flow_WetNatRich,which=1)
xmatSmooth <- extract(Model_covary_flow_WetNatRich,which=2)
yvalues=xmatSmooth[[1]]%*%coef(Model_covary_flow_WetNatRich)$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin[[1]] * coef(Model_covary_flow_WetNatRich)$'bols(TSLW, intercept = FALSE)'
lines(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l")

#################

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary_flow_WetNatRich,which="d3Mon_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_WetNatRich, which=7)[[1]]
plot(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l",xlab='log10(d3Mon_meandepth+1)', ylab='f(d3Mon_meandepth)')
rug(sort(envdata$d3Mon_meandepth+md3Mon_meandepth))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.rich.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.d3Mon_meandepth.RData")

for(i in 1:N) {
mycoefs=model.coefs.rich[[i]]
datenSmall=datenSmall.rich[[i]]
xmatLin <- extracts.rich.d3Mon_meandepth[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3Mon_meandepth, intercept = FALSE)`
lines(sort(datenSmall$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(datenSmall$d3Mon_meandepth+md3Mon_meandepth)], type="l", col=alpha("grey", 0.2))
}

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary_flow_WetNatRich,which="d3Mon_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_WetNatRich, which=7)[[1]]
lines(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l")

#################

md3Mon_wet<-mean(log10(sorteddata$d3Mon_wet+1))
xmatLin <- extract(Model_covary_flow_WetNatRich,which="d3Mon_wet")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_WetNatRich, which=5)[[1]]
plot(sort(envdata$d3Mon_wet+md3Mon_wet),yvalues[order(envdata$d3Mon_wet+md3Mon_wet)], type="l",xlab='log10(d3Mon_wet+1)', ylab='f(d3Mon_wet)')
rug(sort(envdata$d3Mon_wet+md3Mon_wet))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.rich.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.rich.d3Mon_meandepth.RData")

for(i in 1:N) {
mycoefs=model.coefs.rich[[i]]
datenSmall=datenSmall.rich[[i]]
xmatLin <- extracts.rich.d3Mon_wet[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3Mon_wet, intercept = FALSE)`
lines(sort(datenSmall$d3Mon_wet+md3Mon_wet),yvalues[order(datenSmall$d3Mon_wet+md3Mon_wet)], type="l", col=alpha("grey", 0.2))
}

md3Mon_wet<-mean(log10(sorteddata$d3Mon_wet+1))
xmatLin <- extract(Model_covary_flow_WetNatRich,which="d3Mon_wet")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_WetNatRich, which=5)[[1]]
lines(sort(envdata$d3Mon_wet+md3Mon_wet),yvalues[order(envdata$d3Mon_wet+md3Mon_wet)], type="l")

#################


dev.off()








############################################################################################################################################################
# Explore interactions using boosted regression trees. Note that I had to assume poisson as there is no nbinom for this method. 

library(dismo)
library(gbm)

mypreds=daten[,c("WetNatRich","TSLW","FF","d3Mon_wet","d3Mon_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth","Freq_d1","Freq_d3","d90","d365","MeanTemp90","MinTemp90","MaxTemp90")]

Wet_Rich_dismo <- gbm.step(data=mypreds, gbm.x=2:16, gbm.y=1, family="poisson", tree.complexity=2,learning.rate=0.005, bag.fraction=0.5)
#Wet_Nats_simple <-gbm.simplify(Wet_Nats_dismo,n.drop=5)

gbm.plot.fits(Wet_Rich_dismo)

Wet_Nats_dismo.fitvalues <- Wet_Nats_dismo$fitted # extracts the fitted values
Wet_Nats_dismo.residuals <- Wet_Nats_dismo$residuals # extracts the residuals

plot(Wet_Nats_dismo.fitvalues, Wet_Nats_dismo.residuals)

png(paste(image.dir,'Wetland_Rich_marginal_plots_BRT Feb2019.png',sep=''), width=2000, height=2000, units="px", res=300)

gbm.plot(Wet_Rich_dismo,n.plot=6,write.title=FALSE)
dev.off()

find.int<-gbm.interactions(Wet_Rich_dismo)
find.int$interactions
find.int$rank.list

png(paste(image.dir,'Wetland_Native_marginal_plots_BRT interactions Feb2019.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mar = c(0.5, 0.5, 0.5, 0.5),mfrow=c(3,3),cex=0.5)

gbm.perspec(Wet_Nats_dismo, 13,2,col="grey")
gbm.perspec(Wet_Nats_dismo, 11,10)
gbm.perspec(Wet_Nats_dismo, 5,2)
gbm.perspec(Wet_Nats_dismo, 6,2)
gbm.perspec(Wet_Nats_dismo, 7,2)
gbm.perspec(Wet_Nats_dismo, 15,2)
gbm.perspec(Wet_Nats_dismo, 6,1)
gbm.perspec(Wet_Nats_dismo, 4,2)
gbm.perspec(Wet_Nats_dismo, 9,2)

dev.off()

############################################################################################################################################################
# MGCV

library(mgcv)

Variables identified as useful from mboost procedure
# Freq_d3, d3Mon_meandepth, TSLW, d1yrs_wet,d3yrs_meandepth

newdata=sorteddata[!sorteddata$Wet_Natives==0,] # remove sites where there is no abundance

mytry <- gam(Wet_Natives~s(d3Mon_meandepth)+s(log10(TSLW+1))+s(d1yrs_wet)+s(d3yrs_meandepth), family = nb(), data=newdata)

par(mfrow=c(2,2))
plot(mytry, scale=0)
