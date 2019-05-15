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

fitdist(envdata$Wet_Natives, "nbinom")
fitD <- dnbinom(0:100, size=0.4847895, mu=13.8976039)
hist(envdata$Wet_Natives,prob=TRUE)
lines(fitD, lwd="3", col="blue")

###########################################################################################################################
# I explored including previous vegetation abundance (T-1) as potential predictor in current year (T)
# This step could be removed if it was felt that the spatial and spatio-temporal terms delt with any potential autocorrelation 

sitelist=unique(envdata$Site.ID.x)
envdata$Wet_Natives_T1 <-NA

for(s in sitelist) { # loop through each site - because some sites don't have year 4 I have had to create two sets of rules

mydata=envdata[which(envdata$Site.ID.x==s),]

yoi=unique(mydata$WaterYr)
yoi=yoi[yoi != "2016"]

for(y in yoi){

if(y==yoi[1]) next # if y is the first year of the list then skip 
envdata[which(envdata$WaterYr==y & envdata$Site.ID.x==s),c("Wet_Natives_T1")]<-envdata[which(envdata$WaterYr==y-1& envdata$Site.ID.x==s),c("Wet_Natives")]
}

# Do 2016 separately as will have to use 2014 result as predictor as no 2015 data
envdata[which(envdata$WaterYr==2016 & envdata$Site.ID.x==s),c("Wet_Natives_T1")]<-envdata[which(envdata$WaterYr==2014 & envdata$Site.ID.x==s),c("Wet_Natives")]
}

envdata=envdata[!is.na(envdata$Wet_Natives_T1),] # remove sites where there is no T-1 data

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


envdata=envdata[!envdata$Wet_Natives==0,] # remove sites where there is no abundance
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
envdata$Wet_Natives_T1=as.numeric(scale((envdata$Wet_Natives_T1),center=TRUE, scale=FALSE))


###########################################################################################################################				
# Candidate model formulas

dfd=1 # set degrees of freedom to 1

form_legacy <- Wet_Natives~bols(Wet_Natives_T1, intercept=FALSE)+bbs(Wet_Natives_T1, center=TRUE, df=dfd)+bols(interc,intercept=FALSE)


form_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+bols(interc,intercept=FALSE)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) 

			
form_covary <- 	Wet_Natives~
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

				
form_covary_climate <- 	Wet_Natives~
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=dfd)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=dfd)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=dfd)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=dfd)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)
				
				
form_covary_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
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

				
form_tree_spatial <- Wet_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) +
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree       <- Wet_Natives~
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				
### interaction model				
				
form_covary_interact <- 	Wet_Natives~
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

Model_legacy_Wet_Natives <-gamboost(form_legacy,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_spatial_Wet_Natives <-mboost(form_spatial,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_Wet_Natives <-gamboost(form_covary,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_Wet_Natives<-gamboost(form_covary_climate,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow_Wet_Natives<-gamboost(form_covary_flow,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial_Wet_Natives <-gamboost(form_covary_spatial,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_Wet_Natives <-mboost(form_tree_spatial,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
Model_tree_Wet_Natives <-mboost(form_tree,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
Model_covary_interact_Wet_Natives <-gamboost(form_covary_interact,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) 

cv5f <- cv(model.weights(Model_legacy_Wet_Natives), type='subsampling', B=25)
cv_legacy_Wet_Natives <- cvrisk(Model_legacy_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_spatial_Wet_Natives), type='subsampling', B=25)
cv_spatial_Wet_Natives <- cvrisk(Model_spatial_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_Wet_Natives), type='subsampling', B=25)
cv_covar_Wet_Natives <- cvrisk(Model_covary_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_climate_Wet_Natives), type='subsampling', B=25)
cv_covar_climate_Wet_Natives <- cvrisk(Model_covary_climate_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_flow_Wet_Natives), type='subsampling', B=25)
cv_covar_flow_Wet_Natives<- cvrisk(Model_covary_flow_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_spatial_Wet_Natives), type='subsampling', B=25)
cv_covarspatial_Wet_Natives<- cvrisk(Model_covary_spatial_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial_Wet_Natives), type='subsampling', B=25)
cv_tree_spatial_Wet_Natives<- cvrisk(Model_tree_spatial_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_Wet_Natives), type='subsampling', B=25)
cv_tree_Wet_Natives<- cvrisk(Model_tree_Wet_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_interact_Wet_Natives), type='subsampling', B=25)
cv_Model_covary_interact_Wet_Natives<- cvrisk(Model_covary_interact_Wet_Natives, folds=cv5f)

st<-(mstop(cv_legacy_Wet_Natives)) # 343
Model_legacy_Wet_Natives[st]

st<-(mstop(cv_spatial_Wet_Natives)) # 352
Model_spatial_Wet_Natives[st]

st<-(mstop(cv_covar_Wet_Natives)) # 10000
Model_covary_Wet_Natives[st]

st<-(mstop(cv_covar_climate_Wet_Natives)) # 10000
Model_covary_climate_Wet_Natives[st]

st<-(mstop(cv_covar_flow_Wet_Natives)) #10000
Model_covary_flow_Wet_Natives[st]

st<-(mstop(cv_covarspatial_Wet_Natives)) #10000
Model_covary_spatial_Wet_Natives[st]

st<-(mstop(cv_tree_spatial_Wet_Natives )) # 1195
Model_tree_spatial_Wet_Natives[st]

st<-(mstop(cv_tree_Wet_Natives )) # 687
Model_tree_Wet_Natives[st]

st<-(mstop(cv_Model_covary_interact_Wet_Natives)) #10000
Model_covary_interact_Wet_Natives[st]


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

Model_covary_flow_Wet_Natives_final<-gamboost(form_covary_flow,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig

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
tmp.abund.wet.hurdle <- data.frame(cv = c(extrb(cv_legacy_Wet_Natives),extrb(cv_spatial_Wet_Natives), extrb(cv_covar_Wet_Natives), extrb(cv_covar_climate_Wet_Natives),extrb(cv_covar_flow_Wet_Natives),extrb(cv_covarspatial_Wet_Natives),
				 extrb(cv_tree_Wet_Natives),extrb(cv_tree_spatial_Wet_Natives)),
				model = gl(8, length(extrb(cv_spatial_Wet_Natives))),b = factor(rep(1:length(extrb(cv_spatial_Wet_Natives)), 8)))
				  

levels(tmp.abund.wet.hurdle$model) <- nm

save(tmp.abund.wet.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.wet.hurdle.RData") 
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.wet.hurdle.RData")

png("Wet_Natives NegLL comparing different models df1.hurdle.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.abund.wet.hurdle, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:7
		, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.wet.hurdle$model)), par("usr")[3] - 0.008,srt = 30, adj = 1, label = levels(tmp.abund.wet.hurdle$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.wet.hurdle), tmp.abund.wet.hurdle$b, function(x) lines(1:8, tmp.abund.wet.hurdle[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

# Formal model comparison
sgmod_Wet_Natives.hurdle <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.abund.wet.hurdle), mcp(model = "Tukey")))

save(sgmod_Wet_Natives.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_Wet_Natives.hurdle.RData")
load( "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_Wet_Natives.hurdle.RData")

# stability selection of important parameters in additive model
sel_Wet_Natives.hurdle <- stabsel(Model_covary_flow_Wet_Natives,cutoff=0.75, q = 5)

save(sel_Wet_Natives.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Wet_Natives.hurdle.RData") 
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Wet_Natives.hurdle.RData")

# Pseudo r2 but see below for resampled estimates

		m1 <- glm.nb(Wet_Natives ~ 1, data=daten) # null model	
		y <- daten$Wet_Natives
		n <- length(y)
		r2 <- ( 1 - exp( - 2/n * (logLik(Model_covary_flow_Wet_Natives_final) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n})


############################################################################################################################################################
# Evaluate fit of the 'best' model through boot strapping

datenSmall=daten # create copy of data
n=nrow(datenSmall)
set.seed(806)
predicted.in.wet.hurdle <- list()
model.coefs.wet.hurdle <-list()
datenSmall.wet.hurdle <-list()
extracts.wet.TSLW.hurdle <-list()
extracts.wet.d3Mon_meandepth.hurdle <-list()
extracts.wet.d1yrs_wet.hurdle <-list()
extracts.wet.Freq_d3.hurdle <-list()

N = 100

for(i in 1:N) {

cat("iteration:", i, " \n") 

# Randomly permute the observations

indvecL <- sample(1:n, n, replace=FALSE) # random reordering
datenSmall <- daten[indvecL,][1:(n/2),] # subset to half data for test run

#### Run model
		Model_covary_flow_Wet_Natives <-mboost(form_covary_flow,family = NBinomial(),data = datenSmall, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
		cv5f <- cv(model.weights(Model_covary_flow_Wet_Natives), type='subsampling', B=25)
		cv_covar_flow_Wet_Natives<- cvrisk(Model_covary_flow_Wet_Natives, folds=cv5f)
		st<-(mstop(cv_covar_flow_Wet_Natives))
		Model_covary_flow_Wet_Natives[st]
		mycoefs= coef(Model_covary_flow_Wet_Natives) # store each model run coefficients incase we want to plot partial dependency plots for each test run to estimate ci
		extracts.wet.TSLW.hurdle[[i]]=extract(Model_covary_flow_Wet_Natives,which="TSLW")
		extracts.wet.d3Mon_meandepth.hurdle[[i]]=extract(Model_covary_flow_Wet_Natives,which="d3Mon_meandepth")
		extracts.wet.d1yrs_wet.hurdle[[i]]=extract(Model_covary_flow_Wet_Natives,which="d1yrs_wet")
		extracts.wet.Freq_d3.hurdle[[i]]=extract(Model_covary_flow_Wet_Natives,which="Freq_d3")
		model.coefs.wet.hurdle[[i]] = mycoefs
		datenSmall.wet.hurdle[[i]] = datenSmall
#### goodness of fit in-bootstrap
		m1 <- glm.nb(Wet_Natives ~ 1, data=datenSmall) # null model	
		y <- datenSmall$Wet_Natives
		n.y <- length(y)
		r2 <- ( 1 - exp( - 2/n.y * (logLik(Model_covary_flow_Wet_Natives) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n.y})
		predicted.in.wet.hurdle[[i]] = r2
}

# pseudo R2 values with bootsrapped confidence interval		
r2.data=list.rbind(predicted.in.wet.hurdle)
mean(r2.data)
quantile(r2.data, c(.1, .5, .9)) 

# Save out lists in order to plot validation runs on partial plots
save(extracts.wet.TSLW.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.TSLW.hurdle.RData")
save(extracts.wet.d3Mon_meandepth.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d3Mon_meandepth.hurdle.RData")
save(extracts.wet.d1yrs_wet.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d1yrs_wet.hurdle.RData")
save(extracts.wet.d3yrs_meandepth.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d3yrs_meandepth.hurdle.RData")
save(extracts.wet.Freq_d3.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.Freq_d3.hurdle.RData")
save(predicted.in.wet.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.wet.hurdle.RData")
save(model.coefs.wet.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.hurdle.RData") 
save(datenSmall.wet.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.wet.hurdle.RData") 

# load R data files later if required
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.TSLW.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d3Mon_meandepth.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d1yrs_wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d3yrs_meandepth.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.Freq_d3.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.wet.hurdle.RData")

############################################################################################################################################################
# Plot observed versus fitted values

png("Wet_Natives predict versus observed df=1 hurdle.png",width=12, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))

# Predictions for out-of-bootstrap data
predictions<-predict(Model_covary_flow_Wet_Natives_final,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
plot(predictions[,1],daten$Wet_Natives)
abline(0,1)

dev.off()

############################################################################################################################################################
# Plot partial dependency plots from best additive model
N= 100
png(paste(image.dir,'Wetland_Native_marginal_plots_Feb2019.hurdle.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

#################
# plot using TSLW
mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatSmooth <- extract(Model_covary_flow_Wet_Natives_final,which="TSLW")
xmatLin <- extract(Model_covary_flow_Wet_Natives_final,which=1)
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_flow_Wet_Natives_final)$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin[[1]] * coef(Model_covary_flow_Wet_Natives_final)$'bols(TSLW, intercept = FALSE)'
plot(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l",xlab='log10(TSLW+1)', ylab='f(TSLW)', ylim=c(-0.25, 0.25))
rug(sort(envdata$TSLW+mTSLW))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.TSLW.hurdle.RData")

for(i in 1:N) {
mycoefs=model.coefs.wet.hurdle[[i]]
datenSmall=datenSmall.wet.hurdle[[i]]
xmatSmooth <- extracts.wet.TSLW.hurdle[[i]]
yvalues=xmatSmooth[[2]]%*%mycoefs$`bbs(TSLW, df = dfd, center = TRUE)`
lines(sort(datenSmall$TSLW+mTSLW),yvalues[order(datenSmall$TSLW+mTSLW)], type="l", col=alpha("grey", 0.2))
}

mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatSmooth <- extract(Model_covary_flow_Wet_Natives_final,which="TSLW")
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_flow_Wet_Natives_final)$`bbs(TSLW, df = dfd, center = TRUE)`
lines(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l")

#################

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary_flow_Wet_Natives_final,which="d3Mon_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Wet_Natives_final, which=7)[[1]]
plot(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l",xlab='log10(d3Mon_meandepth+1)', ylab='f(d3Mon_meandepth)')
rug(sort(envdata$d3Mon_meandepth+md3Mon_meandepth))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d3Mon_meandepth.hurdle.RData")

for(i in 1:N) {
mycoefs=model.coefs.wet.hurdle[[i]]
datenSmall=datenSmall.wet.hurdle[[i]]
xmatLin <- extracts.wet.d3Mon_meandepth.hurdle[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3Mon_meandepth, intercept = FALSE)`
lines(sort(datenSmall$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(datenSmall$d3Mon_meandepth+md3Mon_meandepth)], type="l", col=alpha("grey", 0.2))
}

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary_flow_Wet_Natives_final,which="d3Mon_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Wet_Natives_final, which=7)[[1]]
lines(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l")

#################

md3yrs_meandepth<-mean((sorteddata$d3yrs_meandepth))
xmatLin<- extract(Model_covary_flow_Wet_Natives_final,which="d3yrs_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Wet_Natives_final, which=15)[[1]]
plot(sort(envdata$d3yrs_meandepth+md3yrs_meandepth),yvalues[order(envdata$d3yrs_meandepth+md3yrs_meandepth)], type="l",xlab='d3yrs_meandepth', ylab='f(d3yrs_meandepth)')
rug(sort(envdata$d3yrs_meandepth+md3yrs_meandepth))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d3yrs_meandepth.hurdle.RData")

for(i in 1:N) {
mycoefs=model.coefs.wet.hurdle[[i]]
datenSmall=datenSmall.wet.hurdle[[i]]
xmatLin <- extracts.wet.d3yrs_meandepth.hurdle[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3yrs_meandepth, intercept = FALSE)`
lines(sort(datenSmall$d3yrs_meandepth+md3yrs_meandepth),yvalues[order(datenSmall$d3yrs_meandepth+md3yrs_meandepth)], type="l", col=alpha("grey", 0.2))
}

md3yrs_meandepth<-mean((sorteddata$d3yrs_meandepth))
xmatLin<- extract(Model_covary_flow_Wet_Natives_final,which="d3yrs_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Wet_Natives_final, which=15)[[1]]
lines(sort(envdata$d3yrs_meandepth+md3yrs_meandepth),yvalues[order(envdata$d3yrs_meandepth+md3yrs_meandepth)], type="l")

#################
# plot using d1yrs_wet
md1yrs_wet<-mean((sorteddata$d1yrs_wet))
xmatSmooth <- extract(Model_covary_flow_Wet_Natives_final,which="d1yrs_wet")
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_flow_Wet_Natives_final)$`bbs(d1yrs_wet, df = dfd, center = TRUE)`
plot(sort(envdata$d1yrs_wet+md1yrs_wet),yvalues[order(envdata$d1yrs_wet+md1yrs_wet)], type="l",xlab='d1yrs_wet', ylab='f(d1yrs_wet)', ylim=c(-0.25, 0.25))
rug(sort(envdata$d1yrs_wet+md1yrs_wet))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d1yrs_wet.hurdle.RData")

for(i in 1:N) {
mycoefs=model.coefs.wet.hurdle[[i]]
datenSmall=datenSmall.wet.hurdle[[i]]
xmatSmooth <- extracts.wet.d1yrs_wet.hurdle[[i]]
yvalues=xmatSmooth[[2]]%*%mycoefs$`bbs(d1yrs_wet, df = dfd, center = TRUE)`
lines(sort(datenSmall$d1yrs_wet+md1yrs_wet),yvalues[order(datenSmall$d1yrs_wet+md1yrs_wet)], type="l", col=alpha("grey", 0.2))
}

md1yrs_wet<-mean((sorteddata$d1yrs_wet))
xmatSmooth <- extract(Model_covary_flow_Wet_Natives_final,which="d1yrs_wet")
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_flow_Wet_Natives_final)$`bbs(d1yrs_wet, df = dfd, center = TRUE)`
lines(sort(envdata$d1yrs_wet+md1yrs_wet),yvalues[order(envdata$d1yrs_wet+md1yrs_wet)], type="l")


dev.off()


############################################################################################################################################################
# MGCV to check model fit

library(mgcv)

Variables identified as useful from mboost procedure
# Freq_d3, d3Mon_meandepth, TSLW, d1yrs_wet,d3yrs_meandepth

newdata=sorteddata[!sorteddata$Wet_Natives==0,] # remove sites where there is no abundance
mytry <- gam(Wet_Natives~s(d3Mon_meandepth)+s(log10(TSLW+1))+s(d1yrs_wet)+s(d3yrs_meandepth), family = nb(), data=newdata)

gam.check(mytry)

par(mfrow=c(2,2))
plot(mytry, scale=0)

resid.ssq <- sum(residuals(mytry,type="pearson")^2)  ## sum of squares of Pearson resids
resid.df <- nrow(sorteddata)-length(coef(mytry))        ## estimated resid df (N-p)
resid.ssq/resid.df   
