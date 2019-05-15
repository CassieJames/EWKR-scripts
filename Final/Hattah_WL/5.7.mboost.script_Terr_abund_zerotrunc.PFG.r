###########################################################################################################################
#### Model abundance of dryland species at Hattah Lakes wetlands
#### C James December 2018 TropWATER
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
# Have a look at response metric. Note that Dry_Natives is occurrence of species in 1:15 sub-quadrats but this value is then summed over all the 
# wetland species in that 1X15m transect so it tends to have a strange peak at 15 because there were lots of occasions with a single species that 
# occurred across all the 1m quads. This is an artifact of the method and might be a good reason to argue against the field approach?!
# Therefore the fit is generally good but underestimates in that interval if you look at the histogram.
# Envdata=subset(envdata, envdata$Terr_Natives>0) # I explored just modelling positive counts but the nbiom model looked okay with all the data
# negative binomial looks like a sensible distribution

fitdist(envdata$Terr_Natives, "nbinom")
fitD <- dnbinom(0:100, size=0.1274003, mu=4.7684495)
hist(envdata$Terr_Natives,prob=TRUE)
lines(fitD, lwd="3", col="blue")

###########################################################################################################################
# I explored including previous vegetation abundance (T-1) as potential predictor in current year (T)
sitelist=unique(envdata$Site.ID.x)
envdata$Terr_Natives_T1 <-NA

for(s in sitelist) { # loop through each site - because some sites don't have year 4 I have had to create two sets of rules

mydata=envdata[which(envdata$Site.ID.x==s),]

yoi=unique(mydata$WaterYr)
yoi=yoi[yoi != "2016"]

for(y in yoi){

if(y==yoi[1]) next # if y is the first year of the list then skip 
envdata[which(envdata$WaterYr==y & envdata$Site.ID.x==s),c("Terr_Natives_T1")]<-envdata[which(envdata$WaterYr==y-1& envdata$Site.ID.x==s),c("Terr_Natives")]
}

# Do 2016 separately as will have to use 2014 result as predictor as no 2015 data
envdata[which(envdata$WaterYr==2016 & envdata$Site.ID.x==s),c("Terr_Natives_T1")]<-envdata[which(envdata$WaterYr==2014 & envdata$Site.ID.x==s),c("Terr_Natives")]
}

envdata=envdata[!is.na(envdata$Terr_Natives_T1),] # remove sites where there is no T-1 data

###########################################################################################################################
# Sort out environmental data - centering is strongly recommended for this method

# Create flood frequency by adding the pre 2005 to the post 2005 estimates from the different models
envdata$FF=envdata$Flood_frequency+envdata$Freq_ALL

envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond record
envdata$TSLW[envdata$Innundated==TRUE]<-1 # if site is recorded as inundated change TSLW to 1 day
envdata$d3Mon_wet[is.na(envdata$d3Mon_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d1yrs_wet[is.na(envdata$d1yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3yrs_wet[is.na(envdata$d3yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3Mon_meandepth[is.na(envdata$d3Mon_meandepth)]<-0 # if site is never wet during time period replace with 0
envdata$d1yrs_meandepth[is.na(envdata$d1yrs_meandepth)]<-0 # if site is never wet during time period replace with 0
envdata$d3yrs_meandepth[is.na(envdata$d3yrs_meandepth)]<-0 # if site is never wet during time period replace with 0


sorteddata = envdata # take copy prior to centering data (this is for partial plots later on
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
envdata$Terr_Natives_T1=as.numeric(scale((envdata$Terr_Natives_T1),center=TRUE, scale=FALSE))


###########################################################################################################################				
# Candidate model formulas

dfd=1 # set degrees of freedom to 3

form_legacy <- Terr_Natives~bols(Terr_Natives_T1, intercept=FALSE)+bbs(Terr_Natives_T1, center=TRUE, df=dfd)+bols(interc,intercept=FALSE)


form_spatial <- Terr_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+bols(interc,intercept=FALSE)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) 

			
form_covary <- 	Terr_Natives~
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
				
form_covary_flow <- 	Terr_Natives~
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

				
form_covary_climate <- 	Terr_Natives~
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=dfd)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=dfd)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=dfd)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=dfd)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)
				
				
form_covary_spatial <- Terr_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
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

				
form_tree_spatial <- Terr_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) +
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree       <- Terr_Natives~
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				
form_covary_interact <- 	Terr_Natives~
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
				bbs(Freq_d3, by = FF,center=TRUE, df=dfd)+
				bbs(d90, by = d3yrs_meandepth,center=TRUE, df=dfd)+
				bbs(Freq_d3, by = d3Mon_meandepth,center=TRUE, df=dfd)+
				bbs(MeanTemp90, by = d3Mon_meandepth,center=TRUE, df=dfd)+
				bbs(Freq_d3, by = d3Mon_wet,center=TRUE, df=dfd)+
				bbs(d3Mon_meandepth, by = FF,center=TRUE, df=dfd)+
				bbs(MeanTemp90, by = TSLW,center=TRUE, df=dfd)+
				bbs(d3yrs_meandepth, by = TSLW,center=TRUE, df=dfd)+
				bbs(d3yrs_wet, by = d3Mon_meandepth,center=TRUE, df=dfd)+
				bbs(d3Mon_wet, by = FF,center=TRUE, df=dfd)
				
				
###########################################################################################################################
# Run models

daten=envdata

Model_legacy_Terr_Natives <-mboost(form_legacy,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.005)) # orig
Model_spatial_Terr_Natives <-mboost(form_spatial,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.005)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_Terr_Natives <-gamboost(form_covary,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.005)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_Terr_Natives<-gamboost(form_covary_climate,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.1)) # orig
Model_covary_flow_Terr_Natives<-gamboost(form_covary_flow,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.1)) # orig
Model_covary_spatial_Terr_Natives <-gamboost(form_covary_spatial,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.1)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_Terr_Natives <-mboost(form_tree_spatial,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.1)) # 
Model_tree_Terr_Natives <-mboost(form_tree,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.1)) # 
Model_covary_interact_Terr_Natives <-gamboost(form_covary_interact,family = Hurdle(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.1)) # 



cv5f <- cv(model.weights(Model_legacy_Terr_Natives), type='subsampling', B=25)
cv_legacy_Terr_Natives <- cvrisk(Model_legacy_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_spatial_Terr_Natives), type='subsampling', B=25)
cv_spatial_Terr_Natives <- cvrisk(Model_spatial_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_Terr_Natives), type='subsampling', B=25)
cv_covar_Terr_Natives <- cvrisk(Model_covary_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_climate_Terr_Natives), type='subsampling', B=25)
cv_covar_climate_Terr_Natives <- cvrisk(Model_covary_climate_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_flow_Terr_Natives), type='subsampling', B=25)
cv_covar_flow_Terr_Natives<- cvrisk(Model_covary_flow_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_spatial_Terr_Natives), type='subsampling', B=25)
cv_covarspatial_Terr_Natives<- cvrisk(Model_covary_spatial_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial_Terr_Natives), type='subsampling', B=25)
cv_tree_spatial_Terr_Natives<- cvrisk(Model_tree_spatial_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_Terr_Natives), type='subsampling', B=25)
cv_tree_Terr_Natives<- cvrisk(Model_tree_Terr_Natives, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_interact_Terr_Natives), type='subsampling', B=25)
cv_interact_Terr_Natives<- cvrisk(Model_covary_interact_Terr_Natives, folds=cv5f)

st<-(mstop(cv_legacy_Terr_Natives))
Model_legacy_Terr_Natives[st]

st<-(mstop(cv_spatial_Terr_Natives))
Model_spatial_Terr_Natives[st]

st<-(mstop(cv_covar_Terr_Natives))
Model_covary_Terr_Natives[st]

st<-(mstop(cv_covar_climate_Terr_Natives))
Model_covary_climate_Terr_Natives[st]

st<-(mstop(cv_covar_flow_Terr_Natives))
Model_covary_flow_Terr_Natives[st]

st<-(mstop(cv_covarspatial_Terr_Natives))
Model_covary_spatial_Terr_Natives[st]

st<-(mstop(cv_tree_spatial_Terr_Natives ))
Model_tree_spatial_Terr_Natives[st]

st<-(mstop(cv_tree_Terr_Natives ))
Model_tree_Terr_Natives[st]

st<-(mstop(cv_interact_Terr_Natives))
Model_covary_interact_Terr_Natives[st]


###########################################################################################################################				
# best model

dfd=1

form_covary_flow <- 	Terr_Natives~
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

Model_covary_flow_Terr_Natives_final<-gamboost(form_covary_flow,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig

cv5f <- cv(model.weights(Model_covary_flow_Terr_Natives_final), type='subsampling', B=25)
cv_covar_flow_Terr_Natives_final<- cvrisk(Model_covary_flow_Terr_Natives_final, folds=cv5f)

st<-(mstop(cv_covar_flow_Terr_Natives_final))
Model_covary_flow_Terr_Natives_final[st]


###########################################################################################################################
# Evaluation of different models using multiplicity adjusted all-pairwise comparisons

library(multcomp)
library(lme4)
library(rlist)

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(Legacy)","(Spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro+Spatial)","(Tree)","(Tree+Spatial)")
tmp.abund.terr.hurdle <- data.frame(cv = c(extrb(cv_legacy_Terr_Natives),extrb(cv_spatial_Terr_Natives), extrb(cv_covar_Terr_Natives), extrb(cv_covar_climate_Terr_Natives),extrb(cv_covar_flow_Terr_Natives),extrb(cv_covarspatial_Terr_Natives),
				 extrb(cv_tree_Terr_Natives),extrb(cv_tree_spatial_Terr_Natives)),
				model = gl(8, length(extrb(cv_spatial_Terr_Natives))),b = factor(rep(1:length(extrb(cv_spatial_Terr_Natives)), 8)))
				  

levels(tmp.abund.terr.hurdle$model) <- nm

save(tmp.abund.terr.hurdle,file="c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.terr.hurdle.RData")

load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.terr.hurdle.RData")

png("Terr_Natives NegLL comparing different models df1 hurdle.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.abund.terr.hurdle, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:7
		, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.terr.hurdle$model)), par("usr")[3] - 0.008,srt = 30, adj = 1, label = levels(tmp.abund.terr.hurdle$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.terr.hurdle), tmp.abund.terr.hurdle$b, function(x) lines(1:8, tmp.abund.terr.hurdle[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

# Formal model comparison
sgmod_Terr_Natives.hurdle <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.abund.terr.hurdle), mcp(model = "Tukey")))
save(sgmod_Terr_Natives.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_Terr_Natives.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_Terr_Natives.hurdle.RData")

# stability selection of important parameters in best model
sel_Terr_Natives.hurdle <- stabsel(Model_covary_flow_Terr_Natives,cutoff=0.75, q = 5)
save(sel_Terr_Natives.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Terr_Natives.hurdle.RData")

load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Terr_Natives.hurdle.RData")

# Pseudo r2 but see below for resampled estimates

		m1 <- glm.nb(Terr_Natives ~ 1, data=daten) # null model	
		y <- daten$Terr_Natives
		n <- length(y)
		r2 <- ( 1 - exp( - 2/n * (logLik(Model_covary_flow_Terr_Natives_final) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n})
	
		
############################################################################################################################################################
# Evaluate fit of the 'best' model

datenSmall=daten # create copy of data
n=nrow(datenSmall)
set.seed(806)
predicted.in.terr.hurdle <- list()
model.coefs.terr.hurdle <-list()
datenSmall.terr.hurdle <-list()
extracts.terr.TSLW.hurdle <-list()
extracts.terr.d3Mon_meandepth.hurdle <-list()
extracts.terr.d3Mon_wet.hurdle <-list()

N = 100 # this will be raised but kept to a minimum atm

for(i in 1:N) {

cat("iteration:", i, " \n") 

# Randomly permute the observations

indvecL <- sample(1:n, n, replace=FALSE) # random reordering
datenSmall <- daten[indvecL,][1:(n/2),] # subset to half data for test run

#### Run model
		Model_covary_flow_Terr_Natives.test <-mboost(form_covary_flow,family = NBinomial(),data = datenSmall, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
		cv5f <- cv(model.weights(Model_covary_flow_Terr_Natives.test), type='subsampling', B=25)
		cv_covary_flow_Terr_Natives.test<- cvrisk(Model_covary_flow_Terr_Natives.test, folds=cv5f)
		st<-(mstop(cv_covary_flow_Terr_Natives.test))
		Model_covary_flow_Terr_Natives.test[st]
		mycoefs= coef(Model_covary_flow_Terr_Natives.test) # store each model run coefficients incase we want to plot partial dependency plots for each test run to estimate ci
		model.coefs.terr.hurdle[[i]] = mycoefs # store coeffcients of each validation run
		datenSmall.terr.hurdle[[i]]= datenSmall # store actual data from each validation run
		extracts.terr.TSLW.hurdle[[i]]=extract(Model_covary_flow_Terr_Natives.test,which="TSLW")
		extracts.terr.d3Mon_meandepth.hurdle[[i]]=extract(Model_covary_flow_Terr_Natives.test,which="d3Mon_meandepth")
		extracts.terr.d3Mon_wet.hurdle[[i]]=extract(Model_covary_flow_Terr_Natives.test, which="d3Mon_wet")
#### goodness of fit in-bootstrap
		m1 <- glm.nb(Terr_Natives ~ 1, data=datenSmall) # null model	
		y <- datenSmall$Terr_Natives
		n.y <- length(y)
		r2 <- ( 1 - exp( - 2/n.y * (logLik(Model_covary_flow_Terr_Natives.test) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n.y})
		predicted.in.terr.hurdle[[i]] = r2
}
		
r2.data=list.rbind(predicted.in.terr.hurdle)
mean(r2.data)
quantile(r2.data, c(.1, .5, .9)) 

save(predicted.in.terr.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.terr.RData")
save(model.coefs.terr.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/“model.coefs.terr.RData")
save(extracts.terr.TSLW.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.TSLW.RData")
save(extracts.terr.d3Mon_meandepth.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_meandepth.RData")
save(extracts.terr.d3Mon_wet.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_wet.RData")
save(datenSmall.terr.hurdle, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.terr.RData") 

load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.terr.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/“model.coefs.terr.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.TSLW.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_meandepth.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_wet.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.terr.hurdle.RData") 

############################################################################################################################################################
# Modelled versus observed

png("Terr_Natives predict versus observed df=1 hurdle.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))

# Predictions for out-of-bootstrap data
predictions<-predict(Model_covary_flow_Terr_Natives_final, type='response')
rownames(predictions)=rownames(daten) 
predictions=as.data.frame(predictions)
plot((predictions[,1]),daten$Terr_Natives)
abline(0,1)

dev.off()

############################################################################################################################################################
# Plot partial dependency plots from best additive model
N=100
png(paste(image.dir,'Dryland_Native_marginal_plots_Feb2019.hurdle.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

#################
# plot using TSLW
mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives_final,which=1)
xmatSmooth <- extract(Model_covary_flow_Terr_Natives_final,which=2)
yvalues=xmatSmooth[[1]]%*%coef(Model_covary_flow_Terr_Natives_final)$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin[[1]] * coef(Model_covary_flow_Terr_Natives_final)$'bols(TSLW, intercept = FALSE)'
plot(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l",xlab='log10(TSLW+1)', ylab='f(TSLW)', ylim=c(-1.5,1.5))
rug(sort(envdata$TSLW+mTSLW))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.terr.hurdle.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.TSLW.hurdle.RData")

for(i in 1:N) {
mycoefs=model.coefs.terr.hurdle[[i]]
datenSmall=datenSmall.terr.hurdle[[i]]
xmats <- extracts.terr.TSLW.hurdle[[i]]
xmatLin=xmats$'bols(TSLW, intercept = FALSE)'
xmatSmooth=xmats$`bbs(TSLW, df = dfd, center = TRUE)`
yvalues=xmatSmooth%*%mycoefs$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin * mycoefs$'bols(TSLW, intercept = FALSE)'
lines(sort(datenSmall$TSLW+mTSLW),yvalues[order(datenSmall$TSLW+mTSLW)], type="l", col=alpha("grey", 0.2))
}

mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives_final,which=1)
xmatSmooth <- extract(Model_covary_flow_Terr_Natives_final,which=2)
yvalues=xmatSmooth[[1]]%*%coef(Model_covary_flow_Terr_Natives_final)$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin[[1]] * coef(Model_covary_flow_Terr_Natives_final)$'bols(TSLW, intercept = FALSE)'
lines(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l")

#################

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives_final,which=7)
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Terr_Natives_final, which=7)[[1]]
plot(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l",xlab='log10(d3Mon_meandepth+1)', ylab='f(d3Mon_meandepth)')
rug(sort(envdata$d3Mon_meandepth+md3Mon_meandepth))

for(i in 1:N) {
mycoefs=model.coefs.terr.hurdle[[i]]
datenSmall=datenSmall.terr.hurdle[[i]]
xmatLin <- extracts.terr.d3Mon_meandepth.hurdle[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3Mon_meandepth, intercept = FALSE)`
lines(sort(datenSmall$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(datenSmall$d3Mon_meandepth+md3Mon_meandepth)], type="l", col=alpha("grey", 0.2))
}

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives_final,which=7)
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Terr_Natives_final, which=7)[[1]]
lines(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l",xlab='log10(d3Mon_meandepth+1)', ylab='f(d3Mon_meandepth)')

#################

md3Mon_wet<-mean(log10(sorteddata$d3Mon_wet+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives_final,which=5)
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Terr_Natives_final, which=5)[[1]]
plot(sort(envdata$d3Mon_wet+md3Mon_wet),yvalues[order(envdata$d3Mon_wet+md3Mon_wet)], type="l",xlab='log10(d3Mon_wet+1)', ylab='f(d3Mon_wet)')
rug(sort(envdata$d3Mon_wet+md3Mon_wet))

for(i in 1:N) {

cat("iteration:", i, " \n") 
mycoefs=model.coefs.terr.hurdle[[i]]
datenSmall=datenSmall.terr.hurdle[[i]]
xmatLin <- extracts.terr.d3Mon_wet.hurdle[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3Mon_wet, intercept = FALSE)`
lines(sort(datenSmall$d3Mon_wet+md3Mon_wet),yvalues[order(datenSmall$d3Mon_wet+md3Mon_wet)], type="l", col=alpha("grey", 0.2))
}

md3Mon_wet<-mean(log10(sorteddata$d3Mon_wet+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives_final,which=5)
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Terr_Natives_final, which=5)[[1]]
lines(sort(envdata$d3Mon_wet+md3Mon_wet),yvalues[order(envdata$d3Mon_wet+md3Mon_wet)], type="l")

dev.off()

############################################################################################################################################################
# MGCV to check model fit

library(mgcv)

Variables identified as useful from mboost procedure
# Freq_d3, d3Mon_meandepth, TSLW, d1yrs_wet,d3yrs_meandepth

newdata=sorteddata[!sorteddata$Terr_Natives==0,] # remove sites where there is no abundance
mytry <- gam(Terr_Natives~s(d3Mon_meandepth)+s(log10(TSLW+1))+s(d3Mon_wet), family = nb(), data=newdata)

gam.check(mytry)

par(mfrow=c(2,2))
plot(mytry, scale=0)

resid.ssq <- sum(residuals(mytry,type="pearson")^2)  ## sum of squares of Pearson resids
resid.df <- nrow(sorteddata)-length(coef(mytry))        ## estimated resid df (N-p)
resid.ssq/resid.df   
