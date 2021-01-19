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
data.matrix.env.data=read.csv("Hattah wetlands response by metrics_aquatics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"
envdata=data.matrix.env.data # copy data

###########################################################################################################################
# Have a look at response metric. Note that Wet_Natives is occurrence of species in 1:15 sub-quadrats but this value is then summed over all the 
# wetland species in that 1X15m transect so it tends to have a strange peak at 15 because there were lots of occasions with a single species that 
# occurred across all the 1m quads. This is an artifact of the method and might be a good reason to argue against the field approach?!
# Therefore the fit is generally good but underestimates in that interval if you look at the histogram.
# Envdata=subset(envdata, envdata$Wet_Natives>0) # I explored just modelling positive counts but the nbiom model looked okay with all the data
# negative binomial looks like a sensible distribution

envdata$Aqua_Natives[envdata$Aqua_Natives>0] <-1

envdata$Aqua_Natives=as.factor(envdata$Aqua_Natives)
fitdist(envdata$Aqua_Natives, "logis")
fitD=dlogis(envdata$Aqua_Natives, location = 0.08426389, scale = 0.183577, log = FALSE)

hist(envdata$Aqua_Natives,prob=TRUE)
lines(fitD, lwd="3", col="blue")

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
envdata$MeanTemp365=as.numeric(scale((envdata$MeanTemp365),center=TRUE, scale=FALSE))
envdata$MinTemp365=as.numeric(scale((envdata$MinTemp365),center=TRUE, scale=FALSE))
envdata$MaxTemp365=as.numeric(scale((envdata$MaxTemp365),center=TRUE, scale=FALSE))
envdata$Freq_d1=as.numeric(scale((envdata$Freq_d1),center=TRUE, scale=FALSE))
envdata$Freq_d3=as.numeric(scale((envdata$Freq_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d1=as.numeric(scale((envdata$CTF_prop_d1),center=TRUE, scale=FALSE))
envdata$CTF_prop_d3=as.numeric(scale((envdata$CTF_prop_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d5=as.numeric(scale((envdata$CTF_prop_d5),center=TRUE, scale=FALSE))
envdata$FF=as.numeric(scale((envdata$FF),center=TRUE, scale=FALSE))
envdata$WaterYr=as.numeric(envdata$WaterYr)
envdata <- as.data.frame(cbind(interc=1, envdata)) # and a column vector of 1's for the intercept



###########################################################################################################################				
# Candidate model formulas

dfd=1 # set degrees of freedom to 1

form_legacy <- Aqua_Natives~bols(Wet_Natives_T1, intercept=FALSE)+bbs(Wet_Natives_T1, center=TRUE, df=dfd)+bols(interc,intercept=FALSE)


form_spatial <- Aqua_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+bols(interc,intercept=FALSE)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) 

			
form_covary <- 	Aqua_Natives~
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
				bols(MeanTemp365, intercept=FALSE)+bbs(MeanTemp365, center=TRUE, df=dfd)+
				bols(MinTemp365, intercept=FALSE)+bbs(MinTemp365, center=TRUE, df=dfd)+
				bols(MaxTemp365, intercept=FALSE)+bbs(MaxTemp365, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)
				
form_covary_flow <- 	Aqua_Natives~
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

				
form_covary_climate <- 	Aqua_Natives~
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=dfd)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=dfd)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=dfd)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=dfd)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=dfd)+
				bols(MeanTemp365, intercept=FALSE)+bbs(MeanTemp365, center=TRUE, df=dfd)+
				bols(MinTemp365, intercept=FALSE)+bbs(MinTemp365, center=TRUE, df=dfd)+
				bols(MaxTemp365, intercept=FALSE)+bbs(MaxTemp365, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)
				
				
form_covary_spatial <- Aqua_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
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
				bols(MeanTemp365, intercept=FALSE)+bbs(MeanTemp365, center=TRUE, df=dfd)+
				bols(MinTemp365, intercept=FALSE)+bbs(MinTemp365, center=TRUE, df=dfd)+
				bols(MaxTemp365, intercept=FALSE)+bbs(MaxTemp365, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)			
							

tctrl = partykit::ctree_control(stump = FALSE) 

				
form_tree_spatial <- Aqua_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) +
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,MeanTemp365,MinTemp365,MaxTemp365,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree       <- Aqua_Natives~
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,MeanTemp365,MinTemp365,MaxTemp365,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				

###########################################################################################################################
# Run models

daten=envdata

Model_spatial_Wet_Natives <-mboost(form_spatial,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_Wet_Natives <-gamboost(form_covary,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_Wet_Natives<-gamboost(form_covary_climate,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow_Wet_Natives<-gamboost(form_covary_flow,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial_Wet_Natives <-gamboost(form_covary_spatial,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_Wet_Natives <-mboost(form_tree_spatial,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
Model_tree_Wet_Natives <-mboost(form_tree,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 


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



###########################################################################################################################
### Best model

dfd=1

form_covary_spatial <- Aqua_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
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

Model_covary_spatial_Aqua_Natives_final<-gamboost(form_covary_spatial,family = Binomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig

cv5f <- cv(model.weights(Model_covary_spatial_Aqua_Natives_final), type='subsampling', B=25)
cv_covar_spatial_Aqua_Natives_final<- cvrisk(Model_covary_spatial_Aqua_Natives_final, folds=cv5f)

st<-(mstop(cv_covar_spatial_Aqua_Natives_final)) #10000
Model_covary_spatial_Aqua_Natives_final[st]

###########################################################################################################################
# Evaluation of different models using multiplicity adjusted all-pairwise comparisons

library(multcomp)
library(lme4)
library(rlist)

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(Spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro+spatial)","(Tree)","(Tree+spatial)")
tmp.abund.aqua <- data.frame(cv = c(extrb(cv_spatial_Wet_Natives), extrb(cv_covar_Wet_Natives), extrb(cv_covar_climate_Wet_Natives),extrb(cv_covar_flow_Wet_Natives),extrb(cv_covarspatial_Wet_Natives),
				 extrb(cv_tree_Wet_Natives),extrb(cv_tree_spatial_Wet_Natives)),
				model = gl(7, length(extrb(cv_spatial_Wet_Natives))),b = factor(rep(1:length(extrb(cv_spatial_Wet_Natives)), 7)))
				  

levels(tmp.abund.aqua$model) <- nm

save(tmp.abund.aqua, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.aqua.RData") 
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.aqua.RData")

png("Aqua_Natives NegLL comparing different models May 2019.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.abund.aqua, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:7
		, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.aqua$model)), par("usr")[3] - 0.008,srt = 30, adj = 1, label = levels(tmp.abund.aqua$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.aqua), tmp.abund.aqua$b, function(x) lines(1:7, tmp.abund.aqua[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

# Formal model comparison
sgmod_aqua_Natives <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.abund.aqua), mcp(model = "Tukey")))

save(sgmod_aqua_Natives, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_Wet_Natives.RData")
load( "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_aqua_Natives.RData")

# stability selection of important parameters in additive model
sel_aqua_Natives <- stabsel(Model_covary_Wet_Natives,cutoff=0.75, q = 5)

save(sel_aqua_Natives, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_aqua_Natives.RData") 
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_aqua_Natives.RData")

# Pseudo r2 but see below for resampled estimates

		m1 <- glm(Aqua_Natives ~ 1, data=daten, family=binomial) # null model	
		y <- daten$Aqua_Natives
		n <- length(y)
		r2 <- ( 1 - exp( - 2/n * (logLik(Model_covary_spatial_Aqua_Natives_final) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n})

############################################################################################################################################################
# Evaluate fit of the 'best' model through boot strapping

datenSmall=daten # create copy of data
n=nrow(datenSmall)
set.seed(806)
predicted.in.aqua <- list()
model.coefs.aqua <-list()
datenSmall.aqua <-list()
extracts.aqua.TSLW <-list()
extracts.aqua.d3Mon_meandepth <-list()
extracts.aqua.d3yrs_wet <-list()
extracts.aqua.Freq_d1 <-list()

N = 100

for(i in 1:N) {

cat("iteration:", i, " \n") 

# Randomly permute the observations

indvecL <- sample(1:n, n, replace=FALSE) # random reordering
datenSmall <- daten[indvecL,][1:(n/2),] # subset to half data for test run


#### Run model
		Model_covary_spatial_Aqua_Natives <-mboost(form_covary_spatial,family = Binomial(),data = datenSmall, control=boost_control(mstop=10000,trace=FALSE,nu=0.01)) # 
		cv5f <- cv(model.weights(Model_covary_spatial_Aqua_Natives), type='subsampling', B=25)
		cv_covar_spatial_Aqua_Natives<- cvrisk(Model_covary_spatial_Aqua_Natives, folds=cv5f)
		st<-(mstop(cv_covar_spatial_Aqua_Natives))
		Model_covary_spatial_Aqua_Natives[st]
		mycoefs= coef(Model_covary_spatial_Aqua_Natives) # store each model run coefficients incase we want to plot partial dependency plots for each test run to estimate ci
		extracts.aqua.TSLW[[i]]=extract(Model_covary_spatial_Aqua_Natives,which="TSLW")
		extracts.aqua.d3Mon_meandepth[[i]]=extract(Model_covary_spatial_Aqua_Natives,which="d3Mon_meandepth")
		extracts.aqua.d3yrs_wet[[i]]=extract(Model_covary_spatial_Aqua_Natives,which="d3yrs_wet")
		extracts.aqua.Freq_d1[[i]]=extract(Model_covary_spatial_Aqua_Natives,which="Freq_d1")
		model.coefs.aqua[[i]] = mycoefs
		datenSmall.aqua[[i]] = datenSmall
#### goodness of fit in-bootstrap
		m1 <- glm(Aqua_Natives ~ 1, family=binomial,data=datenSmall) # null model	
		y <- datenSmall$Aqua_Natives
		n.y <- length(y)
		r2 <- ( 1 - exp( - 2/n.y * (logLik(Model_covary_spatial_Aqua_Natives) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n.y})
		predicted.in.aqua[[i]] = r2
}

# pseudo R2 values with bootsrapped confidence interval		
r2.data=list.rbind(predicted.in.aqua)
mean(r2.data)
quantile(r2.data, c(.1, .5, .9)) 

# Save out lists in order to plot validation runs on partial plots
save(extracts.aqua.TSLW, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.TSLW.RData")
save(extracts.aqua.d3Mon_meandepth, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.d3Mon_meandepth.RData")
save(extracts.aqua.d3yrs_wet, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.d3yrs_wet.RData")
save(extracts.aqua.Freq_d1, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.Freq_d1.RData")
save(predicted.in.aqua, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.aqua.RData")
save(model.coefs.aqua, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.aqua.RData") 
save(datenSmall.aqua, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.aqua.RData") 

# load R data files later if required
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.aqua.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.aqua.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.TSLW.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.d3Mon_meandepth.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.d1yrs_wet.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.d3yrs_meandepth.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.Freq_d1.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.aqua.RData")

############################################################################################################################################################
# Plot observed versus fitted values

png("Aqua_Natives predict versus observed May 2019.png",width=12, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))

# Predictions for out-of-bootstrap data
predictions<-predict(Model_covary_spatial_Aqua_Natives_final,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
plot(daten$Aqua_Natives,predictions[,1],xlab="Presence/Absence of aquatic flora", ylab="Probability of presence of aquatic flora")


dev.off()

Model_covary_spatial_Aqua_Natives_final

############################################################################################################################################################
# Plot partial dependency plots from best additive model
N= 100
library(scales)

png(paste(image.dir,'Aquatic_Native_marginal_plots_May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

is.not.null <- function(x) !is.null(x)

#################
# plot using TSLW
mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatLin <- extract(Model_covary_spatial_Aqua_Natives_final,which="TSLW")
yvalues=xmatLin[[1]]%*%coef(Model_covary_spatial_Aqua_Natives_final, which=8)[[1]]
plot(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l",xlab='log10(TSLW+1)', ylab='f(TSLW)', ylim=c(-0.25, 0.25))
rug(sort(envdata$TSLW+mTSLW))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.aqua.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.TSLW.RData")

for(i in 1:N) {
mycoefs=model.coefs.aqua[[i]]
datenSmall=datenSmall.aqua[[i]]
xmatLin <- extracts.aqua.TSLW[[i]]
if (is.not.null(mycoefs$`bols(TSLW, intercept = FALSE)`)){
yvalues=xmatLin[[1]]%*%mycoefs$`bols(TSLW, intercept = FALSE)`
lines(sort(datenSmall$TSLW+mTSLW),yvalues[order(datenSmall$TSLW+mTSLW)], type="l", col=alpha("grey", 0.2))}
else {
print("missing coef")
}}



mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatLin <- extract(Model_covary_spatial_Aqua_Natives_final,which="TSLW")
yvalues=xmatLin[[1]]%*%coef(Model_covary_spatial_Aqua_Natives_final, which=8)[[1]]
lines(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l")

#################

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatSmooth <- extract(Model_covary_spatial_Aqua_Natives_final,which="d3Mon_meandepth")
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_spatial_Aqua_Natives_final, which=15)[[1]]
plot(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l",xlab='log10(d3Mon_meandepth+1)', ylab='f(d3Mon_meandepth)', ylim=c(-0.4, 0.4))
rug(sort(envdata$d3Mon_meandepth+md3Mon_meandepth))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.wet.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.wet.d3Mon_meandepth.RData")

for(i in 1:N) {
mycoefs=model.coefs.aqua[[i]]
datenSmall=datenSmall.aqua[[i]]
xmatSmooth <- extracts.aqua.d3Mon_meandepth[[i]]
yvalues=xmatSmooth[[2]]%*%mycoefs$`bbs(d3Mon_meandepth, df = dfd, center = TRUE)`
lines(sort(datenSmall$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(datenSmall$d3Mon_meandepth+md3Mon_meandepth)], type="l", col=alpha("grey", 0.2))
}

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatSmooth <- extract(Model_covary_spatial_Aqua_Natives_final,which="d3Mon_meandepth")
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_spatial_Aqua_Natives_final, which=15)[[1]]
lines(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l")

#################

md3yrs_wet<-mean((sorteddata$d3yrs_wet))
xmatLin<- extract(Model_covary_spatial_Aqua_Natives_final,which="d3yrs_wet")
yvalues=xmatLin[[1]]%*%coef(Model_covary_spatial_Aqua_Natives_final, which=20)[[1]]
plot(sort(envdata$d3yrs_wet+md3yrs_wet),yvalues[order(envdata$d3yrs_wet+md3yrs_wet)], type="l",xlab='d3yrs_wet', ylab='f(d3yrs_wet)')
rug(sort(envdata$d3yrs_wet+md3yrs_wet))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.aqua.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.d3yrs_wet.RData")

for(i in 1:N) {
mycoefs=model.coefs.aqua[[i]]
datenSmall=datenSmall.aqua[[i]]
xmatLin <- extracts.aqua.d3yrs_wet[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3yrs_wet, intercept = FALSE)`
lines(sort(datenSmall$d3yrs_wet+md3yrs_wet),yvalues[order(datenSmall$d3yrs_wet+md3yrs_wet)], type="l", col=alpha("grey", 0.2))
}

md3yrs_meandepth<-mean((sorteddata$d3yrs_wet))
xmatLin<- extract(Model_covary_spatial_Aqua_Natives_final,which="d3yrs_wet")
yvalues=xmatLin[[1]]%*%coef(Model_covary_spatial_Aqua_Natives_final, which=20)[[1]]
lines(sort(envdata$d3yrs_wet+md3yrs_wet),yvalues[order(envdata$d3yrs_wet+md3yrs_wet)], type="l")

#################
# plot using Freq_d1
mFreq_d1<-mean((sorteddata$Freq_d1))
xmatLin <- extract(Model_covary_spatial_Aqua_Natives_final,which="Freq_d1")
yvalues=xmatLin[[1]]%*%coef(Model_covary_spatial_Aqua_Natives_final)$`bols(Freq_d1, intercept = FALSE)`
plot(sort(envdata$Freq_d1+mFreq_d1),yvalues[order(envdata$Freq_d1+mFreq_d1)], type="l",xlab='Freq_d1', ylab='f(Freq_d1)')
rug(sort(envdata$Freq_d1+mFreq_d1))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.aqua.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.aqua.Freq_d1.RData")

for(i in 1:N) {
mycoefs=model.coefs.aqua[[i]]
datenSmall=datenSmall.aqua[[i]]
xmatLin <- extracts.aqua.Freq_d1[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(Freq_d1, intercept = FALSE)`
lines(sort(datenSmall$Freq_d1+mFreq_d1),yvalues[order(datenSmall$Freq_d1+mFreq_d1)], type="l", col=alpha("grey", 0.2))
}

mFreq_d1<-mean((sorteddata$Freq_d1))
xmatLin <- extract(Model_covary_spatial_Aqua_Natives_final,which="Freq_d1")
yvalues=xmatLin[[1]]%*%coef(Model_covary_spatial_Aqua_Natives_final)$`bols(Freq_d1, intercept = FALSE)`
lines(sort(envdata$Freq_d1+mFreq_d1),yvalues[order(envdata$Freq_d1+mFreq_d1)], type="l")


dev.off()


############################################################################################################################################################
# Explore interactions using boosted regression trees. Note that I had to assume poisson as there is no nbinom for this method. 

library(dismo)
library(Metrics)
library(groupdata2) # has fold function that allows user to specify group

hehehe=fold(daten,k=10,id_col='Site.ID.x') # creates folds of data whilst keeping site.id in same group
set.seed=806
mypreds=daten[,c("Aqua_Natives","TSLW","FF","d3Mon_wet","d3Mon_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth","Freq_d1","Freq_d3","d90","d365","MeanTemp90","MinTemp90","MaxTemp90","MeanTemp365","MinTemp365","MaxTemp365")]
mypredsv2=daten[,c("TSLW","FF","d3Mon_wet","d3Mon_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth","Freq_d1","Freq_d3","d90","d365","MeanTemp90","MinTemp90","MaxTemp90","MeanTemp365","MinTemp365","MaxTemp365")]

mypreds$Aqua_Natives <-as.numeric(as.character(mypreds[,1]))

Aqua_Nats_dismo=gbm.step(data=mypreds, gbm.x = 2:19, gbm.y = 1,family = "bernoulli", tree.complexity = 2,learning.rate = 0.01, bag.fraction = 0.5,fold.vector=hehehe$.folds)

int.null.deviance = Aqua_Nats_dismo$self.statistics$mean.null 
int.residual.deviance = Aqua_Nats_dismo$cv.statistics$deviance.mean 
int.dev = (int.null.deviance-int.residual.deviance)/int.null.deviance # percent deviance explained in the training dataset
myrmse = rmse(mypreds[,1], Aqua_Nats_dismo$fitted)
int.dev
myrmse

gbm.plot.fits(Aqua_Nats_dismo)

Aqua_Nats_dismo.fitvalues <- Aqua_Nats_dismo$fitted # extracts the fitted values
Aqua_Nats_dismo.residuals <- Aqua_Nats_dismo$residuals # extracts the residuals

plot(Aqua_Nats_dismo.fitvalues, Aqua_Nats_dismo.residuals)

png(paste(image.dir,'Aquatic_Native_marginal_plots_BRT May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

gbm.plot(Aqua_Nats_dismo,n.plot=3,write.title=FALSE)
dev.off()

find.int<-gbm.interactions(Aqua_Nats_dismo)
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

library(gbm)


mypreds=daten[,c("Wet_Natives","TSLW","FF","d3Mon_wet","d3Mon_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth","Freq_d1","Freq_d3","d90","d365","MeanTemp90","MinTemp90","MaxTemp90")]
mydat=daten[,c("TSLW","FF","d3Mon_wet","d3Mon_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth","Freq_d1","Freq_d3","d90","d365","MeanTemp90","MinTemp90","MaxTemp90")]

test.gbm <- gbm.step(mypreds$Wet_Natives ~ .,
                  data = mydat, 
                  distribution = "poisson", 
                  bag.fraction = 0.5, 
                  shrinkage = 0.05, 
                  n.minobsinnode = 50, 
				  interaction.depth = 2,
                  n.trees = 100,
                  cv.folds = 5, 
                  keep.data = FALSE)
				  

############################################################################################################################################################
# MGCV

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
