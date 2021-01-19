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

descdist(envdata$Terr_Natives, boot = 1000)
fitdist(envdata$Terr_Natives, "nbinom")
fitD <- dnbinom(0:100, size=0.1233853, mu=4.3215481)
hist(envdata$Terr_Natives,prob=TRUE)
lines(fitD, lwd="3", col="blue")

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
				bols(MeanTemp365, intercept=FALSE)+bbs(MeanTemp365, center=TRUE, df=dfd)+
				bols(MinTemp365, intercept=FALSE)+bbs(MinTemp365, center=TRUE, df=dfd)+
				bols(MaxTemp365, intercept=FALSE)+bbs(MaxTemp365, center=TRUE, df=dfd)+
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
				bols(MeanTemp365, intercept=FALSE)+bbs(MeanTemp365, center=TRUE, df=dfd)+
				bols(MinTemp365, intercept=FALSE)+bbs(MinTemp365, center=TRUE, df=dfd)+
				bols(MaxTemp365, intercept=FALSE)+bbs(MaxTemp365, center=TRUE, df=dfd)+
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
				bols(MeanTemp365, intercept=FALSE)+bbs(MeanTemp365, center=TRUE, df=dfd)+
				bols(MinTemp365, intercept=FALSE)+bbs(MinTemp365, center=TRUE, df=dfd)+
				bols(MaxTemp365, intercept=FALSE)+bbs(MaxTemp365, center=TRUE, df=dfd)+
				bols(interc,intercept=FALSE)			
							

tctrl = partykit::ctree_control(stump = FALSE) 

				
form_tree_spatial <- Terr_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,intercept=FALSE)+
				bols(Northing, intercept=FALSE)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=dfd,center=TRUE,differences=1)+
				bols(Easting,by=WaterYr,intercept=FALSE)+
				bols(Northing,by=WaterYr,intercept=FALSE)+
				bols(Easting, by = Northing, intercept = FALSE) %X% bols(WaterYr, intercept = FALSE) +
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,MeanTemp365,MinTemp365,MaxTemp365,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree     <- Terr_Natives~
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,MeanTemp365,MinTemp365,MaxTemp365,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

				
				
###########################################################################################################################
# Run models

daten=envdata

Model_spatial_Terr_Natives <-mboost(form_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_Terr_Natives <-mboost(form_covary,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_Terr_Natives<-mboost(form_covary_climate,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow_Terr_Natives<-mboost(form_covary_flow,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial_Terr_Natives <-mboost(form_covary_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_Terr_Natives <-mboost(form_tree_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 
Model_tree_Terr_Natives <-mboost(form_tree,family = NBinomial(),data = daten, control=boost_control(mstop=10000,trace=TRUE,nu=0.01)) # 



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

Model_covary_flow_Terr_Natives_final<-mboost(form_covary_flow,family =NBinomial() ,data = daten, control=boost_control(mstop=10000,trace=TRUE, nu=0.01)) # orig


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
nm <- c("(Spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro+Spatial)","(Tree)","(Tree+Spatial)")
tmp.abund.terr <- data.frame(cv = c(extrb(cv_spatial_Terr_Natives), extrb(cv_covar_Terr_Natives), extrb(cv_covar_climate_Terr_Natives),extrb(cv_covar_flow_Terr_Natives),extrb(cv_covarspatial_Terr_Natives),
				 extrb(cv_tree_Terr_Natives),extrb(cv_tree_spatial_Terr_Natives)),
				model = gl(7, length(extrb(cv_spatial_Terr_Natives))),b = factor(rep(1:length(extrb(cv_spatial_Terr_Natives)), 7)))
				  

levels(tmp.abund.terr$model) <- nm

save(tmp.abund.terr,file="c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.terr.RData")

load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.terr.RData")

png("Terr_Natives NegLL comparing different models May 2019.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.abund.terr, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:7
		, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.terr$model)), par("usr")[3] - 0.008,srt = 30, adj = 1, label = levels(tmp.abund.terr$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.terr), tmp.abund.terr$b, function(x) lines(1:7, tmp.abund.terr[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

# Formal model comparison
sgmod_Terr_Natives <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.abund.terr), mcp(model = "Tukey")))
save(sgmod_Terr_Natives, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_Terr_Natives.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sgmod_Terr_Natives.RData")

# stability selection of important parameters in best model
sel_Terr_Natives <- stabsel(Model_covary_flow_Terr_Natives,cutoff=0.75, q = 5)
save(sel_Terr_Natives, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Terr_Natives.RData")

load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/sel_Terr_Natives.RData")

# Pseudo r2 but see below for resampled estimates

		m1 <- glm.nb(Terr_Natives ~ 1, data=daten) # null model	
		y <- daten$Terr_Natives
		n <- length(y)
		r2 <- ( 1 - exp( - 2/n * (logLik(Model_covary_flow_Terr_Natives_final) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n})
	
###########################################################################################################################
# Evaluation of different models using multiplicity adjusted all-pairwise comparisons - INTERACTIONS included

library(multcomp)
library(lme4)
library(rlist)

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(Spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro+Spatial)","(Tree)","(Tree+Spatial)")
tmp.abund.terr.interact <- data.frame(cv = c(extrb(cv_legacy_Terr_Natives),extrb(cv_spatial_Terr_Natives), extrb(cv_covar_Terr_Natives), extrb(cv_covar_climate_Terr_Natives),extrb(cv_covar_flow_Terr_Natives),extrb(cv_covarspatial_Terr_Natives),
				 extrb(cv_tree_Terr_Natives),extrb(cv_tree_spatial_Terr_Natives),extrb(cv_interact_Terr_Natives)),
				model = gl(7, length(extrb(cv_spatial_Terr_Natives))),b = factor(rep(1:length(extrb(cv_spatial_Terr_Natives)), 7)))
				  

levels(tmp.abund.terr.interact$model) <- nm

save(tmp.abund.terr.interact,file="c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.terr.interact.RData")

load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/tmp.abund.terr.interact.RData")

png("Terr_Natives NegLL comparing different models interact.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.abund.terr.interact, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:7
		, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.terr.interact$model)), par("usr")[3] - 0.008,srt = 30, adj = 1, label = levels(tmp.abund.terr.interact$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.terr.interact), tmp.abund.terr.interact$b, function(x) lines(1:9, tmp.abund.terr.interact[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()	

		
############################################################################################################################################################
# Evaluate fit of the 'best' model

datenSmall=daten # create copy of data
n=nrow(datenSmall)
set.seed(806)
predicted.in.terr <- list()
model.coefs.terr <-list()
datenSmall.terr <-list()
extracts.terr.TSLW <-list()
extracts.terr.d3Mon_meandepth <-list()
extracts.terr.d3Mon_wet <-list()

N = 100 # this will be raised but kept to a minimum atm

for(i in 1:N) {

cat("iteration:", i, " \n") 

# Randomly permute the observations

indvecL <- sample(1:n, n, replace=FALSE) # random reordering
datenSmall <- daten[indvecL,][1:(n/2),] # subset to half data for test run

#### Run model
		Model_covary_flow_Terr_Natives.test <-mboost(form_covary_flow,family = NBinomial(),data = datenSmall, control=boost_control(mstop=10000,trace=FALSE,nu=0.01)) # 
		cv5f <- cv(model.weights(Model_covary_flow_Terr_Natives.test), type='subsampling', B=25)
		cv_covary_flow_Terr_Natives.test<- cvrisk(Model_covary_flow_Terr_Natives.test, folds=cv5f)
		st<-(mstop(cv_covary_flow_Terr_Natives.test))
		Model_covary_flow_Terr_Natives.test[st]
		mycoefs= coef(Model_covary_flow_Terr_Natives.test) # store each model run coefficients incase we want to plot partial dependency plots for each test run to estimate ci
		model.coefs.terr[[i]] = mycoefs # store coeffcients of each validation run
		datenSmall.terr[[i]]= datenSmall # store actual data from each validation run
		extracts.terr.TSLW[[i]]=extract(Model_covary_flow_Terr_Natives.test,which="TSLW")
		extracts.terr.d3Mon_meandepth[[i]]=extract(Model_covary_flow_Terr_Natives.test,which="d3Mon_meandepth")
		extracts.terr.d3Mon_wet[[i]]=extract(Model_covary_flow_Terr_Natives.test, which="d3Mon_wet")
#### goodness of fit in-bootstrap
		m1 <- glm.nb(Terr_Natives ~ 1, data=datenSmall) # null model	
		y <- datenSmall$Terr_Natives
		n.y <- length(y)
		r2 <- ( 1 - exp( - 2/n.y * (logLik(Model_covary_flow_Terr_Natives.test) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n.y})
		predicted.in.terr[[i]] = r2
}
		
r2.data=list.rbind(predicted.in.terr)
mean(r2.data)
quantile(r2.data, c(.1, .5, .9)) 

save(predicted.in.terr, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.terr.RData")
save(model.coefs.terr, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/“model.coefs.terr.RData")
save(extracts.terr.TSLW, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.TSLW.RData")
save(extracts.terr.d3Mon_meandepth, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_meandepth.RData")
save(extracts.terr.d3Mon_wet, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_wet.RData")
save(datenSmall.terr, file = "c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.terr.RData") 
#
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/predicted.in.terr.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/“model.coefs.terr.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.TSLW.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_meandepth.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.d3Mon_wet.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/datenSmall.terr.RData") 

############################################################################################################################################################
# Modelled versus observed

png("Terr_Natives predict versus observed df=1.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))

# Predictions for out-of-bootstrap data
predictions<-predict(Model_covary_flow_Terr_Natives, type='response')
rownames(predictions[,1])=rownames(daten) 
predictions=as.data.frame(predictions)
plot((predictions[,1]),daten$Terr_Natives)
abline(0,1)

dev.off()

############################################################################################################################################################
# Plot partial dependency plots from best additive model

is.not.null <- function(x) !is.null(x)
N=100
png(paste(image.dir,'Dryland_Native_marginal_plots_May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

#################
# plot using TSLW
mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives,which=1)
xmatSmooth <- extract(Model_covary_flow_Terr_Natives,which=2)
yvalues=xmatSmooth[[1]]%*%coef(Model_covary_flow_Terr_Natives)$`bbs(TSLW, df = dfd, center = TRUE)`+ xmatLin[[1]] * coef(Model_covary_flow_Terr_Natives)$'bols(TSLW, intercept = FALSE)'
plot(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l",xlab='log10(TSLW+1)', ylab='f(TSLW)', ylim=c(-1.5,1.5))
rug(sort(envdata$TSLW+mTSLW))

# run through boot strapped samples for plotting on partial effects model
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/model.coefs.terr.RData")
load("c:/Users/jc246980/Documents/Current projects/MD Vegetation/Rdata files/extracts.terr.TSLW.RData")

for(i in 1:N) {
mycoefs=model.coefs.terr[[i]]
datenSmall=datenSmall.terr[[i]]
xmats <- extracts.terr.TSLW[[i]]
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
mycoefs=model.coefs.terr[[i]]
datenSmall=datenSmall.terr[[i]]
xmatLin <- extracts.terr.d3Mon_meandepth[[i]]
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
plot(sort(envdata$d3Mon_wet+md3Mon_wet),yvalues[order(envdata$d3Mon_wet+md3Mon_wet)], type="l",xlab='(d3Mon_wet+1)', ylab='f(d3Mon_wet)')
rug(sort(envdata$d3Mon_wet+md3Mon_wet))

for(i in 1:N) {

cat("iteration:", i, " \n") 
mycoefs=model.coefs.terr[[i]]
datenSmall=datenSmall.terr[[i]]
xmatLin <- extracts.terr.d3Mon_wet[[i]]
if (is.not.null(mycoefs$`bols(d3Mon_wet, intercept = FALSE)`)){
xmatLin <- extracts.terr.d3Mon_wet[[i]]
yvalues=xmatLin[[1]]%*%mycoefs$`bols(d3Mon_wet, intercept = FALSE)`
lines(sort(datenSmall$d3Mon_wet+md3Mon_wet),yvalues[order(datenSmall$d3Mon_wet+md3Mon_wet)], type="l", col=alpha("grey", 0.2))}
else {
print("missing coef")
}}


md3Mon_wet<-mean(log10(sorteddata$d3Mon_wet+1))
xmatLin <- extract(Model_covary_flow_Terr_Natives_final,which=5)
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Terr_Natives_final, which=5)[[1]]
lines(sort(envdata$d3Mon_wet+md3Mon_wet),yvalues[order(envdata$d3Mon_wet+md3Mon_wet)], type="l")

dev.off()

############################################################################################################################################################
# Explore interactions using boosted regression trees. Note that I had to assume poisson as there is no nbinom for this method. 

library(dismo)
library(gbm)

mypreds=daten[,c("Terr_Natives","TSLW","FF","d3Mon_wet","d3Mon_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth","Freq_d1","Freq_d3","d90","d365","MeanTemp90","MinTemp90","MaxTemp90")]

Terr_Nats_dismo <- gbm.step(data=mypreds, gbm.x=2:16, gbm.y=1, family="poisson", tree.complexity=2,learning.rate=0.01, bag.fraction=0.5)
#Terr_Nats_simple <-gbm.simplify(Terr_Nats_dismo,n.drop=5)

Terr_Nats_dismo.fitvalues <- Terr_Nats_dismo$fitted # extracts the fitted values
Terr_Nats_dismo.residuals <- Terr_Nats_dismo$residuals # extracts the residuals

plot(Terr_Nats_dismo.fitvalues, Terr_Nats_dismo.residuals)

png(paste(image.dir,'Dryland_Native_marginal_plots_BRT Dec2018.png',sep=''), width=2000, height=2000, units="px", res=300)
gbm.plot(Terr_Nats_dismo,n.plot=12,write.title=FALSE)
dev.off()

find.int<-gbm.interactions(Terr_Nats_dismo)
find.int$interactions
find.int$rank.list

png(paste(image.dir,'Dryland_Native_marginal_plots_BRT interactions Feb2019.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mar = c(0.5, 0.5, 0.5, 0.5),mfrow=c(3,3),cex=0.5)

gbm.perspec(Terr_Nats_dismo, 10,2)
gbm.perspec(Terr_Nats_dismo, 13,4)
gbm.perspec(Terr_Nats_dismo, 10,4)
gbm.perspec(Terr_Nats_dismo, 4,2)
gbm.perspec(Terr_Nats_dismo, 10,8)
gbm.perspec(Terr_Nats_dismo, 7,4)
gbm.perspec(Terr_Nats_dismo, 13,1)
gbm.perspec(Terr_Nats_dismo, 11,8)
gbm.perspec(Terr_Nats_dismo, 13,2)
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
