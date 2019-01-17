###########################################################################################################################
#### Model abundance of wetland species at Hattah Lakes wetlands
#### C James December 2018 TropWATER
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

# Because the legacy/carry over model was not good at explaining any variation I have removed the below step as we loose some data due to not being able
# to include the first year of monitoring in the model
envdata=envdata[!is.na(envdata$Wet_Natives_T1),] # remove sites where there is no T-1 data

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
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable
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
# Candidate models
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
				
###########################################################################################################################

daten=envdata

# Run models
Model_legacy_Wet_Natives <-mboost(form_legacy,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_spatial_Wet_Natives <-mboost(form_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_Wet_Natives <-gamboost(form_covary,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_Wet_Natives<-gamboost(form_covary_climate,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow_Wet_Natives<-gamboost(form_covary_flow,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial_Wet_Natives <-gamboost(form_covary_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_Wet_Natives <-mboost(form_tree_spatial,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # 
Model_tree_Wet_Natives <-mboost(form_tree,family = NBinomial(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # 

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

#st<-(mstop(cv_legacy_Wet_Natives))
#Model_legacy_Wet_Natives[st]

st<-(mstop(cv_spatial_Wet_Natives))
Model_spatial_Wet_Natives[st]

st<-(mstop(cv_covar_Wet_Natives))
Model_covary_Wet_Natives[st]

st<-(mstop(cv_covar_climate_Wet_Natives))
Model_covary_climate_Wet_Natives[st]

st<-(mstop(cv_covar_flow_Wet_Natives))
Model_covary_flow_Wet_Natives[st]

st<-(mstop(cv_covarspatial_Wet_Natives))
Model_covary_spatial_Wet_Natives[st]

st<-(mstop(cv_tree_spatial_Wet_Natives ))
Model_tree_spatial_Wet_Natives[st]

st<-(mstop(cv_tree_Wet_Natives ))
Model_tree_Wet_Natives[st]


###########################################################################################################################
# Evaluation of different models using multiplicity adjusted all-pairwise comparisons

library(multcomp)
library(lme4)
library(rlist)

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(Legacy)","(Spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro/spatial)","(Tree)","(Tree/spatial)")
tmp.abund.wet <- data.frame(cv = c(extrb(cv_legacy_Wet_Natives),extrb(cv_spatial_Wet_Natives), extrb(cv_covar_Wet_Natives), extrb(cv_covar_climate_Wet_Natives),extrb(cv_covar_flow_Wet_Natives),extrb(cv_covarspatial_Wet_Natives),
				 extrb(cv_tree_Wet_Natives),extrb(cv_tree_spatial_Wet_Natives)),
				model = gl(8, length(extrb(cv_spatial_Wet_Natives))),b = factor(rep(1:length(extrb(cv_spatial_Wet_Natives)), 8)))
				  

levels(tmp.abund.wet$model) <- nm


png("Wet_Natives NegLL comparing different models df1.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.abund.wet, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:7
		, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.abund.wet$model)), par("usr")[3] - 0.008,srt = 30, adj = 1, label = levels(tmp.abund.wet$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.abund.wet), tmp.abund.wet$b, function(x) lines(1:8, tmp.abund.wet[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

# Formal model comparison
sgmod_Wet_Natives <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.abund.wet), mcp(model = "Tukey")))


# stability selection of important parameters in additive model
sel_Wet_Natives <- stabsel(Model_covary_flow_Wet_Natives,cutoff=0.75, q = 5)

# Pseudo r2 but see below for resampled estimates

		m1 <- glm.nb(Wet_Natives ~ 1, data=daten) # null model	
		y <- daten$Wet_Natives
		n <- length(y)
		r2 <- ( 1 - exp( - 2/n * (logLik(Model_tree_Wet_Natives) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n})
		
		
Model_covary_flow_Wet_Natives$fitted
		
############################################################################################################################################################
# Evaluate fit of the 'best' model

datenSmall=daten # create copy of data
n=nrow(datenSmall)
set.seed(806)
predicted.in <- list()
model.coefs <-list()

N = 100 # this will be raised but kept to a minimum atm

for(i in 1:N) {

# Randomly permute the observations

indvecL <- sample(1:n, n, replace=FALSE) # random reordering
datenSmall <- daten[indvecL,][1:(n/2),] # subset to half data for test run


#### Run model
		Model_tree_Wet_Natives.test <-mboost(form_tree,family = NBinomial(),data = datenSmall, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # 
		cv5f <- cv(model.weights(Model_tree_Wet_Natives.test), type='subsampling', B=5)
		cv_tree_Wet_Natives.test<- cvrisk(Model_tree_Wet_Natives.test, folds=cv5f)
		st<-(mstop(cv_tree_Wet_Natives.test))
		Model_tree_Wet_Natives.test[st]
		mycoefs= coef(Model_tree_Wet_Natives.test) # store each model run coefficients incase we want to plot partial dependency plots for each test run to estimate ci
		model.coefs[[i]] = mycoefs
#### goodness of fit in-bootstrap
		m1 <- glm.nb(Wet_Natives ~ 1, data=datenSmall) # null model	
		y <- datenSmall$Wet_Natives
		n.y <- length(y)
		r2 <- ( 1 - exp( - 2/n.y * (logLik(Model_tree_Wet_Natives.test) - logLik(m1)) ) ) / ( 1 - exp(logLik(m1))^{2/n.y})
		predicted.in[[i]] = r2
}
		
r2.data=list.rbind(predicted.in)
mean(r2.data)
quantile(r2.data, c(.1, .5, .9)) 
	
#plot the predicted versus the observed. Its not wonderful but Jim tells me he has seen worse!

png("Wet_Natives predict versus observed df=1.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))

# Predictions for out-of-bootstrap data
predictions<-predict(Model_tree_Wet_Natives,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
plot(predictions[,1],daten$Wet_Natives)
abline(0,1)

dev.off()

############################################################################################################################################################

png(paste(image.dir,'Wetland_Native_marginal_plots_Dec2018.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

# plot using TSLW
mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatSmooth <- extract(Model_covary_flow_Wet_Natives,which="TSLW")
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_flow_Wet_Natives)$`bbs(TSLW, df = dfd, center = TRUE)`
plot(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l",xlab='log10(TSLW+1)', ylab='f(TSLW)')
rug(sort(envdata$TSLW+mTSLW))

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary_flow_Wet_Natives,which="d3Mon_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Wet_Natives, which=7)[[1]]
plot(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l",xlab='log10(d3Mon_meandepth+1)', ylab='f(d3Mon_meandepth)')
rug(sort(envdata$d3Mon_meandepth+md3Mon_meandepth))

md3yrs_meandepth<-mean((sorteddata$d3yrs_meandepth))
xmatLin<- extract(Model_covary_flow_Wet_Natives,which="d3yrs_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary_flow_Wet_Natives, which=15)[[1]]
plot(sort(envdata$d3yrs_meandepth+md3yrs_meandepth),yvalues[order(envdata$d3yrs_meandepth+md3yrs_meandepth)], type="l",xlab='d3yrs_meandepth', ylab='f(d3yrs_meandepth)')
rug(sort(envdata$d3yrs_meandepth+md3yrs_meandepth))

# plot using d1yrs_wet
md1yrs_wet<-mean((sorteddata$d1yrs_wet))
xmatSmooth <- extract(Model_covary_flow_Wet_Natives,which="d1yrs_wet")
yvalues=xmatSmooth[[2]]%*%coef(Model_covary_flow_Wet_Natives)$`bbs(d1yrs_wet, df = dfd, center = TRUE)`
plot(sort(envdata$d1yrs_wet+md1yrs_wet),yvalues[order(envdata$d1yrs_wet+md1yrs_wet)], type="l",xlab='d1yrs_wet', ylab='f(d1yrs_wet)')
rug(sort(envdata$d1yrs_wet+md1yrs_wet))

dev.off()

############################################################################################################################################################
# Explore interactions using boosted regression trees. Note that I had to assume poisson as there is no nbinom for this method. 

library(dismo)
library(gbm)

mypreds=daten[,c("Wet_Natives","TSLW","FF","d3Mon_wet","d3Mon_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth","Freq_d1","Freq_d3","d90","d365","MeanTemp90","MinTemp90","MaxTemp90")]

Wet_Nats_dismo <- gbm.step(data=mypreds, gbm.x=2:16, gbm.y=1, family="poisson", tree.complexity=5,learning.rate=0.01, bag.fraction=0.5)
#Wet_Nats_simple <-gbm.simplify(Wet_Nats_dismo,n.drop=5)

Wet_Nats_dismo.fitvalues <- Wet_Nats_dismo$fitted # extracts the fitted values
Wet_Nats_dismo.residuals <- Wet_Nats_dismo$residuals # extracts the residuals

plot(Wet_Nats_dismo.fitvalues, Wet_Nats_dismo.residuals)

png(paste(image.dir,'Wetland_Native_marginal_plots_BRT Dec2018.png',sep=''), width=2000, height=2000, units="px", res=300)
gbm.plot(Wet_Nats_dismo,n.plot=12,write.title=FALSE)
dev.off()

find.int<-gbm.interactions(Wet_Nats_dismo)
find.int$interactions
find.int$rank.list

png(paste(image.dir,'Wetland_Native_marginal_plots_BRT interactions Dec2018.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mar = c(0.5, 0.5, 0.5, 0.5),mfrow=c(3,3),cex=0.5)

gbm.perspec(Wet_Nats_dismo, 4,2)
gbm.perspec(Wet_Nats_dismo, 13,1)
gbm.perspec(Wet_Nats_dismo, 7,4)
gbm.perspec(Wet_Nats_dismo, 11,8)
gbm.perspec(Wet_Nats_dismo, 15,2)
gbm.perspec(Wet_Nats_dismo, 13,2)
gbm.perspec(Wet_Nats_dismo, 7,2)
gbm.perspec(Wet_Nats_dismo, 6,1)
gbm.perspec(Wet_Nats_dismo, 8,2)

dev.off()