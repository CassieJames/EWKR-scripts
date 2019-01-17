##############################################################################
#### Model richness and abundance of wetland species at Hattah Lakes wetlands
#### C James October 2018 

library(corrplot)
library(mboost)
library(countreg)
library(partykit)


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"
envdata=data.matrix.env.data # copy data
envdata$FF=envdata$Flood_frequency+envdata$Freq_ALL
 
envdata=subset(envdata, envdata$Terr_Natives>0) # just model positive counts

# sort out environmental data - centre and log
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond record
envdata$TSLW[envdata$Innundated==TRUE]<-1 # if site is recorded as inundated change TSLW to 1 day
envdata$d3Mon_wet[is.na(envdata$d3Mon_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d1yrs_wet[is.na(envdata$d1yrs_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d3yrs_wet[is.na(envdata$d3yrs_wet)]<-0 # if site us never wet during time period replace with 0
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
envdata$Freq_d1=as.numeric(scale((envdata$Freq_d1),center=TRUE, scale=FALSE))
envdata$Freq_d3=as.numeric(scale((envdata$Freq_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d1=as.numeric(scale((envdata$CTF_prop_d1),center=TRUE, scale=FALSE))
envdata$CTF_prop_d3=as.numeric(scale((envdata$CTF_prop_d3),center=TRUE, scale=FALSE))
envdata$CTF_prop_d5=as.numeric(scale((envdata$CTF_prop_d5),center=TRUE, scale=FALSE))
envdata$FF=as.numeric(scale((envdata$FF),center=TRUE, scale=FALSE))

envdata$WaterYr=as.numeric(envdata$WaterYr)


envdata <- as.data.frame(cbind(interc=1, envdata)) # and a column vector of 1's for the intercept


####################################################################################################################################################				


df=1

form_spatial <- Terr_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=df,center=TRUE,differences=1)+bols(interc,intercept=FALSE,df=df)

			
form_covary <- 	Terr_Natives~
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=df)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=df)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=df)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=df)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=df)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=df)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=df)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=df)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=df)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=df)+
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=df)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=df)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=df)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=df)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=df)+
				bols(interc,intercept=FALSE,df=df)
				
form_covary_flow <- 	Terr_Natives~
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=df)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=df)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=df)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=df)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=df)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=df)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=df)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=df)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=df)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=df)+bols(interc,intercept=FALSE,df=df)

				
form_covary_climate <- 	Terr_Natives~
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=df)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=df)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=df)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=df)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=df)+
				bols(interc,intercept=FALSE,df=df)
				
				
form_covary_spatial <- Terr_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=df)+
				bols(FF, intercept=FALSE)+bbs(FF, center=TRUE, df=df)+
				bols(d3Mon_wet, intercept=FALSE)+bbs(d3Mon_wet, center=TRUE, df=df)+
				bols(d3Mon_meandepth, intercept=FALSE)+bbs(d3Mon_meandepth, center=TRUE, df=df)+
				bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=df)+
				bols(d1yrs_meandepth, intercept=FALSE)+bbs(d1yrs_meandepth, center=TRUE, df=df)+
				bols(d3yrs_wet, intercept=FALSE)+bbs(d3yrs_wet, center=TRUE, df=df)+
				bols(d3yrs_meandepth, intercept=FALSE)+bbs(d3yrs_meandepth, center=TRUE, df=df)+
				bols(Freq_d1, intercept=FALSE)+bbs(Freq_d1, center=TRUE, df=df)+
				bols(Freq_d3, intercept=FALSE)+bbs(Freq_d3, center=TRUE, df=df)+
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=df)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=df)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=df)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=df)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=df)+
				bols(interc,intercept=FALSE,df=df)			
							

tctrl = partykit::ctree_control(stump = FALSE) 

				
form_tree_spatial <- Terr_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE,df=df)

				
daten=envdata

# Run models
Model_spatial_Terr_Natives <-mboost(form_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_Terr_Natives <-gamboost(form_covary,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_Terr_Natives<-gamboost(form_covary_climate,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow_Terr_Natives<-gamboost(form_covary_flow,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial_Terr_Natives <-gamboost(form_covary_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_Terr_Natives <-mboost(form_tree_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # 

###

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


####################################################################################
##### Evaluation of different models using multiplicity adjusted all-pairwise comparisons
library("multcomp")
library("lme4")

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro/spatial)","(tree/spatial)")
tmp.ter <- data.frame(cv = c(extrb(cv_spatial_Terr_Natives), extrb(cv_covar_Terr_Natives), extrb(cv_covar_climate_Terr_Natives),extrb(cv_covar_flow_Terr_Natives),extrb(cv_covarspatial_Terr_Natives),
				 extrb(cv_tree_spatial_Terr_Natives)),
				model = gl(6, length(extrb(cv_spatial_Terr_Natives))),b = factor(rep(1:length(extrb(cv_spatial_Terr_Natives)), 6)))
				  

levels(tmp.ter$model) <- nm


png("Terr_Natives NegLL comparing different models.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp.ter, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.ter$model)), par("usr")[3] - 0.09,srt = 30, adj = 1, label = levels(tmp.ter$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.ter), tmp.ter$b, function(x) lines(1:6, tmp.ter[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

sgmod_Terr_Natives <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.ter), mcp(model = "Tukey")))

### conclusion is that spatial model is best model with models including covarites performing worse (significantly so in the case of the hydrological variables)
### Spatial effects had higher selection probability than spatio-temporal terms or intercept
### stability selection
sel_Terr_Natives <- stabsel(Model_spatial_Terr_Natives,cutoff=0.75, q = 1)



