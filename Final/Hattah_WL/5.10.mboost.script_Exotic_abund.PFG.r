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
envdata$ExoticAbund=envdata$Terr_Exotics+envdata$Wet_Exotics

envdata$ExoticAbund[envdata$ExoticAbund > 0] <- 1 


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

form_spatial <- ExoticAbund~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=df,center=TRUE,differences=1)+bols(interc,intercept=FALSE,df=df)

			
form_covary <- 	ExoticAbund~
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
				
form_covary_flow <- 	ExoticAbund~
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

				
form_covary_climate <- 	ExoticAbund~
				bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=df)+
				bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=df)+
				bols(MeanTemp90, intercept=FALSE)+bbs(MeanTemp90, center=TRUE, df=df)+
				bols(MinTemp90, intercept=FALSE)+bbs(MinTemp90, center=TRUE, df=df)+
				bols(MaxTemp90, intercept=FALSE)+bbs(MaxTemp90, center=TRUE, df=df)+
				bols(interc,intercept=FALSE,df=df)
				
				
form_covary_spatial <- ExoticAbund~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=df,center=TRUE,differences=1)+
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

				
form_tree_spatial <- ExoticAbund~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=df,center=TRUE,differences=1)+
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE,df=df)

				
daten=envdata

# Run models
Model_spatial_ExoticAbund <-mboost(form_spatial,family = Poisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_ExoticAbund<-gamboost(form_covary,family = Poisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate_ExoticAbund<-gamboost(form_covary_climate,family = Poisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow_ExoticAbund<-gamboost(form_covary_flow,family = Poisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial_ExoticAbund <-gamboost(form_covary_spatial,family = Poisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial_ExoticAbund <-mboost(form_tree_spatial,family = Poisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # 

###

cv5f <- cv(model.weights(Model_spatial_ExoticAbund), type='subsampling', B=25)
cv_spatial_ExoticAbund <- cvrisk(Model_spatial_ExoticAbund, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_ExoticAbund), type='subsampling', B=25)
cv_covar_ExoticAbund <- cvrisk(Model_covary_ExoticAbund, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_climate_ExoticAbund), type='subsampling', B=25)
cv_covar_climate_ExoticAbund <- cvrisk(Model_covary_climate_ExoticAbund, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_flow_ExoticAbund), type='subsampling', B=25)
cv_covar_flow_ExoticAbund<- cvrisk(Model_covary_flow_ExoticAbund, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_spatial_ExoticAbund), type='subsampling', B=25)
cv_covarspatial_ExoticAbund<- cvrisk(Model_covary_spatial_ExoticAbund, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial_ExoticAbund), type='subsampling', B=25)
cv_tree_spatial_ExoticAbund<- cvrisk(Model_tree_spatial_ExoticAbund, folds=cv5f)


st<-(mstop(cv_spatial_ExoticAbund))
Model_spatial_ExoticAbund[st]

st<-(mstop(cv_covar_ExoticAbund))
Model_covary_ExoticAbund[st]

st<-(mstop(cv_covar_climate_ExoticAbund))
Model_covary_climate_ExoticAbund[st]

st<-(mstop(cv_covar_flow_ExoticAbund))
Model_covary_flow_ExoticAbund[st]

st<-(mstop(cv_covarspatial_ExoticAbund))
Model_covary_spatial_ExoticAbund[st]

st<-(mstop(cv_tree_spatial_ExoticAbund ))
Model_tree_spatial_ExoticAbund[st]


####################################################################################
##### Evaluation of different models using multiplicity adjusted all-pairwise comparisons
library("multcomp")
library("lme4")

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro/spatial)","(tree/spatial)")
tmp.ExoticAbund <- data.frame(cv = c(extrb(cv_spatial_ExoticAbund ), extrb(cv_covar_ExoticAbund ), extrb(cv_covar_climate_ExoticAbund ),extrb(cv_covar_flow_ExoticAbund ),extrb(cv_covarspatial_ExoticAbund ),
				 extrb(cv_tree_spatial_ExoticAbund )),
				model = gl(6, length(extrb(cv_spatial_ExoticAbund ))),b = factor(rep(1:length(extrb(cv_spatial_ExoticAbund )), 6)))
				  

levels(tmp.ExoticAbund$model) <- nm
		
		
png("Exotic Abundance NegLL comparing different models.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
        boxplot(cv ~  model, data = tmp.ExoticAbund, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp.ExoticAbund$model)), par("usr")[3] - 0.09,srt = 30, adj = 1, label = levels(tmp.ExoticAbund$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp.ExoticAbund), tmp.ExoticAbund$b, function(x) lines(1:6, tmp.ExoticAbund[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()

sgmod_ExoticAbund <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp.ExoticAbund), mcp(model = "Tukey")))


### stability selection
sel.ExoticAbund <- stabsel(Model_covary_flow_ExoticAbund,cutoff=0.75, q = 10)

###

labels=c("bols(d3yrs_wet)","bols(d3Mon_meandepth)","bbs(TSLW)","bbs(d3Mon_meandepth)", "bols(TSLW)", "bbs(d1yrs_meandepth)", "bbs(d1yrs_wet)",
"bbs(Freq_d1)","bbs(MinTemp90)", "bols(MaxTemp90)","bols(d90)","bols(Freq_d1)", "bols(MeanTemp90)","bols(Freq_d3)","bols(FF)","bols(MinTemp90)",
"bbs(Freq_d3)","bols(d1yrs_wet)","bbs(d3yrs_wet)","bols(d365)","bols(d3yrs_meandepth)","bbs(FF)","bbs(MaxTemp90)","bbs(d90)","bbs(d3yrs_meandepth)",
"bols(d3Mon_wet)","bols(interc)","bbs(MeanTemp90)","bbs(d365)","bols(d1yrs_meandepth)","bbs(d3Mon_wet)")

selmax=(as.data.frame(sel$max))
selmax=cbind(selmax,selmax)
colnames(selmax)[1]="Stable.selection"
selmax=selmax[order(selmax$Stable.selection),]

ggplot(selmax, aes(Stable.selection,rownames(selmax))) +
geom_segment(aes(x = 0, y=rownames(selmax), xend = Stable.selection, yend=rownames(selmax)), color = "grey50") +
        geom_point()+ ylab(NULL)
    
		
png("Wet_Natives stability selection plot.png",width=35, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
ggplot(selmax, aes(Stable.selection,rownames(selmax))) +
geom_segment(aes(x = 0, y=rownames(selmax), xend = Stable.selection, yend=rownames(selmax)), color = "grey50") +
        geom_point()+ ylab(NULL)
		
		
dev.off()

##################################################################################################################
# removed poorer predictor

form_covary <- 	Wet_Natives~
				bbs(TSLW, center=TRUE, df=df)+
				bols(FF, intercept=FALSE)+
				bols(d3Mon_wet, intercept=FALSE)+
				bols(d3Mon_meandepth, intercept=FALSE)+
				bbs(d1yrs_wet, center=TRUE, df=df)+
				bbs(d1yrs_meandepth, center=TRUE, df=df)+
				bols(d3yrs_wet, intercept=FALSE)+
				bols(d3yrs_meandepth, intercept=FALSE)+
				bbs(Freq_d1, center=TRUE, df=df)+
				bbs(Freq_d3, center=TRUE, df=df)+
				bols(d90, intercept=FALSE)+
				bols(d365, intercept=FALSE)+
				bols(MeanTemp90, intercept=FALSE)+
				bbs(MinTemp90, center=TRUE, df=df)+
				bols(MaxTemp90, intercept=FALSE)+
				bols(interc,intercept=FALSE,df=df)
				
Model_covary <-gamboost(form_covary,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) 

cv5f <- cv(model.weights(Model_covary), type='subsampling', B=25)
cv_covar <- cvrisk(Model_covary, folds=cv5f)

st<-(mstop(cv_covar))
Model_covary[st]

############################################################################################################################################################

png(paste(image.dir,'Wetland_Native_Abund_marginal_plots_Dec2018.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

# plot using TSLW
mTSLW<-mean(log10(sorteddata$TSLW+1))
xmatSmooth <- extract(Model_covary,which="TSLW")
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Model_covary)$`bbs(TSLW, df = df, center = TRUE)`
plot(sort(envdata$TSLW+mTSLW),yvalues[order(envdata$TSLW+mTSLW)], type="l",xlab='log10(TSLW+1)', ylab='f(TSLW)')
rug(sort(envdata$TSLW+mTSLW))

md3Mon_meandepth<-mean(log10(sorteddata$d3Mon_meandepth+1))
xmatLin <- extract(Model_covary,which="d3Mon_meandepth")
yvalues=xmatLin[[1]]%*%coef(Model_covary, which=4)[[1]]
plot(sort(envdata$d3Mon_meandepth+md3Mon_meandepth),yvalues[order(envdata$d3Mon_meandepth+md3Mon_meandepth)], type="l",xlab='log10(d3Mon_meandepth+1)', ylab='f(d3Mon_meandepth)')
rug(sort(envdata$d3Mon_meandepth+md3Mon_meandepth))

md3yrs_wet<-mean((sorteddata$d3yrs_wet))
xmatLin<- extract(Model_covary,which="d3yrs_wet")
yvalues=xmatLin[[1]]%*%coef(Model_covary, which=7)[[1]]
plot(sort(envdata$d3yrs_wet+md3yrs_wet),yvalues[order(envdata$d3yrs_wet+md3yrs_wet)], type="l",xlab='d3yrs_wet', ylab='f(d3yrs_wet)')
rug(sort(envdata$d3yrs_wet+md3yrs_wet))

# plot using d1yrs_wet
md1yrs_wet<-mean((sorteddata$d1yrs_wet))
xmatSmooth <- extract(Model_covary,which="d1yrs_wet")
yvalues=xmatSmooth[[1]]%*%coef(Model_covary)$`bbs(d1yrs_wet, df = df, center = TRUE)`
plot(sort(envdata$d1yrs_wet+md1yrs_wet),yvalues[order(envdata$d1yrs_wet+md1yrs_wet)], type="l",xlab='d1yrs_wet', ylab='f(d1yrs_wet)')
rug(sort(envdata$d1yrs_wet+md1yrs_wet))

dev.off()
