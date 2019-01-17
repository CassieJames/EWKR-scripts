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

################################################################################################
#specify model formulas with df=3

form_spatial <- Wet_Natives~ brandom(Site.ID.x,df=3)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)

form_covariates <- 	Wet_Natives~brandom(Site.ID.x,df=3)+
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
				bols(interc,intercept=FALSE)
				
				
				
form_covary_flow <- Wet_Natives~brandom(Site.ID.x,df=3)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(d3Mon_wet, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, center=TRUE, df=3)+
				bbs(d1yrs_wet, center=TRUE, df=3)+
				bbs(d1yrs_meandepth, center=TRUE, df=3)+
				bbs(d3yrs_wet, center=TRUE, df=3)+
				bbs(d3yrs_meandepth, center=TRUE, df=3)+
				bbs(Freq_d1, center=TRUE, df=3)+
				bbs(Freq_d3, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)
				
form_covary_climate <- 	Wet_Natives~brandom(Site.ID.x,df=3)+
				bbs(d90, center=TRUE, df=3)+
				bbs(d365, center=TRUE, df=3)+
				bbs(MeanTemp90, center=TRUE, df=3)+
				bbs(MinTemp90, center=TRUE, df=3)+
				bbs(MaxTemp90, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)
				

				
form_covary_flow_d90 <- Wet_Natives~brandom(Site.ID.x,df=3)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(d3Mon_wet, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, center=TRUE, df=3)+
				bbs(d1yrs_wet, center=TRUE, df=3)+
				bbs(d1yrs_meandepth, center=TRUE, df=3)+
				bbs(d3yrs_wet, center=TRUE, df=3)+
				bbs(d3yrs_meandepth, center=TRUE, df=3)+
				bbs(Freq_d1, center=TRUE, df=3)+
				bbs(Freq_d3, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)+
				bbs(TSLW,by=d90, center=TRUE, df=3)+
				bbs(d3Mon_wet,by=d90, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, by=d90,center=TRUE, df=3)+
				bbs(d1yrs_wet, by=d90,center=TRUE, df=3)+
				bbs(d1yrs_meandepth, by=d90,center=TRUE, df=3)+
				bbs(d3yrs_wet, by=d90,center=TRUE, df=3)+
				bbs(d3yrs_meandepth, by=d90,center=TRUE, df=3)+
				bbs(Freq_d1, by=d90,center=TRUE, df=3)+
				bbs(Freq_d3, by=d90,center=TRUE, df=3)+
				bols(interc,intercept=FALSE)
				
form_covary_flow_d90_spatial <- Wet_Natives~brandom(Site.ID.x,df=3)+
				bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(d3Mon_wet, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, center=TRUE, df=3)+
				bbs(d1yrs_wet, center=TRUE, df=3)+
				bbs(d1yrs_meandepth, center=TRUE, df=3)+
				bbs(d3yrs_wet, center=TRUE, df=3)+
				bbs(d3yrs_meandepth, center=TRUE, df=3)+
				bbs(Freq_d1, center=TRUE, df=3)+
				bbs(Freq_d3, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)+
				bbs(TSLW,by=d90, center=TRUE, df=3)+
				bbs(d3Mon_wet,by=d90, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, by=d90,center=TRUE, df=3)+
				bbs(d1yrs_wet, by=d90,center=TRUE, df=3)+
				bbs(d1yrs_meandepth, by=d90,center=TRUE, df=3)+
				bbs(d3yrs_wet, by=d90,center=TRUE, df=3)+
				bbs(d3yrs_meandepth, by=d90,center=TRUE, df=3)+
				bbs(Freq_d1, by=d90,center=TRUE, df=3)+
				bbs(Freq_d3, by=d90,center=TRUE, df=3)+
				bols(interc,intercept=FALSE)				

form_covary_flow_d90_vary <- 	Wet_Natives~brandom(Site.ID.x,df=3)+
				bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=FF,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(d3Mon_wet, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, center=TRUE, df=3)+
				bbs(d1yrs_wet, center=TRUE, df=3)+
				bbs(d1yrs_meandepth, center=TRUE, df=3)+
				bbs(d3yrs_wet, center=TRUE, df=3)+
				bbs(d3yrs_meandepth, center=TRUE, df=3)+
				bbs(Freq_d1, center=TRUE, df=3)+
				bbs(Freq_d3, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)+
				bbs(TSLW,by=d90, center=TRUE, df=3)+
				bbs(d3Mon_wet,by=d90, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, by=d90,center=TRUE, df=3)+
				bbs(d1yrs_wet, by=d90,center=TRUE, df=3)+
				bbs(d1yrs_meandepth, by=d90,center=TRUE, df=3)+
				bbs(d3yrs_wet, by=d90,center=TRUE, df=3)+
				bbs(d3yrs_meandepth, by=d90,center=TRUE, df=3)+
				bbs(Freq_d1, by=d90,center=TRUE, df=3)+
				bbs(Freq_d3, by=d90,center=TRUE, df=3)+
				bols(interc,intercept=FALSE)
				

form_covary_spatial <- Wet_Natives~ brandom(Site.ID.x,df=3)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
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
				bols(interc,intercept=FALSE)	

form_hydro_spatial <- Wet_Natives~ brandom(Site.ID.x,df=3)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(d3Mon_wet, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, center=TRUE, df=3)+
				bbs(d1yrs_wet, center=TRUE, df=3)+
				bbs(d1yrs_meandepth, center=TRUE, df=3)+
				bbs(d3yrs_wet, center=TRUE, df=3)+
				bbs(d3yrs_meandepth, center=TRUE, df=3)+
				bbs(Freq_d1, center=TRUE, df=3)+
				bbs(Freq_d3, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)	
				

tctrl = partykit::ctree_control(stump = FALSE) 
				
form_tree <- Wet_Natives~brandom(Site.ID.x,df=3)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

				
form_tree_spatial <- Wet_Natives~brandom(Site.ID.x,df=3)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

											

form_vary<- Wet_Natives~ brandom(Site.ID.x,df=3)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
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
				bols(interc,intercept=FALSE)


form_tree_spatial_vary <- Wet_Natives~brandom(Site.ID.x,df=3)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=FF,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				
####################################################################################################################################################				
#specify model formulas with df=1

form_spatial <- Wet_Natives~ brandom(Site.ID.x,df=1)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)

form_covary <- 	Wet_Natives~brandom(Site.ID.x,df=1)+
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
				bols(interc,intercept=FALSE)
				
form_covary_flow <- 	Wet_Natives~brandom(Site.ID.x,df=1)+
				bbs(TSLW, center=TRUE, df=1)+
				bbs(d3Mon_wet, center=TRUE, df=1)+
				bbs(d3Mon_meandepth, center=TRUE, df=1)+
				bbs(d1yrs_wet, center=TRUE, df=1)+
				bbs(d1yrs_meandepth, center=TRUE, df=1)+
				bbs(d3yrs_wet, center=TRUE, df=1)+
				bbs(d3yrs_meandepth, center=TRUE, df=1)+
				bbs(Freq_d1, center=TRUE, df=1)+
				bbs(Freq_d3, center=TRUE, df=1)+
				bols(interc,intercept=FALSE)
				
form_covary_climate <- 	Wet_Natives~brandom(Site.ID.x,df=1)+
				bbs(d90, center=TRUE, df=1)+
				bbs(d365, center=TRUE, df=1)+
				bbs(MeanTemp90, center=TRUE, df=1)+
				bbs(MinTemp90, center=TRUE, df=1)+
				bbs(MaxTemp90, center=TRUE, df=1)+
				bols(interc,intercept=FALSE)
				
				
form_covary_spatial <- Wet_Natives~ brandom(Site.ID.x,df=1)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
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
				bols(interc,intercept=FALSE)				

form_hydro_spatial <- Wet_Natives~ brandom(Site.ID.x,df=1)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=1)+
				bbs(d3Mon_wet, center=TRUE, df=1)+
				bbs(d3Mon_meandepth, center=TRUE, df=1)+
				bbs(d1yrs_wet, center=TRUE, df=1)+
				bbs(d1yrs_meandepth, center=TRUE, df=1)+
				bbs(d3yrs_wet, center=TRUE, df=1)+
				bbs(d3yrs_meandepth, center=TRUE, df=1)+
				bbs(Freq_d1, center=TRUE, df=1)+
				bbs(Freq_d3, center=TRUE, df=1)
							

tctrl = partykit::ctree_control(stump = FALSE) 
				
form_tree <- Wet_Natives~brandom(Site.ID.x,df=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

				
form_tree_spatial <- Wet_Natives~brandom(Site.ID.x,df=1)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree_spatial_hydro <- Wet_Natives~brandom(Site.ID.x,df=1)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,tree_controls=tctrl)+bols(interc,intercept=FALSE)
											
				
form_vary<- Wet_Natives~ brandom(Site.ID.x,df=1)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=FF,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
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
				bols(interc,intercept=FALSE)


form_tree_spatial_vary <- Wet_Natives~brandom(Site.ID.x,df=1)+bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=FF,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				
daten=envdata

# Run models
Model_spatial <-mboost(form_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary <-gamboost(form_covary,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate<-gamboost(form_covary_climate,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow<-gamboost(form_covary_flow,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial <-gamboost(form_covary_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_hydro_spatial <-gamboost(form_hydro_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 100
Model_tree <-mboost(form_tree,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) 
Model_tree_spatial <-mboost(form_tree_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # 
Model_tree_spatial_hydro <-mboost(form_tree_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # reduced down the mstop value for trial runs as it takes a while
Model_vary <- mboost(form_vary,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01))
Model_tree_spatial_vary <- mboost(form_vary,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01))

###

cv5f <- cv(model.weights(Model_spatial), type='subsampling', B=25)
cv_spatial <- cvrisk(Model_spatial, folds=cv5f)

cv5f <- cv(model.weights(Model_covary), type='subsampling', B=25)
cv_covar <- cvrisk(Model_covary, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_climate), type='subsampling', B=25)
cv_covar_climate <- cvrisk(Model_covary_climate, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_flow), type='subsampling', B=25)
cv_covar_flow<- cvrisk(Model_covary_flow, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_spatial), type='subsampling', B=25)
cv_covarspatial<- cvrisk(Model_covary_spatial, folds=cv5f)

cv5f <- cv(model.weights(Model_hydro_spatial), type='subsampling', B=25)
cv_hydrospatial<- cvrisk(Model_hydro_spatial, folds=cv5f)

cv5f <- cv(model.weights(Model_tree), type='subsampling', B=25)
cv_tree<- cvrisk(Model_tree, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial), type='subsampling', B=25)
cv_tree_spatial<- cvrisk(Model_tree_spatial, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial_hydro), type='subsampling', B=25)
cv_tree_spatial_hydro<- cvrisk(Model_tree_spatial_hydro, folds=cv5f)

cv5f <- cv(model.weights(Model_vary), type='subsampling', B=25)
cv_vary<- cvrisk(Model_vary, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial_vary), type='subsampling', B=25)
cv_varytree<- cvrisk(Model_tree_spatial_vary, folds=cv5f)



st<-(mstop(cv_spatial))
Model_spatial[st]

st<-(mstop(cv_covar))
Model_covary[st]

st<-(mstop(cv_covar_climate))
Model_covary_climate[st]

st<-(mstop(cv_covar_flow))
Model_covary_flow[st]

st<-(mstop(cv_covarspatial))
Model_covary_spatial[st]

st<-(mstop(cv_hydrospatial))
Model_hydro_spatial[st]

st<-(mstop(cv_tree))
Model_tree[st]

st<-(mstop(cv_tree_spatial ))
Model_tree_spatial[st]

st<-(mstop(cv_tree_spatial_hydro))
Model_tree_spatial_hydro[st]

st<-(mstop(cv_vary))
Model_vary[st]

st<-(mstop(cv_varytree))
Model_tree_spatial_vary[st]





####################################################################################
##### Evaluation of different models using multiplicity adjusted all-pairwise comparisons
library("multcomp")
library("lme4")

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro/spatial)","(hydro/spatial)","(Climate+Hydro/vary)","(tree)","(tree/spatial)", 
        "tree/hydro/spatial", "(tree/vary)")
tmp <- data.frame(cv = c(extrb(cv_spatial), extrb(cv_covar), extrb(cv_covar_climate),extrb(cv_covar_flow),extrb(cv_covarspatial), extrb(cv_hydrospatial),extrb(cv_vary),
				extrb(cv_tree), extrb(cv_tree_spatial),extrb(cv_tree_spatial_hydro),extrb(cv_varytree)),
				model = gl(11, length(extrb(cv_spatial))),b = factor(rep(1:length(extrb(cv_spatial)), 11)))
				  

levels(tmp$model) <- nm



png("Wet_Natives NegLL comparing different models.png",width=25, height=25, units='cm', res=300, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),cex=1,oma=c(2,0,1,0.5))
		
boxplot(cv ~  model, data = tmp, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
out <- tapply(1:nrow(tmp), tmp$b, function(x) lines(1:11, tmp[x,"cv"], col = rgb(0,0,0,0.1)))

axis(1, at = 1:11, label = levels(tmp$model), tick = FALSE, las = 3, cex.axis = 0.85)
axis(2)

sgmod <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp), mcp(model = "Tukey")))

dev.off()

####################################################################################################################################################				
#specify model formulas with df=1 with random intercept excluded

form_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)

form_covary <- 	Wet_Natives~
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
				bols(interc,intercept=FALSE)
				
form_covary_flow <- 	Wet_Natives~
				bbs(TSLW, center=TRUE, df=1)+
				bbs(d3Mon_wet, center=TRUE, df=1)+
				bbs(d3Mon_meandepth, center=TRUE, df=1)+
				bbs(d1yrs_wet, center=TRUE, df=1)+
				bbs(d1yrs_meandepth, center=TRUE, df=1)+
				bbs(d3yrs_wet, center=TRUE, df=1)+
				bbs(d3yrs_meandepth, center=TRUE, df=1)+
				bbs(Freq_d1, center=TRUE, df=1)+
				bbs(Freq_d3, center=TRUE, df=1)+
				bols(interc,intercept=FALSE)
				
form_covary_climate <- 	Wet_Natives~
				bbs(d90, center=TRUE, df=1)+
				bbs(d365, center=TRUE, df=1)+
				bbs(MeanTemp90, center=TRUE, df=1)+
				bbs(MinTemp90, center=TRUE, df=1)+
				bbs(MaxTemp90, center=TRUE, df=1)+
				bols(interc,intercept=FALSE)
				
				
form_covary_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
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
				bols(interc,intercept=FALSE)				

form_hydro_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=1)+
				bbs(d3Mon_wet, center=TRUE, df=1)+
				bbs(d3Mon_meandepth, center=TRUE, df=1)+
				bbs(d1yrs_wet, center=TRUE, df=1)+
				bbs(d1yrs_meandepth, center=TRUE, df=1)+
				bbs(d3yrs_wet, center=TRUE, df=1)+
				bbs(d3yrs_meandepth, center=TRUE, df=1)+
				bbs(Freq_d1, center=TRUE, df=1)+
				bbs(Freq_d3, center=TRUE, df=1)
							

tctrl = partykit::ctree_control(stump = FALSE) 
				
form_tree <- Wet_Natives~
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

				
form_tree_spatial <- Wet_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree_spatial_hydro <- Wet_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,tree_controls=tctrl)+bols(interc,intercept=FALSE)
											
				
form_vary<- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=FF,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
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
				bols(interc,intercept=FALSE)


form_tree_spatial_vary <- Wet_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=FF,boundary.knots=NULL,df=1,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				
daten=envdata

####################################################################################################################################################				
#specify model formulas with df=3 with random intercept excluded

form_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)

form_covary <- 	Wet_Natives~
				bbs(TSLW, center=TRUE, df=3)+
				bbs(FF, center=TRUE, df=3)+
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
				bols(interc,intercept=FALSE)
				
form_covary_flow <- 	Wet_Natives~
				bbs(TSLW, center=TRUE, df=3)+
				bbs(FF, center=TRUE, df=3)+
				bbs(d3Mon_wet, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, center=TRUE, df=3)+
				bbs(d1yrs_wet, center=TRUE, df=3)+
				bbs(d1yrs_meandepth, center=TRUE, df=3)+
				bbs(d3yrs_wet, center=TRUE, df=3)+
				bbs(d3yrs_meandepth, center=TRUE, df=3)+
				bbs(Freq_d1, center=TRUE, df=3)+
				bbs(Freq_d3, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)
				
form_covary_climate <- 	Wet_Natives~
				bbs(d90, center=TRUE, df=3)+
				bbs(d365, center=TRUE, df=3)+
				bbs(MeanTemp90, center=TRUE, df=3)+
				bbs(MinTemp90, center=TRUE, df=3)+
				bbs(MaxTemp90, center=TRUE, df=3)+
				bols(interc,intercept=FALSE)
				
				
form_covary_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(FF, center=TRUE, df=3)+
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
				bols(interc,intercept=FALSE)				

form_hydro_spatial <- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(FF, center=TRUE, df=3)+
				bbs(d3Mon_wet, center=TRUE, df=3)+
				bbs(d3Mon_meandepth, center=TRUE, df=3)+
				bbs(d1yrs_wet, center=TRUE, df=3)+
				bbs(d1yrs_meandepth, center=TRUE, df=3)+
				bbs(d3yrs_wet, center=TRUE, df=3)+
				bbs(d3yrs_meandepth, center=TRUE, df=3)+
				bbs(Freq_d1, center=TRUE, df=3)+
				bbs(Freq_d3, center=TRUE, df=3)
							

tctrl = partykit::ctree_control(stump = FALSE) 
				
form_tree <- Wet_Natives~
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

				
form_tree_spatial <- Wet_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)

form_tree_spatial_hydro <- Wet_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by=WaterYr,boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				btree(TSLW,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,tree_controls=tctrl)+bols(interc,intercept=FALSE)
											
				
form_vary<- Wet_Natives~ bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bbs(TSLW, center=TRUE, df=3)+
				bbs(FF, center=TRUE, df=3)+
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
				bols(interc,intercept=FALSE)


form_tree_spatial_vary <- Wet_Natives~bspatial(Easting, Northing, knots = 6, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				bspatial(Easting, Northing, knots = 6, by= WaterYr, boundary.knots=NULL,df=3,center=TRUE,differences=1)+
				btree(TSLW,FF,d3Mon_wet,d3Mon_meandepth,d1yrs_wet,d1yrs_meandepth,d3yrs_wet,d3yrs_meandepth,Freq_d1,Freq_d3,d90,d365,MeanTemp90,MinTemp90,MaxTemp90,tree_controls=tctrl)+
				bols(interc,intercept=FALSE)
				
daten=envdata

# Run models
Model_spatial <-mboost(form_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary <-gamboost(form_covary,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_covary_climate<-gamboost(form_covary_climate,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_flow<-gamboost(form_covary_flow,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # orig
Model_covary_spatial <-gamboost(form_covary_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE, nu=0.01)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
Model_tree_spatial <-mboost(form_tree_spatial,family = MBztpoisson(),data = daten, control=boost_control(mstop=1000,trace=TRUE,nu=0.01)) # 

###

cv5f <- cv(model.weights(Model_spatial), type='subsampling', B=25)
cv_spatial <- cvrisk(Model_spatial, folds=cv5f)

cv5f <- cv(model.weights(Model_covary), type='subsampling', B=25)
cv_covar <- cvrisk(Model_covary, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_climate), type='subsampling', B=25)
cv_covar_climate <- cvrisk(Model_covary_climate, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_flow), type='subsampling', B=25)
cv_covar_flow<- cvrisk(Model_covary_flow, folds=cv5f)

cv5f <- cv(model.weights(Model_covary_spatial), type='subsampling', B=25)
cv_covarspatial<- cvrisk(Model_covary_spatial, folds=cv5f)

cv5f <- cv(model.weights(Model_tree_spatial), type='subsampling', B=25)
cv_tree_spatial<- cvrisk(Model_tree_spatial, folds=cv5f)


st<-(mstop(cv_spatial))
Model_spatial[st]

st<-(mstop(cv_covar))
Model_covary[st]

st<-(mstop(cv_covar_climate))
Model_covary_climate[st]

st<-(mstop(cv_covar_flow))
Model_covary_flow[st]

st<-(mstop(cv_covarspatial))
Model_covary_spatial[st]


st<-(mstop(cv_tree_spatial ))
Model_tree_spatial[st]


####################################################################################
##### Evaluation of different models using multiplicity adjusted all-pairwise comparisons
library("multcomp")
library("lme4")

extrb <- function(obj, m = mstop(obj))
    obj[, attr(obj, "mstop") == m]
nm <- c("(spatial)", "(Climate+Hydro)","(Climate)","(Hydro)", "(Climate+Hydro/spatial)","(tree/spatial)")
tmp <- data.frame(cv = c(extrb(cv_spatial), extrb(cv_covar), extrb(cv_covar_climate),extrb(cv_covar_flow),extrb(cv_covarspatial),
				 extrb(cv_tree_spatial)),
				model = gl(6, length(extrb(cv_spatial))),b = factor(rep(1:length(extrb(cv_spatial)), 6)))
				  

levels(tmp$model) <- nm


png("Wet_Natives NegLL comparing different models.png",width=25, height=25, units='cm', res=300, pointsize=20, bg='white')
        par(mar=c(5,4,1,1),cex=1,oma=c(3,2,1,1))
		
boxplot(cv ~  model, data = tmp, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
		axis(1, at = 1:6, label = FALSE, tick = FALSE, las = 3, cex.axis = 0.75)
		text(seq_along(levels(tmp$model)), par("usr")[3] - 0.09,srt = 30, adj = 1, label = levels(tmp$model), xpd = TRUE, font = 1, cex=0.7)
		axis(2)
		out <- tapply(1:nrow(tmp), tmp$b, function(x) lines(1:6, tmp[x,"cv"], col = rgb(0,0,0,0.1)))
		box(which = "plot", lty = "solid")

dev.off()


sgmod <- summary(glht(lmer(cv ~ model + (1 | b), data = tmp), mcp(model = "Tukey")))
png("Wet_Natives NegLL comparing different models.png",width=25, height=25, units='cm', res=300, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),cex=1,oma=c(2,0,1,0.5))
		
boxplot(cv ~  model, data = tmp, axes = FALSE,
        ylab = "Out-of bootstrap neg. log-lik", xlab = "")
out <- tapply(1:nrow(tmp), tmp$b, function(x) lines(1:13, tmp[x,"cv"], col = rgb(0,0,0,0.1)))

axis(1, at = 1:13, label = levels(tmp$model), tick = FALSE, las = 3, cex.axis = 0.85)
axis(2)

### stability selection
sel <- stabsel(Model_covary_spatial,cv_covarspatial, q = 5)


par(mfrow=c(4, 4)) 


predicted<-list() # creates a whole bunch of empty lists
predicted.insample <-list()
nuvec<-list()
null.model.coef<-list()
null.model.nu<-list()


# Predictions for out-of-bootstrap data
predictions<-predict(Model_covary_spatial,type='response')
rownames(predictions)=rownames(daten) # I removed rows where the TSLW was beyond record but the predicted dataframe has consec row numbers
predictions=as.data.frame(predictions)
plot(predictions[,1],daten$Wet_Natives)
abline(0,1)


############################################################################################################################################################

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