##############################################################################
#### Model richness and abundance of wetland species at Hattah Lakes wetlands
#### C James October 2018 

library(corrplot)
library(mboost)
library(countreg)
library(party)
library(mgcv)


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"
envdata=data.matrix.env.data # copy data

# sort out environmental data - centre and log
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond record
envdata$TSLW[envdata$Innundated==TRUE]<-1 # if site is recorded as inundated change TSLW to 1 day
envdata$d3Mon_wet[is.na(envdata$d3Mon_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d1yrs_wet[is.na(envdata$d1yrs_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d3yrs_wet[is.na(envdata$d3yrs_wet)]<-0 # if site us never wet during time period replace with 0
envdata$d3Mon_meandepth[is.na(envdata$d3Mon_meandepth)]<-0 # if site us never wet during time period replace with 0
envdata$d1yrs_meandepth[is.na(envdata$d1yrs_meandepth)]<-0 # if site us never wet during time period replace with 0
envdata$d3yrs_meandepth[is.na(envdata$d3yrs_meandepth)]<-0 # if site us never wet during time period replace with 0
envdata$TSLW=as.numeric(scale(log10(envdata$TSLW+1),center=TRUE, scale=FALSE)) #This is a very strongly skewed variable 


mod1 <- gam( Wet_Natives ~ s(TSLW)+ s( Site.ID , bs="re" )+ s( Easting , Northing )  , family = ziP , data=envdata )