# GDM modelling trials on Hattah Lakes Floodplain dataset
# C. James 26 August 2016
library(gdm)

# Set up data into correct format
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Final_Metrics_Hattah_FP.csv") # load data

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))
mycodes=data.frame(read.csv("HL_FP_sp_codes.csv"))
mydata=merge(mydata,mycodes,by.x="Scientific_name", by.y="Scientific.name")
takeout=c("Inundate", "Lea.litt", "Bar.grou")
mydata=mydata[!mydata$sp_code %in% takeout,]

env.data=data.matrix[,c("Row.names", "d30", "d365", "Inundated", "Flood_frequency", "TSLF", "Easting", "Northing")]
env.data$TSLF[is.na(env.data$TSLF)] <- 10000
colnames(env.data)[1]=c("Unique_site_year")

data.trial<- merge(mydata[,c("Unique_site_year","sp_code","Abundance")],env.data,by="Unique_site_year")
sppTab=data.trial[,c("sp_code","Unique_site_year", "Abundance")]
envTab<-data.trial[,c("Unique_site_year","d30","d365","Inundated", "Flood_frequency", "TSLF", "Easting", "Northing")]
envTab$Inundated <-as.factor(envTab$Inundated)

gdm.data <- formatsitepair(sppTab, 2, abundance=TRUE,sppColumn="sp_code",siteColumn="Unique_site_year",abundColumn="Abundance", XColumn="Easting", YColumn="Northing",predData=envTab)




gdm.hattah.t1 <-gdm(gdm.data,geo=TRUE)
length(gdm.hattah.t1$predictors)
plot(gdm.hattah.t1, plot.layout=c(3,3))