# GDM modelling trials on Hattah Lakes Floodplain dataset
# C. James 26 August 2016
library(gdm)

# Set up data into correct format
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Final_Metrics_Hattah_FP.csv") # load data


spp.data=mydata[,c("Unique_site_year", "Scientific_name", "Abundance")]
spp.data=merge(spp.data,data.matrix[,c("Row.names", "Easting", "Northing")], by.x="Unique_site_year", by.y="Row.names")
spp.data=spp.data[,c("Easting", "Northing", "Unique_site_year", "Scientific_name", "Abundance")]
env.data=data.matrix[,c("Row.names", "d30", "d365", "Inundated", "Flood_frequency", "TSLF", "Easting", "Northing")]
env.data$TSLF[is.na(env.data$TSLF)] <- 10000

colnames(env.data)[1]=c("Unique_site_year")




gdm.data <- formatsitepair(spp.data, 2,  sppColumn="Scientific_name", abundance=TRUE, abundColumn="Abundance",XColumn="Easting",YColumn="Northing", siteColumn="Unique_site_year", predData=env.data)

gdm.hattah.t1 <-gdm(gdm.data,geo=F)
length(gdm.hattah.t1$predictors)
plot(gdm.hattah.t1, plot.layout=c(3,3))