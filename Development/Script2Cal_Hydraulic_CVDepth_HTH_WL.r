##########################################################################################################
### Script to determine hydrological variability for each site
### Cassie James
### 16 May 2018
library(hydrostats)
library(tidyr)
library(ggplot2)
library(zoo) # sorts out date issues with origin needing to be supplied error

# step 1 - upload all required data 

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Depth.dat=read.csv("Hattah_sites_time_series_data.csv", check.names=FALSE) # load hydrological information but prevent r form converting column names to a more 'friendly' format
Depth_long <- gather(Depth.dat, site, depth, gather_cols=2:406, factor_key=TRUE) # had to use numbers for columns as function would not recognise names

Depth_long$site <-as.character(Depth_long$site)
temp <- strsplit(Depth_long$site, split="_")
SITE_only <- sapply(temp, "[", 1)
Depth_long$SITE <-SITE_only # this dataset has both ends of the transect - will need to use grep to match with vegetation dataset which does not include a or b
Depth_long$Date=as.Date(as.character(Depth_long$DATE), format = "%m/%d/%Y") # R needs to recognise column as dates

Hydrology=Depth_long # take copy of data

# Bring in wetland data
data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL_FGcorrections.csv"))
Sites_sampled <-unique(mydata[,c("Site.ID","Unique_site_year","Date.of.collection")]) # creates unique list of sites and dates sampled
Sites_sampled$Date.of.collection<-as.Date(as.character(Sites_sampled$Date.of.collection), format = "%d/%m/%Y")
Sites_sampled$Site.ID <- as.character(Sites_sampled$Site.ID)


# Create empty data frame to receive results
Output= matrix(NA,nrow=nrow(Sites_sampled), ncol=6)
colnames(Output)=c("Site.ID","Unique_site_year","Date.of.collection","HydroYear", "CV_depth","TSLW")
Output=as.data.frame(Output)
Output[,1:3] <-Sites_sampled[,1:3]
Output$Site.ID <-gsub("CCNT", "NCT", Output$Site.ID)
Output$Site.ID <-gsub("KRT", "KT", Output$Site.ID)
Output$Unique_site_year <-gsub("CCNT", "NCT", Output$Unique_site_year)
Output$Unique_site_year <-gsub("KRT", "KT", Output$Unique_site_year)

# Run through each row of the site and date list and apply function to extract depth and duration info over the different time frames
# This function identifies the water year (July 1st to 30th June using the year going into as the water year attribute)
wtr_yr <- function(dates, start_month=7) {
  # Convert dates into POSIXlt
  dates.posix = as.POSIXlt(dates)
  # Year offset
  offset = ifelse(dates.posix$mon >= start_month - 1, 1, 0)
  # Water year
  adj.year = dates.posix$year + 1900 + offset
  # Return the water year
  adj.year
}

for (i in 1:nrow(Sites_sampled)){ # 
soi=Sites_sampled[i,1]
soi=ifelse(grepl("CCNT",soi),gsub("CCNT","NCT",soi),soi) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie
soi=ifelse(grepl("KRT",soi),gsub("KRT","KT",soi),soi) 
Hydrology = as.data.frame(Depth.dat[,grepl(soi, colnames(Depth.dat), fixed=TRUE)]) # extracts timeseries data for both ends of transect
Hydrology=cbind(Depth.dat$DATE,Hydrology)
Hydrology$Depthmean <-rowMeans(Hydrology[2:ncol(Hydrology)]) # as there is information to both ends of the transect I have taken the mean of this
colnames(Hydrology)[1]="Date"
Hydrology=Hydrology[,c("Date", "Depthmean")]
colnames(Hydrology)[2]="Q"
Hydrology <- within(Hydrology, Date <- as.Date(as.character(Date), format = "%m/%d/%Y"))
doi=Sites_sampled[i,3]
doi=as.Date(paste(doi,sep=""), format="%Y-%m-%d")  

#Script to calculate days since last wet
Startdate=as.Date(paste("2005-01-01",sep=""), format="%Y-%m-%d")  
TSLWdata=Hydrology[Hydrology$Date %in% Startdate:doi,]
TSLWdata$hydro.year=wtr_yr(TSLWdata$Date) 
tdata=TSLWdata[which(TSLWdata$Q >0),]
TSLW=tdata[nrow(tdata),"Date"]

d1095=doi-1095	# twelve months ago
TSLWdata3=TSLWdata[TSLWdata$Date %in% as.Date(doi:d1095),]

if(nrow(tdata)>0){
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),6] <- as.numeric(doi-TSLW, units = "days")}else{
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),6] <- "NA"
}
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),4] <- TSLWdata[nrow(TSLWdata),"hydro.year"]
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),5] <- (sd(TSLWdata3$Q,na.rm=TRUE)/mean(TSLWdata3$Q))
}

Output$sites<- Output$Site.ID
Output$sites[grepl("BIT",Output$sites)] <- "BIT"
Output$sites[grepl("BLT",Output$sites)] <- "BLT"
Output$sites[grepl("BOT",Output$sites)] <- "BOT"
Output$sites[grepl("BRT",Output$sites)] <- "BRT"
Output$sites[grepl("CCS",Output$sites)] <- "CCS"
Output$sites[grepl("HT",Output$sites)] <- "HT"
Output$sites[grepl("LHAT",Output$sites)] <- "LHAT"
Output$sites[grepl("KT",Output$sites)] <- "KT"
Output$sites[grepl("YT",Output$sites)] <- "YT"
Output$sites[grepl("MOT",Output$sites)] <- "MOT"
Output$sites[grepl("NCT",Output$sites)] <- "NCT"
Output$sites[grepl("NN",Output$sites)] <- "NN"

Output <-replace(Output,is.na(Output),0)

mydata_mean <- aggregate(Output, by = list(Output$sites,Output$HydroYear), FUN = mean)
mydata_max <- aggregate(Output, by = list(Output$sites,Output$HydroYear), FUN = max)

CV_Depth <-cbind(mydata_mean,mydata_max)

CV_Depth=CV_Depth[,c(1,2,7,15,17)]
colnames(CV_Depth)=c("Sites","HydroYear","MeanDepthCV","MaxDepthCV","TSLW")


rownames(CV_Depth)=paste(CV_Depth$Sites,"_",right(CV_Depth$HydroYear,2),sep="")
groups=as.factor(groups)

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
write.csv(CV_Depth , file = "CVDepth_HTH_WL.csv")



######################################################################################################################################################################
#### relating CV_depth to alpha diversity

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
data.matrix.HTWL=read.csv("Spp_site_matrix summarised to wetland_HTH_WL_June 2018.csv",row.names = 1 ) # load data - object name is 'tada' ... :)

mydataPA=data.matrix.HTWL
mydataPA[mydataPA>0] <-1


data.matrix.HTWL$SpeciesRich <-rowSums(mydataPA)
Output$sites<- Output$Site.ID
Output$sites[grepl("BIT",Output$sites)] <- "BIT"
Output$sites[grepl("BLT",Output$sites)] <- "BLT"
Output$sites[grepl("BOT",Output$sites)] <- "BOT"
Output$sites[grepl("BRT",Output$sites)] <- "BRT"
Output$sites[grepl("CCS",Output$sites)] <- "CCS"
Output$sites[grepl("HT",Output$sites)] <- "HT"
Output$sites[grepl("LHAT",Output$sites)] <- "LHAT"
Output$sites[grepl("KT",Output$sites)] <- "KT"
Output$sites[grepl("YT",Output$sites)] <- "YT"
Output$sites[grepl("MOT",Output$sites)] <- "MOT"
Output$sites[grepl("NCT",Output$sites)] <- "NCT"
Output$sites[grepl("NN",Output$sites)] <- "NN"

meeee=merge(data.matrix.HTWL,CV_Depth, by="row.names")

meeee$TSLW[meeee$TSLW=="NA"] <-2000


read.csv("HattahWL_betawithinsites.csv",row.names=1) 
