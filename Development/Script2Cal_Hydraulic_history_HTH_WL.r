##########################################################################################################
### Script to determine flood history for Hattah Lakes wetland sites
### Cassie James
### 16 May 2018
library(hydrostats)
library(tidyr)
library(ggplot2)
library(zoo) # sorts out date issues with origin needing to be supplied error

# step 1 - upload all required data 

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Depth.dat=read.csv("Hattah_sites_time_series_data.csv", check.names=FALSE) # load hydrological information but prevent r form converting column names to a more 'friendly' format

# Bring in wetland data and new dates

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL.csv"))
mydata$Site.ID=gsub("CCNT","NCT",mydata$Site.ID) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie - my attempt to make everything consistent
mydata$Site.ID=gsub("KRT","KT",mydata$Site.ID)
mydata$Unique_site_year=gsub("CCNT","NCT",mydata$Unique_site_year) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie
mydata$Unique_site_year=gsub("KRT","KT",mydata$Unique_site_year)

date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)

mydata2=merge(mydata, newdates, by="Unique_site_year" )
Sites_sampled <-unique(mydata2[,c("Site.ID.x","Date.of.collection.y")]) # creates unique list of sites and dates

# Create empty data frame to receive results
Output= matrix(NA,nrow=nrow(Sites_sampled), ncol=7)
colnames(Output)=c("Site.ID", "Date.of.collection", "TSLW", "Dry_Annual", "Shallow_Annual", "Medium_Annual", "Deep_Annual")
Output=as.data.frame(Output)
Output[,1] <-Sites_sampled$Site.ID.x 
Output[,2] <-Sites_sampled$Date.of.collection.y


# Loop through sites

for (i in 1:nrow(Sites_sampled)){ # 
soi=Sites_sampled[i,1]
DepOI = as.data.frame(Depth.dat[,grepl(soi, colnames(Depth.dat), fixed=TRUE)]) # import modelled data for site of interest (soi)
DepOI=cbind(Depth.dat$DATE,DepOI) # bind date and depth
colnames(DepOI)[1]="Date"
DepOI$Date=as.Date(as.character(DepOI$Date), format = "%m/%d/%Y") 

# This bit just means across the two quadrat ends where they exist
mycols=ncol(DepOI)
if(mycols>=3){
DepOImean=as.data.frame(rowMeans(DepOI[,2:mycols])); DepOImean=as.data.frame(DepOImean); DepOI=cbind(DepOI$Date,DepOImean); colnames(DepOI)=c("Date", "Depth") } else {DepOI=DepOI} 

doi=Sites_sampled[i,2] # extracts date if interest
Depth_on_day <-DepOI[which(DepOI$Date==doi),2]

# Determination of metrics for annual time period

d365 =doi-364 
loi365=DepOI[DepOI$Date %in% as.Date(doi:d365),] # extracts year of data prior to sampling date

Dry<-nrow(loi365[loi365$Depth == 0,])
Shallow<-nrow(loi365[loi365$Depth <= 10 & loi365$Depth > 0 ,])
Medium<-nrow(loi365[loi365$Depth <= 50 & loi365$Depth > 10 ,])
Deep<-nrow(loi365[loi365$Depth > 50,])

#Script to calculate days since last wet
Startdate=as.Date(paste("2005-01-01",sep=""), format="%Y-%m-%d")  
TSLWdata=DepOI[DepOI$Date %in% Startdate:doi,]

tdata=TSLWdata[which(TSLWdata$Depth >0),]

tdatarows=nrow(tdata)
if(tdatarows>0){
TSLW=tdata[nrow(tdata),"Date"]-1}else{
TSLW="NA"}



if(tdatarows>0){
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),3] <- as.numeric(doi-TSLW, units = "days")}else{
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),3] <- "NA"
}
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),4] <- Dry/365
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),5] <- Shallow/365
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),6] <- Medium/365
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),7] <- Deep/365

}








# Create empty data frame to receive results
Output= matrix(NA,nrow=length(sitelist), ncol=2)
colnames(Output)=c("Date.of.collection","TSLF")
rownames(Output)=sitelist
Output=as.data.frame(Output)

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

for (s in sitelist){ # 

tdata=Depth_long[grep(s,Depth_long$site),]

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
Hydrology.dat=as.data.frame(t(sapply(doi, myfunction)))

#Script to calculate days since last wet
Startdate=as.Date(paste("2005-01-01",sep=""), format="%Y-%m-%d")  
TSLWdata=Hydrology[Hydrology$Date %in% Startdate:doi,]
tdata=TSLWdata[which(TSLWdata$Q >0),]
TSLW=tdata[nrow(tdata),"Date"]

#script to determine number of years inundated (frequency)
TSLWdata$hydro.year=wtr_yr(TSLWdata$Date) 
mydata=aggregate(TSLWdata$Q, by=list(TSLWdata$hydro.year),sum)
mydatawet=mydata[which(mydata$x>0),]

#CV = sd(TSLWdata$Q, na.rm=TRUE)/mean(TSLWdata$Q, na.rm=TRUE)*100 leave for now as too strongly dependent upon the magnitude
   


Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),4:11] <- Hydrology.dat
if(nrow(tdata)>0){
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),12] <- as.numeric(doi-TSLW, units = "days")}else{
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),12] <- "NA"
}
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),13] <- nrow(mydatawet)/nrow(mydata)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),14] <- mydata[nrow(mydata),"x"]
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),15] <- TSLWdata[nrow(TSLWdata),"hydro.year"]

}

### Codes altered to avoid confusion and allow merging with other datasets
Output$Unique_site_year <-gsub("LHT", "LHAT", Output$Unique_site_year) # Changed so that scripts were not confusing HT and LHT
Output$Unique_site_year <-gsub("CHT", "CCS", Output$Unique_site_year) # Changed so that scripts were not confusing CHT and HT
Output$Site.ID <-gsub("LHT", "LHAT", Output$Site.ID) # Changed so that scripts were not confusing HT and LHT
Output$Site.ID <-gsub("CHT", "CCS", Output$Site.ID) # Changed so that scripts were not confusing CHT and HT
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
write.csv(Output , file = "Hydraulics_HTH_WL.csv")