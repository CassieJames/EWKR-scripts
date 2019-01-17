##########################################################################################################
### Script to determine flood history for Hattah Lakes wetland sites
### Cassie James
### 4 September 2018
library(hydrostats)
library(tidyr)
library(ggplot2)
library(zoo) # sorts out date issues with origin needing to be supplied error

# step 1 - upload all required data 

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Depth.dat=read.csv("Hattah_sites_time_series_data.csv", check.names=FALSE) # load hydrological information but prevent r form converting column names to a more 'friendly' format

# Bring in wetland data and new/corrected dates

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
Output= matrix(NA,nrow=nrow(Sites_sampled), ncol=17)
colnames(Output)=c("Site.ID", "Date.of.collection", "TSLW", "MRI_length", "MRI_shallow", "MRI_medium", "MRI_deep","d1yrs_dry", "d1yrs_shallow", "d1yrs_medium", "d1yrs_deep", "d3yrs_dry", "d3yrs_shallow", "d3yrs_medium", "d3yrs_deep", "Freq_d1", "Freq_d3")
Output=as.data.frame(Output)
Output[,1] <-Sites_sampled$Site.ID.x 
Output[,2] <-Sites_sampled$Date.of.collection.y


# Loop through sites

for (i in 1:nrow(Sites_sampled)){ # 
soi=Sites_sampled[i,1]

DepOI = as.data.frame(Depth.dat[,grepl(soi, colnames(Depth.dat), fixed=TRUE)]) # import modelled data for site of interest (soi)
DepOI=cbind(Depth.dat$DATE,DepOI) # bind date and depth together
colnames(DepOI)[1]="Date"
DepOI$Date=as.Date(as.character(DepOI$Date), format = "%m/%d/%Y") 

# This bit just means across the two quadrat ends where they exist
mycols=ncol(DepOI)
if(mycols>=3){
DepOImean=as.data.frame(rowMeans(DepOI[,2:mycols])); DepOImean=as.data.frame(DepOImean); DepOI=cbind(DepOI$Date,DepOImean); colnames(DepOI)=c("Date", "Depth") } else {DepOI=DepOI} 

doi=Sites_sampled[i,2] # extracts date if interest
Depth_on_day <-DepOI[which(DepOI$Date==doi),2]

# Script to calculate days since last wet

Startdate=as.Date(paste("2005-01-01",sep=""), format="%Y-%m-%d")  
TSLWdata=DepOI[DepOI$Date %in% Startdate:doi,]
tdata=TSLWdata[which(TSLWdata$Depth >0),]
tdatarows=nrow(tdata)
if(tdatarows>0){
TSLW=tdata[nrow(tdata),"Date"]-1}else{
TSLW="NA"}

rtdates =doi-Startdate
rtd=DepOI[DepOI$Date %in% as.Date(doi:Startdate),] # extracts full record prior to sampling
colnames(rtd)[2] ="Q"


HSpell=high.spell.lengths(rtd,threshold=0.001) 

if(!is.na(HSpell[1])) {
MRI_startdate=HSpell[nrow(HSpell),1] # most recent inundation start date MRI is Most Recent Inundation
MRI_length=HSpell[nrow(HSpell),2] # most recent inundation length in days
MRI_enddate=MRI_startdate+MRI_length
MRI_record =DepOI[DepOI$Date %in% as.Date(MRI_startdate:MRI_enddate),] # extracts last inundation event from record
colnames(MRI_record)[2] ="Q"

MRI_shallow<-nrow(MRI_record[MRI_record$Q <= .10 & MRI_record$Q > 0 ,])
MRI_medium<-nrow(MRI_record[MRI_record$Q <= .50 & MRI_record$Q > .10 ,])
MRI_deep<-nrow(MRI_record[MRI_record$Q > .50,]) }else{

MRI_length=0
MRI_shallow=0
MRI_medium=0
MRI_deep=0}

d1yrs = doi-(365*1)
d1yrs_record =DepOI[DepOI$Date %in% as.Date(doi:d1yrs),] # extracts last inundation event from record
colnames(d1yrs_record)[2] ="Q"

d1yrs_dry<-nrow(d1yrs_record[d1yrs_record$Q <= 0 ,])
d1yrs_shallow<-nrow(d1yrs_record[d1yrs_record$Q <= .10 & d1yrs_record$Q > 0 ,])
d1yrs_medium<-nrow(d1yrs_record[d1yrs_record$Q <= .50 & d1yrs_record$Q > .10 ,])
d1yrs_deep<-nrow(d1yrs_record[d1yrs_record$Q > .50,])

d3yrs = doi-(365*3)
d3yrs_record =DepOI[DepOI$Date %in% as.Date(doi:d3yrs),] # extracts last inundation event from record
colnames(d3yrs_record)[2] ="Q"

d3yrs_dry<-nrow(d3yrs_record[d3yrs_record$Q <= 0 ,])
d3yrs_shallow<-nrow(d3yrs_record[d3yrs_record$Q <= .10 & d3yrs_record$Q > 0 ,])
d3yrs_medium<-nrow(d3yrs_record[d3yrs_record$Q <= .50 & d3yrs_record$Q > .10 ,])
d3yrs_deep<-nrow(d3yrs_record[d3yrs_record$Q > .50,])

if(!is.na(high.spell.lengths(d1yrs_record, thresh=0.00001)[1])){
Freq_d1 <-nrow(high.spell.lengths(d1yrs_record, thresh=0.00001))}else{Freq_d1=0}

if(!is.na(high.spell.lengths(d3yrs_record, thresh=0.00001)[1])){
Freq_d3 <-nrow(high.spell.lengths(d3yrs_record, thresh=0.00001))}else{Freq_d3=0}


if(tdatarows>0){
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),3] <- as.numeric(doi-TSLW, units = "days")}else{
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),3] <- "NA"
}

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),4] <- MRI_length
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),5] <- MRI_shallow
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),6] <- MRI_medium
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),7] <- MRI_deep

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),8] <- d1yrs_dry/365
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),9] <- d1yrs_shallow/365
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),10] <- d1yrs_medium/365
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),11] <- d1yrs_deep/365

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),12] <- d3yrs_dry/(365+365+365)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),13] <- d3yrs_shallow/(365+365+365)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),14] <- d3yrs_medium/(365+365+365)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),15] <- d3yrs_deep/(365+365+365)

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),16] <- Freq_d1
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),17] <- Freq_d3


}



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

Output$WaterYr <-wtr_yr(Output$Date.of.collection)


### Codes altered to avoid confusion and allow merging with other datasets
Output$Site.ID <-gsub("LHT", "LHAT", Output$Site.ID) # Changed so that scripts were not confusing HT and LHT
Output$Site.ID <-gsub("CHT", "CCS", Output$Site.ID) # Changed so that scripts were not confusing CHT and HT
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
write.csv(Output , file = "Hydraulics_HTH_WL.csv")