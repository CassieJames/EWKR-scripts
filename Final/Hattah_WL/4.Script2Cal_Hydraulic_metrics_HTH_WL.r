###################################################################################################################################
### Script to determine short term flood frequency for Hattah Lakes wetland sites
### Written by  C.S.James (JCU, TropWATER)
### GNU General Public License .. feel free to use / distribute ... no warranties
### 27th August 2019
###################################################################################################################################
### Notes: This script determines some hydraulic variables based on the BIGMOD modelled waterdepths - bigmod data is only available from 1/1/2005
###################################################################################################################################
library(hydrostats)
library(tidyr)
library(ggplot2)
library(zoo) # sorts out date issues with origin needing to be supplied error

# step 1 - upload all required data and sort out site names

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Depth.dat=read.csv("Hattah_sites_time_series_dataV3.csv", check.names=FALSE) # load hydrological information but prevent r from converting column names to a more 'friendly' format
colnames(Depth.dat)=gsub("KRT","KT",colnames(Depth.dat))
colnames(Depth.dat)=gsub("CCNT","NCT",colnames(Depth.dat))
colnames(Depth.dat)=gsub("CHT","CCS",colnames(Depth.dat))
colnames(Depth.dat)=gsub("LHT","LHAT",colnames(Depth.dat))

# REMEMBER to sort out issue with names CHT4N and CHTN4 in hydrology file

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL.csv"))
mydata$Site.ID=gsub("CCNT","NCT",mydata$Site.ID) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie - my attempt to make everything consistent
mydata$Site.ID=gsub("KRT","KT",mydata$Site.ID)
mydata$Site.ID=gsub("CHT","CCS",mydata$Site.ID)
mydata$Site.ID=gsub("LHT","LHAT",mydata$Site.ID)
mydata$Unique_site_year=gsub("CCNT","NCT",mydata$Unique_site_year) 
mydata$Unique_site_year=gsub("KRT","KT",mydata$Unique_site_year)
mydata$Unique_site_year=gsub("CHT","CCS",mydata$Unique_site_year)
mydata$Unique_site_year=gsub("LHT","LHAT",mydata$Unique_site_year)

date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CHT","CCS",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("LHT","LHAT",newdates$Unique_site_year)

mydata2=merge(mydata, newdates, by="Unique_site_year")
Sites_sampled <-unique(mydata2[,c("Site.ID.x","Date.of.collection.y", "Unique_site_year")]) # creates unique list of sites and dates

# step 2 - Create empty data frame to receive results

Output= matrix(NA,nrow=nrow(Sites_sampled), ncol=29)
colnames(Output)=c("Site.ID", "Date.of.collection", "TSLW", "MRI_length", "MRI_meandepth","MRI_shallow","MRI_deep",
"d3Mon_wet","d3Mon_meandepth","d3Mon_shallow", "d3Mon_deep","d3Mon_cv","d3Mon_cvmonthall","d3Mon_drylength",
"d1yrs_wet","d1yrs_meandepth","d1yrs_shallow", "d1yrs_deep", "d1yrs_cv","d1yrs_cvmonthall","d1yrs_drylength",
"d3yrs_wet","d3yrs_meandepth", "d3yrs_shallow", "d3yrs_deep","d3yrs_cv","d3yrs_cvmonthall","d3yrs_drylength","Unique_site_year")
Output=as.data.frame(Output)
Output[,1] <-Sites_sampled$Site.ID.x 
Output[,2] <-Sites_sampled$Date.of.collection.y
Output[,29] <-Sites_sampled$Unique_site_year


# step 3 - Loop through sites

for (i in 1:nrow(Sites_sampled)){ # 
soi=Sites_sampled[i,1]

DepOI = as.data.frame(Depth.dat[,grepl(soi, colnames(Depth.dat), fixed=TRUE)]) # import modelled data for site of interest (soi)
DepOI=cbind(Depth.dat$DATE,DepOI) # bind date and depth together
colnames(DepOI)[1]="Date"
DepOI$Date=as.Date(as.character(DepOI$Date), format = "%d/%m/%Y") 

# This bit just means across the two quadrat ends where they exist
mycols=ncol(DepOI)

if(mycols>=3){ # note its three columns because of Date column as well as start and end of transect - note that a few transects don't have start and end info
DepOImean=as.data.frame(rowMeans(DepOI[,2:mycols])); DepOImean=as.data.frame(DepOImean); DepOI=cbind(DepOI$Date,DepOImean); colnames(DepOI)=c("Date", "Depth") } else {DepOI=DepOI} 

doi=Sites_sampled[i,2] # extracts date of interest
Depth_on_day <-DepOI[which(DepOI$Date==doi),2]

# step 3.1 - Script to calculate days since last wet

Startdate=as.Date(paste("2005-01-01",sep=""), format="%Y-%m-%d")  
TSLWdata=DepOI[DepOI$Date %in% Startdate:doi,] # extracts data up until date of interest
tdata=TSLWdata[which(TSLWdata$Depth >0),] # subsets data to only dates with water
tdatarows=nrow(tdata)
if(tdatarows>0){
TSLW=tdata[nrow(tdata),"Date"]-1
TSLW=doi-TSLW
TSLW=TSLW[[1]]
}else{
TSLW="NA"}


# step 3.2 - Script to determine character of most recent inundation (MRI)

rtdates =doi-Startdate # provides number of dates between doi and start of record
rtd=DepOI[DepOI$Date %in% as.Date(doi:Startdate),] # extracts full record prior to sampling
colnames(rtd)[2] ="Q"

HSpell=high.spell.lengths(rtd,threshold=0.0001) 

if(!is.na(HSpell[[1]][1])) {
MRI_startdate=HSpell[nrow(HSpell),1] # most recent inundation start date MRI is Most Recent Inundation
MRI_length=HSpell[nrow(HSpell),2] # most recent inundation length in days
MRI_enddate=MRI_startdate+MRI_length
MRI_record =DepOI[DepOI$Date %in% as.Date(MRI_startdate:MRI_enddate),] # extracts last inundation event from record
colnames(MRI_record)[2] ="Q"

MRI_shallow<-nrow(MRI_record[MRI_record$Q <= .10 & MRI_record$Q > 0 ,])
MRI_deep<-nrow(MRI_record[MRI_record$Q > .10 ,])
MRI_meandepth<-mean(MRI_record[MRI_record$Q > 0, 2])}else{
MRI_meandepth=0
MRI_length=0
MRI_shallow=0
MRI_deep=0}

##### step 3.3 - step determines conditional mean depth and proportion of time wet metrics for non nested time series 3 months, 1 year and 3 years

d3Mon = doi-(90)
d3Mon_record =DepOI[DepOI$Date %in% as.Date(doi:d3Mon),] # extracts last inundation event from record
colnames(d3Mon_record)[2] ="Q"

d3Mon_meandepth<-mean(d3Mon_record[d3Mon_record$Q >0,2 ])
d3Mon_dry<-nrow(d3Mon_record[d3Mon_record$Q <= 0 ,])
d3Mon_wet<-nrow(d3Mon_record[d3Mon_record$Q >0 ,])
d3Mon_shallow<-nrow(d3Mon_record[d3Mon_record$Q <= .10 & d3Mon_record$Q > 0 ,])
d3Mon_deep<-nrow(d3Mon_record[d3Mon_record$Q > .10,])
d3Mon_cv<-daily.cv(d3Mon_record[d3Mon_record$Q > 0 ,])
d3Mon_cvmonthall<-monthly.cv(d3Mon_record)
d3Mon_drylength<-low.spell.lengths(d3Mon_record,threshold=0.0001)
d3Mon_drylength=max(d3Mon_drylength[,2])

d1yrs = d3Mon -(365*1)
d1yrs_record =DepOI[DepOI$Date %in% as.Date(d3Mon:d1yrs),] # extracts last inundation event from record
colnames(d1yrs_record)[2] ="Q"

d1yrs_meandepth<-mean(d1yrs_record[d1yrs_record$Q >0 ,2])
d1yrs_dry<-nrow(d1yrs_record[d1yrs_record$Q <= 0 ,])
d1yrs_wet<-nrow(d1yrs_record[d1yrs_record$Q >0 ,])
d1yrs_shallow<-nrow(d1yrs_record[d1yrs_record$Q <= .10 & d1yrs_record$Q > 0 ,])
d1yrs_deep<-nrow(d1yrs_record[d1yrs_record$Q > .10,])
d1yrs_cv<-daily.cv(d1yrs_record[d1yrs_record$Q > 0 ,])
d1yrs_cvmonthall<-monthly.cv(d1yrs_record)
d1yrs_drylength<-low.spell.lengths(d1yrs_record,threshold=0.0001)
d1yrs_drylength=max(d1yrs_drylength[,2])

if(doi <= as.Date(paste("2008-04-01",sep=""), format="%Y-%m-%d")){  
d3yrs = d1yrs-(round(365.242189*2,0))
d3yrs_record =DepOI[DepOI$Date %in% as.Date(d1yrs:d3yrs),] # extracts last inundation event from record
colnames(d3yrs_record)[2] ="Q"}

if(doi > as.Date(paste("2008-04-01",sep=""), format="%Y-%m-%d")){  
d3yrs = d1yrs-(round(2*365.242189,0))
d3yrs_record =DepOI[DepOI$Date %in% as.Date(d1yrs:d3yrs),] # extracts last inundation event from record
colnames(d3yrs_record)[2] ="Q"}

d3yrs_meandepth<-mean(d3yrs_record[d3yrs_record$Q >0 ,2])
d3yrs_dry<-nrow(d3yrs_record[d3yrs_record$Q <= 0 ,])
d3yrs_wet<-nrow(d3yrs_record[d3yrs_record$Q >0 ,])
d3yrs_shallow<-nrow(d3yrs_record[d3yrs_record$Q <= .10 & d3yrs_record$Q > 0 ,])
d3yrs_deep<-nrow(d3yrs_record[d3yrs_record$Q > .10,])
d3yrs_cv<-monthly.cv(d3yrs_record)
d3yrs_cv<-daily.cv(d3yrs_record[d3yrs_record$Q > 0 ,])
d3yrs_cvmonthall<-monthly.cv(d3yrs_record)
d3yrs_drylength<-low.spell.lengths(d3yrs_record,threshold=0.0001)
d3yrs_drylength=max(d3yrs_drylength[,2])

##### step 3.4 - Determines time periods with nested time series 1 year, 3 years, 5 and 10 years for short to medium term FF calcs

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),3] <- TSLW
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),4] <- MRI_length
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),5] <- MRI_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),6] <- MRI_shallow
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),7] <- MRI_deep

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),8] <- d3Mon_wet/91
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),9] <- d3Mon_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),10] <- d3Mon_shallow/91
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),11] <- d3Mon_deep/91
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),12] <- d3Mon_cv
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),13] <- d3Mon_cvmonthall
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),14] <- d3Mon_drylength

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),15] <- d1yrs_wet/366
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),16] <- d1yrs_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),17] <- d1yrs_shallow/366
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),18] <- d1yrs_deep/366
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),19] <- d1yrs_cv
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),20] <- d1yrs_cvmonthall
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),21] <- d1yrs_drylength

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),22] <- d3yrs_wet/(365+365+366)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),23] <- d3yrs_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),24] <- d3yrs_shallow/(365+365+366)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),25] <- d3yrs_deep/(365+365+366)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),26] <- d3yrs_cv
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),27] <- d3yrs_cvmonthall
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),28] <- d3yrs_drylength


Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),29] <- Sites_sampled[i,3]

if(i==1) {

}}



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
Output$DOI_MRI_start=Output$Date.of.collection-Output$MRI_length # tells us start of the MRI - to check whether it is being truncated by the record


### Codes altered to avoid confusion and allow merging with other datasets
Output$Site.ID <-gsub("LHT", "LHAT", Output$Site.ID) # Changed so that scripts were not confusing HT and LHT
Output$Site.ID <-gsub("CHT", "CCS", Output$Site.ID) # Changed so that scripts were not confusing CHT and HT
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
write.csv(Output , file = "Hydraulics_HTH_WL_medium and deep merged December 2020.csv")