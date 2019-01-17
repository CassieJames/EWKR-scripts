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
Depth.dat=read.csv("Hattah_sites_time_series_dataV3.csv", check.names=FALSE) # load hydrological information but prevent r form converting column names to a more 'friendly' format

# REMEMBER to sort out issue with names CHT4N and CHTN4 in hydrology file

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

mydata2=merge(mydata, newdates, by="Unique_site_year")
Sites_sampled <-unique(mydata2[,c("Site.ID.x","Date.of.collection.y", "Unique_site_year")]) # creates unique list of sites and dates

# Create empty data frame to receive results
Output= matrix(NA,nrow=nrow(Sites_sampled), ncol=33)
colnames(Output)=c("Site.ID", "Date.of.collection", "TSLW", "MRI_length", "MRI_meandepth","MRI_shallow","MRI_deep","d3Mon_wet","d3Mon_meandepth","d3Mon_shallow", "d3Mon_deep","d1yrs_wet","d1yrs_meandepth","d1yrs_shallow", 
"d1yrs_deep", "d3yrs_wet","d3yrs_meandepth", "d3yrs_shallow", "d3yrs_deep", "Freq_d1", "Freq_d3", "Freq_d5","Freq_ALL","Length_d1", "Length_d3", "Length_d5","CTF_prop_d1","CTF_average_d1","CTF_prop_d3","CTF_average_d3","CTF_prop_d5","CTF_average_d5", "Unique_site_year")
Output=as.data.frame(Output)
Output[,1] <-Sites_sampled$Site.ID.x 
Output[,2] <-Sites_sampled$Date.of.collection.y
Output[,33] <-Sites_sampled$Unique_site_year


# Loop through sites

for (i in 1:nrow(Sites_sampled)){ # 
soi=Sites_sampled[i,1]

DepOI = as.data.frame(Depth.dat[,grepl(soi, colnames(Depth.dat), fixed=TRUE)]) # import modelled data for site of interest (soi)
DepOI=cbind(Depth.dat$DATE,DepOI) # bind date and depth together
colnames(DepOI)[1]="Date"
DepOI$Date=as.Date(as.character(DepOI$Date), format = "%d/%m/%Y") 

# This bit just means across the two quadrat ends where they exist
mycols=ncol(DepOI)
if(mycols>=3){
DepOImean=as.data.frame(rowMeans(DepOI[,2:mycols])); DepOImean=as.data.frame(DepOImean); DepOI=cbind(DepOI$Date,DepOImean); colnames(DepOI)=c("Date", "Depth") } else {DepOI=DepOI} 

doi=Sites_sampled[i,2] # extracts date of interest
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

##### Section determines conditional mean depth and proportion of time wet metrics for non nested time series 3 months, 1 year and 3 years

d3Mon = doi-(90)
d3Mon_record =DepOI[DepOI$Date %in% as.Date(doi:d3Mon),] # extracts last inundation event from record
colnames(d3Mon_record)[2] ="Q"

d3Mon_meandepth<-mean(d3Mon_record[d3Mon_record$Q >0,2 ])
d3Mon_dry<-nrow(d3Mon_record[d3Mon_record$Q <= 0 ,])
d3Mon_wet<-nrow(d3Mon_record[d3Mon_record$Q >0 ,])
d3Mon_shallow<-nrow(d3Mon_record[d3Mon_record$Q <= .10 & d3Mon_record$Q > 0 ,])
d3Mon_deep<-nrow(d3Mon_record[d3Mon_record$Q > .10,])

d1yrs = d3Mon -(365*1)
d1yrs_record =DepOI[DepOI$Date %in% as.Date(d3Mon:d1yrs),] # extracts last inundation event from record
colnames(d1yrs_record)[2] ="Q"

d1yrs_meandepth<-mean(d1yrs_record[d1yrs_record$Q >0 ,2])
d1yrs_dry<-nrow(d1yrs_record[d1yrs_record$Q <= 0 ,])
d1yrs_wet<-nrow(d1yrs_record[d1yrs_record$Q >0 ,])
d1yrs_shallow<-nrow(d1yrs_record[d1yrs_record$Q <= .10 & d1yrs_record$Q > 0 ,])
d1yrs_deep<-nrow(d1yrs_record[d1yrs_record$Q > .10,])

if(doi <= as.Date(paste("2008-04-01",sep=""), format="%Y-%m-%d")){  
d3yrs = d1yrs-(658)
d3yrs_record =DepOI[DepOI$Date %in% as.Date(d1yrs:d3yrs),] # extracts last inundation event from record
colnames(d3yrs_record)[2] ="Q"}

if(doi > as.Date(paste("2008-04-01",sep=""), format="%Y-%m-%d")){  
d3yrs = d1yrs-(2*365)
d3yrs_record =DepOI[DepOI$Date %in% as.Date(d1yrs:d3yrs),] # extracts last inundation event from record
colnames(d3yrs_record)[2] ="Q"}

d3yrs_meandepth<-mean(d3yrs_record[d3yrs_record$Q >0 ,2])
d3yrs_dry<-nrow(d3yrs_record[d3yrs_record$Q <= 0 ,])
d3yrs_wet<-nrow(d3yrs_record[d3yrs_record$Q >0 ,])
d3yrs_shallow<-nrow(d3yrs_record[d3yrs_record$Q <= .10 & d3yrs_record$Q > 0 ,])
d3yrs_deep<-nrow(d3yrs_record[d3yrs_record$Q > .10,])

##### Determines time periods with nested time series 1 year, 3 years and 5 years

d1yrs = doi -(365*1)
d1yrs_record =DepOI[DepOI$Date %in% as.Date(doi:d1yrs),] # extracts last inundation event from record
colnames(d1yrs_record)[2] ="Q"

if(is.na(low.spell.lengths(d1yrs_record, thresh=0.00001)[[1]][1])){ # this code checks to see whether there are any low spell events and if site wet for whole period it returns 0
Freq_d1 <-0
}else{
if(low.spell.lengths(d1yrs_record, thresh=0.00001)[[2]][1]>=366){
Freq_d1<-0
}else{
Freq_d1 <-nrow(low.spell.lengths(d1yrs_record, thresh=0.00001))
}
}

HSpell=high.spell.lengths(d1yrs_record,threshold=0.0001)

if(!is.na(HSpell[[1]][1])) {
d1yrs_length=mean(HSpell$spell.length) # most recent inundation start date MRI is Most Recent Inundation
d1yrs_lengthsd=sd(HSpell$spell.length)}else{
d1yrs_length=0
d1yrs_lengthsd=0}

CTF_1yrs=CTF(d1yrs_record,threshold=0.0001)

d3yrs = doi-(365*3)
d3yrs_record =DepOI[DepOI$Date %in% as.Date(doi:d3yrs),] # extracts full three years prior for flood frequency estimates
colnames(d3yrs_record)[2] ="Q"

if(is.na(low.spell.lengths(d3yrs_record, thresh=0.00001)[[1]][1])){ # this code checks to see whether there are any low spell events and if site wet for whole period it returns 0
Freq_d3 <-0
}else{
if(low.spell.lengths(d3yrs_record, thresh=0.00001)[[2]][1]>=1096){
Freq_d3<-0
}else{
Freq_d3 <-nrow(low.spell.lengths(d3yrs_record, thresh=0.00001))
}
}

HSpell=high.spell.lengths(d3yrs_record,threshold=0.0001) 
if(!is.na(HSpell[[1]][1])) {
d3yrs_length=mean(HSpell$spell.length) # most recent inundation start date MRI is Most Recent Inundation
d3yrs_lengthsd=sd(HSpell$spell.length)}else{
d3yrs_length=0
d3yrs_lengthsd=0}

CTF_3yrs=CTF(d3yrs_record,threshold=0.0001)


d5yrs = doi-(365*5)
d5yrs_record =DepOI[DepOI$Date %in% as.Date(doi:d5yrs),] # extracts last inundation event from record
colnames(d5yrs_record)[2] ="Q"

if(is.na(low.spell.lengths(d5yrs_record, thresh=0.00001)[[1]][1])){ # this code checks to see whether there are any low spell events and if site wet for whole period it returns 0
Freq_d5 <-0
}else{
if(low.spell.lengths(d5yrs_record, thresh=0.00001)[[2]][1]>=1826){
Freq_d5<-0
}else{
Freq_d5 <-nrow(low.spell.lengths(d5yrs_record, thresh=0.00001))
}
}

HSpell=high.spell.lengths(d5yrs_record,threshold=0.0001) 
if(!is.na(HSpell[[1]][1])) {
d5yrs_length=mean(HSpell$spell.length) # most recent inundation start date MRI is Most Recent Inundation
d5yrs_lengthsd=sd(HSpell$spell.length)}else{
d5yrs_length=0
d5yrs_lengthsd=0}

CTF_5yrs=CTF(d5yrs_record,threshold=0.0001)

##### Extract all record to start of 2005 for flood frequency calculation

Record.start="01/01/2005"
Record.start=as.Date(as.character(Record.start), format = "%d/%m/%Y") 
dALL = doi-Record.start

dALL_record =DepOI[DepOI$Date %in% as.Date(Record.start:doi),] # extracts last inundation event from record
colnames(dALL_record)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
dALL_lengths=high.spell.lengths(dALL_record, threshold=0.00001, ind.days=5) # de
if(!is.na(dALL_lengths[[1]][1])) {
Freq_ALL=nrow(dALL_lengths) }else{
Freq_ALL=0}


if(tdatarows>0){
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),3] <- as.numeric(doi-TSLW, units = "days")}else{
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),3] <- "NA"
}

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),4] <- MRI_length
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),5] <- MRI_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),6] <- MRI_shallow
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),7] <- MRI_deep

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),8] <- d3Mon_wet/91
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),9] <- d3Mon_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),10] <- d3Mon_shallow/91
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),11] <- d3Mon_deep/91

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),12] <- d1yrs_wet/366
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),13] <- d1yrs_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),14] <- d1yrs_shallow/366
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),15] <- d1yrs_deep/366

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),16] <- d3yrs_wet/(365+365+366)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),17] <- d3yrs_meandepth
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),18] <- d3yrs_shallow/(365+365+366)
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),19] <- d3yrs_deep/(365+365+366)

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),20] <- Freq_d1
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),21] <- Freq_d3
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),22] <- Freq_d5
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),23] <- Freq_ALL

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),24] <- d1yrs_length
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),25] <- d3yrs_length
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),26] <- d5yrs_length

Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),27] <- CTF_1yrs[1][[1]]
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),28] <- CTF_1yrs[2][[1]]
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),29] <- CTF_3yrs[1][[1]]
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),30] <- CTF_3yrs[2][[1]]
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),31] <- CTF_5yrs[1][[1]]
Output[grepl(soi, Output$Site.ID, fixed=TRUE)& grepl(doi, Output$Date.of.collection),32] <- CTF_5yrs[2][[1]]



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
Output$DOI_MRI_start=Output$Date.of.collection-Output$MRI_length # tells us start of the MRI - to check whether it is being truncated by the record

### Codes altered to avoid confusion and allow merging with other datasets
Output$Site.ID <-gsub("LHT", "LHAT", Output$Site.ID) # Changed so that scripts were not confusing HT and LHT
Output$Site.ID <-gsub("CHT", "CCS", Output$Site.ID) # Changed so that scripts were not confusing CHT and HT
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
write.csv(Output , file = "Hydraulics_HTH_WL_medium and deep merged.csv")