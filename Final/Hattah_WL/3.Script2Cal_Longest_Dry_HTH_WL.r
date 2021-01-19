###################################################################################################################################
### Script to calculate long term cease to flow metrics
### Written by  C.S.James (JCU, TropWATER)
### GNU General Public License .. feel free to use / distribute ... no warranties
### 27th August 2019
###################################################################################################################################
### Notes: Script involves cobbling together the hydrodynamic model outputs (starts in 1/1/2005) with the RIMFIM CTF data based on the Euston gauged data 
### 	   There are particular issues where the time frames overlap the bounadry between the two hydrology datasets. Also needed to account for the 
###        different durations that are needed to exceed before wetlands begin to fill (used estimates for Butcher and Hale, 2011). Therefore dry periods 
###        were estimated as CTF + period of time needed to dry (estimated as 365 days)
###
###################################################################################################################################
# Load libraries

library(hydrostats)
library(lubridate)
library(tidyr)
library(ggplot2)
library(zoo) # sorts out date issues with origin needing to be supplied error

# step 1 - Bring in wetland data and tidy up for dates and site names as some site labels are different between different datasets.

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

tada=merge(mydata, newdates, by="Unique_site_year")

# step 2 - upload Bigmod hydrodynamic data (available from 1/1/2005)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Depth.dat=read.csv("Hattah_sites_time_series_dataV3.csv", check.names=FALSE) # load hydrological information but prevent r from converting column names to a more 'friendly' format
colnames(Depth.dat)=gsub("KRT","KT",colnames(Depth.dat))
colnames(Depth.dat)=gsub("CCNT","NCT",colnames(Depth.dat))
colnames(Depth.dat)=gsub("CHT","CCS",colnames(Depth.dat))
colnames(Depth.dat)=gsub("LHT","LHAT",colnames(Depth.dat))

# step 3 - upload Euston gauge data, merge with CTF data and tidy up names

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Euston flow data/"; setwd(data.dir)
Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Flood.dat=data.frame(read.csv("HWL_FLOOD_DAT_Dec2017.csv")) # load hydrological information regarding CTF
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 

Flood.dat$Site.ID=gsub("KRT","KT",Flood.dat$Site.ID)
Flood.dat$Site.ID=gsub("CCNT","NCT",Flood.dat$Site.ID)
Flood.dat$Site.ID=gsub("CHT","CCS",Flood.dat$Site.ID)
Flood.dat$Site.ID=gsub("LHT","LHAT",Flood.dat$Site.ID)
mydataEust=merge(tada, Flood.dat, by.x="Site.ID.x",by.y="Site.ID",all.x=TRUE)

Sites_sampled <-unique(tada[,c("Site.ID.x","Date.of.collection.y", "Unique_site_year","Inundated")]) # creates unique list of sites and dates

# step 4 - Create empty data frame to receive results

Output= matrix(NA,nrow=nrow(Sites_sampled), ncol=12)
colnames(Output)=c("Site.ID", "Date.of.collection",  "Unique_site_year", "TSLW_LT","Max.CTF.30yrs", "P.CTF.30yrs","Max.CTF.10yrs", "P.CTF.10yrs", "n.events.30yrs", "FF.30","n.events.10yrs", "FF.10")
Output=as.data.frame(Output)
Output[,1] <-Sites_sampled$Site.ID.x 
Output[,2] <-Sites_sampled$Date.of.collection.y
Output[,3] <-Sites_sampled$Unique_site_year

sitelist=Sites_sampled$Unique_site_year


Hattah =sitelist[grepl("HT|CCS",sitelist)]
Lhat =sitelist[grepl("LHAT",sitelist)]
Mournpall =sitelist[grepl("MOT",sitelist)] 
Brockie = sitelist[grepl("BRT",sitelist)]
Bitterang = sitelist[grepl("BIT|NCT",sitelist)]
Kramen = sitelist[grepl("KT",sitelist)]
Bulla = sitelist[grepl("BLT",sitelist)]
Yerang =sitelist[grepl("YT",sitelist)] 
NipNip = sitelist[grepl("BOT|NN",sitelist)]



# step 5 - Loop through sites

for (i in 1:nrow(Sites_sampled)){ # 
soi=Sites_sampled[i,1]

# extract from Bigmod
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

s=Sites_sampled[i,3]

# extract from RIMFIM
CTF=mydataEust[which(mydataEust$Unique_site_year==s),c("RIMFIM_CTF")] # extract CTF info for the site
CTF=CTF[1]  # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes - another data issue to sort out)
Date=unique(mydataEust[which(mydataEust$Unique_site_year==s),c("Date.of.collection.y")])
Date=as.Date(as.character(Date), format = "%Y-%m-%d") # make sure r recognises it as a date!
Date=Date[1] # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes)
if(Date>as.Date(as.character("01/01/2012"), format = "%d/%m/%Y")) { # this part of the script finds Chalka Creek and if the sample date is after 2012 it changes the CTF to 25 000
if(grepl("CCS",s,fixed=TRUE)) {CTF=25000}
}

Record.start="31/12/2004"
Record.start=as.Date(as.character(Record.start), format = "%d/%m/%Y") 

# 30 years

d30 =Date-round((365.242189*30),0) # determine 30 year span from sampling date
d30=d30[1] # Problem with duplicates again

TSLWdata=DepOI[DepOI$Date %in% d30:doi,] # extracts data up until date of interest

loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Record.start),] # subset flow data to 30 years of interest but remove record post 1/1/2005
colnames(loi30)=c("Date", "Q") # rename columns for ease
DOI.30=ts.format(loi30, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI.30, threshold=CTF, ind.days=30) 

if(s %in% Hattah){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) 
eof=Frequency[j,1]+days(Frequency[j,2]+730) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F"

}}}

if(s %in% Lhat){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) 
eof=Frequency[j,1]+days(Frequency[j,2]+365) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F"

}}}

if(s %in% Mournpall){
Frequency=(Duration[which(Duration$spell.length>=30),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(30) 
eof=Frequency[j,1]+days(Frequency[j,2]+730) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F"

}}}

if(s %in% Yerang){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) 
eof=Frequency[j,1]+days(Frequency[j,2]+450) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F"

}}}

if(s %in% Brockie){
Frequency=(Duration[which(Duration$spell.length>=30),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(30) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+365) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% NipNip){
Frequency=(Duration[which(Duration$spell.length>=42),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(42) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+365) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% Bitterang){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+(730)) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% Bulla){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+(550)) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% Kramen){
Frequency=(Duration[which(Duration$spell.length>=15),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(15) 
eof=Frequency[j,1]+days(Frequency[j,2]+425) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.30[which(DOI.30$Date %within% int),2] <- "F"
}}}

DOI.30[which(DOI.30$Q != "F"),2] <-0
DOI.30[which(DOI.30$Q == "F"),2] <-1 # trun flow data and 'F's to 1s and 0s for ease of calculation

colnames(TSLWdata)=c("Date", "Q") # rename columns for ease
full.record=rbind(DOI.30,TSLWdata) # bind recent hydrodynmic data to estimates of inundation from rimfim modelling
full.record$Q=as.numeric(full.record$Q)

full.record[which(full.record$Q >1),2] <-1 
high_spells_30=high.spells(full.record, thresh=0.00001,plot = FALSE)#
agg=aggregate(Q ~ cut(Date, "1 year"), full.record, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
FF.30=nrow(tbl[(tbl$'Median'>=0.5),])

CTF.30.yrs=CTF(full.record, threshold = 0.001)


tdata=full.record[which(full.record$Q>0),] # subsets data to only dates with water
tdatarows=nrow(tdata)
if(tdatarows>0){
TSLW_LT=tdata[nrow(tdata),"Date"]+1
TSLW_LT=as.Date(as.character(TSLW_LT), format = "%Y-%m-%d") 
TSLW_LT=(Date-TSLW_LT)
TSLW_LT=TSLW_LT[[1]]}else{
TSLW_LT="NA"}

Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("TSLW_LT", colnames(Output), fixed=TRUE)] <- TSLW_LT


#### 

d10=Date-round((365.242189*10),0) # determine 10 year span from sampling date
d10=d10[1] # Problem with duplicates again

TSLWdata=DepOI[DepOI$Date %in% d10:doi,] # extracts data up until date of interest

loi10=Euston_discharge[Euston_discharge$Date %in% (d10:Record.start),] # subset flow data to 30 years of interest but remove record post 1/1/2005
colnames(loi10)=c("Date", "Q") # rename columns for ease
DOI.10=ts.format(loi10, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI.10, threshold=CTF, ind.days=30) 

if(s %in% Hattah){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) 
eof=Frequency[j,1]+days(Frequency[j,2]+730) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F"

}}}

if(s %in% Lhat){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) 
eof=Frequency[j,1]+days(Frequency[j,2]+365) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F"

}}}

if(s %in% Mournpall){
Frequency=(Duration[which(Duration$spell.length>=30),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(30) 
eof=Frequency[j,1]+days(Frequency[j,2]+730) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F"

}}}

if(s %in% Yerang){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) 
eof=Frequency[j,1]+days(Frequency[j,2]+450) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F"

}}}

if(s %in% Brockie){
Frequency=(Duration[which(Duration$spell.length>=30),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(30) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+365) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% NipNip){
Frequency=(Duration[which(Duration$spell.length>=42),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(42) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+365) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% Bitterang){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+(730)) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% Bulla){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(25) # start of inundation period
eof=Frequency[j,1]+days(Frequency[j,2]+(550)) # estimate end of inundation period
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F" # replace flow data with F where  periods estimated
}}}

if(s %in% Kramen){
Frequency=(Duration[which(Duration$spell.length>=15),])
if (nrow(Frequency) >0){
for (j in 1:nrow(Frequency)){ 
sof=Frequency[j,1]+days(15) 
eof=Frequency[j,1]+days(Frequency[j,2]+425) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI.10[which(DOI.10$Date %within% int),2] <- "F"
}}}


DOI.10[which(DOI.10$Q != "F"),2] <-0
DOI.10[which(DOI.10$Q == "F"),2] <-1 # trun flow data and 'F's to 1s and 0s for ease of calculation

colnames(TSLWdata)=c("Date", "Q") # rename columns for ease
full.record=rbind(DOI.10,TSLWdata) # bind recent hydrodynmic data to estimates of inundation from rimfim modelling
full.record$Q=as.numeric(full.record$Q)

full.record[which(full.record$Q >1),2] <-1 
high_spells_10=high.spells(full.record, thresh=0.00001,plot = FALSE)#
agg=aggregate(Q ~ cut(Date, "1 year"), full.record, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
FF.10=nrow(tbl[(tbl$'Median'>=0.5),])

CTF.10.yrs=CTF(full.record, threshold = 0.001)


#### 
Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("FF.30", colnames(Output), fixed=TRUE)] <- FF.30
Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("n.events.30yrs", colnames(Output), fixed=TRUE)] <- high_spells_30$n.events
Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("Max.CTF.30yrs", colnames(Output), fixed=TRUE)] <- CTF.30.yrs[5]
Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("P.CTF.30yrs", colnames(Output), fixed=TRUE)] <- CTF.30.yrs[1]

Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("FF.10", colnames(Output), fixed=TRUE)] <- FF.10
Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("n.events.10yrs", colnames(Output), fixed=TRUE)] <- high_spells_10$n.events
Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("Max.CTF.10yrs", colnames(Output), fixed=TRUE)] <- CTF.10.yrs[5]
Output[grepl(s, Output$Unique_site_year, fixed=TRUE),grepl("P.CTF.10yrs", colnames(Output), fixed=TRUE)] <- CTF.10.yrs[1]

}


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
write.csv(Output , file = "Flood_HTH_WL_Cease2flowDec2020.csv") # save data out