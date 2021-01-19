###################################################################################################################################
### Script to determine long term flood frequency for Hattah Lakes wetland sites
### Written by  C.S.James (JCU, TropWATER)
### GNU General Public License .. feel free to use / distribute ... no warranties
### 27th August 2019
###################################################################################################################################
### Notes: This script calculates FF up until 1/1/2005 when more accurate modelling including pumping activities can be added using 
### bigmod output. I have used the RIMFIM CTF (reduced for Chalka in around 2012/2013) and the flood duration required to flood 
### according to Butcher and Hale (2011)
###################################################################################################################################
# Load libraries
library(hydrostats)
library(lubridate)

# step 1 - upload all required data (hydrology record from Euston, CTF data and samples dates)
data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Euston flow data/"; setwd(data.dir)

Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Flood.dat=data.frame(read.csv("HWL_FLOOD_DAT_Dec2017.csv")) # load hydrological information 

Flood.dat$Site.ID=gsub("KRT","KT",Flood.dat$Site.ID)
Flood.dat$Site.ID=gsub("CCNT","NCT",Flood.dat$Site.ID)
Flood.dat$Site.ID=gsub("CHT","CCS",Flood.dat$Site.ID)
Flood.dat$Site.ID=gsub("LHT","LHAT",Flood.dat$Site.ID)

Flood.dat$Wetland=gsub("KRT","KT",Flood.dat$Wetland)
Flood.dat$Wetland=gsub("CCNT","NCT",Flood.dat$Wetland)
Flood.dat$Wetland=gsub("CHT","CCS",Flood.dat$Wetland)
Flood.dat$Wetland=gsub("LHT","LHAT",Flood.dat$Wetland)

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
newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv"))
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CHT","CCS",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("LHT","LHAT",newdates$Unique_site_year)

mydata2=merge(mydata, newdates, by="Unique_site_year") # this replaces the incorrect date of collection in the original veg database with new dates and adds longs and lats

##############################################################################################
# Time frames 30 years
# create a data frame to receive results

sitelist=unique(mydata2$Unique_site_year)
tdata=unique(mydata2[,c("Unique_site_year","Site.ID.x","Date.of.collection.y", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID.x",by.y="Site.ID",all.x=TRUE)
tdata$Flood_frequency <- 0

Hattah =sitelist[grepl("HT|CCS",sitelist)]
Lhat =sitelist[grepl("LHAT",sitelist)]
Mournpall =sitelist[grepl("MOT",sitelist)] 
Brockie = sitelist[grepl("BRT",sitelist)]
Bitterang = sitelist[grepl("BIT|NCT",sitelist)]
Kramen = sitelist[grepl("KT",sitelist)]
Bulla = sitelist[grepl("BLT",sitelist)]
Yerang =sitelist[grepl("YT",sitelist)] 
NipNip = sitelist[grepl("BOT|NN",sitelist)]

for(s in sitelist) {


CTF=tdata[which(tdata$Unique_site_year==s),c("RIMFIM_CTF")] # extract CTF info for the site
CTF=CTF[1]  # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes - another data issue to sort out)
Site.ID=tdata[which(tdata$Unique_site_year==s),c("Site.ID.x")] # extract site id
Date=tdata[which(tdata$Unique_site_year==s),c("Date.of.collection.y")]
Date=as.Date(as.character(Date), format = "%d/%m/%Y") # make sure r recognises it as a date!
Date=Date[1] # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes)
if(Date>as.Date(as.character("01/01/2012"), format = "%d/%m/%Y")) { # this part of the script finds Chalka Creek and if the sample date is after 2012 it changes the CTF to 25 000
if(grepl("CCS",s,fixed=TRUE)) {CTF=25000}
}
Record.start="01/01/2005"
Record.start=as.Date(as.character(Record.start), format = "%d/%m/%Y") 
d30 =Date-round((365.242189*30),0); d30=d30[1] # determine 30 year span and deal with problem with duplicates again

loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Record.start),] # subset flow data to 30 years of interest but remove record post 1/1/2005
colnames(loi30)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
DOI=ts.format(loi30, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF, ind.days=30) # determine high spells for specified 30 year record based on CTF threshold - 5 day separation of peaks required to determine separate events

# Water takes several days to reach Lake Lockie and then a further 3 weeks to reach other lakes
if(s %in% Hattah){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(25) 
eof=Frequency[i,1]+days(Frequency[i,2]+730) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])

}

if(s %in% Lhat){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(25) 
eof=Frequency[i,1]+days(Frequency[i,2]+365) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}

if(s %in% Mournpall){
Frequency=(Duration[which(Duration$spell.length>=30),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(30) 
eof=Frequency[i,1]+days(Frequency[i,2]+730) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}

if(s %in% Yerang){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(25) 
eof=Frequency[i,1]+days(Frequency[i,2]+450) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}

if(s %in% Brockie){
Frequency=(Duration[which(Duration$spell.length>=30),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(30) 
eof=Frequency[i,1]+days(Frequency[i,2]+365) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}

# Info sheet on RAMSAR wetlands states that water reaches furthest lakes around 1 month (several days to reach Lake Lockie and a further 3 weeks to reach further
if(s %in% NipNip){
Frequency=(Duration[which(Duration$spell.length>=42),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(42) 
eof=Frequency[i,1]+days(Frequency[i,2]+360) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}

# Moxham and kenny 2016 Lake Bitterang report states that water resides in lakes for an average of 2 years according to SKM report from 2004
if(s %in% Bitterang){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(25) 
eof=Frequency[i,1]+days(Frequency[i,2]+(365*2)) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}

if(s %in% Bulla){
Frequency=(Duration[which(Duration$spell.length>=25),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(25) 
eof=Frequency[i,1]+days(Frequency[i,2]+(550)) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}

if(s %in% Kramen){
Frequency=(Duration[which(Duration$spell.length>=15),])
if (nrow(Frequency) >0){
for (i in 1:nrow(Frequency)){ 
sof=Frequency[i,1]+days(15) 
eof=Frequency[i,1]+days(Frequency[i,2]+(365+60)) # estimate 
int <- interval(ymd(sof), ymd(eof))
DOI[which(DOI$Date %within% int),2] <- "F"
}}
DOI[which(DOI$Q != "F"),2] <-0
DOI[which(DOI$Q == "F"),2] <-1 # turn flow data and 'F's to 1s and 0s for ease of calculation
DOI$Q=as.numeric(DOI$Q)
agg=aggregate(Q ~ cut(Date, "1 year"), DOI, summary) # aggregates data by year
cols <- lapply(c(1, 4) , function(x) agg)
tbl <- do.call(data.frame, cols)
tbl=as.data.frame((tbl[,2]))
Frequency=nrow(tbl[(tbl$'3rd Qu.'==1),])
}



tdata[grepl(s, tdata$Unique_site_year, fixed=TRUE),grepl("Flood_frequency", colnames(tdata), fixed=TRUE)] <- Frequency # place in data frame
}
HTH_WL_Flood_data_d30=tdata # copy data so the file name is more informative



write.csv(HTH_WL_Flood_data_d30 , file = "Flood_HTH_WL_flood frequency prior to 2005 December 2020.csv") # save data out