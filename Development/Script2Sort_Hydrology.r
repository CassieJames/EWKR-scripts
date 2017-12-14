# Script to plot hydrology data and to determine basic hydrology metrics (long term flood frequency based on 30 year record and time since last flood)
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016
# Note: probably very long winded way of getting a reasonable approximation for flood frequency that takes pumping activities into account....

image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Euston flow data/"
setwd(image.dir) # set image directory for output


# Load hydrology files 

Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv")) # read in full record
Euston_discharge_samp=Euston_discharge[11324:14775,] # subset record to record from 2007 (must sort this code out so that it picks up the start of 2007 in a better way!)
colnames(Euston_discharge) = c("Date","Flow")
windowsFonts(A=windowsFont("Calibri"))
# set font styles

# Create hydrology and rainfall plot 

png(paste(image.dir,'Euston_hydrograph.png',sep=''), width=2000, height=3000, units="px", res=300)
par(mfrow=c(3,1))
par(mar=c(3,4,2,2)+0.5, mgp = c(0, 0.5, 0))

d <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
Flow<-Euston_discharge$Flow
Flow[is.na(Flow)] <- 0 
Years=seq(1976,2016,1)

plot(d,Euston_discharge$Flow, pch="", xlab="", ylab="", ylim=c(-10000,180000), tck=0.02, axes=FALSE, family="A")
axis(2, at=seq(0,180000,20000), lwd=0, las=2, tck=0)
axis(1, at=seq(d[1], d[14775], by = "year"), labels=Years, las=0, tck=0.02,family="A")
abline(h=36700, col="red", lwd=2)
lines(d, Flow, lwd=1)
mtext("Murray @ Euston (404203C)",side=3, family="A", cex=0.8,font=2, adj=0)
mtext("Stream Discharge ML/Day", side=4, line=3, cex=0.8, col="Black", family="A") # its loosing this off the side of the png so need to sort out margins

Years=seq(2007,2016,1)
ds <- as.Date(Euston_discharge_samp$Date, format="%d/%m/%Y") 
Flow_s<-Euston_discharge_samp$Flow
Flow_s[is.na(Flow_s)] <- 0 
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
Dates=data.frame(read.csv("HTH_FP_dates.csv"))
Dates2 <- within(Dates, Sample.date <- as.Date((Sample.date), format = "%d/%m/%Y")) # make sure r recognises dates as dates
Dates2$height <-c(-5000)

plot(ds,Euston_discharge_samp$Flow, pch="", xlab="", ylab="", ylim=c(-10000,180000), tck=0.02, axes=FALSE, family="A")
axis(2, at=seq(0,180000,20000), lwd=0, las=2, tck=0)
axis(1, at=seq(d[11324], d[14775], by = "year"), labels=Years, las=0, tck=0.02,family="A")
abline(h=36700, col="red", lwd=2) # puts a red line at the CTF value for the Hattah Lakes system for info.
lines(ds, Flow_s, lwd=1)
points(Dates2$Sample.date, Dates2$height, pch=17)
mtext("Murray @ Euston (404203C)",side=3, family="A", cex=0.8,font=2, adj=0)
mtext("Stream Discharge ML/Day", side=4, line=3, cex=0.8, col="Black", family="A")

Rainfall=data.frame(read.csv("IDCJAC0009_076043_1800_Data.csv"))
Rainfall2 <- within(Rainfall, Date <- as.Date(as.character(Date), format = "%d/%m/%Y"))
Rainfall_samp = Rainfall2[Rainfall2$Date %in% c(as.Date('2007-01-01'):as.Date('2016-06-13')),]

plot(Rainfall_samp$Date, Rainfall_samp$Rainfall, xlab=NA, ylab=NA,pch="", ylim=c(0,100), tck=0.02, axes=FALSE, family="A" )
axis(2, at=seq(0,100,10), lwd=0, las=2, tck=0)
axis(1, at=seq(d[11324], d[14775], by = "year"), labels=Years, las=0, tck=0.02,family="A")
lines(Rainfall_samp$Date, Rainfall_samp$Rainfall, lwd=1)
mtext("Rainfall @ NULKWYNE KIAMAL (76043)",side=3, family="A", cex=0.8,font=2, adj=0)

dev.off()

##########################################################################################################
### script to determine long term flood frequency for Hattah Lakes floodplain sites
library(hydrostats)

# step 1 - upload all required data (hydrology record from Euston, CTF data and samples dates)

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Euston flow data/"; setwd(data.dir)
Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Flood.dat=data.frame(read.csv("HFP_FLOOD_DAT2.csv")) # load hydrological information
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))

# Step 2 - create a data frame to receive results

sitelist=unique(mydata$Unique_site_year)
tdata=unique(mydata[,c("Unique_site_year","Site.ID","Date.of.collection", "Flow", "Inundated")])
tdata=merge(tdata, Flood.dat, by="Site.ID",all.x=TRUE)
tdata$Flood_frequency <- 0
tdata$TSLF <-1900 # place random estimate in to start with


# Step 3 - for each unique site and date combination work out frequency and duration of flooding for previous 30 years

for(s in sitelist) {

CTF=tdata[which(tdata$Unique_site_year==s),c("RIMFIM_GIS")] # extract CTF info for the site
CTF=CTF[1]  # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes - another data issue to sort out)
Site.ID=tdata[which(tdata$Unique_site_year==s),c("Site.ID")] # extract site id
Alt.CTF=Flood.dat[which(Flood.dat$Site.ID==Site.ID),c("Override_RIMFIM")] # Look up site id in flood dat table and extract override RIM FIM value if present
if(!is.na(Alt.CTF)) {CTF=Alt.CTF} # if Alt.CTF is not "NA" then replace CTF with Alt.CTF (this overrides the existing CTF where the CTF is contradicted by field obs)

Date=tdata[which(tdata$Unique_site_year==s),c("Date.of.collection")]
Date=as.Date(as.character(Date), format = "%d/%m/%Y") # make sure r recognises it as a date!
Date=Date[1] # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes)
d30 =Date-10950 # determine approximately 30 year span to extra flow data over - this is a cheat method which I will correct to an exact date extraction when I get time
d30=d30[1] # Problem with duplicates again
loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Date),] # subset flow data to 30 years of interest

# This stage adds the pumping activities
pump.dat=Flood.dat[which(Flood.dat$Site.ID==Site.ID),c("y2006", "y2010", "y2013", "y2014")] # Look up site id in flood.dat table and returns list of pumping dates
pump.dat[pump.dat==""]<-NA # replace empty cells with NA
pump.2006 =  pump.dat[1]
pump.2010 =  pump.dat[2]
pump.2013 =  pump.dat[3]
pump.2014 =  pump.dat[4]

pump.activity= matrix(NA,nrow=4, ncol=2)
rownames(pump.activity)=c('pump.2006', 'pump.2010', 'pump.2013', 'pump.2014')
colnames(pump.activity)=c('Date','Flow' )
pump.activity=as.data.frame(pump.activity) # change to data frame to allow mix of data types in matrix (i.e. date and numerical)
pump.activity$Date <- c(as.Date('2006-12-01'),as.Date('2010-12-01'),as.Date('2013-12-01'), as.Date('2014-09-24'))
if(grepl("TR6OFT", Site.ID)) {pump.activity$Date <- c(as.Date('2006-09-25'),as.Date('2010-12-01'),as.Date('2013-12-01'), as.Date('2014-09-24'))} # change date for TR6 as they would have flooded earlier in 2006 (sep)

if(!is.na(pump.2006)) {pump.activity["pump.2006","Flow"]=CTF+1000} # if pumping is not NA (i.e. it occurred) make flow equal to CTF+1000 (this allows r to see the pumping activity as a flow event)
if(!is.na(pump.2013)) {pump.activity["pump.2013","Flow"]=CTF+1000} # if pumping is not NA (i.e. it occurred) make flow equal to CTF+1000
if(!is.na(pump.2014)) {pump.activity["pump.2014","Flow"]=CTF+1000} # if pumping is not NA (i.e. it occurred) make flow equal to CTF+1000
pump.activity["pump.2010","Flow"]=0 # for the time being don't adjust for 2010 and leave flow as zero

loi30=rbind(loi30,pump.activity) # bind pumping activities to end of flow record
loi30=loi30[order(as.Date(loi30$Date, format="%Y-%m-%d")),] # sort data into time series for high spell analysis
loi30=loi30[loi30$Date<Date,] # remove pumping activities that occur after sampling date!
colnames(loi30)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
DOI=ts.format(loi30, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF, ind.days=60) # determine high spells for specified 30 year record based on CTF threshold - 60 day separation of peaks required to determine separate events
Duration$start.date <- as.Date(Duration$start.date, format = "%Y-%m-%d") # ensure dates are in same format
Duration=Duration[Duration$start.date %in% (d30:Date),] #subset output to only period of interest - because I tagged pump events onto loi30 its included events that happened AFTER the sampling - these need to be removed
Frequency=nrow(Duration) # thirty year flood frequency - number of spell events in output
TSLF=Duration[nrow(Duration),1] # date of last flood prior to sampling will be bottom row of returned output
TSLF=Date-TSLF # determine differences between date of sampling and date of last flood in number of days
TSLF=as.numeric(TSLF, units="days") # retrieve only value and not text
tdata[grep(s, tdata$Unique_site_year),c("Flood_frequency")] <- Frequency # place in data frame
if(length(TSLF)==0) {TSLF='NA'}
tdata[grep(s, tdata$Unique_site_year),c("TSLF")] <- TSLF # TSLF = time since last flood (this is the date of the first day when the flows at Euston gauge exceeded the CTF)
}

HTH_FP_Flood_data=tdata # copy data so the file name is more informative
write.csv(HTH_FP_Flood_data , file = "Flood_HTH_FP_pumps_corrected 60 day interval.csv") # save data out

##########################################################################################################
### script to determine long term flood frequency for Hattah Lakes wetland sites
library(hydrostats)

# step 1 - upload all required data (hydrology record from Euston, CTF data and samples dates)

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Euston flow data/"; setwd(data.dir)
plot.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/WL Hydrology RIMFIM hydrology plots/"

Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Flood.dat=data.frame(read.csv("HWL_FLOOD_DAT_Dec2017.csv")) # load hydrological information
Flood.dat$Y2005_Early <- as.Date(Flood.dat$Y2005_Early, format="%d/%m/%Y") 
Flood.dat$Y2005_Late <- as.Date(Flood.dat$Y2005_Late, format="%d/%m/%Y") 
Flood.dat$Y2006_Early <- as.Date(Flood.dat$Y2006_Early, format="%d/%m/%Y") 
Flood.dat$Y2006_Late <- as.Date(Flood.dat$Y2006_Late, format="%d/%m/%Y") 
Flood.dat$Y2009_Early <- as.Date(Flood.dat$Y2009_Early, format="%d/%m/%Y") 
Flood.dat$Y2009_Late <- as.Date(Flood.dat$Y2009_Late, format="%d/%m/%Y") 
Flood.dat$Y2010_Early <- as.Date(Flood.dat$Y2010_Early, format="%d/%m/%Y") 
Flood.dat$Y2010_Late <- as.Date(Flood.dat$Y2010_Late, format="%d/%m/%Y") 
Flood.dat$Y2013 <- as.Date(Flood.dat$Y2013, format="%d/%m/%Y") 
Flood.dat$Y2014 <- as.Date(Flood.dat$Y2014, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL.csv"))
date.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv"))
mydata2=merge(mydata, newdates, by="Unique_site_year") # this replaces the incorrect date of collection in the original veg database with new dates and adds longs and lats

# Step 2 - create a data frame to receive results

sitelist=unique(mydata2$Unique_site_year)
tdata=unique(mydata2[,c("Unique_site_year","Site.ID.x","Date.of.collection.y", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID.x",by.y="Site.ID",all.x=TRUE)
tdata$Flood_frequency <- 0
tdata$TSLF <-1900 # place random estimate in to start with


# Step 3 - for each unique site and date combination work out frequency and duration of flooding for previous 30 years

for(s in sitelist) {

CTF=tdata[which(tdata$Unique_site_year==s),c("RIMFIM_CTF")] # extract CTF info for the site
CTF=CTF[1]  # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes - another data issue to sort out)
Site.ID=tdata[which(tdata$Unique_site_year==s),c("Site.ID.x")] # extract site id
Date=tdata[which(tdata$Unique_site_year==s),c("Date.of.collection.y")]
Date=as.Date(as.character(Date), format = "%d/%m/%Y") # make sure r recognises it as a date!
Date=Date[1] # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes)
if(Date>as.Date(as.character("01/01/2012"), format = "%d/%m/%Y")) { # this part of the script finds Chalka Creek and if the sample date is after 2012 it changes the CTF to 25 000
if(grepl("CHT",s,fixed=TRUE)) {CTF=25000}
}
d30 =Date-10950 # determine approximately 30 year span to extra flow data over - this is a cheat method which I will correct to an exact date extraction when I get time
d30=d30[1] # Problem with duplicates again
loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Date),] # subset flow data to 30 years of interest


Site.ID <- factor(Site.ID, levels=levels(Flood.dat$Site.ID)) # quick fix because r wants the same number of factors for site.ID 
# This stage adds the pumping activities
pump.dat=Flood.dat[which(Flood.dat$Site.ID==Site.ID),] # Look up site id in flood.dat table and returns list of pumping dates

pump.2005.Early =  pump.dat[5]
pump.2005.Late =  pump.dat[6]
pump.2006.Early =  pump.dat[7]
pump.2006.Late =  pump.dat[8]
pump.2009.Early =  pump.dat[9]
pump.2009.Late =  pump.dat[10]
pump.2010.Early =  pump.dat[11]
pump.2010.Late =  pump.dat[12]
pump.2013 =  pump.dat[13]
pump.2014=  pump.dat[14]

pump.activity= matrix(NA,nrow=10, ncol=2)
rownames(pump.activity)=c("pump.2005.Early","pump.2005.Late","pump.2006.Early","pump.2006.Late","pump.2009.Early","pump.2009.Late","pump.2010.Early","pump.2010.Late", "pump.2013","pump.2014")
colnames(pump.activity)=c('Date','Flow' )
pump.activity=as.data.frame(pump.activity) # change to data frame to allow mix of data types in matrix (i.e. date and numerical)
pump.activity$Date <- as.data.frame(t(pump.dat[5:14]))[,1]


if(!is.na(pump.2005.Early )) {pump.activity["pump.2005.Early","Flow"]=CTF+1000} # if pumping is not NA (i.e. it occurred) make flow equal to CTF+1000 (this allows r to see the pumping activity as a flow event)
if(!is.na(pump.2005.Late)) {pump.activity["pump.2005.Late","Flow"]=CTF+1000} # if pumping is not NA (i.e. it occurred) make flow equal to CTF+1000
if(!is.na(pump.2006.Early)) {pump.activity["pump.2006.Early","Flow"]=CTF+1000} 
if(!is.na(pump.2006.Late)) {pump.activity["pump.2006.Late","Flow"]=CTF+1000} 
if(!is.na(pump.2009.Early)) {pump.activity["pump.2009.Early","Flow"]=CTF+1000} 
if(!is.na(pump.2009.Late)) {pump.activity["pump.2009.Late","Flow"]=CTF+1000} 
if(!is.na(pump.2009.Late)) {pump.activity["pump.2009.Late","Flow"]=CTF+1000} 
if(!is.na(pump.2010.Early)) {pump.activity["pump.2010.Early","Flow"]=CTF+1000} 
if(!is.na(pump.2010.Late)) {pump.activity["pump.2010.Late","Flow"]=CTF+1000} 
if(!is.na(pump.2013)) {pump.activity["pump.2013","Flow"]=CTF+1000} 
if(!is.na(pump.2014)) {pump.activity["pump.2014","Flow"]=CTF+1000} 

colnames(loi30)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
colnames(pump.activity)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
loi30=rbind(loi30,pump.activity) # bind pumping activities to end of flow record
loi30=loi30[order(as.Date(loi30$Date, format="%Y-%m-%d")),] # sort data into time series for high spell analysis
loi30=loi30[loi30$Date<Date,] # remove pumping activities that occur after sampling date!
loi30=na.omit(loi30)
DOI=ts.format(loi30, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF, ind.days=60) # determine high spells for specified 30 year record based on CTF threshold - 5 day separation of peaks required to determine separate events
Duration$start.date <- as.Date(Duration$start.date, format = "%Y-%m-%d") # ensure dates are in same format
Duration=Duration[Duration$start.date %in% (d30:Date),] #subset output to only period of interest - because I tagged pump events onto loi30 its included events that happened AFTER the sampling - these need to be removed
Frequency=nrow(Duration) # thirty year flood frequency - number of spell events in output
TSLF=Duration[nrow(Duration),1] # date of last flood prior to sampling will be bottom row of returned output
TSLF=Date-TSLF # determine differences between date of sampling and date of last flood in number of days
TSLF=as.numeric(TSLF, units="days") # retrieve only value and not text
tdata[grepl(s, tdata$Unique_site_year, fixed=TRUE),grepl("Flood_frequency", colnames(tdata), fixed=TRUE)] <- Frequency # place in data frame
if(length(TSLF)==0) {TSLF='NA'}
tdata[grepl(s, tdata$Unique_site_year, fixed=TRUE),grepl("TSLF", colnames(tdata), fixed=TRUE)] <- TSLF # TSLF = time since last flood (this is the date of the first day when the flows at Euston gauge exceeded the CTF)
#write.csv(loi30 , paste(plot.dir,s,".csv",sep="")) # save data out

}

HTH_WL_Flood_data=tdata # copy data so the file name is more informative
write.csv(HTH_WL_Flood_data , file = "Flood_HTH_WL_pumps_corrected_60 day interval using RIMFIM CTF.csv") # save data out