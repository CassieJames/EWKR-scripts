# Script to plot hydrology data and to determine basic hydrology metrics (long term flood frequency based on 30 year record and time since last flood)
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016


image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Euston flow data/"
setwd(image.dir) # set image directory for output

# Load hydrology files 

Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv")) # read in full record
Euston_discharge_samp=Euston_discharge[11324:14775,] # subset record to record from 2007
colnames(Euston_discharge) = c("Date","Flow")
windowsFonts(A=windowsFont("Calibri"))
# set font styles

# Create plot 

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
mtext("Stream Discharge ML/Day", side=4, line=3, cex=0.8, col="Black", family="A")


Years=seq(2007,2016,1)
ds <- as.Date(Euston_discharge_samp$Date, format="%d/%m/%Y") 
Flow_s<-Euston_discharge_samp$Flow
Flow_s[is.na(Flow_s)] <- 0 

plot(ds,Euston_discharge_samp$Flow, pch="", xlab="", ylab="", ylim=c(-10000,180000), tck=0.02, axes=FALSE, family="A")
axis(2, at=seq(0,180000,20000), lwd=0, las=2, tck=0)
axis(1, at=seq(d[11324], d[14775], by = "year"), labels=Years, las=0, tck=0.02,family="A")
abline(h=36700, col="red", lwd=2)
lines(ds, Flow_s, lwd=1)
points(SampDates$Date, SampDates$height, pch=17)
mtext("Murray @ Euston (404203C)",side=3, family="A", cex=0.8,font=2, adj=0)
mtext("Stream Discharge ML/Day", side=4, line=3, cex=0.8, col="Black", family="A")


data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
Dates=data.frame(read.csv("HTH_FP_dates.csv"))
Rainfall=data.frame(read.csv("IDCJAC0009_076043_1800_Data.csv"))
Dates2 <- within(Dates, Sample.date <- as.Date((Sample.date), format = "%d/%m/%Y")) # make sure r recognises dates as dates
Rainfall2 <- within(Rainfall, Date <- as.Date(as.character(Date), format = "%d/%m/%Y"))
Rainfall_samp = Rainfall2[Rainfall2$Date %in% c(as.Date('2007-01-01'):as.Date('2016-06-13')),]

plot(Rainfall_samp$Date, Rainfall_samp$Rainfall, xlab=NA, ylab=NA,pch="", ylim=c(0,100), tck=0.02, axes=FALSE, family="A" )
axis(2, at=seq(0,100,10), lwd=0, las=2, tck=0)
axis(1, at=seq(d[11324], d[14775], by = "year"), labels=Years, las=0, tck=0.02,family="A")
lines(Rainfall_samp$Date, Rainfall_samp$Rainfall, lwd=1)
points(SampDates$Date, SampDates$height, pch=17)
mtext("Rainfall @ NULKWYNE KIAMAL (76043)",side=3, family="A", cex=0.8,font=2, adj=0)

dev.off()



### script to determine long term flood frequency for HFP sites
library(hydrostats)

# step 1 - upload all required data (hydrology record from Euston, CTF data and samples dates)

image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Euston flow data/"; setwd(image.dir)
Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
Flood.dat=data.frame(read.csv("HFP_FLOOD_DAT.csv")) # flood information from field records etc
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))

# Step 2 - create a data frame to receive results

sitelist=unique(mydata$Unique_site_year)
tdata=unique(mydata[,c("Unique_site_year","Site.ID","Date.of.collection", "Flow", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID",by.y="HFP_SITES",all.x=TRUE)
tdata$Flood_frequency <- 0
tdata$TSLF <-1900 # place random estimate in to start with


# Step 3 - for each unique site and date combination work out frequency and duration of flooding for previous 30 years

for(s in sitelist) {

CTF=tdata[which(tdata$Unique_site_year==s),c("RIMFIM_GIS")] # extract CTF info for the site
CTF=CTF[1]  # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes)
Date=tdata[which(tdata$Unique_site_year==s),c("Date.of.collection")]
Date=as.Date(as.character(Date), format = "%d/%m/%Y") # make sure r recognises it as a date!
Date=Date[1] # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes)
d30 =Date-10950 # determine approximately 30 year span to extra flow data over - this is a cheat method which i will correct to an exact date extraction when I get time
d30=d30[1] # DITTO
loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Date),] # subset flow data to 30 years of interest
colnames(loi30)=c("Date", "Q") # rename columns to high.spell.lengths function recognises name
DOI=ts.format(loi30, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF) # determine high spells for specified 30 year record based on CTF threshold - 5 day separation of peaks required to determine separate events
Frequency=nrow(Duration) # thirty year flood frequency - number of spell events in output
TSLF=Duration[nrow(Duration),1] # date of last flood prior to sampling will be bottom row of returned output
TSLF <- as.Date((TSLF), format = "%Y-%m-%d") # ensure dates are in same format
TSLF=Date-TSLF # determine differences between date of sampling and date of last flood as number of days
TSLF=as.numeric(TSLF, units="days") # retrieve only value and not text
tdata[grep(s, tdata$Unique_site_year),c("Flood_frequency")] <- Frequency # place in dataframe
tdata[grep(s, tdata$Unique_site_year),c("TSLF")] <- TSLF
}

HTH_FP_Flood_data=tdata # copy data
write.csv(HTH_FP_Flood_data , file = "Flood_HTH_FP.csv") # save data