
# Script to plot hydrology data
# C James
# 29th July 2016

Dates.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/"; setwd(Dates.dir)
SampDates<-data.frame(read.csv(file="HTH_FP_samplDates.csv")) # Read in sampling dates
SampDates$height = -10000
SampDates$Date<- as.Date(SampDates$Date, format="%d/%m/%Y") # make sure r recognises Date as a date and not as a factor
image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Euston flow data/"
setwd(image.dir) # set image directory for output


# Load hydrology files for MDB

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
#mtext("Date", side=1, line=1.5, cex=0.8, col="Black", family="A")



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
library(hydrostats))

# step 1 - upload all required data (hydrology record from Euston, CTF data and samples dates)

image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Euston flow data/"
setwd(image.dir)
Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
Flood.dat=data.frame(read.csv("HFP_FLOOD_DAT.csv"))

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))

# Step 2 - create a data frame to receive results

sitelist=unique(mydata$Unique_site_year)
tdata=unique(mydata[,c("Unique_site_year","Site.ID","Date.of.collection", "Flow", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID",by.y="HFP_SITES",all.x=TRUE)
tdata$Flood_frequency <- 0
tdata$TSLF <-1900


# Step 3 - for each unique site and date combination work out frequency and duration of flooding for previous 30 years
for(s in sitelist) {

CTF=tdata[which(tdata$Unique_site_year==s),c("RIMFIM_GIS")]
CTF=CTF[1]
Date=tdata[which(tdata$Unique_site_year==s),c("Date.of.collection")]
Date=as.Date(as.character(Date), format = "%d/%m/%Y") # make sure r recognises it as a date!
Date=Date[1]
d30 =Date-10950 	# approximately 30 year span to extra flow data over
d30=d30[1]
loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Date),] # subset flow data to 30 years of interest
colnames(loi30)=c("Date", "Q")
DOI=ts.format(loi30, format="%Y-%m-%d")
Duration=high.spell.lengths(DOI, threshold=CTF)
Frequency=nrow(Duration) # thirty year flood frequency
TSLF=Duration[nrow(Duration),1]
TSLF <- as.Date((TSLF), format = "%Y-%m-%d")
TSLF=Date-TSLF
TSLF=as.numeric(TSLF, units="days")
tdata[grep(s, tdata$Unique_site_year),c("Flood_frequency")] <- Frequency
tdata[grep(s, tdata$Unique_site_year),c("TSLF")] <- TSLF
}

HTH_FP_Flood_data=tdata
write.csv(HTH_FP_Flood_data , file = "Flood_HTH_FP.csv")