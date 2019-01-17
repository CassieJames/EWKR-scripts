##########################################################################################################
### script to determine long term flood frequency for Hattah Lakes wetland sites
library(hydrostats)

# step 1 - upload all required data (hydrology record from Euston, CTF data and samples dates)

data.dir="C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Euston flow data/"; setwd(data.dir)

Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
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
Flood.dat$Site.ID=gsub("KRT","KT",Flood.dat$Site.ID)
Flood.dat$Site.ID=gsub("CCNT","NCT",Flood.dat$Site.ID)


data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL.csv"))
mydata$Site.ID=gsub("CCNT","NCT",mydata$Site.ID) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie - my attempt to make everything consistent
mydata$Site.ID=gsub("KRT","KT",mydata$Site.ID)
mydata$Unique_site_year=gsub("CCNT","NCT",mydata$Unique_site_year) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie
mydata$Unique_site_year=gsub("KRT","KT",mydata$Unique_site_year)

date.dir="C:/Users/jc246980/Documents/MD Vegetation/Environmental data/"; setwd (date.dir) 
newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)

mydata2=merge(mydata, newdates, by="Unique_site_year" ) # this replaces the incorrect date of collection in the original veg database with new dates and adds longs and lats

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
Date=Date[1] # where there are duplicates take only the first date (these are invariably the same site but with contradicting inundation field notes)
if(Date>as.Date(as.character("01/01/2012"), format = "%d/%m/%Y")) { # this part of the script finds Chalka Creek and if the sample date is after 2012 it changes the CTF to 25 000
if(grepl("CHT",s,fixed=TRUE)) {CTF=25000}
}
d30 =Date-10950 # determine approximately 30 year span to extra flow data over - this is a cheat method which I will correct to an exact date extraction when I get time
d30=d30[1] # Problem with duplicates again
loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Date),] # subset flow data to 30 years of interest

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