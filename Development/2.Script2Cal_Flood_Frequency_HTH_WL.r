
##########################################################################################################
### script to determine long term flood frequency for Hattah Lakes wetland sites
### this script calculates FF up until 1/1/2005 when more accurate modelling including pumping activities can be added
library(hydrostats)

# step 1 - upload all required data (hydrology record from Euston, CTF data and samples dates)

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Euston flow data/"; setwd(data.dir)
plot.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/WL Hydrology RIMFIM hydrology plots/"

Euston_discharge<-data.frame(read.csv(file="MeanWaterFlow_Euston.csv"))
Euston_discharge$Date <- as.Date(Euston_discharge$Date, format="%d/%m/%Y") 
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Flood.dat=data.frame(read.csv("HWL_FLOOD_DAT_Dec2017.csv")) # load hydrological information 
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL.csv"))
date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv"))
mydata2=merge(mydata, newdates, by="Unique_site_year") # this replaces the incorrect date of collection in the original veg database with new dates and adds longs and lats

# 30 years FF
# create a data frame to receive results

sitelist=unique(mydata2$Unique_site_year)
tdata=unique(mydata2[,c("Unique_site_year","Site.ID.x","Date.of.collection.y", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID.x",by.y="Site.ID",all.x=TRUE)
tdata$Flood_frequency <- 0


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
Record.start="01/01/2005"
Record.start=as.Date(as.character(Record.start), format = "%d/%m/%Y") 
d30 =Date-round((365.242189*30),0); d30=d30[1] # determine 30 year span and deal with problem with duplicates again

loi30=Euston_discharge[Euston_discharge$Date %in% (d30:Record.start),] # subset flow data to 30 years of interest but remove record post 1/1/2005
colnames(loi30)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
DOI=ts.format(loi30, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF, ind.days=5) # determine high spells for specified 30 year record based on CTF threshold - 5 day separation of peaks required to determine separate events
Frequency=nrow(Duration) # thirty year flood frequency - number of spell events in output
tdata[grepl(s, tdata$Unique_site_year, fixed=TRUE),grepl("Flood_frequency", colnames(tdata), fixed=TRUE)] <- Frequency # place in data frame

}
HTH_WL_Flood_data_d30=tdata # copy data so the file name is more informative

###########
# repeat for shorter time frames 20 years
# create a data frame to receive results

sitelist=unique(mydata2$Unique_site_year)
tdata=unique(mydata2[,c("Unique_site_year","Site.ID.x","Date.of.collection.y", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID.x",by.y="Site.ID",all.x=TRUE)
tdata$Flood_frequency <- 0

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
Record.start="01/01/2005"
Record.start=as.Date(as.character(Record.start), format = "%d/%m/%Y") 
d20 =Date-round((365.242189*20),0) # determine 20 year span
d20=d20[1] # Problem with duplicates again

loi20=Euston_discharge[Euston_discharge$Date %in% (d20:Record.start),] # subset flow data to 30 years of interest but remove record post 1/1/2005
colnames(loi20)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
DOI=ts.format(loi20, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF, ind.days=5) # determine high spells for specified 30 year record based on CTF threshold - 5 day separation of peaks required to determine separate events
Frequency=nrow(Duration) # thirty year flood frequency - number of spell events in output
tdata[grepl(s, tdata$Unique_site_year, fixed=TRUE),grepl("Flood_frequency", colnames(tdata), fixed=TRUE)] <- Frequency # place in data frame


}
HTH_WL_Flood_data_d20=tdata # copy data so the file name is more informative

###########
# repeat for shorter time frames 10 years
# create a data frame to receive results

sitelist=unique(mydata2$Unique_site_year)
tdata=unique(mydata2[,c("Unique_site_year","Site.ID.x","Date.of.collection.y", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID.x",by.y="Site.ID",all.x=TRUE)
tdata$Flood_frequency <- 0
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
Record.start="01/01/2005"
Record.start=as.Date(as.character(Record.start), format = "%d/%m/%Y") 
d10 =Date-round((365.242189*10),0) # determine 10 year span
d10=d10[1] # Problem with duplicates again

if(d10<Record.start){ # make sure that its only collecting data prior to Bigmod data)
loi10=Euston_discharge[Euston_discharge$Date %in% (d10:Record.start),] # subset flow data to 10 years of interest but remove record post 1/1/2005
colnames(loi10)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
DOI=ts.format(loi10, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF, ind.days=5) # determine high spells for specified 30 year record based on CTF threshold - 5 day separation of peaks required to determine separate events
Frequency=nrow(Duration) # thirty year flood frequency - number of spell events in output
} else {Frequency=0}
tdata[grepl(s, tdata$Unique_site_year, fixed=TRUE),grepl("Flood_frequency", colnames(tdata), fixed=TRUE)] <- Frequency # place in data frame
}
HTH_WL_Flood_data_d10=tdata # copy data so the file name is more informative

###########
# repeat for shorter time frames 5 years
# create a data frame to receive results

sitelist=unique(mydata2$Unique_site_year)
tdata=unique(mydata2[,c("Unique_site_year","Site.ID.x","Date.of.collection.y", "Inundated")])
tdata=merge(tdata, Flood.dat, by.x="Site.ID.x",by.y="Site.ID",all.x=TRUE)
tdata$Flood_frequency <- 0

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
Record.start="01/01/2005"
Record.start=as.Date(as.character(Record.start), format = "%d/%m/%Y") 
d5 =Date-round((365.242189*5),0) # determine 5 year span
d5=d5[1] # Problem with duplicates again

if(d5<Record.start){ # make sure that its only collecting data prior to Bigmod data)
loi5=Euston_discharge[Euston_discharge$Date %in% (d5:Record.start),] # subset flow data to 5 years of interest but remove record post 1/1/2005
colnames(loi5)=c("Date", "Q") # rename columns so high.spell.lengths function recognises name
DOI=ts.format(loi5, format="%Y-%m-%d") # change into standard date format so that function recognises date
Duration=high.spell.lengths(DOI, threshold=CTF, ind.days=5) # determine high spells for specified 30 year record based on CTF threshold - 5 day separation of peaks required to determine separate events
Frequency=nrow(Duration) # thirty year flood frequency - number of spell events in output
} else {Frequency=0}
tdata[grepl(s, tdata$Unique_site_year, fixed=TRUE),grepl("Flood_frequency", colnames(tdata), fixed=TRUE)] <- Frequency # place in data frame
}
HTH_WL_Flood_data_d5=tdata # copy data so the file name is more informative

HTH_WL_Flood_data=cbind(HTH_WL_Flood_data_d5, HTH_WL_Flood_data_d10[, ncol(HTH_WL_Flood_data_d10)], HTH_WL_Flood_data_d20[, ncol(HTH_WL_Flood_data_d20)], HTH_WL_Flood_data_d30[, ncol(HTH_WL_Flood_data_d30)])
setwd("C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/")

colnames(HTH_WL_Flood_data)[21] <- "d30yearsFF"
colnames(HTH_WL_Flood_data)[20] <- "d20yearsFF"
colnames(HTH_WL_Flood_data)[19] <- "d10yearsFF"
colnames(HTH_WL_Flood_data)[18] <- "d5yearsFF"

write.csv(HTH_WL_Flood_data , file = "Flood_HTH_WL_flood frequency prior to 2005V3.csv") # save data out