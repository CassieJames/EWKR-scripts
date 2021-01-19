###################################################################################################################################
### Script to determine short term flood frequency for Hattah Lakes wetland sites
### Written by  C.S.James (JCU, TropWATER)
### GNU General Public License .. feel free to use / distribute ... no warranties
### 27th August 2019
###################################################################################################################################
### Notes: Script aggregates environmental data into a single file
###################################################################################################################################
date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 

newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("LHT","LHAT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CHT","CCS",newdates$Unique_site_year)
newdates$Site.ID=gsub("KRT","KT",newdates$Site.ID)
newdates$Site.ID=gsub("CCNT","NCT",newdates$Site.ID)
newdates$Site.ID=gsub("LHT","LHAT",newdates$Site.ID)
newdates$Site.ID=gsub("CHT","CCS",newdates$Site.ID)

#create season variable and sort out unique ID columns to match cleaned species matrix
library(zoo)
      yq <- as.yearqtr(as.yearmon(newdates$Date.of.collection , "%m/%d/%Y") + 1/12)
      newdates$Season <- factor(format(yq, "%q"), levels = 1:4, 
      labels = c("SU", "AU", "WI", "SP"))
newdates$Unique_site_year_season = paste(newdates$Unique_site_year,"_",newdates$Season,sep="")

newdates$Unique_site_year_season=paste("HAT_",newdates$Unique_site_year_season,sep="")
newdates$Unique_site_year_season=gsub("_13","_2013",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_16","_2016",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_15","_2015",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_14","_2014",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_12","_2012",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_11","_2011",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_10","_2010",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_09","_2009",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_08","_2008",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_07","_2007",newdates$Unique_site_year_season)
newdates$Unique_site_year_season=gsub("_06","_2006",newdates$Unique_site_year_season)


#Hydrology
data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd(data.dir)
Hydrodata=data.frame(read.csv("Flood_HTH_WL_flood frequency prior to 2005 December 2020.csv")) # load corrected hydro data with duplicates and date errors removed (V2 is version where contradictions in 'inundation' have been corrected)
Hydrodata$Unique_site_year=gsub("KRT","KT",Hydrodata$Unique_site_year)
Hydrodata$Unique_site_year=gsub("CCNT","NCT",Hydrodata$Unique_site_year)
Hydrodata$Unique_site_year=gsub("LHT","LHAT",Hydrodata$Unique_site_year)
Hydrodata$Unique_site_year=gsub("CHT","CCS",Hydrodata$Unique_site_year)
Hydrodata=Hydrodata[!duplicated(Hydrodata), ] # remove duplicates
Hydrodata=Hydrodata[,c(3,4,5,6,20)] # tidy up and remove extra columns not needed

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
Bigmod=data.frame(read.csv("Hydraulics_HTH_WL_medium and deep merged December 2020.csv")) # import metrics determined from Bigmod model
Bigmod$Unique_site_year=gsub("KRT","KT",Bigmod$Unique_site_year)
Bigmod$Unique_site_year=gsub("CCNT","NCT",Bigmod$Unique_site_year)
Bigmod$Unique_site_year=gsub("LHT","LHAT",Bigmod$Unique_site_year)
Bigmod$Unique_site_year=gsub("CHT","CCS",Bigmod$Unique_site_year)

Hydrodata=merge(Bigmod, Hydrodata, by="Unique_site_year", All.X=TRUE)

#Rainfall
data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes rainfall/"; setwd(data.dir)
Rainfall.dat=data.frame(read.csv("Rainfall_HTH_WL.csv")) # load rainfall data which is already sorted into dates x rainfall metrics
Rainfall.dat <- within(Rainfall.dat, Date <- as.Date(as.character(Date), format = "%Y-%m-%d")) # ensure dates are recognised in rainfall dat
Rainfall.dat=Rainfall.dat[!duplicated(Rainfall.dat), ] # remove duplicates
Rainfall.dat=Rainfall.dat[,c(2:10)]

#Temperature
data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes temperature/"; setwd(data.dir)
Temperature.dat=data.frame(read.csv("Temperature_HTH_WL.csv")) # load temperature data which is already sorted into dates x rainfall metrics
Temperature.dat <- within(Temperature.dat, Date <- as.Date(as.character(Date), format = "%Y-%m-%d")) # ensure dates are recognised in rainfall dat
Temperature.dat=Temperature.dat[,c(2:14)]
Temperature.dat=Temperature.dat[!duplicated(Temperature.dat), ] # remove duplicates

#Veg overstorey structure and location
data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes veg structure/"; setwd(data.dir)
VegStructure.dat=data.frame(read.csv("Hattah WL Veg Structure.csv")) # load veg structure data
VegStructure.dat<-VegStructure.dat[,c(2,14)] 
VegStructure.dat$Site.ID=gsub("CHT","CCS",VegStructure.dat$Site.ID)
VegStructure.dat$Site.ID=gsub("LHT","LHAT",VegStructure.dat$Site.ID)
VegStructure.dat$Site.ID=gsub("CCNT","NCT",VegStructure.dat$Site.ID)
VegStructure.dat$Site.ID=gsub("KRT","KT",VegStructure.dat$Site.ID)


#Merge various datasets
mydata_env=merge(newdates, Rainfall.dat, by.x="Date.of.collection", by.y="Date", all.x=TRUE) # merge species data with rainfall data by date to create a unique site-sample date by rainfall metrics table
mydata_env=merge(mydata_env, Hydrodata, by="Unique_site_year", all.x=TRUE) # I am not sure why but when I merge its duplicating selected lines of the database....
mydata_env=merge(mydata_env, Temperature.dat, by.x="Date.of.collection.x", by.y="Date",all.x=TRUE) # merge 
mydata_env=merge(mydata_env, VegStructure.dat, by.x="Site.ID.x", by.y="Site.ID",all.x=TRUE) # merge 

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
CTF.dat=data.frame(read.csv("Flood_HTH_WL_Cease2flowDec2020.csv")) # load CTF data
CTF.dat$Site.ID=gsub("CHT","CCS",CTF.dat$Site.ID)
CTF.dat$Site.ID=gsub("LHT","LHAT",CTF.dat$Site.ID)
CTF.dat$Unique_site_year=gsub("CHT","CCS",CTF.dat$Unique_site_year)
CTF.dat$Unique_site_year=gsub("LHT","LHAT",CTF.dat$Unique_site_year)

mydata_env=merge(mydata_env, CTF.dat, by="Unique_site_year", all.x=TRUE)

date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
write.csv(mydata_env, file = "Hattah Lakes wetlands transect based env data January 2021.csv") # save data out 

